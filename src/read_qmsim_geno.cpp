/**
 * @file read_qmsim_geno.cpp
 * @brief Fast C++ reader for QMSim p1_mrk_qtl_*.txt genotype files.
 *
 * QMSim format:
 *   - Line 1: header (skipped).
 *   - Subsequent lines: ID  <genotype_string>[>]
 *     where '>' means the genotype continues on the next line.
 *   - Each digit encodes a diploid genotype at one locus:
 *       0 = 0|0   (hom ref)
 *       1 = het   (unphased — phase assigned deterministically)
 *       2 = 1|1   (hom alt)
 *       3 = 0|1   (het, paternal=0, maternal=1)
 *       4 = 1|0   (het, paternal=1, maternal=0)
 *   - Number of loci = length(genotype_string).
 *
 * Two output modes controlled by @c geno_mode:
 *
 *   - "haplotype" : integer matrix, 2*n_ind rows x n_loci cols.
 *                   Row 2*i   = paternal haplotype (0/1).
 *                   Row 2*i+1 = maternal haplotype (0/1).
 *                   Both rows carry the animal ID as rowname.
 *                   Unphased het (code 1) is resolved deterministically
 *                   using a hash of (individual_index, locus_index, seed).
 *
 *   - "allele"    : integer matrix, n_ind rows x n_loci cols.
 *                   Cell values are the raw digits from the file (0-4).
 *                   Each row carries the animal ID as rowname.
 *
 * @author Generated for HAPTRACE R package
 */

#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdint>

/* ====================================================================== */
/*  Mix64 hash for deterministic phasing of code '1'                      */
/* ====================================================================== */

/**
 * @brief 64-bit finaliser (splitmix64 / Stafford variant).
 * @param x Input value.
 * @return Hashed value.
 */
static inline uint64_t Mix64(uint64_t x)
{
	x ^= x >> 30;
	x *= 0xbf58476d1ce4e5b9ULL;
	x ^= x >> 27;
	x *= 0x94d049bb133111ebULL;
	x ^= x >> 31;
	return x;
}

/**
 * @brief Deterministically assign phase for unphased het (code 1).
 * @param ind_idx  Individual index (0-based).
 * @param loc_idx  Locus index (0-based).
 * @param seed     User-supplied seed.
 * @return true if paternal=1, maternal=0 (i.e. code 4);
 *         false if paternal=0, maternal=1 (i.e. code 3).
 */
static inline bool PhaseHet(int ind_idx, int loc_idx, uint64_t seed)
{
	uint64_t h = seed;
	h ^= static_cast<uint64_t>(ind_idx) * 0x9e3779b97f4a7c15ULL;
	h ^= static_cast<uint64_t>(loc_idx) * 0x517cc1b727220a95ULL;
	return (Mix64(h) & 1ULL) != 0;
}

/* ====================================================================== */
/*  Helper: strip continuation marker                                     */
/* ====================================================================== */

/**
 * @brief Strip trailing whitespace and the continuation marker '>'.
 * @param s Input string (modified in place).
 * @return true if the line ended with '>' (continuation expected).
 */
static bool StripContinuation(std::string& s)
{
	while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r'))
	{
		s.pop_back();
	}
	if (!s.empty() && s.back() == '>')
	{
		s.pop_back();
		return true;
	}
	return false;
}

/* ====================================================================== */
/*  Main exported function                                                */
/* ====================================================================== */

/**
 * @brief Read a QMSim p1_mrk_qtl_*.txt file and return a single matrix.
 *
 * @param file_path  Path to the QMSim genotype file.
 * @param geno_mode  "haplotype" (default) or "allele".
 * @param phase_seed Seed for deterministic phasing of unphased het
 *                   (code 1).  Only used in "haplotype" mode.
 *                   Default: 42.
 * @return Rcpp::IntegerMatrix with rownames set to animal IDs.
 */
// [[Rcpp::export]]
Rcpp::IntegerMatrix read_qmsim_geno(const std::string& file_path,
                                    const std::string& geno_mode = "haplotype",
                                    int phase_seed = 42)
{
	/* ------------------------------------------------------------------ */
	/*  Validate mode                                                     */
	/* ------------------------------------------------------------------ */
	const bool mode_hap = (geno_mode == "haplotype");
	const bool mode_allele = (geno_mode == "allele");

	if (!mode_hap && !mode_allele)
	{
		Rcpp::stop(
			"Unknown geno_mode '%s'. "
			"Supported: \"haplotype\", \"allele\".",
			geno_mode.c_str());
	}

	const uint64_t seed = static_cast<uint64_t>(phase_seed);

	/* ------------------------------------------------------------------ */
	/*  Open file                                                         */
	/* ------------------------------------------------------------------ */
	std::ifstream ifs(file_path, std::ios::in);
	if (!ifs.is_open())
	{
		Rcpp::stop("Cannot open file: %s", file_path.c_str());
	}

	/* ------------------------------------------------------------------ */
	/*  Skip header                                                       */
	/* ------------------------------------------------------------------ */
	std::string line;
	if (!std::getline(ifs, line))
	{
		Rcpp::stop("File is empty: %s", file_path.c_str());
	}

	/* ------------------------------------------------------------------ */
	/*  Stream every individual                                           */
	/* ------------------------------------------------------------------ */
	std::vector<int> ids;
	std::vector<std::string> geno_strings;

	ids.reserve(4096);
	geno_strings.reserve(4096);

	while (std::getline(ifs, line))
	{
		/* Skip blank lines */
		if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos)
		{
			continue;
		}

		/* ---- parse ID + first genotype fragment ---------------------- */
		std::istringstream iss(line);
		int id = 0;
		std::string geno_part;
		if (!(iss >> id >> geno_part))
		{
			Rcpp::stop("Malformed line (cannot parse ID + genotype): %s", line.c_str());
		}

		/* ---- handle continuation '>' --------------------------------- */
		bool continues = StripContinuation(geno_part);
		std::string geno_full(geno_part);

		while (continues)
		{
			if (!std::getline(ifs, line))
			{
				Rcpp::stop("Unexpected end of file during continuation for ID %d", id);
			}
			geno_part = line;
			const std::size_t start = geno_part.find_first_not_of(" \t");
			if (start != std::string::npos)
			{
				geno_part = geno_part.substr(start);
			}
			continues = StripContinuation(geno_part);
			geno_full += geno_part;
		}

		ids.push_back(id);
		geno_strings.push_back(std::move(geno_full));
	}
	ifs.close();

	if (ids.empty())
	{
		Rcpp::stop("No individuals found in file: %s", file_path.c_str());
	}

	/* ------------------------------------------------------------------ */
	/*  Determine dimensions                                              */
	/* ------------------------------------------------------------------ */
	const std::size_t n_loci_sz = geno_strings[0].size();
	if (n_loci_sz == 0)
	{
		Rcpp::stop("Genotype string for ID %d is empty.", ids[0]);
	}
	const int n_loci = static_cast<int>(n_loci_sz);
	const int n_ind = static_cast<int>(ids.size());

	/* ------------------------------------------------------------------ */
	/*  MODE: haplotype                                                   */
	/*  2*n_ind rows  x  n_loci cols   (values 0/1)                       */
	/*                                                                    */
	/*  Decoding per digit:                                               */
	/*    0 -> pat=0, mat=0                                               */
	/*    1 -> unphased het, resolved by PhaseHet()                       */
	/*    2 -> pat=1, mat=1                                               */
	/*    3 -> pat=0, mat=1                                               */
	/*    4 -> pat=1, mat=0                                               */
	/* ------------------------------------------------------------------ */
	if (mode_hap)
	{
		const int n_rows = 2 * n_ind;
		Rcpp::IntegerMatrix mat(n_rows, n_loci);

		/* Build rownames: each ID appears twice */
		Rcpp::CharacterVector rnames(n_rows);
		for (int i = 0; i < n_ind; ++i)
		{
			const std::string sid = std::to_string(ids[i]);
			rnames[2 * i] = sid;
			rnames[2 * i + 1] = sid;
		}

		/* Fill */
		for (int i = 0; i < n_ind; ++i)
		{
			const std::string& gs = geno_strings[i];
			if (static_cast<int>(gs.size()) != n_loci)
			{
				Rcpp::stop(
					"ID %d has %zu loci but expected %d.",
					ids[i], gs.size(), n_loci);
			}

			const char* p = gs.data();
			const int row_pat = 2 * i;
			const int row_mat = 2 * i + 1;

			for (int j = 0; j < n_loci; ++j)
			{
				const char ch = p[j];

				switch (ch)
				{
					case '0': /* 0|0 */
						mat(row_pat, j) = 0;
						mat(row_mat, j) = 0;
						break;
					case '2': /* 1|1 */
						mat(row_pat, j) = 1;
						mat(row_mat, j) = 1;
						break;
					case '3': /* 0|1 */
						mat(row_pat, j) = 0;
						mat(row_mat, j) = 1;
						break;
					case '4': /* 1|0 */
						mat(row_pat, j) = 1;
						mat(row_mat, j) = 0;
						break;
					case '1': /* unphased het — deterministic assignment */
						if (PhaseHet(i, j, seed))
						{
							mat(row_pat, j) = 1;
							mat(row_mat, j) = 0;
						}
						else
						{
							mat(row_pat, j) = 0;
							mat(row_mat, j) = 1;
						}
						break;
					default:
						Rcpp::stop(
							"Invalid genotype code '%c' at locus %d for ID %d. "
							"Expected 0, 1, 2, 3, or 4.",
							ch, j + 1, ids[i]);
				}
			}
		}

		Rcpp::rownames(mat) = rnames;
		return mat;
	}

	/* ------------------------------------------------------------------ */
	/*  MODE: allele                                                      */
	/*  n_ind rows  x  n_loci cols   (raw digits from file: 0-4)          */
	/* ------------------------------------------------------------------ */
	{
		Rcpp::IntegerMatrix mat(n_ind, n_loci);

		/* Build rownames */
		Rcpp::CharacterVector rnames(n_ind);
		for (int i = 0; i < n_ind; ++i)
		{
			rnames[i] = std::to_string(ids[i]);
		}

		/* Fill */
		for (int i = 0; i < n_ind; ++i)
		{
			const std::string& gs = geno_strings[i];
			if (static_cast<int>(gs.size()) != n_loci)
			{
				Rcpp::stop(
					"ID %d has %zu loci but expected %d.",
					ids[i], gs.size(), n_loci);
			}

			const char* p = gs.data();
			for (int j = 0; j < n_loci; ++j)
			{
				const char ch = p[j];
				if (ch < '0' || ch > '4')
				{
					Rcpp::stop(
						"Invalid genotype code '%c' at locus %d for ID %d. "
						"Expected 0, 1, 2, 3, or 4.",
						ch, j + 1, ids[i]);
				}
				mat(i, j) = ch - '0';
			}
		}

		Rcpp::rownames(mat) = rnames;
		return mat;
	}
}
