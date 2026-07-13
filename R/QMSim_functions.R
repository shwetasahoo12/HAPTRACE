##### Functions used for reading files from QMSim #####

#' Read QTL file obtained from QMSim
#'
#' @description
#' This function reads the effect file generated from QMSim with the file prefix "effect_qtl".
#' It generates the effects for the defined QTLs in QMSim.
#'
#' @param filename Requires the name of the QTL file with allele effect obtained from QMSim.
#'
#' @returns Returns a readable QTL effect
#' @seealso [HAPTRACE::CalQTLeffect()]
#' @export
#'
#' @examples
#' \donttest{
#' QTLfile <- ReadQTL(filename = "effect_qtl_001.txt")
#' }
#'
ReadQTL = function(filename) {
  if (missing(filename)) {
    stop("QTL file name is not defined")
  }
  RawQTLfile = base::readLines(filename)
  RawQTLfile = base::gsub(" +", " ", gsub("[0,1,2]:", "", RawQTLfile))
  RawQTLfile = base::strsplit(RawQTLfile, " ")
  RawQTLfile = lapply(RawQTLfile, function(x)
    x[1:4])
  CorQTLfile = do.call(rbind, RawQTLfile)
  CorQTLfile = CorQTLfile[-1, ]
  CorQTLfile = as.data.frame(CorQTLfile, stringsAsFactors = FALSE)
  CorQTLfile$CA1_A2 = as.numeric(CorQTLfile$V3) +  as.numeric(CorQTLfile$V4)
  return(CorQTLfile)
}

#' Calculate QTL effect from file generated from ReadQTL (For QMSim generated files)
#'
#' @description
#' This function requires the file with information on markerID, chromosome and marker information with prefix "lm_mrk_qtl",
#' total number of SNP markers ,and QTL effects generated from ReadQTL function. It generates
#' SNP effects of length of total number of markers along with QTLs.
#'
#'
#' @param filename Requires the name of the file with markerID, chr, and marker position.
#' @param nSNP Total number of SNP markers.
#' @param QTLfile QTL effects generated from ReadQTL function.
#'
#' @returns Returns a vector of the SNP effects with length equals to the total number of SNP markers.
#' @seealso [HAPTRACE::ReadQTL()]
#' @export
#'
#' @importFrom utils read.table
#'
#' @examples
#' \donttest{
#' popQM <- read_qmsim_geno("p1_mrk_qtl_001.txt")
#' QTL <- ReadQTL(filename = "effect_qtl_001.txt")
#' QTLeffect <- CalQTLeffect(filename = "lm_mrk_qtl_001.txt",
#' nSNP = ncol(popQM),QTLfile = QTL)
#' }
#'
#'
CalQTLeffect = function(filename, nSNP, QTLfile) {
  if (missing(filename)) {
    stop("QTL file generated from QMSim is not defined")
  }
  QTLposfile = utils::read.table(filename, header = TRUE, stringsAsFactors = FALSE)
  if (missing(nSNP)) {
    stop("Please define the number of markers (SNPs)")
  }

  TBVeffect = vector(length = nSNP)
  TBVeffect[] = 0
  QTLindex = as.matrix(grep("Q", QTLposfile$ID))
  if (missing(QTLfile)) {
    stop("Please define the number of QTL file")
  }
  TBVeffect[QTLindex] = QTLfile$CA1_A2
  TBVeffect = as.vector(as.numeric(TBVeffect))
  return(TBVeffect)
}

