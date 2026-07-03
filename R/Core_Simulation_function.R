#### Core functions of the simulation ####

#' Add Recombination
#'
#' @description
#' Function to add recombination to the haplotype information of the animals.
#'
#' @param map Mapfile of the population comprising of marker ID, chromosome number and marker position.
#' @param Haplotype Haplotype information of the population where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param recL Average number of recombinations to be introduced per animal.
#' @param nChr Number of haploid chromosomes (n) defined in the population.
#' @param recProb A vector of recombination probability between the markers. Its length should be total
#' number of markers minus 1. User can define the recombination probability between markers
#' helping in simulating a desired population.
#'
#' @returns Returns a list of recombination points and number of recombinations per animal.
#' @export
#'
#' @importFrom stats rpois
#'
#' @examples popG <- matrix(data = sample(c(0,1), 40, replace = TRUE), nrow = 4, ncol = 10)
#' popGmap <- generateMAP(nChr = 2, n_markers = 5, len_Chr = 10)
#' popGRec <- RecombinationPoint(map = popGmap, Haplotype = popG, recL = 2, nChr = 2)
RecombinationPoint <- function(map, Haplotype, recL, nChr, recProb) {
  if (missing(map))
    stop("Please define the map file")
  if (missing(Haplotype))
    stop("Please define the Haplotype information of the animals")
  if (!is.matrix(Haplotype))
    stop("Haplotype information is not in matrix form")
  if (missing(recL))
    stop("Please define the recL parameter")
  if (missing(nChr))
    stop("Please define the haploid number of chromosomes(N)")
  ## Chromosome structure
  chr_marker_counts <- as.vector(table(map$Chr))
  chr_starts <- c(1, base::cumsum(chr_marker_counts[-length(chr_marker_counts)]) +
                    1)
  chr_ends <- base::cumsum(chr_marker_counts)
  n_ind <- nrow(Haplotype) / 2
  n_markers <- ncol(Haplotype)
  ## Poisson-distributed intra-chromosomal crossovers
  nrec <- stats::rpois(n_ind, recL)
  ##Exclude boundary positions from crossover candidates
  boundary_set <- c(chr_starts, chr_ends)
  candidates <- base::setdiff(seq_len(n_markers), boundary_set)
  recList <- list()
  for (i in seq_len(n_ind)) {
    ## Step1:Independent assortment
    ## Each chromosome independently chooses strand(0 or 1)
    strand_choice <- base::sample(c(0L, 1L), nChr, replace = TRUE)
    ##Convert to swap points where strand flips
    assortment_swaps <- integer(0)
    ##If Chr1 is strand 1,swap at position 1
    if (strand_choice[1] == 1L) {
      assortment_swaps <- c(assortment_swaps, chr_starts[1])
    }
    ## For subsequent chromosomes :swap at start if
    ## different from previous chromosome
    for (k in 2:nChr) {
      if (strand_choice[k] != strand_choice[k - 1]) {
        assortment_swaps <- c(assortment_swaps, chr_starts[k])
      }
    }
    ## Step2:Intra-chromosomal crossovers
    if (nrec[i] > 0 && length(candidates) > 0) {
      n_draw <- min(nrec[i], length(candidates))
      if (missing(recProb)) {
        additional <- base::sample(candidates, n_draw, replace = FALSE)
      } else {
        cand_prob <- recProb[candidates]
        cand_prob <- cand_prob / sum(cand_prob)
        additional <- base::sample(candidates,
                             n_draw,
                             prob = cand_prob,
                             replace = FALSE)
      }
    } else {
      additional <- base::integer(0)
    }
    recpoints <- base::sort(base::unique(c(assortment_swaps, additional)))
    recList[[i]] <- recpoints
  }
  list(recList = recList, nrec = nrec)
}

# At the TOP of the file
.DISABLE_FUNCTIONS <- list(
  MutationRate = TRUE
)


#' Adding Mutation
#'
#' @description
#' Function to add mutation rate to the haplotype information of the animals.
#'
#' @param Haplotype Haplotype information of the population where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param rate Mutation rate. The default mutation rate is = 2.5*10^-5.
#' However, user can also define mutation rate as per requirement.
#'
#' @returns Returns a haplotype matrix with mutations.
#' @export
#'
#' @examples popH <- matrix(data = sample(c(0,1), 60, replace = TRUE), nrow = 6, ncol = 10)
#' popH_Mutation <- MutationRate(Haplotype = popH, rate = 2.5*10^-5)
MutationRate <- function(Haplotype, rate) {
  if (.DISABLE_FUNCTIONS$MutationRate) {
    return(Haplotype)
  }
  else if (!is.matrix(Haplotype)) {
    stop("Haplotype information is not in matrix format")
  }
  ### Calculate number of mutations according to Haplotype ###
  mutationNumber <- rate * nrow(Haplotype) * ncol(Haplotype)
  if (mutationNumber == 0L)
    return(Haplotype)
  ###Create mutation indices
  columnA <- base::sample(1:nrow(Haplotype), mutationNumber, replace = TRUE)
  columnB <- base::sample(1:ncol(Haplotype), mutationNumber, replace = TRUE)

  # Combine columns into a matrix
  mutation <- as.matrix(base::cbind(columnA, columnB))
  matrix_value <- c(0, 1)
  ## Haplotype
  if (all(Haplotype %in% matrix_value)) {
    Haplotype[mutation[, 1:2]] <- 1 - Haplotype[mutation[, 1:2]]
  }
  else{
    ## Coded haplotype
    Haplotype[mutation[, 1:2]] <- Haplotype[mutation[, 1:2]] * (-1)
  }
  return(Haplotype)
}


#' Select Sires
#'
#' @description
#' Function to select the sire.
#'
#'
#' @param Haplotype Haplotype information of the population where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param nSire Number of sires to be selected as parents for the simulation.
#' @param nProgeny Number of progeny to be simulated in the population.
#'
#' @returns Returns a list of sire IDs and haplotype information of the sire.
#' @export
#'
#' @examples popH <- matrix(data = sample(c(0,1), 60, replace = TRUE), nrow = 6, ncol = 10)
#' popH_Sire <- SireSample(Haplotype = popH, nSire = 1, nProgeny = 4)
SireSample = function(Haplotype, nSire, nProgeny) {
  if (missing(Haplotype) || missing(nSire) || missing(nProgeny)) {
    stop("Please define the all function parameters")
  }
  else if (!is.matrix(Haplotype)) {
    stop("Haplotype information is not in matrix format")
  }
  else if (nrow(Haplotype) %% 2 != 0) {
    stop("Haplotype matrix must have even rows")
  }
  else if (nSire > nProgeny) {
    stop("Please ensure that number of sires should be smaller than number of progeny")
  }
  else {
    ### Randomly select first index
    columnA <- as.matrix(base::sample(seq(1, nrow(Haplotype), by = 2), nSire, replace = FALSE))
    ### Get the second index
    columnB <- columnA + 1
    ## combined indices of the animals selected
    row_indices <- as.matrix(base::cbind(columnA, columnB))
    ### replicate randomly for the number of progeny
    replication_num <- as.integer(base::diff(c(0, base::sort(base::sample(
      1:(round(nProgeny) - 1), (nSire - 1)
    )), round(nProgeny))))
    replication_num <- as.vector(replication_num)
    # replicate the row indices with randomized replication
    rep_sire <- row_indices[rep(1:nrow(row_indices), replication_num), ]
    ## randomise sire
    random_sire_rep <- as.vector(t(rep_sire))
    SireMatrix <- Haplotype[random_sire_rep, ]
    SireID <- rownames(Haplotype)[random_sire_rep]
    return(list(Sire = SireID, Matrix = SireMatrix))
  }
}


#' Select Dams
#'
#' @description
#' Function to select the dams.
#'
#' @param Haplotype Haplotype information of the population where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param nDam Number of Dams to be selected as parents for the simulation.
#' @param nProgeny Number of progeny to be simulated in the population.
#'
#' @returns Return a list of dam IDs and haplotype information of the dam.
#' @export
#'
#' @examples pop_J <- matrix(data = sample(c(0,1), 60, replace = TRUE), nrow = 6, ncol = 10)
#' pop_J_dam <- DamSample(Haplotype = pop_J, nDam = 1, nProgeny = 4)
DamSample = function(Haplotype, nDam, nProgeny) {
  if (missing(Haplotype) || missing(nDam) || missing(nProgeny)) {
    stop("Please define the all function parameters")
  }
  else if (!is.matrix(Haplotype)) {
    stop("Haplotype information is not in matrix format")
  }
  else if (nrow(Haplotype) %% 2 != 0) {
    stop("Haplotype matrix must have even rows")
  }
  else if (nDam > nProgeny) {
    stop("Please ensure that number of dams should be smaller than number of progeny")
  }
  else {
    ### Randomly select first index
    columnA <- as.matrix(base::sample(seq(1, nrow(Haplotype), by = 2), nDam, replace = FALSE))
    ### Get the second index
    columnB <- columnA + 1
    ## combined indices of the animals selected
    row_indices <- as.matrix(base::cbind(columnA, columnB))
    ### replicate randomly for the number of progeny
    replication_num <- as.integer(base::diff(c(0, base::sort(base::sample(
      1:(round(nProgeny) - 1), (nDam - 1)
    )), round(nProgeny))))
    replication_num <- as.vector(replication_num)
    # replicate the row indices with randomized replication
    rep_dam <- row_indices[base::rep(1:nrow(row_indices), replication_num), ]
    ## randomise sire
    random_dam_rep <- as.vector(t(rep_dam))
    DamMatrix <- Haplotype[random_dam_rep, ]
    DamID <- rownames(Haplotype)[random_dam_rep]
    return(list(Dam = DamID, Matrix = DamMatrix))
  }
}


#' Function for mating sire and dam
#'
#' @description
#' Function to generate progeny from selected sires and dams.
#'
#' @param Sire Haplotype information of the Sire in matrix format.
#' @param Dam Haplotype information of the Dam in matrix format.
#'
#' @returns Returns haplotype information of the progeny.
#' @export
#'
#' @examples Sire_K <- matrix(data = sample(c(0,1), 30, replace = TRUE), nrow = 3, ncol = 10)
#' Dam_K <- matrix(data = sample(c(0,1), 30, replace = TRUE), nrow = 3, ncol = 10)
#' popK <- MatingPop(Sire = Sire_K, Dam = Dam_K)
MatingPop = function(Sire, Dam) {
  if (!is.matrix(Sire) || !is.matrix(Dam)) {
    stop("Both sire and dam haplotype information needs to be in matrix form")
  }
  else if (nrow(Sire) != nrow(Dam)) {
    stop("Unequal numbers of sire and dam haplotype information.
         Please check sire and dam information.")
  }
  #### get the next generation #####
  Progeny <- matrix(0, nrow = 2 * nrow(Sire), ncol =  ncol(Sire))
  Progeny[seq(1, nrow(Progeny), by = 2), ] <- Sire
  Progeny[seq(2, nrow(Progeny), by = 2), ] <- Dam
  return(Progeny)
}


#' Random selection of sire and dam gametes
#'
#' @param Gamete_Matrix A matrix comprising haplotype information of sire or dam gametes after recombination
#'
#' @returns Returns matrix comprising random selection of gametes
#' @export
#'
#' @examples mat = matrix(data = sample(c(-1,1,-2,2), 100*10, replace = TRUE)
#' , nrow = 100, ncol = 10)
#' mat = PopulationID(Population = mat, PopCode = "M")
#' rev_mat = Random_gametes(Gamete_Matrix = mat)
#'

Random_gametes = function(Gamete_Matrix) {
  if (!is.matrix(Gamete_Matrix)) {
    stop("Input file should be a matrix")
  }
  else if (nrow(Gamete_Matrix) %% 2 != 0) {
    stop("Haplotype matrix must have even rows")
  }
  n_gametes <- nrow(Gamete_Matrix) / 2
  sel_row_idx <- base::sapply(1:n_gametes, function(p) {
    base::sample(c(2 * p - 1, 2 * p), size = 1)
  })
  selected_gamete <- Gamete_Matrix[sel_row_idx, , drop = FALSE]
  return(selected_gamete)
}
