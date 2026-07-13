##### Functions for extracting data from simulation study ####

#' Recode Haplotype with different breed codes
#'
#' @description
#' Function to recode the haplotype according to different breed codes as specified by users.
#'
#'
#' @param df A matrix consisting of haplotype information of the population where
#' rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param popBase Code specified according to the specific breed.
#'
#' @returns Returns a matrix consisting of coded haplotype information of the animals
#' according to the defined popBase.
#' @export
#'
#' @examples popA <- matrix(data = sample(c(0,1), 60, replace = TRUE), nrow = 6, ncol = 10)
#' PA_Haplo <- Recode_Haplotype(df = popA, popBase = 1)
Recode_Haplotype <- function(df, popBase) {
  matrix_value = c(0, 1)
  if (!all(df %in% matrix_value)) {
    stop(
      "Matrix consisting haplotype information contain values other than 0 and 1.
         Please check your marker file."
    )
  }
  if (!is.matrix(df)) {
    stop("Haplotype information is not in matrix form")
  }
  # Replace all old value with new value #
  df[df == 0] <- -popBase
  df[df == 1] <-  popBase
  return(df)
}

#' Recode Coded haplotypes according to breed of origin into 0 and 1 format
#'
#' @description
#' Function to recode the coded haplotype back to its original form.
#'
#'
#' @param df A matrix consisting of coded haplotype information of the animals
#' according to the breed of origin.
#'
#' @returns Returns a matrix consisting of a haplotype information in (0,1) format.
#' @export
#'
#' @examples popB <- matrix(data = sample(c(-1,1), 60, replace = TRUE), nrow = 6, ncol = 10)
#' coded_haplo <- Coded_to_Haplo(df = popB)
Coded_to_Haplo <- function(df) {
  # Replace all old value with new value #
  if (!is.matrix(df)) {
    stop("Coded Haplotype information is not in matrix form")
  }
  df[df < 0] <- 0
  df[df > 0] <- 1
  return(df)
}


#' Calculate Allele Frequency
#'
#' @description
#' Function to calculate allele frequency at each locus of a population.
#'
#' @param df A matrix consisting of a genotype information of the animals in (0,1,2) format.
#'
#' @returns Returns the value of allele frequency at each locus.
#' @export
#'
#' @examples popC <- matrix(data = sample(c(0,1,2), 60, replace = TRUE), nrow = 6, ncol = 10)
#' allelefreq <- AlleleFreq(df = popC)
AlleleFreq = function(df) {
  matrix_value = c(0, 1, 2)
  if (!all(df %in% matrix_value)) {
    stop(
      "Matrix consisting genotype information contain values other than 0, 1 and 2.
         Please check your marker file."
    )
  }
  if (!is.matrix(df)) {
    stop("Haplotype information is not in matrix form")
  }
  allelefrequency <- colMeans(df) / 2
  return(allelefrequency)
}

#' Make Genotype from Haplotype
#'
#' @description
#' Function to generate genotype information from haplotype (0,1) information.
#'
#' @param df Haplotype (0,1) information of the population where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#'
#' @returns Returns the genotype information of the animal.
#' @export
#'
#' @examples popD <- matrix(data = sample(c(0,1), 60, replace = TRUE), nrow = 6, ncol = 10)
#' genotype <- MakeGeno(df = popD)
MakeGeno = function(df) {
  ### Select odd and even rows if dataframe(df) is a matrix ###
  #G1 is selecting odd rows and G2 is selecting even rows
  matrix_value = c(0, 1)
  if (!all(df %in% matrix_value)) {
    stop(
      "Matrix consisting haplotype information contain values other than 0 and 1.
         Please check your marker file."
    )
  }
  if (!is.matrix(df)) {
    stop("Haplotype information is not in matrix form")
  }
  else if (base::nrow(df) %% 2 != 0) {
    stop("Haplotype matrix must have even rows")
  } else{
    G1 <- df[seq(1, nrow(df), by = 2), ]
    G2 <- df[seq(2, nrow(df), by = 2), ]
    genotype <- G1 + G2
    return(genotype)
  }
}


#' Calculation of True breeding Value (TBV)
#'
#' @description
#' Function to estimate true breeding value.
#'
#' @param Haplotype Haplotype (0,1) information of the population where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param Effect A vector comprising the SNP effects of the population.
#'
#' @returns Returns a haplotype information with TBV estimated from animal's haplotype information.
#' The TBV mentioned in "TBV" column provides final TBV (sum of TBV from haplotype of the animal) for an animal.
#' @seealso [HAPTRACE::Var_PhenoRecords()]
#' @export
#'
#' @examples popE <- matrix(data= sample(c(0,1), 16, replace= TRUE), nrow = 4, ncol = 4)
#' QTLEff <- c(0.2, 0.1, 0.3, 0.1)
#' TBV_haplotype <- TBV_Haplo(Haplotype = popE, Effect = QTLEff)
TBV_Haplo = function(Haplotype, Effect) {
  matrix_value = c(0, 1)
  if (!all(Haplotype %in% matrix_value)) {
    stop(
      "Matrix consisting haplotype information contain values other than 0 and 1.
         Please check your marker file."
    )
  }
  if (base::ncol(Haplotype) != length(Effect)) {
    stop("The number of SNP effects doesn't match the number of markers")
  }
  else if (!is.matrix(Haplotype) || !is.vector(Effect)) {
    stop("Haplotype file must be in matrix form and SNP effects must be in vector form")
  }
  else{
    TBV <- base::tcrossprod(Haplotype, t(Effect))
    TBV <- as.matrix(as.data.frame(TBV[seq(1, length(TBV), 2)] + TBV[seq(2, length(TBV), 2)], stringsAsFactors = FALSE))
    TBV <- as.matrix(base::rep(TBV, each = 2))
    colnames(TBV) <- 'TBV'
    Haplotype <- base::cbind(Haplotype, TBV)
    return(Haplotype)
  }
}

#' Generate variance components and phenotypic records
#'
#' @description
#' Function to generate variance components and phenotype records for the simulated population.
#'
#' @param Haplotype Haplotype information (0,1) of the population where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param h2 Heritability of the trait to be simulated.
#' @param trait_mean Mean of the trait to be simulated.
#' @param VarE Error Variance obtained after rescaling the QTL effects
#'
#' @returns Additive genetic variance, phenotypic variance, error variance and haplotype marker information with phenotypic records.
#' @seealso [HAPTRACE::TBV_Haplo()]
#' @export
#'
#' @importFrom stats rnorm var
#'
#' @examples popF <- matrix(data= sample(c(0,1), 16, replace= TRUE), nrow = 4, ncol = 4)
#' QTLEff <- c(0.2, 0.1, 0.3, 0.1)
#' TBV_haplotype <- TBV_Haplo(Haplotype = popF, Effect = QTLEff)
#' Var_Ph_R <- Var_PhenoRecords(Haplotype = TBV_haplotype, h2 = 0.3,
#' trait_mean = 20, VarE = 0.6)
Var_PhenoRecords = function(Haplotype, h2, trait_mean, VarE) {
  if (!is.matrix(Haplotype)) {
    stop("Provided haplotype data is not in matrix format")
  }
  else if (!("TBV" %in% colnames(Haplotype))) {
    stop("Haplotype data must contain a true breeding value 'TBV' column.")
  }
  else if (missing(h2)) {
    stop("Please define heritability of the trait to be simulated")
  }
  else if (h2 <= 0) {
    stop(
      "The heritability of the trait less than or equal to zero! Check the heritability again."
    )
  }
  else if (h2 == 1) {
    warning("The heritability of the trait is 1! Please proceed with caution")
  }
  else if (h2 > 1) {
    stop("Heritability of the trait cannot exceed 1")
  }
  else if (missing(trait_mean)) {
    stop("Please define the mean of the trait")
  }
  else if (missing(VarE)) {
    stop("Please define the error variance (VarE)")
  }
  else {
    AddVar <- stats::var(Haplotype[, 'TBV'])
    ErrorVar <- VarE
    PhenoVar <- AddVar + ErrorVar
    error <- as.matrix(stats::rnorm(
      nrow(Haplotype) / 2,
      mean = 0,
      sd = base::sqrt(ErrorVar)
    ))
    error <- as.matrix(base::rep(error, each = 2))
    phenotype <- trait_mean + Haplotype[, 'TBV'] + error
    colnames(phenotype) <- 'phenotype'
    Haplotype <- base::cbind(Haplotype, phenotype)
    return(
      list(
        VarAdd = AddVar,
        VarPheno = PhenoVar,
        VarError = ErrorVar,
        Haplo_pheno = Haplotype
      )
    )
  }
}

#' Generate MAP file
#'
#' @description
#' Function to generate MAP file with marker ID, chromosome and marker position columns.
#'
#' @param nChr Number of haploid chromosomes (n) defined in the population.
#' @param n_markers Number of SNP markers to be generated per chromosome.
#' This parameter is independent of QTL effects at this stage.
#' @param len_Chr Length of the chromosome. It can be a single value or vector specifying different length per chromosome in cM.
#' If users wants to use base pairs(bp) as a measure to determine position, it can be obtained by bp = cM *10^6 .The base pair
#' position should be calculated by users as it is not part of this function.
#'
#'
#' @returns Returns a map file of the population with marker ID, chromosome number and position of the markers on the chromosome.
#' @seealso [HAPTRACE::QTLeffects()], [HAPTRACE::mapQTLfile()]
#' @export
#'
#' @examples map <- generateMAP(nChr = 2, n_markers = 10, len_Chr = 5)
generateMAP = function(nChr, n_markers, len_Chr) {
  if (missing(nChr)) {
    stop("Please define the number of chromosomes")
  }
  else if (nChr > 50) {
    stop("Haploid number of chromosomes (N) exceeded 50. Please check and try again.")
  }
  else if (missing(n_markers)) {
    stop("Please define the number of markers")
  }
  else if (missing(len_Chr)) {
    stop("Please define the length of the chromosome")
  }

  MarkerID <- as.matrix(paste("M", 1:(n_markers * nChr), sep = ""))
  Chr <- as.matrix(base::sort(base::rep(1:nChr, n_markers)))
  Positionlist <- list()
  if (length(len_Chr) == 1) {
    for (i in 1:nChr) {
      Positionlist[[i]] <- base::lapply(len_Chr, function(x)
        base::sort(stats::runif(
          n_markers, min = 0, max = x
        )))
    }
  }
  else if (is.vector(len_Chr)) {
    Positionlist <- base::lapply(len_Chr, function(x)
      base::sort(stats::runif(
        n_markers, min = 0, max = x
      )))
  }
  else {
    stop("The length of the chromosomes should be defined in vector form")
  }
  Position <- round(as.matrix(unlist(Positionlist)), 8)
  MAP <- base::cbind(MarkerID, Chr, Position)
  colnames(MAP) <- c("ID", "Chr", "Position")
  MAP <- base::as.data.frame(MAP, stringsAsFactors = FALSE)
  MAP$Chr <- as.numeric(MAP$Chr)
  MAP$Position <- as.numeric(MAP$Position)
  return(MAP)
}


#' Calculation QTL effects
#'
#' @description
#' Function to estimate QTL effects with output of QTL ID, chromosome, QTL position and QTL effects
#'
#' @param n_QTL_Chr Number of QTL to be simulated per chromosome
#' @param nChr Number of haploid chromosomes (n) defined in the population.
#' @param len_Chr Length of the chromosome.It can be either provided single constant value or vector specifying varied length of the chromosomes.
#' @param shape_gamma Shape of the QTL effects if simulated using gamma distribution.
#' @param effect_distribution Define the distribution of QTL effects. It can be either gamma, normal or uniform distribution.
#'
#' @returns file comprising of QTL ID, chromosomes, position of the QTL and QTL effects.
#' @seealso [HAPTRACE::generateMAP()], [HAPTRACE::mapQTLfile()]
#' @export
#'
#' @importFrom stats runif rgamma rnorm var
#'
#' @examples QTL <- QTLeffects(n_QTL_Chr = 1, nChr = 2, len_Chr = 10,
#' shape_gamma = 0.4, effect_distribution = "gamma")
QTLeffects = function(n_QTL_Chr,
                      nChr,
                      len_Chr,
                      shape_gamma = 0.4,
                      effect_distribution = "gamma") {
  if (missing(n_QTL_Chr)) {
    stop("Please define the number of QTLs to be simulated per chromosome")
  }
  else if (missing(nChr)) {
    stop(
      "Define the number of chromosomes. Please ensure that haploid number of chromosomes (N) are provided."
    )
  }
  else if (nChr > 50) {
    stop("Haploid number of chromosomes (N) exceeded 50. Please check and try again.")
  }
  else if (missing(len_Chr)) {
    stop("Please define the length of the chromosome in cM")
  }
  distribution <- c("gamma", "normal", "uniform")
  if (!effect_distribution %in% distribution) {
    stop("The effect distribution must be one of : ",
         paste(distribution, collapse = ", "))
  }
  QTLID <- as.matrix(paste("Q", 1:(n_QTL_Chr * nChr), sep = ""))
  Chr <- as.matrix(base::sort(base::rep(1:nChr, n_QTL_Chr)))
  QTLpos <- list()
  if (length(len_Chr) == 1) {
    for (i in 1:nChr) {
      QTLpos[[i]] <- base::lapply(len_Chr, function(x)
        base::sort(stats::runif(
          n_QTL_Chr, min = 0, max = x
        )))
    }
  }
  else if (is.vector(len_Chr)) {
    QTLpos <- base::lapply(len_Chr, function(x)
      base::sort(stats::runif(
        n_QTL_Chr, min = 0, max = x
      )))
  }
  else {
    stop("The length of the chromosomes should be defined in vector form")
  }
  QTLpos_Chr <- round(as.matrix(unlist(QTLpos)), 8)
  ### QTL effects
  if (effect_distribution == "gamma") {
    QTLeffI <- as.matrix(stats::rgamma(n_QTL_Chr * nChr, shape = shape_gamma))
    Neg <- matrix(
      data = base::sample(c(-1, 1), size = n_QTL_Chr * nChr, replace = TRUE),
      nrow = n_QTL_Chr * nChr,
      ncol = 1
    )
    QTLEffR <- Neg * (QTLeffI)
  }
  else if (effect_distribution == "normal") {
    QTLEffR <- as.matrix(stats::rnorm(n_QTL_Chr * nChr, mean = 0, sd = 1))
  }
  else if (effect_distribution == "uniform") {
    QTLeffI <- as.matrix(stats::runif(n_QTL_Chr * nChr, min = 0, max = 1))
    Neg <- matrix(
      data = base::sample(c(-1, 1), size = n_QTL_Chr * nChr, replace = TRUE),
      nrow = n_QTL_Chr * nChr,
      ncol = 1
    )
    QTLEffR <- Neg * (QTLeffI)
  }
  QTL <- cbind(QTLID, Chr , QTLpos_Chr, QTLEffR)
  colnames(QTL) <- c("ID", "Chr", "Position", "Effect")
  QTL <- as.matrix(QTL)
  return(QTL)
}



#' Calculate True Breed Proportion with coded haplotype
#'
#' @description
#' Function to calculate true breed proportion in an animal with coded haplotype information.
#'
#' @param Haplotype A matrix with coded haplotype information of the population where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#'
#' @returns Returns true breed proportion of population with coded haplotype.
#' @seealso [HAPTRACE::Admixplot()]
#' @export
#'
#' @examples P1 <- matrix(data = sample(c(-1, 1, 2, -2), 20, replace = TRUE), nrow = 4, ncol = 5)
#' TBP <- trueBreedComp(Haplotype = P1)
trueBreedComp <- function(Haplotype)
{
  if (!is.matrix(Haplotype)) {
    stop(
      "Haplotype is not in matrix. Please ensure that the haplotype is coded specific to the population"
    )
  } else if (nrow(Haplotype) %% 2 != 0) {
    stop("Haplotype matrix must have even rows")
  }
  else {
    BP <- base::abs(Haplotype)
    levels <- base::unique(as.vector(BP))
    levels <- levels[levels != 0]
    BP1 <- base::apply(BP, 1, function(x)
      base::table(factor(x, levels)))
    BP2 <- BP1 / colSums(BP1, na.rm = TRUE)
    TrueBPM <- (BP2[, seq(1, ncol(BP2), 2)] + BP2[, seq(2, ncol(BP2), 2)]) /
      2
    rownames(TrueBPM) <- paste("P", rownames(TrueBPM), sep = "")
    return(TrueBPM)
  }
}

#' Create combined map and qtl file
#'
#' @description
#' This function combines marker and qtl file together. SNP markers are often observed to be near or away from the QTL in real life scenario.
#' However, QTL position is included along with map file to get combined marker file for simulation purpose in this package.
#' This step is crucial to calculate the true breeding value and simulate the phenotype with desired population parameters.
#'
#' @param map map file with marker ID, chromosome, and marker position
#' @param QTL qtl file with QTL ID, chromosome, marker position and effect
#'
#' @returns combined file with marker and QTL ID, chromosome, position and QTL effects.
#' @seealso [HAPTRACE::generateMAP()], [HAPTRACE::QTLeffects()]
#' @export
#'
#' @examples popmap <- generateMAP(nChr = 5, n_markers = 10,
#' len_Chr = c(10, 20, 40, 50, 60))
#' popQTL <- QTLeffects(n_QTL_Chr = 2, nChr = 5, len_Chr = c(10, 20, 40, 50, 60),
#' shape_gamma = 0.4, effect_distribution = "gamma")
#' pop_mapqtl <- mapQTLfile(map = popmap, QTL = popQTL)
mapQTLfile <- function(map, QTL) {
  if (missing(map)) {
    stop("Please define the map file")
  }
  else if (missing(QTL)) {
    stop("Please define the QTL effect")
  }
  else {
    map$Effect <- 0
    map <- map[, c("ID", "Chr", "Position", "Effect")]
    map$Chr <- as.numeric(map$Chr)
    map$Position <- as.numeric(map$Position)
    map$Effect <- as.numeric(map$Effect)
    QTL <- base::as.data.frame(QTL, stringsAsFactors = FALSE)
    QTL$Chr <- as.numeric(QTL$Chr)
    QTL$Position <- as.numeric(QTL$Position)
    QTL$Effect <- as.numeric(QTL$Effect)
    if (length(unique(map$Chr)) != length(unique(QTL$Chr))) {
      stop("Number of chromosome in map and QTL file should be same")
    }
    map_qtl <- rbind(map, QTL)
    combined_map_qtl <- map_qtl[base::order(map_qtl$Chr, map_qtl$Position), ]
    rownames(combined_map_qtl) <- NULL
    combined_map_qtl$Effect <- as.numeric(combined_map_qtl$Effect)
    return(combined_map_qtl)
  }
}


#' Generate Historical Population
#'
#' @description
#' Function to simulate historical population.
#'
#' @param n_generations Number of generations to be simulated.
#' @param n_animals Number of animals to be simulated.
#' @param map Map file for the population.
#' @param nSire Number of sires to be selected as parents for the simulation.
#' @param nDam Number of Dams to be selected as parents for the simulation.
#' @param nProgeny Number of progeny to be simulated in the population.
#' @param mutationRate Mutation rate to added to the simulated population.
#' @param SelType Selection Type. It can be based on phenotype, random or true breeding value. The default selection type is random.
#' For phenotypic selection, use "Pheno" and for true breeding value, use "TBV".
#' @param Effect SNP effect for calculation of TBV.
#' @param h2 Heritability of the trait to be simulated.
#' @param trait_mean Mean of the trait to be simulated.
#' @param VarE Error Variance obtained after rescaling the QTL effects.
#' @param recL Average number of recombinations to be introduced per animal.
#' @param nChr Number of haploid chromosomes (n) defined in the population.
#' @param recProb A vector of recombination probability between the markers. Its length should be total
#' number of markers minus 1. User can define the recombination probability between markers
#' helping in simulating a desired population.
#' @param minLength Minimum length between recombination points beyond which recombination is allowed. The default value is 0.
#' If a user define minimum length, it will ensure that difference between two recombination points
#' less than minimum length won't allow recombination.
#' @param IndStartVal Starting value for the Animal ID of the historical population.
#' @param prefixID Prefix to be added to the ID of the animal of the historical population
#'
#' @returns Returns haplotype information of the historical population.
#' @seealso [HAPTRACE::MutationRate()], [HAPTRACE::mapQTLfile()], [HAPTRACE::generateMAP()], [HAPTRACE::QTLeffects()]
#'
#' @export
#'
#' @examples HP_map <- generateMAP(nChr = 2, n_markers = 5, len_Chr = 5)
#' QTL_HP <- QTLeffects(n_QTL_Chr = 2, nChr = 2, len_Chr = 5, shape_gamma = 0.4,
#' effect_distribution = "gamma")
#' HP_mapqtl <- mapQTLfile(map = HP_map, QTL = QTL_HP)
#' HP_effect <- as.numeric(HP_mapqtl$Effect)
#' H1 <- generateHP(n_generations = 5, n_animals = 20, map = HP_mapqtl,
#' nSire = 1, nDam = 1, nProgeny = 6, mutationRate = 2.5 * 10^-5, SelType = "random",
#' Effect = HP_effect, h2 = 0.3, trait_mean = 10, VarE = 0.6,recL = 2,nChr = 2, minLength = 0)
#'
generateHP <- function(n_generations,
                       n_animals,
                       map,
                       nSire,
                       nDam,
                       nProgeny,
                       mutationRate = 2.5 * 10^-5,
                       SelType = "random",
                       Effect,
                       h2 = 0.3,
                       trait_mean = 10,
                       VarE,
                       recL = 5,
                       nChr,
                       recProb,
                       minLength = 0,
                       IndStartVal = 1,
                       prefixID = "H") {
  if (missing(n_generations)) {
    stop("Please define the number of generations to be simulated")
  }
  else if (missing(n_animals)) {
    stop("Please define the number of animals to be simulated per generation")
  }
  else if (missing(map)) {
    stop("Please define the map file")
  }
  else if (missing(nSire)) {
    stop("Please define the number of sires to be used for simulation")
  }
  else if (missing(nDam)) {
    stop("Please define the number of dams to be used for simulation")
  }
  else if (missing(nProgeny)) {
    stop("Please define the number of progeny to be simulated per generation")
  }
  else if (nSire > nProgeny) {
    stop("Please ensure that number of sires should be smaller than number of progeny")
  }
  else if (nDam > nProgeny) {
    stop("Please ensure that number of dams should be smaller than number of progeny")
  }
  else if (missing(Effect)) {
    stop("Please define the SNP effect")
  }
  else if (missing(h2)) {
    stop("Please define the heritability of the trait to be simulated")
  }
  else if (h2 == 1) {
    warning("The heritability of the trait is 1! Proceed with caution!")
  }
  else if (h2 > 1) {
    stop("The heritability of the trait is more than 1! Check the heritability again.")
  }
  else if (missing(trait_mean)) {
    stop("Please define the mean of the trait")
  }
  else if (missing(VarE)) {
    stop("Please define the error variance of the trait")
  }
  else if (missing(recL)) {
    stop("Please define this parameter")
  }
  else if (missing(nChr)) {
    stop(
      "Define the number of chromosomes. Please ensure that haploid number of chromosomes (N) are provided"
    )
  }
  else if (nChr > 50) {
    stop("Haploid number of chromosomes (N) exceeded 50. Please check and try again.")
  }
  else {
    generation <- 1
    n_markers <- length(map$Chr)
    if (length(Effect) != n_markers) {
      stop(
        "The length of the effect ",
        length(Effect),
        " must be equal to the number of markers ",
        n_markers,
        "."
      )
    }
    HP <- matrix(
      data = base::sample(
        c(0, 1),
        2 * n_animals * n_markers,
        replace = TRUE,
        prob = c(0.5, 0.5)
      ),
      ncol = n_markers,
      nrow = 2 * n_animals
    )
    colnames(HP) <- paste("M", 1:ncol(HP), sep = "")
    rownames(HP) <- paste("H", rep(1:(nrow(HP) / 2), each = 2), sep = "_")
    PopInfo <- list()
    recSire <- list()
    recDam <- list()
    OnlyF1 <- list()
    F1Haplotype <- list()
    Parents_list <- list()
    for (generation in 1:n_generations) {
      message("Working on Generation: ", generation)
      Newpop <- PopSelect(
        map,
        HP,
        nSire,
        nDam,
        mutationRate,
        nProgeny,
        SelType,
        Effect,
        h2,
        trait_mean,
        VarE,
        recL,
        nChr,
        recProb,
        minLength
      )
      PopInfodata <- Newpop$PopInfo
      recpointSire <- Newpop$rec_sire
      recpointDam <- Newpop$rec_dam

      ##New population
      OnlyF1 <- MatingPop(Newpop$Sire, Newpop$Dam)

      ###Get Sire and Dam
      #### Parents in next generation ####
      Parents <- data.frame(
        AnimID = paste(prefixID, IndStartVal:(IndStartVal + nrow(OnlyF1) / 2 - 1), sep =
                         "_"),
        SireID = Newpop$SireID ,
        DamID = Newpop$DamID,
        Sex = as.numeric(c(rep(
          1, nrow(OnlyF1) / 4
        ), rep(
          2, nrow(OnlyF1) / 4
        )))
      )

      ##### Recent population is foundation for next generation
      rownames(OnlyF1) <- rep(Parents$AnimID, each = 2)
      populationSim <- OnlyF1
      F1Haplotype[[generation]] <- OnlyF1
      PopInfo[[generation]] <- PopInfodata
      recSire[[generation]] <- recpointSire
      recDam[[generation]] <- recpointDam
      Parents_list[[generation]] <- Parents
      IndStartVal <- (IndStartVal + nrow(OnlyF1) / 2 - 1) + 1
    }
    return(
      list(
        Haplotype = do.call(rbind, F1Haplotype),
        parents = do.call(rbind, Parents_list),
        population = do.call(rbind, PopInfo),
        recpointS = do.call(rbind, recSire),
        recpointD = do.call(rbind, recDam)
      )
    )
  }
}


#' Generate PLINK pedigree file
#'
#' @description
#' Function to generate ped file for PLINK analysis.
#'
#' @param Genotype A matrix containing the genotype information of the animals where rows
#' @param AnimalID Animal ID for the animals in the population
#' @param filename Naming a ped file
#'
#' @returns A ped file for PLINK analysis.
#' As of now, it only saves animal ID and genotypes.
#' Saving pedigree in ped file is in our plan for the update.
#' @export
#'
#' @importFrom utils write.table
#' @importFrom data.table fwrite as.data.table
#'
#' @examples
#' \donttest{
#' pop_plink <- matrix(data = sample(c(0,1,2), 60, replace = TRUE), nrow = 6, ncol = 10)
#' rownames(pop_plink) <- paste("PP", 1:nrow(pop_plink), sep = "_")
#' output_path <- file.path(tempdir(), "pop")
#' plinkped(Genotype = pop_plink, AnimalID = rownames(pop_plink), filename = output_path)
#' }
#'
#'
plinkped <- function (Genotype, AnimalID, filename)
{
  if (missing(Genotype)) {
    stop("Genotype information is not defined")
  }
  else if (!is.matrix(Genotype)) {
    stop("Please provide genotype information in matrix form")
  }
  else if (missing(AnimalID)) {
    stop("Please define AnimalID")
  }
  else if (missing(filename)) {
    stop("Please specify the filename for ped file")
  }
  else {
    G1 <- Genotype
    G2 <- Genotype
    G1[G1 == 0] <- 1
    G2[G2 == 0] <- 5
    G2[G2 == 1] <- 6
    G2[G2 == 2] <- 7
    G2[G2 == 5] <- 1
    G2[G2 == 6] <- 2
    G2[G2 == 7] <- 2
    T1 <- matrix(0, nrow = nrow(G1), ncol = 2 * ncol(G1))
    T1[, seq(1, ncol(T1), by = 2)] <- G1
    T1[, seq(2, ncol(T1), by = 2)] <- G2
    rownames(T1) <- AnimalID
    T1 <- as.data.table(T1, keep.rownames = TRUE)
    fwrite(
      T1,
      file = paste0(filename, ".ped"),
      sep = " ",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
}


#' Format PLINK MAP file from QMSim
#'
#' @description
#' Function to format map file from QMSim for PLINK analysis.
#'
#'
#' @param filepath Map file generated from the QMSim.
#' @param outputpath Path directory for the output file
#' @param filename Name of the map file to be extracted
#'
#' @returns Returns a map file that can be used for PLINK analysis.
#' @export
#'
#' @importFrom data.table fread fwrite
#'
#' @examples
#' \dontrun{
#' output_folder <- file.path(tempdir(), "map")
#' dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
#' QMSim_PlinkMAP(filepath = "lm_mrk_qtl_001.txt",
#' outputpath = output_folder, filename = "pop")
#' }
#'
QMSim_PlinkMAP = function(filepath, outputpath, filename) {
  if (missing(filepath)) {
    stop("Please specify the file path where map file is saved")
  }
  MAP <- data.table::fread(filepath)
  MAP$PositionR <- 0
  MAP$Physical_pos <- (MAP$Position) * 1000000

  # Combine output_path and filename with .txt extension
  output_path_file <- file.path(outputpath, paste0(filename, ".map"))

  ###Reposition the columns now in map file ###
  MAP_REVISED <- MAP[, c('Chr', 'ID', 'PositionR', 'Physical_pos')]
  data.table::fwrite(
    MAP_REVISED,
    file = output_path_file,
    sep = " ",
    row.names = F,
    quote = F,
    col.names = F
  )
}


#' Adding Population Code before simulating the population
#'
#' @description
#' Function to add Population code to the population which will help in generating parent IDs.
#'
#' @param Population Haplotype information of the base population to be used for simulating the population.
#' @param PopCode Code to be assigned for the base population which will appear in IDs.
#'
#' @returns Returns a population with provided PopCode, which will appear as PopulationID.
#' @export
#'
#' @examples PopN_Check <- matrix(data = sample(c(0,1), 200, replace = TRUE), nrow = 10, ncol = 20)
#' PopN_CheckID <- PopulationID(Population = PopN_Check, PopCode = "P1")
PopulationID = function(Population, PopCode) {
  if (missing(Population)) {
    stop("Population information is missing")
  }
  else if (!is.matrix(Population)) {
    stop("Please specify the Haplotype information of the population")
  }
  else if (missing(PopCode)) {
    stop("Please specify the PopCode")
  }
  rownames(Population) <- paste(PopCode, rep(1:(nrow(Population) / 2), each =
                                               2), sep = "_")
  return(Population)
}


#' Reformat genotype data for BLUPF90 analysis
#'
#' @description
#' Function to format genotype file for BLUPF90 analysis.
#'
#' @param Genodata A matrix comprising genotype information (0,1,2) of the animals where rownames depicts the animal ID.
#' @param savefile File name in which the genotype data will be saved.
#'
#' @returns Returns a text file comprising of genotype data for its use in BLUPF90.
#' @importFrom gdata write.fwf
#'
#' @export
#'
#' @examples
#' \donttest{
#' Geno <- matrix(data = sample(c(0,1,2), 100, replace = TRUE), nrow = 10, ncol = 10)
#' rownames(Geno) <- paste("G", 1:nrow(Geno), sep = "_")
#' output_path <- file.path(tempdir(), "genotype")
#' Reformat_BLUPF90(Genodata = Geno, savefile = output_path)
#' }
#'
Reformat_BLUPF90 = function(Genodata, savefile) {
  if (missing(Genodata)) {
    stop("Genotype information is missing")
  } else if (!is.matrix(Genodata)) {
    stop("Genotype information is not in matrix form")
  }
  else if (missing(savefile)) {
    stop("Please specify the file name in which genotype file needs to be saved")
  }
  else {
    GenoID <-  base::rownames(Genodata)
    Geno <- data.frame(GenoID = GenoID, Genodata)
    # Combine
    Geno$Genodata <- base::apply(Geno[, -1], 1, function(x)
      paste(x, collapse = ""))
    final_Geno <- Geno[, c("GenoID", "Genodata")]
    # Write fixed width file
    gdata::write.fwf(
      final_Geno,
      file = paste0(savefile, ".txt"),
      sep = " ",
      rownames = FALSE,
      colnames = FALSE
    )
  }
}
