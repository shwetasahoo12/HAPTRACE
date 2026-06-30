#' True Breeding Value (TBV) estimation according to breed of origin
#'
#' @description
#' Function to estimate TBV for animals having different SNP effects according to breed of origin.
#'
#' @param Haplotype A matrix comprising a coded haplotype information of the population where rows are number of animals
#' and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param Breed_Effect A matrix comprising the breed effect defined according to the breed of origin
#' such that row depicts the breeds and columns depicts the SNP effects.
#'
#' @returns Returns a matrix of coded Haplotype with TBV estimated according to the breed of origin.
#' @export
#'
#' @examples pop_BH <- matrix(sample(c(-1,1,-2,2,-3,3), 16, replace = TRUE), nrow = 4, byrow = TRUE)
#' Breed_Eff <- matrix(runif(4 * 4, 0, 1), nrow = 4, byrow = TRUE)
#' Breed_Hap <- BreedEff_Haplo(Haplotype = pop_BH, Breed_Effect = Breed_Eff)
#'
BreedEff_Haplo = function(Haplotype, Breed_Effect) {
  if (!is.matrix(Haplotype)) {
    stop("Provided haplotype data is not in matrix format")
  } else if (!is.matrix(Breed_Effect)) {
    stop("Provided Breed-specific effect is not in matrix format. Check the structure of the breed effect.")
  }
  Haplotype_R <- Haplotype
  Pos_Substitution <- (Haplotype > 0 &
                         Haplotype <= nrow(Breed_Effect))
  Index = base::cbind(Haplotype[Pos_Substitution], col(Haplotype)[Pos_Substitution])
  Haplotype_R[Pos_Substitution] <- Breed_Effect[Index]
  Haplotype_R[!Pos_Substitution] <- 0
  TBV = base::rowSums(Haplotype_R)
  TBV = as.matrix(as.data.frame(TBV[seq(1, length(TBV), 2)] + TBV[seq(2, length(TBV), 2)], stringsAsFactors = FALSE))
  TBV = as.matrix(base::rep(TBV, each = 2))
  colnames(TBV) = 'TBV'
  Haplotype = base::cbind(Haplotype, TBV)
  return(Haplotype)
}


#' Selecting Population with breed effect according to breed of origin
#'
#' @description
#' Function to simulate a population with different breed effect according to breed of origin.
#'
#' @param map Map file for the population comprising of marker ID, chromosome number and marker position.
#' @param Haplotype A matrix with coded haplotype information of the population where rows are number of animals
#' and columns are number of markers.
#' In case of haplotype, two rows provide information on genomic material of an animal.
#' @param nSire Number of sires to be used for simulation.
#' @param nDam Number of dams to be used for simulation.
#' @param rate Mutation rate to introduce mutation in the simulated population.
#' @param nProgeny Number of progeny to be simulated in each generation.
#' @param SelType A character value specifying the selection criteria of the simulated population.
#' Selection criteria can be based on phenotype, random or true breeding value. The default selection type is "random".
#' For phenotypic selection, use "Pheno" and for true breeding value, use "TBV".
#' @param Breed_Effect A matrix comprising breed effect according to the breed of origin where rows contain information
#' on SNP effects.
#' @param h2 A numeric value describing the heritability of the trait to be simulated ranging from 0 to 1.
#' @param trait_mean A positive integer value describing mean of the phenotype of the trait to be simulated.
#' @param VarE A positive integer value depicting error variance obtained after rescaling the QTL effects.
#' @param recL A numeric value describing the average number of recombinations to be introduced per animal.
#' @param nChr A numeric value describing the number of haploid chromosomes (N) defined in the simulated population.
#' @param recProb A vector of recombination probability between the markers. Its length should be total
#' number of markers minus 1. User can define the recombination probability between markers
#' helping in simulating a desired population.
#' @param minLength Minimum length between recombination points beyond which recombination is allowed. The default value is 0.
#' If a user define minimum length, it will ensure that difference between two recombination points
#' less than minimum length won't allow recombination.
#'
#'
#' @returns Returns selected sire and dams along with with their ID, recombination points of sire and dam, and population coded haplotype with phenotype and TBV.
#' @seealso [HAPTRACE::BreedEff_Haplo()]
#'
#' @export
#'
#' @importFrom hsphase addSwitch
#'
#' @examples PopS <- matrix(data = sample(c(-1,1, -2,2), 200, replace = TRUE),
#' nrow = 20, ncol = 10)
#' PopS_map <- generateMAP(2, 5, 5)
#' PopS_QTL <-  matrix(c(-0.001, 0.0003, 0, 0, 0, 0.0032, 0.0023, -0.0012,0, 0,
#'  -0.002, 0.004, 0, 0, 0.002, 0, 0.002, 0,0, -0.0012), byrow = 2, nrow = 2, ncol = 10)
#' Selpop <- BreedPopSelect(map = PopS_map, Haplotype = PopS, nSire = 2, nDam = 2, rate = 2.5 * 10^-5,
#' nProgeny = 6, SelType = "random", Breed_Effect = PopS_QTL, trait_mean = 5,VarE = 0.6,
#' h2 = 0.2, recL = 1, nChr = 2, minLength = 0)
#'
BreedPopSelect = function(map,
                          Haplotype,
                          nSire,
                          nDam,
                          rate = 2.5 * 10^-5,
                          nProgeny,
                          SelType = "random",
                          Breed_Effect,
                          h2,
                          trait_mean,
                          VarE,
                          recL,
                          nChr,
                          recProb,
                          minLength) {
  if (missing(map)) {
    stop("Please define the map file")
  }
  else if (missing(Haplotype)) {
    stop("Please define the haplotype information of the animals")
  }
  else if (!is.matrix(Haplotype)) {
    stop("Haplotype information is not in matrix form")
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
  else if (nDam > nProgeny) {
    stop("Please ensure that number of dams should be smaller than number of progeny")
  }
  else if (missing(Breed_Effect)) {
    stop("Please define the SNP effect")
  }
  else if (missing(h2)) {
    stop("Please define the heritability of the trait to be simulated")
  }
  else if (h2 <= 0) {
    stop("The heritability of the trait is 0 or negative! Check the heritaility and try again.")
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
    stop("Define the number of chromosomes. Please ensure that haploid number of chromosomes (N) are provided")
  }
  else if (nChr > 50) {
    stop("Haploid number of chromosomes (N) exceeded 50. Please check and try again.")
  }

  ### Calculate TBV ###
  pop = BreedEff_Haplo(Haplotype, Breed_Effect)
  ## Calculate phenotype
  popTBVpheno = Var_PhenoRecords(pop, h2, trait_mean, VarE)

  ### Sire and Dam with TBV and phenotype
  popinfo = popTBVpheno$Haplo_pheno

  #### Divide the population in two parts
  ## Sires
  Sirepop = popinfo[1:(nrow(popinfo) / 2), ]
  ## Dams
  Dampop = popinfo[(nrow(popinfo) / 2 + 1):nrow(popinfo), ]


  ##### Sampling of random sires and dams for cross ####
  #################################################
  if (SelType == "TBV") {
    #### Order the haplotype based on TBV
    Sirepop = Sirepop[rev(order(Sirepop[, 'TBV'])), ]
    Dampop = Dampop[rev(order(Dampop[, 'TBV'])), ]

    ## Select Sire and Dam based on no. of sires and dam
    Sirepop = Sirepop[1:(2 * nSire), -c(ncol(Sirepop), ncol(Sirepop) - 1)]
    Dampop = Dampop[1:(2 * nDam), -c(ncol(Dampop), ncol(Dampop) - 1)]

    ### Select sires randomly (few sires are selected)
    SireR = SireSample(Sirepop, nSire, nProgeny)
    ### Select dams randomly (Usually all dams are selected)
    DamR = DamSample(Dampop, nDam, nProgeny)

  } else if (SelType == "Pheno") {
    ## Order the haplotype based on phenotype
    Sirepop = Sirepop[rev(order(Sirepop[, 'phenotype'])), ]
    Dampop = Dampop[rev(order(Dampop[, 'phenotype'])), ]

    ## Select Sire and Dam based on no. of sires and dam
    Sirepop = Sirepop[1:(2 * nSire), -c(ncol(Sirepop), ncol(Sirepop) - 1)]
    Dampop = Dampop[1:(2 * nDam), -c(ncol(Dampop), ncol(Dampop) - 1)]
    ### Select sires randomly (few sires are selected)
    SireR = SireSample(Sirepop, nSire, nProgeny)
    ### Select dams randomly (Usually all dams are selected)
    DamR = DamSample(Dampop, nDam, nProgeny)
  } else if (SelType == "random") {
    ## Select Sire and Dam based on no. of sires and dam
    Sirepop = Sirepop[, -c(ncol(Sirepop), ncol(Sirepop) - 1)]
    Dampop = Dampop[, -c(ncol(Dampop), ncol(Dampop) - 1)]
    ### Select sires randomly (few sires are selected)
    SireR = SireSample(Sirepop, nSire, nProgeny)
    ### Select dams randomly (Usually all dams are selected)
    DamR = DamSample(Dampop, nDam, nProgeny)
  }

  ## Get Sire and Dam haplotype details
  SireSelect = SireR$Matrix
  DamSelect = DamR$Matrix

  ## Add recombination points for Sire and Dams
  ## Sire
  RecPointListS = RecombinationPoint(map, SireSelect, recL, nChr, recProb)
  RecSirelist = RecPointListS$recList
  SireSelectFinal <- addSwitch(SireSelect, RecSirelist, minLength)
  rownames(SireSelectFinal) = rownames(SireSelect)

  ##Dam
  RecPointListD = RecombinationPoint(map, DamSelect, recL, nChr, recProb)
  RecDamlist = RecPointListD$recList
  DamSelectFinal <- addSwitch(DamSelect, RecDamlist, minLength)
  rownames(DamSelectFinal) = rownames(DamSelect)

  ###### Selection of random gametes of sire and dam
  sel_sire = Random_gametes(Gamete_Matrix = SireSelectFinal)
  sel_dam = Random_gametes(Gamete_Matrix = DamSelectFinal)

  ## Sire
  ## Add mutation points
  sel_sire_final = MutationRate(sel_sire, rate)

  ## Sire
  ## Add mutation points
  sel_dam_final = MutationRate(sel_dam, rate)

  ####Get parents ID###
  SireID = rownames(sel_sire_final)
  DamID = rownames(sel_dam_final)

  return(
    list(
      Sire = sel_sire_final,
      Dam = sel_dam_final,
      SireID = SireID,
      DamID = DamID,
      PopInfo = popinfo,
      rec_sire = RecPointListS,
      rec_dam = RecPointListD
    )
  )
}



#' Generate Population with specific breed effect
#'
#' @description
#' Generate Population of a crossbred with breed effect specific to the breed of origin.
#'
#' @param ngenerations Number of generations to be simulated.
#' @param map Map file for the population comprising of marker ID, chromosome number and marker position.
#' @param populationSim Coded haplotype information of the population where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param nSire Number of sires to be selected as parents for the simulation.
#' @param nDam Number of Dams to be selected as parents for the simulation.
#' @param mutationRate Mutation rate to added to the simulated population.
#' @param SelType Selection Type. It can be based on phenotype, random or true breeding value. The default selection type is "random".
#' For phenotypic selection, use "Pheno" and for true breeding value, use "TBV".
#' @param Breed_Effect A matrix comprising breed effect according to the breed of origin where rows contain information
#' on SNP effects.
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
#' @param IndStartVal Starting value for the Animal ID of the simulated population.
#' @param prefixID Prefix to be added to the ID of the animal.
#' @param nProgeny Number of progeny to be simulated in the population.
#'
#'
#' @returns Returns a list of Haplotype information, pedigree of the simulated animals, population with TBV and phenotype, recombination point list of sire and dam.
#' @seealso [HAPTRACE::BreedEff_Haplo()], [HAPTRACE::BreedPopSelect()]
#'
#' @export
#'
#' @importFrom hsphase addSwitch
#'
#' @examples PopGen <- matrix(data = sample(c(-1, 1), 200, replace = TRUE), nrow = 20, ncol = 10)
#' PopGen <- PopulationID(Population = PopGen, PopCode = "PG1")
#' PopGen_map <- generateMAP(nChr = 2, n_markers = 5, len_Chr = 5)
#' QTL_PopGen <- matrix(data = c(-0.001, 0.0003, 0, 0, 0,
#' 0.0032, 0.0023, -0.0012,0, 0), nrow = 1, ncol = 10)
#' PopF1 <- GenBreedPop(ngenerations = 2, map = PopGen_map, populationSim = PopGen,
#' nSire = 1, nDam = 1,mutationRate = 2.5 * 10^-5, SelType = "TBV", h2 = 0.3, trait_mean = 5,
#' VarE = 0.6, Breed_Effect = QTL_PopGen, recL = 1, nChr = 2, minLength = 0,
#' IndStartVal = 1,prefixID = "A", nProgeny = 10)
#'
GenBreedPop <- function(ngenerations = 5,
                        map = map,
                        populationSim,
                        nSire = 54,
                        nDam = 1000,
                        mutationRate = 2.5 * 10^-5,
                        SelType = "random",
                        h2 = 0.3,
                        trait_mean,
                        VarE,
                        Breed_Effect,
                        recL,
                        nChr,
                        recProb,
                        minLength,
                        IndStartVal = 1,
                        prefixID = "A",
                        nProgeny = 1000)
{

  if (missing(map)) {
    stop("Please define the map file")
  }
  else if (missing(populationSim)) {
    stop("Please define the haplotype information of the population")
  }
  if (!is.matrix(populationSim)) {
    stop("Haplotype information is not in matrix form")
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
  else if (missing(Breed_Effect)) {
    stop("Please define the SNP effect")
  }
  else if (missing(h2)) {
    stop("Please define the heritability of the trait to be simulated")
  }
  else if (h2 <= 0) {
    stop("The heritability of the trait is 0 or negative! Check the heritaility and try again.")
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
    stop("Define the number of chromosomes. Please ensure that haploid number of chromosomes (N) are provided")
  }
  else if (nChr > 50) {
    stop("Haploid number of chromosomes (N) exceeded 50. Please check and try again.")
  }

  PopInfo = list()
  recSire = list()
  recDam = list()
  OnlyF1 = list()
  F1Haplotype = list()
  Parents_list = list()
  generation <- 1
  for (generation in 1:ngenerations) {
    message("Working on Generation: ", generation)
    ### New population ##
    Newpop <- BreedPopSelect(
      map,
      populationSim,
      nSire,
      nDam,
      mutationRate,
      nProgeny,
      SelType,
      Breed_Effect,
      h2,
      trait_mean,
      VarE,
      recL,
      nChr,
      recProb,
      minLength
    )

    PopInfodata = Newpop$PopInfo
    recpointSire = Newpop$rec_sire
    recpointDam = Newpop$rec_dam

    ##New population
    OnlyF1 = MatingPop(Newpop$Sire, Newpop$Dam)

    ###Get Sire and Dam
    #### Parents in next generation ####
    Parents = data.frame(
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
    populationSim = OnlyF1
    F1Haplotype[[generation]] = OnlyF1
    PopInfo[[generation]] = PopInfodata
    recSire[[generation]] = recpointSire
    recDam[[generation]] = recpointDam
    Parents_list[[generation]] = Parents
    IndStartVal = (IndStartVal + nrow(OnlyF1) / 2 - 1) + 1
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


#' Function for crossing two population with breed effects according to breed of origin
#'
#' @description
#' Function to cross two population and estimation of TBV with breed effect according to the breed of origin.
#'
#'
#' @param map Map file for the population comprising of marker ID, chromosome number and marker position.
#' @param populationSim_A Coded haplotype information of the population A where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param populationSim_B Coded haplotype information of the population B where rows are number of animals and columns are number of markers in matrix form.
#' In case of haplotype, 2 rows provide information on genomic material of an animal.
#' @param nSireA Number of sires to be selected as parents for the simulation from population A.
#' @param nDamA Number of dams to be selected as parents for the simulation from population A.
#' @param nSireB Number of sires to be selected as parents for the simulation from population B.
#' @param nDamB Number of dams to be selected as parents for the simulation from population B.
#' @param Sire_From The population from which sires should being selected for cross population.
#' If sire should be selected from population A, then use "A" and if sire should be selected from population B,
#' then use "B".
#' @param mutationRate Mutation rate to added to the simulated population.
#' @param SelType Selection Type. It can be based on phenotype, random or true breeding value. The default selection type is "random".
#' For phenotypic selection, use "Pheno" and for true breeding value, use "TBV".
#' @param h2 Heritability of the trait to be simulated.
#' @param trait_mean Mean of the trait to be simulated.
#' @param VarE_PopA Error Variance obtained after rescaling the QTL effects for PopA.
#' @param VarE_PopB Error Variance obtained after rescaling the QTL effects for PopB.
#' @param Breed_Effect_A A matrix comprising breed effect according to the breed of origin for population A.
#' @param Breed_Effect_B A matrix comprising breed effect according to the breed of origin for population B.
#' @param IndStartVal Starting value for the animal ID.
#' @param prefixID Prefix to be added to the ID of the animal.
#' @param nProgeny Number of progeny to be simulated in the population.
#' @param recL Average number of recombinations to be introduced per animal.
#' @param nChr Number of haploid chromosomes (n) defined in the population.
#' @param recProb A vector of recombination probability between the markers. Its length should be total
#' number of markers minus 1. User can define the recombination probability between markers
#' helping in simulating a desired population.
#' @param minLength Minimum length between recombination points beyond which recombination is allowed. The default value is 0.
#' If a user define minimum length, it will ensure that difference between two recombination points
#' less than minimum length won't allow recombination.
#'
#'
#'
#' @returns Returns a list of Haplotype information, ID of the simulated animals, population with TBV and phenotype, recombination point list of sire and dam.
#' @seealso [HAPTRACE::BreedEff_Haplo()], [HAPTRACE::BreedPopSelect()]
#' @export
#'
#' @importFrom hsphase addSwitch
#'
#' @examples crossmap <- generateMAP(nChr = 5, n_markers = 2, len_Chr = 20)
#' QTL_cross <- matrix(data = c(-0.001, 0.0003, 0, 0, 0, 0.0032, 0.0023, -0.0012,0, 0)
#' , nrow = 1, ncol = 10)
#' P1cross <- matrix(data = sample(c(-1, 1), 40, replace = TRUE), nrow = 4, ncol = 10)
#' P2cross <- matrix(data = sample(c(-2, 2), 40, replace = TRUE), nrow = 4, ncol = 10)
#' P1cross <- PopulationID(Population = P1cross, PopCode = "P1")
#' P2cross <- PopulationID(Population = P2cross, PopCode = "P2")
#' F2cross <- BreedPopCross(map = crossmap, populationSim_A = P1cross,
#' populationSim_B = P2cross, nSireA = 1, nDamA = 1, nSireB = 1, nDamB = 1,
#' Sire_From = "A",mutationRate = 2.5 * 10^-5, SelType = "TBV", h2 = 0.3,
#' Breed_Effect_A = QTL_cross, Breed_Effect_B = QTL_cross, IndStartVal = 1,
#' prefixID = "Cr", nProgeny = 10, recL = 1, nChr = 2, minLength = 0,
#' trait_mean = 5, VarE_PopA = 0.6, VarE_PopB = 0.6)
#'
BreedPopCross <- function(map = map,
                          populationSim_A,
                          populationSim_B,
                          nSireA = 54,
                          nDamA = 1000,
                          nSireB = 26,
                          nDamB = 1000,
                          Sire_From,
                          mutationRate = 2.5 * 10^-5,
                          SelType,
                          h2 = 0.3,
                          trait_mean,
                          VarE_PopA,
                          VarE_PopB,
                          Breed_Effect_A,
                          Breed_Effect_B,
                          IndStartVal = 1,
                          prefixID = "Cr",
                          nProgeny = 1000,
                          recL,
                          nChr,
                          recProb,
                          minLength)
{
  if (missing(map)) {
    stop("Please define the map file")
  }
  else if (missing(populationSim_A)) {
    stop("Please define the haplotype information of the population A")
  }
  else if (missing(populationSim_B)) {
    stop("Please define the haplotype information of the population B")
  }
  if (!is.matrix(populationSim_A)) {
    stop("Haplotype information of population A is not in matrix form")
  }
  if (!is.matrix(populationSim_B)) {
    stop("Haplotype information of population B is not in matrix form")
  }
  else if (missing(nSireA)) {
    stop("Please define the number of sires to be used for simulation for population A")
  }
  else if (missing(nDamA)) {
    stop("Please define the number of dams to be used for simulation for population A")
  }
  else if (missing(nSireB)) {
    stop("Please define the number of sires to be used for simulation for population B")
  }
  else if (missing(nDamB)) {
    stop("Please define the number of dams to be used for simulation for population B")
  }
  else if (missing(Sire_From)) {
    stop("Please define from which population sire should be selected")
  }
  else if (missing(h2)) {
    stop("Please define the heritability of the trait to be simulated")
  }
  else if (h2 <= 0) {
    stop("The heritability of the trait is 0 or negative! Check the heritaility and try again.")
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
  else if (missing(VarE_PopA)) {
    stop("Please define the error variance of the trait for population A")
  }
  else if (missing(VarE_PopB)) {
    stop("Please define the error variance of the trait for population B")
  }
  else if (missing(Breed_Effect_A)) {
    stop("Please define the SNP effect for Population A")
  }
  else if (missing(Breed_Effect_B)) {
    stop("Please define the SNP effect for Population B")
  }
  else if (missing(nProgeny)) {
    stop("Please define the number of progeny to be simulated per generation")
  }
  else if (nSireA > nProgeny) {
    stop("Please ensure that number of sires should be smaller than number of progeny")
  }
  else if (nDamA > nProgeny) {
    stop("Please ensure that number of dams should be smaller than number of progeny")
  }
  else if (nSireB > nProgeny) {
    stop("Please ensure that number of sires should be smaller than number of progeny")
  }
  else if (nDamB > nProgeny) {
    stop("Please ensure that number of dams should be smaller than number of progeny")
  }
  else if (missing(recL)) {
    stop("Please define this parameter")
  }
  else if (missing(nChr)) {
    stop("Define the number of chromosomes. Please ensure that haploid number of chromosomes (N) are provided")
  }
  else if (nChr > 50) {
    stop("Haploid number of chromosomes (N) exceeded 50. Please check and try again.")
  }

  message("Working on two populations cross")
  OnlyF3 = list()
  Parents_list3 = list()
  #### get the next generation #####
  ### New population ##
  PopBaseA = BreedPopSelect(
    map,
    populationSim_A,
    nSireA,
    nDamA,
    mutationRate,
    nProgeny,
    SelType,
    Breed_Effect_A,
    h2,
    trait_mean,
    VarE_PopA,
    recL,
    nChr,
    recProb,
    minLength
  )
  PopBaseB = BreedPopSelect(
    map,
    populationSim_B,
    nSireB,
    nDamB,
    mutationRate,
    nProgeny,
    SelType,
    Breed_Effect_B,
    h2,
    trait_mean,
    VarE_PopB,
    recL,
    nChr,
    recProb,
    minLength
  )
  if(Sire_From == "A"){
    #### Sire and Dam ###
    SireG = PopBaseA$Sire
    DamG = PopBaseB$Dam
    ##New population
    OnlyF3 = MatingPop(SireG, DamG)
    #########Select parents for 2 population
    Parents_F3 = data.frame(
      AnimID = paste(prefixID, IndStartVal:(IndStartVal + nrow(OnlyF3) / 2 - 1), sep =
                       "_"),
      SireID = PopBaseA$SireID,
      DamID = PopBaseB$DamID,
      Sex = as.numeric(c(rep(1, nrow(
        OnlyF3
      ) / 4), rep(2, nrow(
        OnlyF3
      ) / 4)))
    )
  }
  else if(Sire_From == "B"){
    #### Sire and Dam ###
    SireG = PopBaseB$Sire
    DamG = PopBaseA$Dam
    ##New population
    OnlyF3 = MatingPop(SireG, DamG)
    #########Select parents for 2 population
    Parents_F3 = data.frame(
      AnimID = paste(prefixID, IndStartVal:(IndStartVal + nrow(OnlyF3) / 2 - 1), sep =
                       "_"),
      SireID = PopBaseB$SireID,
      DamID = PopBaseA$DamID,
      Sex = as.numeric(c(rep(1, nrow(
        OnlyF3
      ) / 4), rep(2, nrow(
        OnlyF3
      ) / 4)))
    )
  }

  rownames(OnlyF3) <- rep(Parents_F3$AnimID, each = 2)
  return(
    list(
      Haplotype = OnlyF3,
      parents = Parents_F3,
      popA = PopBaseA$PopInfo,
      recpointS_popA = PopBaseA$rec_sire,
      recpointD_popA = PopBaseA$rec_dam,
      popB = PopBaseB$PopInfo,
      recpointS_popB = PopBaseB$rec_sire,
      recpointD_popB = PopBaseB$rec_dam
    )
  )
}
