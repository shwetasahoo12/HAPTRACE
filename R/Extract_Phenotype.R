#' Extract TBV and Phenotype
#'
#' @param df A matrix comprising the haplotype information along with TBV and Phenotype.
#' Ensure that rownames of the matrix is Animal ID
#'
#' @returns A datafile comprising Animal ID, TBV and phenotype.
#' @export
#'
#' @examples popEP <- matrix(data= sample(c(0,1), 16, replace= TRUE), nrow = 4, ncol = 4)
#' popEP <- PopulationID(Population = popEP, PopCode = "EP")
#' QTLEff <- c(0.2, 0.1, 0.3, 0.1)
#' TBV_haplotype2 <- TBV_Haplo(Haplotype = popEP, Effect = QTLEff)
#' Var_Ph_R <- Var_PhenoRecords(Haplotype = TBV_haplotype2, h2 = 0.3,
#' trait_mean = 20, VarE = 0.6)
#' Haplo_TBV_Ph <- Var_Ph_R$Haplo_pheno
#' Datafile <- Ext_Pheno(Haplo_TBV_Ph)
#'
Ext_Pheno = function(df){
  if (missing(df)) {
    stop("Please define the file")
  }
  ## Select TBV and phenotype
  df2 = df[, (ncol(df) - 1):ncol(df)]
  ## Select alternate rows
  df2Pheno = df2[seq(1, nrow(df2), by = 2), ]
  ID = rownames(df2Pheno)
  df2Pheno = as.data.frame(df2Pheno)
  df2Pheno$AnimID = ID
  ## Arrange the df
  df3 = df2Pheno[, c(3,1,2)]
  return(df3)
}


#' Combine Pedigree and Phenotype data
#'
#' @param Pedigree Pedigree file comprising AnimalID, Sire ID, Dam ID and Sex of the animal.
#' @param Phenotype Phenotype file comprising Animal ID, TBV and phenotype of the animal
#'
#' @returns A datafile comprising AnimalID, Sire ID, Dam ID, Sex, TBV and phenotype of the animal.
#' @export
#'
#' @seealso [HAPTRACE::Ext_Pheno()]
#'
#' @examples PopGen1 <-  matrix(data = sample(c(-1, 1), 200, replace = TRUE), nrow = 20, ncol = 10)
#' PopGen1 <- PopulationID(Population = PopGen1, PopCode = "PG")
#' PopGen1_map <-  generateMAP(2, 5, 5)
#' QTL_PopGen1  <-  c(-0.001, 0.0003, 0, 0, 0, 0.0032, 0.0023, -0.0012,0, 0)
#' PopPP <- GenOnePopulation(ngenerations=2, map=PopGen1_map,
#' PopGen1,nSire=1,nDam=1,mutationRate=2.5 * 10^-5,
#' SelType = "TBV", h2 = 0.3, trait_mean = 5 ,VarE = 0.6, Effect = QTL_PopGen1,
#' recL = 1, nChr = 2, minLength = 0, IndStartVal=1,prefixID="A",nProgeny=10)
#' Phenofile <- Ext_Pheno(df = PopPP$population)
#' Ped_Phenofile <- Ped_Pheno(Pedigree = PopPP$parents, Phenotype = Phenofile)
#'
#'
Ped_Pheno = function(Pedigree, Phenotype){
  if (missing(Pedigree)) {
    stop("Please define the Pedigree file")
  }
  else if (missing(Phenotype)) {
    stop("Please define the Phenotype file")
  }
  if(!("AnimID" %in% names(Phenotype))){
    stop("AnimID is not found in the Phenotype file.
         Make sure column name for animal IDs is AnimID")
  }
  ## Combined datafile for Pedigree and Phenotype
  Combined = merge(Pedigree, Phenotype[, ], by = "AnimID", all.x = TRUE)
  Pref = gsub("[0-9_]", "", Combined$AnimID)
  Number = as.numeric(gsub("[^0-9]", "", Combined$AnimID))
  CombinedData = Combined[order(Pref, Number), ]
  return(CombinedData)
}

