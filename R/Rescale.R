#' Rescaling SNP Effect
#'
#' @description
#' This function allows scaling of the SNP effects based on the scaling factor
#' derived from genotype of the base population and standard deviation of the phenotype.
#' This allows users to scale the SNP effects to get desirable variance components.
#'
#' @param Haplotype A matrix where each row represents a haplotype (0,1)
#' and each column represents a marker.
#' @param Effect A vector with SNP effect of a population to be scaled before the simulation of a population
#' @param Phenosd A value of phenotypic standard deviation which will help in estimating desirable variance components
#' @param h2 Heritability of the trait to be simulated for the population.
#'
#' @returns This function returns a scaled SNP effect, scaled BV along with
#' phenotypic, additive and residual variance.
#' @export
#'
#' @importFrom stats sd
#'
#' @examples popRescale <- matrix(c(0,1,1,0,
#' 1,0,0,0,
#' 1,1,1,1,
#' 1,1,1,0), byrow = TRUE, nrow = 4, ncol = 4)
#' QTLEff = c(0.1, 0.2, 0, 0)
#' RescaleEff <- Rescale(Haplotype = popRescale, Effect = QTLEff,
#' Phenosd = 4, h2 = 0.2)
Rescale = function(Haplotype, Effect, Phenosd, h2) {
  if (missing(Haplotype)) {
    stop("Please define the Haplotype information")
  }
  else if (!is.matrix(Haplotype)) {
    stop("Haplotype information is not in matrix form")
  }
  else if (missing(Effect)) {
    stop("Please define the SNP effect")
  }
  else if (missing(Phenosd)) {
    stop("Please define the standard deviation of the phenotype (trait to be simulated)")
  }
  else if (missing(h2)) {
    stop("Please define the heritability of the trait to be simulated")
  }
  else if (h2 <= 0) {
    stop("The heritability of the trait less than or equal to zero! Check the heritability again.")
  }
  else if (h2 == 1) {
    warning("The heritability of the trait is 1! Proceed with caution!")
  }
  else if (h2 > 1) {
    stop("The heritability of the trait is more than 1! Check the heritability again.")
  }
  else {
    PhenoVar = Phenosd * Phenosd
    AddVar = h2 * PhenoVar
    ErrorVar = (1 - h2) * PhenoVar
    Geno = MakeGeno(Coded_to_Haplo(Haplotype))
    BV = base::tcrossprod(Geno, t(Effect))
    sd_bv = stats::sd(BV)
    if (is.na(sd_bv) || !is.finite(sd_bv) || sd_bv == 0){
      stop("Rescaling cannot be done: check your input parameter again")
    }
    scaleBV = base::sqrt(AddVar) / sd_bv
    BV_R1 = BV * scaleBV
    ScaledQTLEff = Effect * scaleBV
    return(
      list(
        QTLeffect = ScaledQTLEff,
        BreedingValue = BV_R1,
        scale = scaleBV,
        PhenoVar = PhenoVar,
        AddVar = AddVar,
        ErrorVar = ErrorVar
      )
    )
  }
}

