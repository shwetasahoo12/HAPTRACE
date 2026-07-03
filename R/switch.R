#' Adding recombination to the haplotype matrix
#'
#' @param haplotypeMatrix A matrix with haplotype information where rows are number of animals and columns are number of markers.
#'
#' @param switchPoints A list of recombination points per animal.
#' @param minLength Minimum length between recombination points beyond which recombination is allowed. The default value is 0.
#' If a user define minimum length, it will ensure that difference between two recombination points
#' less than minimum length won't allow recombination.
#'
#' @return Haplotype matrix with recombination
#' @export
#'
#' @import Rcpp
#' @import RcppArmadillo
#' @useDynLib HAPTRACE, .registration = TRUE
#'
#' @examples haplotype <- matrix(c(0, 0, 0, 0,
#'1, 1, 1, 1,
#'0, 0, 1, 1,
#'1, 1, 0, 0,
#'1, 1, 1, 1,
#'0, 0, 0, 0), byrow = TRUE, nrow = 6)
#'switchPoints <- list(firstInd = c(2), secondInd = c(1, 3), latsInd = 0)
#'addSwitch(haplotype, switchPoints, 0)

addSwitch <- function(haplotypeMatrix, switchPoints, minLength)
{

  ID = rownames(haplotypeMatrix)
  if(nrow(haplotypeMatrix)/2!=length(switchPoints))
  {
    stop("The number of elements in the switchPoints must be equal to the number of individuals")
  }

  for (i in 1:length(switchPoints))
  {
    switchPoints[[i]] <- c(0, switchPoints[[i]], ncol(haplotypeMatrix))
  }
  result <- .Call("switchAdd", haplotypeMatrix, switchPoints, minLength, PACKAGE = "HAPTRACE")
  rownames(result) = ID
  result
}

