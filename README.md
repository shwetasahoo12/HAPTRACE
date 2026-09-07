
# HAPTRACE

Haplotype-Based Tracking of Admixed Population for Breed Composition
Estimation (HAPTRACE) is an open-source R package which simulates and
tracks haplotype in an admixed population for accurate estimation of
breed composition. This package offers multiple functionalities
including estimation of true breed composition, visualization of
haplotypes of admixed animals, simulation of purebred and admixed
animals with tracked haplotype, genotype and phenotype.

## Installation

The github repository contains the recent version of HAPTRACE (0.1.2).
You can install the development version of HAPTRACE from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("shwetasahoo12/HAPTRACE")
```

## Example

To see the example simulation of the package, refer to the vignette
below:

``` r
## Example simulation
library(HAPTRACE)
vignette("Introduction_HAPTRACE", package = "HAPTRACE")
#> Warning: vignette 'Introduction_HAPTRACE' not found
```

## Citation

If you use this package in your research, please cite it as:

Sahoo S, Ferdosi M, van der Werf J, de las Heras-Saldana S (2025).
“Proc. Assoc. Advmt. Anim. Breed. Genet.” Proc. Assoc. Advmt. Anim.
Breed. Genet., 26, 323.
