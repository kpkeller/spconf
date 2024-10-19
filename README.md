
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spconf

R Package for Computing Scales of Spatial Smoothing for Confounding
Adjustment

This package is designed to calculate the *effective bandwidth* of a
spatial smoothing matrix, following two different procedures, described
by Keller and Szpiro (2020) and Rainey and Keller (2024). This package
also contains a wrapper function to create a thin-plate regression
spline basis using the mgcv package (Wood, 2011).

## Installation

The `spconf` package can be installed by running

    remotes::install_github("kpkeller/spconf")

## Computing Effective Range

The primary function is `compute_effective_range()`. The function takes
a matrix of spline values `X`, with the assumption that the splines are
nested to that adding additional columns increases the flexibility of
forms that can be fit. For each choice of degrees of freedom, the
function computes the effective range of the smoothing matrix. Using
`smoothedCurve = TRUE`, the effective range is computed using the
procedure introduced by Keller and Szpiro (2020) and requires the `span`
input; and using `smoothedCurve = FALSE`, the effective range is
computed using the procedure introduced by Rainey and Keller (2024)
which no longer requires the `span` input.

``` r
# Example using metric of Rainey and Keller (2024)
M <- 16
tprs_df <- 10
si <- seq(0, 1, length=M+1)[-(M+1)]
gridcoords <- expand.grid(x=si, y=si)
tprsX <- computeTPRS(coords = gridcoords, maxdf = tprs_df+1)
compute_effective_range(X=tprsX, coords=gridcoords, namestem = "tprs",
                        df=3:10, smoothedCurve = FALSE)
#>         3         4         5         6         7         8         9        10 
#> 0.3801727 0.3365728 0.3125000 0.3125000 0.3125000 0.3125000 0.3125000 0.3186887
```

``` r

# Example using metric of Keller and Szpiro (2020)
xloc <- runif(n=100, min=0, max=10)
X <- splines::ns(x=xloc, df=4, intercept=TRUE)
colnames(X) <- paste0("s", 1:ncol(X))
xplot <- 0:10
compute_effective_range(X=X, coords=as.matrix(xloc), df=2:4, 
                        smoothedCurve = TRUE, newd=xplot, namestem="s")
#>        2        3        4 
#> 4.512474 4.472711 4.215862
```

## References

Keller and Szpiro (2020). Selecting a scale for spatial confounding
adjustment. *Journal of the Royal Statistical Society, Series A*
<https://doi.org/10.1111/rssa.12556>.

Rainey and Keller (2024). spconfShiny: An R Shiny application for
calculating the spatial scale of smoothing splines for point data. *PLOS
ONE* <https://doi.org/10.1371/journal.pone.0311440>

Wood (2011). Fast Stable Restricted Maximum Likelihood and Marginal
Likelihood Estimation of Semiparametric Generalized Linear Models.
*Journal of the Royal Statistical Society Series B: Statistical
Methodology* <https://doi.org/10.1111/j.1467-9868.2010.00749.x>
