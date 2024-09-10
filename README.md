
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spconf

R Package for Computing Scales of Spatial Smoothing for Confounding
Adjustment

This package is designed to calculate the *effective bandwidth* of a
spatial smoothing matrix, following two different procedures, one
described by Keller and Szpiro (2020). This package also contains a
wrapper function to create a thin-plate regression spline basis using
the mgcv package (Wood, 2011).

## Installation

The `spconf` package can be installed by running

    devtools::install_github("kpkeller/spconf")

## Computing Effective Range

The primary function is `compute_effective_range()`. The function takes
a matrix of spline values `X`, with the assumption that the splines are
nested to that adding additional columns increases the flexibility of
forms that can be fit. For each choice of degrees of freedom, the
function computes the effective range of the smoothing matrix. Using
`smoothedCurve = TRUE`, the effective range is determined by the
procedure described by Keller and Szpiro (2020) and require the `span`
input; and using `smoothedCurve = FALSE`, the effective range is
determined by finding smallest distance that has a negative weight and
no longer requires the `span` input.

``` r
xloc <- runif(n=100, min=0, max=10)
X <- splines::ns(x=xloc, df=4, intercept=TRUE)
colnames(X) <- paste0("s", 1:ncol(X))
xplot <- 0:10
compute_effective_range(X=X, coords=as.matrix(xloc), df=2:4, 
                        smoothedCurve = TRUE, newd=xplot, namestem="s")
#> Df =  2 
#> Df =  3 
#> Df =  4
#>        2        3        4 
#> 4.397386 4.230755 3.972437
```

``` r

M <- 16
tprs_df <- 10
si <- seq(0, 1, length=M+1)[-(M+1)]
gridcoords <- expand.grid(x=si, y=si)
tprsX <- computeTPRS(coords = gridcoords, maxdf = tprs_df+1)
compute_effective_range(X=tprsX$tprsX, coords=gridcoords, namestem = "tprs",
                        df=3:10, smoothedCurve = FALSE)
#> Df =  3 
#> Df =  4 
#> Df =  5 
#> Df =  6 
#> Df =  7 
#> Df =  8 
#> Df =  9 
#> Df =  10
#>         3         4         5         6         7         8         9        10 
#> 0.3801727 0.3365728 0.3125000 0.3125000 0.3125000 0.3125000 0.3125000 0.3186887
```

## References

Keller and Szpiro (2020). Selecting a scale for spatial confounding
adjustment. *Journal of the Royal Statistical Society, Series A*
<https://doi.org/10.1111/rssa.12556>.

Wood (2011). Fast Stable Restricted Maximum Likelihood and Marginal
Likelihood Estimation of Semiparametric Generalized Linear Models.
*Journal of the Royal Statistical Society Series B: Statistical
Methodology* <https://doi.org/10.1111/j.1467-9868.2010.00749.x>
