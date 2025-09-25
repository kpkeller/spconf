# Functions in File #
# computeS
# fitLoess
# compute_lowCurve
# find_first_zero_cross
# find_zeros_cross
# compute_effective_range
# compute_effective_range_nochecks



#' @title Compute Smoothing Matrix
#' @description Calculates the smoothing (or "hat") matrix from a design matrix.
#' @param x Matrix of spline values, assumed to have full rank. A data frame is coerced into a matrix.
#' @param inds Column indices of smoothing matrix to return (corresponding to rows in \code{x}).
#' @details Given a matrix \code{X} of spline values, this computes  \code{S=X(X'X)^(-1)X'}. When \code{x} has many rows, this can be quite large. The \code{inds} argument can be used to return a subset of columns from \code{S}.
#'
#' @return  An \eqn{N}-by-\eqn{n} matrix, where \eqn{n} is the length of \code{inds} and \eqn{N} is the number of rows in \code{x}.
#' @export
#' @seealso \code{\link{compute_effective_range}}
#' @examples
#' # Simple design matrix case
#' X <- cbind(1, rep(c(0, 1), each=4))
#' S <- computeS(X)
#' # More complex example
#' xloc <- runif(n=100, min=0, max=10)
#' X <- splines::ns(x=xloc, df=4, intercept=TRUE)
#' S <- computeS(X)
#' S2 <- computeS(X, inds=1:4)
computeS <- function(x, inds=1:nrow(x)){
    if (!inherits(x, "matrix") && !inherits(x, "data.frame")) stop("'x' must be a matrix or data frame.")
    x <- as.matrix(x)
    S <- x %*% solve(crossprod(x), t(x[inds,]))
 S
}

#' @title Fit a loess curve
#' @description Wrapper function for fitting and predicting from \code{loess()}.
#' @param y Dependent variable values
#' @param x Independent variable values
#' @param newx Values of x to use for prediction.
#' @param span Controls the amount of smoothing. Passed to \code{\link[stats]{loess}}; see that function for details.
#' @param ... Additional arguments passed to \code{\link[stats]{loess}}
#' @return A vector of the same length of \code{newx} providing the predictions from a loess smooth.
#' @export
#' @importFrom stats loess predict
#' @examples
#' x <- seq(0, 5, length=50)
#' y <- cos(4*x) + rnorm(50, sd=0.5)
#' xplot <- seq(0, 5, length=200)
#' lfit <- fitLoess(y=y, x=x, newx=xplot, span=0.2)
#' plot(x, y)
#' points(xplot, lfit, type="l")
fitLoess <- function(y, x, newx=x, span=0.5,...){
    lowcurve <- stats::loess(y~x, span=span,...)
    lowpred <- stats::predict(lowcurve, newdata=newx)
    lowpred
}


#' @title Compute loess curves for smoothing matrix
#' @description Calculates a loess curve for the smoothing matrix entries, as a function of distance between points.
#' @param S Smoothing matrix, or a subset of columns from a smoothing matrix.
#' @param D Distance matrix, or a subset of columns from a distance matrix.
#' @param newd Distances to use for loess prediction.
#' @param cl Cluster object, or number of cluster instances to create. Defaults to no parallelization.
#' @param span Passed to \code{\link{fitLoess}}
#' @return List with three elements: \eqn{n}-by-\eqn{N} matrix, where \eqn{n} is the length of \code{newd} and \eqn{N} is the number of columns in \code{S}; a vector of length \eqn{n} giving the median curve value; a vector of length \eqn{n} giving the mean curve value.
#' @details For each column in \code{S}, a loess curve is fit to the values as a function of the distances between points, which are taken from the columns of \code{D}. Thus, the order of rows and columns in \code{S} should match the order of rows and columns in \code{D}.
#' For a large number of locations, this procedure may be somewhat slow. The \code{cl} argument can be used to parallelize the operation using \code{\link[parallel]{clusterMap}}.
#' @seealso \code{\link{computeS}} \code{\link{fitLoess}}
#' @export
#' @importFrom parallel makeCluster clusterExport clusterMap stopCluster
#' @importFrom stats median
#' @examples
#'
#' xloc <- runif(n=100, min=0, max=10)
#' X <- splines::ns(x=xloc, df=4, intercept=TRUE)
#' S <- computeS(X)
#' d <- as.matrix(dist(xloc))
#' xplot <- 0:10
#' lC <- compute_lowCurve(S, D=d, newd=xplot)
#' matplot(xplot, lC$SCurve, type="l", col="black")
#' points(xplot, lC$SCurveMedian, type="l", col="red")
compute_lowCurve <- function(S, D, newd, cl=NULL, span=0.1){
    if (nrow(S)!=nrow(D)) stop("Rows in S must correspond to rows in D.")
    if (ncol(S)!=ncol(D)) stop("Columns in S must correspond to columns in D.")
    if (is.null(cl)){
        n <- ncol(S)
        lowCurveM <- matrix(NA, nrow=length(newd), ncol=n)
        for (i in 1:n){#nrow(S)) {
            lowCurveM[,i] <- fitLoess(y=S[,i], x=D[, i], newx=newd, span=span)
        }
    } else  {
        stop_cluster <- FALSE
        if (is.numeric(cl) && cl>0) {
            cl <- makeCluster(cl)
            stop_cluster <- TRUE
        }  else  if (!inherits(cl, "cluster")){
            stop("'cl' should be NULL, a cluster object, or a non-negative number indicating number of cores.")
        }
        S_list <- as.list(as.data.frame(S))
        D_list <- as.list(as.data.frame(D))
        newfitLoess <- function(s, d) fitLoess(s, d, newx=newd, span=span)
        environment(newfitLoess) <- .GlobalEnv
        clusterExport(cl=cl, varlist=c("S_list", "D_list", "fitLoess", "newd", "span"), envir=environment())
        lowCurveM <- clusterMap(cl=cl,newfitLoess , s=S_list, d=D_list, SIMPLIFY=TRUE)
        if (stop_cluster) stopCluster(cl) # Stop process, if function started it.
    }
    lowCurveMedian <- apply(lowCurveM, 1, stats::median, na.rm = TRUE)
    lowCurveMean <- apply(lowCurveM, 1, mean, na.rm = TRUE)
    out <- list(SCurve=lowCurveM,
                SCurveMedian=lowCurveMedian,
                SCurveMean=lowCurveMean)
    out
}


# Function for getting zero via linear interpolation between
# the two points straddling zero.
#' @title Find zero
#' @description
#' Calculates the zero of a function by linear interpolation between the first two points either side of zero.
#' @param x Function values, assumed to be ordered
#' @return Index of first value of \code{x} that lies below 0. Decimal values will be returned using a simple interpolation of the two values straddling 0.
#' @export
#' @seealso \code{\link{find_zeros_cross}}, \code{\link{compute_effective_range}}
find_first_zero_cross <- function(x){
    upper_ind <- min(which(x<0))
    lower_ind <- upper_ind- 1
    frac_ind <- x[lower_ind]/(x[lower_ind] - x[upper_ind])
    lower_ind + frac_ind
}

# Function for finding the smallest distance that has a
# negative smoothing weight for each column of the smoothing matrix
#' @title Find distance to first zero
#' @description For a set of distance and smoothing matrix values, determines the smallest distance that corresponds with negative value for each column of the smoothing matrix.
#' @param D Distance matrix, or a subset of columns from a distance matrix.
#' @param S Smoothing matrix, or a subset of columns from a smoothing matrix.
#' @export
#' @return Vector of length equal to the number of columns in \code{D} and \code{S}. Each value is the smallest observed distance (from a column of \code{D}) that has a negative value in the corresponding column of \code{S}.
find_zeros_cross <- function(D, S){
    if (nrow(S)!=nrow(D)) stop("Rows in S must correspond to rows in D.")
    if (ncol(S)!=ncol(D)) stop("Columns in S must correspond to columns in D.")

    out <- rep(0, ncol(D))

    for(i in 1:ncol(D)){
        x <- cbind(D[,i], S[,i])
        x <- x[order(x[,1]),]

        if(min(x[,2]) >= 0){
            out[i] <- NA
        }else{
            upper_ind <- min(which(x[,2] < 0))
            # lower_ind <- upper_ind - 1
            # frac_dist <- x[lower_ind,2]/(x[lower_ind,2] - x[upper_ind,2])
            # cross_vals[i] <- x[lower_ind,1] + (x[upper_ind,1] - x[lower_ind,1])*frac_dist
            out[i] <- x[upper_ind,1]
        }
    }
    out
}



#' @title Compute effective range
#' @description Calculates the effective range for a spline basis matrix.
#' @param X Matrix of spline values. See \code{namestem} for expected column names.
#' @param coords Matrix of point coordinates. Defaults to the \code{x} and \code{y} columns of \code{X}, but can have a different number of columns for settings with different dimensions.
#' @param df Degrees of freedom for which effective range should be computed.
#' @param nsamp Number of observations from \code{X} from which to sample. Defaults to minimum of 1,000 and \code{nrow(X)}.
#' @param newd Distance values at which to make loess predictions. Should correspond to distances in the same units as \code{coords}. Only needed when \code{smoothedCurve} is \code{TRUE}.
#' @param scale_factor Factor by which range should be scaled. Often physical distance corresponding to resolution of grid. Defaults to 1, so that range is reported on the same scale as distance in \code{coords}. Only needed when \code{smoothedCurve} is \code{TRUE}.
#' @param smoothedCurve Should the effective range be computed using the procedure introduced by Keller and Szpiro, 2020, (\code{TRUE}) or the procedure introduced by Rainey and Keller, 2024, (\code{FALSE}). See Details.
#' @param returnFull Should the mean and median curves be returned (\code{TRUE}), or just the range value of where they first cross zero (\code{FALSE}).
#' @param cl Cluster object, or number of cluster instances to create. Defaults to no parallelization.
#' @param namestem Stem of names of columns of X corresponding to evaluated splines.
#' Defaults to \code{""}, meaning names of the form \code{1}, \code{2}, ...
#' @param inds Optional vector of indices to use as subset. If provided, \code{nsamp} is not used.
#' @param verbose Control printing of a df counter to console.
#' @param span Passed to \code{\link{fitLoess}}. If too small, then can lead to unstable loess estimates. Only needed when \code{smoothedCurve} is \code{TRUE}.
#' @return The effective bandwidth for each value of \code{df}. If \code{returnFull = FALSE}, then this is a vector of the same length as \code{df}. If \code{returnFull = TRUE} and \code{smoothedCurve = TRUE}, this is a list that additionally contains values of the pointwise median and mean of the smoothed curves.
#' @details Using the given spline basis and the inputted coordinates, the effective bandwidth is computed for the given degrees of freedom. This is accomplished by computing a distance matrix from the coordinates and a smoothing matrix from the basis.
#' Setting \code{smoothedCurve = TRUE} (see Keller and Szpiro, 2020, for details), for each column of smoothing weights, a LOESS curve is fit to the smoothing weights as a function of the distances, and the distance where the curve first crosses zero is obtained.
#' Setting \code{smoothedCurve = FALSE} (see Rainey and Keller, 2024, for details), for each column of smoothing weights, the smallest distance that corresponds with the first negative smoothing weight is obtained.
#' Then, for both procedures, the median of the obtained distances is reported as the effective bandwidth.
#' @details The columns of \code{X} are selected by name, and so are assumed to have a numeric value in the column name that indicates the spline number. For example, the columns containing the first three splines should be "1", "2", and "3". IF there is a fixed character prefix, that can be supplied via \code{namestem}. For example, if the columns are "s1", "s2", "s3", then set \code{namestem="s"}.
#' @seealso \code{\link{compute_lowCurve}}
#' @references Keller and Szpiro (2020). Selecting a scale for spatial confounding adjustment. Journal of the Royal Statistical Society, Series A  https://doi.org/10.1111/rssa.12556.
#' @references Rainey and Keller (2024). spconfShiny: An R Shiny application for calculating the spatial scale of smoothing splines for point data. PLOS ONE https://doi.org/10.1371/journal.pone.0311440
#' @export
#' @importFrom flexclust dist2
#' @examples
#' M <- 16
#' tprs_df <- 10
#' si <- seq(0, 1, length=M+1)[-(M+1)]
#' gridcoords <- expand.grid(x=si, y=si)
#' tprsX <- computeTPRS(coords = gridcoords, maxdf = tprs_df+1)
#' compute_effective_range(X=tprsX, coords=gridcoords, df=3:10, smoothedCurve=FALSE)
#'
#'
#' xloc <- runif(n=100, min=0, max=10)
#' X <- splines::ns(x=xloc, df=4, intercept=TRUE)
#' colnames(X) <- paste0("s", 1:ncol(X))
#' xplot <- 0:10
#' compute_effective_range(X=X, coords=as.matrix(xloc), df=2:4, newd=xplot,
#'                         namestem="s", smoothedCurve = TRUE)
#'
compute_effective_range <- function(X, coords=X[, c("x", "y")], df=3, nsamp=min(1000, nrow(X)), smoothedCurve = FALSE, newd=seq(0, 1, 100), scale_factor=1, returnFull=FALSE, cl=NULL,namestem="tprs", inds=NULL,verbose=FALSE, span=0.1){
    ngrid <- nrow(X)
    if (is.null(inds)){
        inds <- sample(ngrid, size=nsamp)
    }
    D <- flexclust::dist2(coords, coords[inds,])
    if(returnFull & smoothedCurve){
        out <- vector("list", length(df))
    } else {
        out <- numeric(length(df))
    }
    names(out) <- df
    if(!all(paste0(namestem, 1:max(df)) %in% colnames(X))) stop(paste0("Column names of X must take the form ", namestem, max(df)))
    for (k in seq_along(df)){
        if(verbose){
            cat("Df = ", df[k], "\n")
        }
        out[[k]] <- compute_effective_range_nochecks(X=X[, paste0(namestem, 1:df[k]), drop=FALSE], inds=inds, smoothedCurve = smoothedCurve, newd=newd, D=D, scale_factor=scale_factor, returnFull=returnFull, cl=cl, span=span)
    }
    out
}

#' @rdname compute_effective_range
#' @param inds Indices of observations to use for computation. Passed to \code{\link{computeS}}.
#' @param D Distance matrix for coordinates.
#' @importFrom stats median
#' @export
compute_effective_range_nochecks <- function(X, inds, newd, D, smoothedCurve = FALSE, scale_factor=1, returnFull=FALSE, cl=NULL, span=0.1){
    S <- computeS(X, inds=inds)
    if(smoothedCurve){
        SCurve <- compute_lowCurve(S=S, D=D, newd=newd, cl=cl, span=span)
        out <-  find_first_zero_cross(x = SCurve$SCurveMedian)*scale_factor
        if (returnFull){
            out <- list(range=out,
                        curve_median=SCurve$SCurveMedian,
                        curve_mean=SCurve$SCurveMean)
        }
    }else{
        zeros <-  find_zeros_cross(D = D, S = S)
        out <- stats::median(zeros, na.rm =T)*scale_factor
    }

    return(out)
}


