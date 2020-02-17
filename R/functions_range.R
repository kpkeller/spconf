
#' @title Compute Smoothing Matrix
#' @description Calculates the matrix \code{X(X'X)^(-1)X'}.
#' @param x Matrix of spline values. A data frame is coerced into a matrix.
#' @param inds Indices to compute smoothing matrix for.
#' @details When \code{x} has many rows, this can be quite large.
#' @return  An \eqn{n}-by-\eqn{N}, where \eqn{n} is the length of \code{inds} and \eqn{N} is the number of rows in \code{x}.
#' @export
computeS <- function(x, inds=1:nrow(x)){
    x <- as.matrix(x)
    S <- x[inds,] %*% solve(crossprod(x), t(x))
 S
}


#' @title Fit a loess curve
#' @description Wrapper function for fitting and predicting from \code{loess()}.
#' @param y Dependent variable values
#' @param x Independent variable values
#' @param newx Values of x to use for prediction.
#' @param span Passed to \code{\link[stats]{loess}}.
#' @param ... Additional arguments passed to \code{\link[stats]{loess}}
#' @export
#' @importFrom stats loess predict
fitLoess <- function(y, x, newx=x, span=0.1,...){
    lowcurve <- stats::loess(y~x, span=span,...)
    lowpred <- stats::predict(lowcurve, newdata=newx)
    lowpred
}


#' @title Compute loess curves for smoothing matrix
#' @description Calculates a loess curve for the smoothing matrix entries, as a function of distance between points.
#' @param S Smoothing matrix.
#' @param dgrid Distance matrix.
#' @param newd Distances to use for loess prediction.
#' @param cl Cluster object, or number of cluster instances to create. Defaults to no parallelization.
#' @seealso \code{\link{computeS}}
#' @export
#' @importFrom parallel makeCluster clusterExport clusterMap stopCluster
#' @importFrom stats median
compute_lowCurve <- function(S, dgrid, newd, cl=NULL){
    if (is.null(cl)){
        n <- nrow(S)
        lowCurveM <- matrix(NA, nrow=length(newd), ncol=n)
        for (i in 1:n){#nrow(S)) {
            lowCurveM[,i] <- fitLoess(y=S[i,], x=dgrid[, i], newx=newd)
        }
    } else  {
        stop_cluster <- FALSE
        if (is.numeric(cl) && cl>0) {
            cl <- makeCluster(cl)
            stop_cluster <- TRUE
        }  else  if (!inherits(cl, "cluster")){
            stop("'cl' should be NULL, a cluster object, or a non-negative number indicating number of cores.")
        }
        S_list <- as.list(as.data.frame(t(S)))
        dgrid_list <- as.list(as.data.frame(dgrid))
        newfitLoess <- function(s, d) fitLoess(s, d, newx=newd)
        environment(newfitLoess) <- .GlobalEnv
        clusterExport(cl=cl, varlist=c("S_list", "dgrid_list", "fitLoess", "newd"), envir=environment())
        lowCurveM <- clusterMap(cl=cl,newfitLoess , s=S_list, d=dgrid_list, SIMPLIFY=TRUE)
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
#' @title Find first zero
#' @description Rudimentary root finding function to calculate first cross of y-axis
#' @param x Function values, assumed to be ordered
#' @return Index of first value of \code{x} that lies below 0. Decimal values will be returned using a basic interpolation of the two values straddling 0.
#' @export
find_first_zero_cross <- function(x){
    upper_ind <- min(which(x<0))
    lower_ind <- upper_ind- 1
    frac_ind <- x[lower_ind]/(x[lower_ind] - x[upper_ind])
    lower_ind + frac_ind
}

#
#
# Add Details here.....
#
#
#' @title Compute effective range
#' @description Calculates the effective range for a spline basis matrix.
#' @param X Matrix of spline values. Currently assumed to have columns \code{x}, \code{y}, \code{s1},\code{s2}, ...
#' @param df Degrees of freedom to include.
#' @param nsamp Number of observations from \code{X} from which to sample. Defaults to minimum of 1,000 and \code{nrow(X)}.
#' @param newd Distance values at which to make loess predictions.
#' @param scale_factor Factor by which range should be scaled. Usually physical distance corresponding to resolution of grid.
#' @param returnFull Should the mean and median curves be returned (TRUE), or just the range value of where they first cross zero (FALSE).
#' @param cl Cluster object, or number of cluster instances to create. Defaults to no parallelization.
#' @param namestem Stem of names of columns of X corresponding to evaluated splines. Defaults to \code{"s"}, meaning
#' names of the form \code{s1}, \code{s2}, ...
#' @param inds Optional vector of indices to use as subset. If provided, \code{nsamp} is not used.
#' @param verbose Control message printing.
#' @export
#' @importFrom flexclust dist2
compute_effective_range <- function(X, df=3, nsamp=min(1000, nrow(X)), newd=seq(0, 100, 1), scale_factor=1, returnFull=FALSE, cl=NULL,namestem="s", inds=NULL,verbose=TRUE){
    ngrid <- nrow(X)
    if (is.null(inds)){
        inds <- sample(ngrid, size=nsamp)
    }
    dgrid <- flexclust::dist2(X[, c("x", "y")], X[inds, c("x", "y")])
    if(returnFull){
        out <- vector("list", length(df))
    } else {
        out <- numeric(length(df))
    }
    names(out) <- df
    if(!all(paste0(namestem, 1:max(df)) %in% names(X))) stop(paste0("Names of X must take the form ", namestem, max(df)))
    for (k in seq_along(df)){
        cat("Df = ", df[k], "\n")
        out[[k]] <- compute_effective_range_nochecks(X=X[, paste0(namestem, 1:df[k]), drop=FALSE], inds=inds, newd=newd, dgrid=dgrid, scale_factor=scale_factor, returnFull=returnFull, cl=cl)
    }
    out
}

#' @rdname compute_effective_range
#' @param inds Indices of observations to use for computation. Passed to \code{\link{computeS}}.
#' @param dgrid Distance matrix.
#' @inheritParams compute_effective_range
#' @export
compute_effective_range_nochecks <- function(X, inds, newd, dgrid, scale_factor=1, returnFull=FALSE, cl=NULL){
    S <- computeS(X, inds=inds)
    SCurve <- compute_lowCurve(S=S, dgrid=dgrid, newd=newd, cl=cl)
    out <-  find_first_zero_cross(SCurve$SCurveMedian)*scale_factor
    if (returnFull){
        out <- list(range=out,
                    curve_median=SCurve$SCurveMedian,
                    curve_mean=SCurve$SCurveMean)
    }
    return(out)
}


