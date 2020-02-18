
#' @title Compute Smoothing Matrix
#' @description Calculates the smothing (or "hat") matrix from a design matrix.
#' @param x Matrix of spline values. A data frame is coerced into a matrix.
#' @param inds Indices to compute smoothing matrix for.
#' @details Given a matrix \code{X} of spline values, this computes  \code{S=X(X'X)^(-1)X'}. When \code{x} has many rows, this can be quite large. Note that \code{x} must be full-rank.
#'
#' @return  An \eqn{n}-by-\eqn{n} matrix, where \eqn{n} is the number of rows in \code{x}.
#' @export
#' @examples
#' # Simple design matrix case
#' X <- cbind(1, rep(c(0, 1), each=4))
#' S <- computeS(X)
#' # More complex example
#' xloc <- runif(n=100, min=0, max=10)
#' X <- splines::ns(x=xloc, df=4, intercept=TRUE)
#' S <- computeS(X)
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
#' @param span Controls the amount of smoothing. Passed to \code{\link[stats]{loess}}; see that function for details.
#' @param ... Additional arguments passed to \code{\link[stats]{loess}}
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
#' @param S Smoothing matrix.
#' @param dgrid Distance matrix.
#' @param newd Distances to use for loess prediction.
#' @param cl Cluster object, or number of cluster instances to create. Defaults to no parallelization.
#' @param span Passed to \code{\link{fitLoess}}
#' @seealso \code{\link{computeS}}
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
#' lC <- compute_lowCurve(S, dgrid=d, newd=xplot)
#' matplot(xplot, lC$SCurve, type="l", col="black")
#' points(xplot, lC$SCurveMedian, type="l", col="red")
compute_lowCurve <- function(S, dgrid, newd, cl=NULL, span=0.1){
    if (is.null(cl)){
        n <- nrow(S)
        lowCurveM <- matrix(NA, nrow=length(newd), ncol=n)
        for (i in 1:n){#nrow(S)) {
            lowCurveM[,i] <- fitLoess(y=S[i,], x=dgrid[, i], newx=newd, span=span)
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
        newfitLoess <- function(s, d) fitLoess(s, d, newx=newd, span=span)
        environment(newfitLoess) <- .GlobalEnv
        clusterExport(cl=cl, varlist=c("S_list", "dgrid_list", "fitLoess", "newd", "span"), envir=environment())
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



#' @title Compute effective range
#' @description Calculates the effective range for a spline basis matrix.
#' @param X Matrix of spline values. Currently assumed to have columns \code{x}, \code{y}, \code{s1},\code{s2}, ...
#' @param df Degrees of freedom to include.
#' @param nsamp Number of observations from \code{X} from which to sample. Defaults to minimum of 1,000 and \code{nrow(X)}.
#' @param newd Distance values at which to make loess predictions.
#' @param scale_factor Factor by which range should be scaled. Usually physical distance corresponding to resolution of grid.
#' @param returnFull Should the mean and median curves be returned (TRUE), or just the range value of where they first cross zero (FALSE).
#' @param cl Cluster object, or number of cluster instances to create. Defaults to no parallelization.
#' @param namestem Stem of names of columns of X corresponding to evaluated splines.
#' Defaults to \code{""}, meaning names of the form \code{1}, \code{2}, ...
#' @param inds Optional vector of indices to use as subset. If provided, \code{nsamp} is not used.
#' @param verbose Control message printing.
#' @param span Passed to \code{\link{fitLoess}}
#' @export
#' @importFrom flexclust dist2
#' @examples
#'
#' xloc <- runif(n=100, min=0, max=10)
#' X <- splines::ns(x=xloc, df=4, intercept=TRUE)
#' colnames(X) <- paste0("s", 1:ncol(X))
#' xplot <- 0:10
#' #er <- compute_effective_range(X=X, df=4, newd=xplot, namestem="s")
#'
#' # Simulation settings
#' M <- 16
#' tprs_df <- 10
#' si <- seq(0, 1, length=M+1)[-(M+1)]
#' gridcoords <- expand.grid(x=si, y=si)
#' tprsSC <- mgcv::smoothCon(mgcv::s(x, y, fx=TRUE, k=tprs_df + 1), data= gridcoords)
#' tprsX <- mgcv::PredictMat(tprsSC[[1]], data= gridcoords)
#' # Re-order the TPRS
#' tprsX <- tprsX[, c(ncol(tprsX) + -1:0, 1:(ncol(tprsX)-2))]
#' tprsX <- cbind(gridcoords ,tprsX)
#' compute_effective_range(tprsX)
compute_effective_range <- function(X, df=3, nsamp=min(1000, nrow(X)), newd=seq(0, 1, 100), scale_factor=1, returnFull=FALSE, cl=NULL,namestem="", inds=NULL,verbose=TRUE, span=0.1){
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
        out[[k]] <- compute_effective_range_nochecks(X=X[, paste0(namestem, 1:df[k]), drop=FALSE], inds=inds, newd=newd, dgrid=dgrid, scale_factor=scale_factor, returnFull=returnFull, cl=cl, span=span)
    }
    out
}

#' @rdname compute_effective_range
#' @param inds Indices of observations to use for computation. Passed to \code{\link{computeS}}.
#' @param dgrid Distance matrix.
#' @export
compute_effective_range_nochecks <- function(X, inds, newd, dgrid, scale_factor=1, returnFull=FALSE, cl=NULL, span=span){
    S <- computeS(X, inds=inds)
    SCurve <- compute_lowCurve(S=S, dgrid=dgrid, newd=newd, cl=cl, span=span)
    out <-  find_first_zero_cross(SCurve$SCurveMedian)*scale_factor
    if (returnFull){
        out <- list(range=out,
                    curve_median=SCurve$SCurveMedian,
                    curve_mean=SCurve$SCurveMean)
    }
    return(out)
}


