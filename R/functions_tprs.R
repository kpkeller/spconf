# Functions in File #
#
# arrangeTPRS
# computeTPRS


#' @title Create TPRS basis
#' @description Compute TPRS basis for given spatial coordinates
#' @param coords Data frame containing the coordinates.
#' @param maxdf Largest number of splines to include in TPRS basis
#' @param rearrange Logical indicator of whether to rearrange the columns of TPRS basis.
#' @param intercept Logical indicator of whether or not to remove the intercept column from the basis when \code{rearrange} is \code{TRUE}.
#' @return An \eqn{n}-by-\eqn{k} matrix of spline basis functions where \eqn{n} is the number of rows in \code{coords} and \eqn{k} is equal to \code{maxdf}
#' @details \code{computeTPRS} creates a thin-plate regression spline (TPRS) basis from a two-dimensional set of coordinate locations using the \code{mgcv} package.
#'
#'  The output from \code{mgcv} is structured to have the linear terms as the last columns of the matrix. Use \code{arrangeTPRS()} to arrange the matrix columns to be in order of increasing resolution. Specifically, it moves the last two columns to the left of the matrix and the third-from last column, which corresponds to the intercept, is optionally removed.
#'
#' @importFrom mgcv smoothCon s PredictMat
#' @export
#' @examples
#' x <- runif(100)
#' y <- runif(100)
#' mat <- computeTPRS(data.frame(x, y), maxdf=4)
computeTPRS <- function(coords, maxdf, rearrange=TRUE, intercept = FALSE){
    x1 <- x2 <- NULL # to avoid warnings
    colnames(coords) <- c('x1', 'x2')
    tprsSC <- mgcv::smoothCon(mgcv::s(x1, x2, fx=T, k=maxdf+1), data=coords)
    tprsX <- mgcv::PredictMat(tprsSC[[1]], data=coords)
    if (rearrange){
        tprsX <- arrangeTPRS(tprs = tprsX, intercept = intercept)
    }

    if (!intercept & rearrange){
        colnames(tprsX) <- paste0("tprs", 1:maxdf)
    } else {
        colnames(tprsX) <- paste0("tprs", 1:(maxdf+1))
    }
    tprsX
}


#' @rdname computeTPRS
#' @param tprs Matrix of TPRS basis values (from \code{computeTPRS}).
#' @export
arrangeTPRS <- function(tprs, intercept=FALSE){
    if (intercept){
        inds <- c(ncol(tprs) + -2:0, 1:(ncol(tprs)-3))
    } else {
        inds <- c(ncol(tprs) + -1:0, 1:(ncol(tprs)-3))
    }
    tprs <- tprs[, inds]
    colnames(tprs) <- 1:ncol(tprs)
    tprs
}

