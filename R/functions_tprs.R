# Functions in File #
#
# arrangeTPRS
# computeTPRS

#' @title Rearrange TPRS basis
#' @description Rearrange the columns of the TPRS basis
#' @param tprs Matrix of TPRS basis values
#' @param intercept Logical indicator of whether or not to remove the intercept column
#' @details This function takes TPRS basis values, as created by the \code{mgcv} package, and rearranges them. The last two columns are moved to the left of the matrix and the third-from last column, which corresponds to the intercept, is optionally removed.
#' @seealso \code{\link{computeTPRS}}
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

#' @title Create TPRS basis
#' @description Compute TPRS basis for given spatial coordinates
#' @param coords Data frame containing the coordinates.
#' @param maxdf Largest number of splines to include in TPRS basis
#' @param rearrange Logical indicator of whether to rearrange the columns of TPRS basis
#' @param intercept Logical indicator of whether or not to remove the intercept column
#' @details This function creates a TPRS basis using the \code{mgcv} package from the given coordinates with the option to rearrange the columns such that last two columns are moved to the left of the matrix and the third-from last column, which corresponds to the intercept, is optionally removed.
#' @importFrom mgcv smoothCon s PredictMat
#' @seealso \code{\link{arrangeTPRS}}
#' @export
#' @examples
#' x <- runif(100)
#' y <- runif(100)
#' mat <- computeTPRS(data.frame(x, y), maxdf=4)
computeTPRS <- function(coords, maxdf, rearrange=TRUE, intercept = FALSE){
    x1 = x2 = NULL
    max_tprsdf <- min(ceiling(0.75*nrow(coords)), maxdf)
    colnames(coords) <- c('x1', 'x2')
    tprsSC <- mgcv::smoothCon(mgcv::s(x1, x2, fx=T, k=max_tprsdf+1), data=coords)
    tprsX <- mgcv::PredictMat(tprsSC[[1]], data=coords)
    if (rearrange){
        tprsX <- arrangeTPRS(tprs = tprsX, intercept = intercept)
    }
    return(list(tprsX = tprsX, df = max_tprsdf))
}
