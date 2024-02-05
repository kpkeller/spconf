
## Added LOESS flag in the function call and when calling function compute_effective_range_nochecks
compute_effective_range <- function(X, coords=X[, c("x", "y")], df=3, nsamp=min(1000, nrow(X)), LOESS = FALSE, newd=seq(0, 1, 100), scale_factor=1, returnFull=FALSE, cl=NULL,namestem="", inds=NULL,verbose=TRUE, span=0.1){
    ngrid <- nrow(X)
    if (is.null(inds)){
        inds <- sample(ngrid, size=nsamp)
    }
    dgrid <- flexclust::dist2(coords, coords[inds,])
    if(returnFull & LOESS){
        out <- vector("list", length(df))
    } else {
        out <- numeric(length(df))
    }
    names(out) <- df
    if(!all(paste0(namestem, 1:max(df)) %in% colnames(X))) stop(paste0("Column names of X must take the form ", namestem, max(df)))
    for (k in seq_along(df)){
        cat("Df = ", df[k], "\n")
        out[[k]] <- compute_effective_range_nochecks(X=X[, paste0(namestem, 1:df[k]), drop=FALSE], inds=inds, newd=newd, dgrid=dgrid, LOESS = LOESS, scale_factor=scale_factor, returnFull=returnFull, cl=cl, span=span)
    }
    out
}

## Added new computing option
compute_effective_range_nochecks <- function(X, inds, newd, dgrid, LOESS = FALSE, scale_factor=1, returnFull=FALSE, cl=NULL, span=0.1){
    S <- computeS(X, inds=inds)
    if(LOESS){
        SCurve <- compute_lowCurve(S=S, dgrid=dgrid, newd=newd, cl=cl, span=span)
        out <-  find_first_zero_cross(x = SCurve$SCurveMedian, LOESS = LOESS)*scale_factor
    }else{
        out <-  find_first_zero_cross(dgrid = dgrid, S = S, LOESS = LOESS)
    }
    if (returnFull & LOESS){
        out <- list(range=out,
                    curve_median=SCurve$SCurveMedian,
                    curve_mean=SCurve$SCurveMean)
    }
    return(out)
}

## Added new computing option

find_first_zero_cross <- function(x=NA, dgrid=NA, S=NA, LOESS=FALSE){

    if(!LOESS){
        cross_vals <- rep(0, ncol(dgrid))

        for(i in 1:ncol(dgrid)){
            x <- cbind(dgrid[,i], S[,i])
            x <- x[order(x[,1]),]

            if(min(x[,2]) >= 0){
                cross_vals[i] <- NA
            }else{
                upper_ind <- min(which(x[,2] < 0))
                # lower_ind <- upper_ind - 1
                # frac_dist <- x[lower_ind,2]/(x[lower_ind,2] - x[upper_ind,2])
                # cross_vals[i] <- x[lower_ind,1] + (x[upper_ind,1] - x[lower_ind,1])*frac_dist
                cross_vals[i] <- x[upper_ind,1]
            }
            out <- median(cross_vals, na.rm =T)
        }
    }else{
        upper_ind <- min(which(x<0))
        lower_ind <- upper_ind- 1
        frac_ind <- x[lower_ind]/(x[lower_ind] - x[upper_ind])
        out <- lower_ind + frac_ind
    }

    return(out)
}
