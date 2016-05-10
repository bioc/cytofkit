#' Dimension reduction of cytof expression data 
#' 
#' Apply dimension reduction on the cytof expression data, 
#' with method \code{isomap}, \code{pca}, or \code{tsne}. 
#' 
#' @param data An expression data matrix.
#' @param method Method chosed for dimensition reduction, must be one of \code{isomap}, \code{pca} or \code{tsne}. 
#' @param distMethod Method for distance calcualtion.
#' @param out_dim The dimensionality of the output.
#' @param isomap_k Number of shortest dissimilarities retained for a point, parameter for \code{isomap} method.
#' @param isomap_ndim Number of axes in metric scaling, parameter for \code{isomap} method.
#' @param isomapFragmentOK What to do if dissimilarity matrix is fragmented, parameter for \code{isomap} method.
#' @return a matrix of the dimension reducted data, with colnames and rownames(if have, same as the input).
#' @author Chen Jinmiao
#' @importFrom vegan vegdist spantree isomap
#' @importFrom Rtsne Rtsne
#' @import stats
#' @export
#' @examples
#' data(iris)
#' in_data <- iris[, 1:4]
#' out_data <- cytof_dimReduction(in_data)
cytof_dimReduction <- function(data, method = "tsne", distMethod = "euclidean", out_dim = 2, isomap_k = 5, isomap_ndim = NULL, isomapFragmentOK = TRUE) {
    
    rnames <- row.names(data)
    data <- as.matrix(data)
    
    if (method == "pca") {
        cat("  PCA...")
        mapped <- prcomp(data, scale = TRUE)$x
    } else if (method == "isomap") {
        cat("  ISOMAP...")
        if (is.null(isomap_ndim)) {
            isomap_ndim <- ncol(data)
        }
        
        ord <- tryCatch(
            {
                dis <- vegdist(data, method = distMethod)
                isomap(dis, ndim = isomap_ndim, k = isomap_k, fragmentedOK = isomapFragmentOK)
            }, error=function(cond) {
                message("Runing isomap failed")
                message("Here's the original error message:")
                message(cond)
                return(NULL)
            }
        )    
        
        if(is.null(ord)){
            mapped <- NULL
        }else{
            mapped <- ord$points
        }
        
    } else if (method == "tsne") {
        cat("  t-SNE...")
        tsne_out <- Rtsne(as.matrix(data), initial_dims = dim(as.matrix(data))[2], 
            dims = 2, perplexity = 30, theta = 0.5, check_duplicates = FALSE, 
            pca = TRUE)
        mapped <- tsne_out$Y
    } else if (is.null(method)){
        return(NULL)
    } else {
        stop("dimReductionMethod [", method, "] doesn't exit for cytofkit!")
    }
    
    ## organize output
    if(ncol(mapped) < out_dim){
        out_dim <- ncol(mapped)
        message("Run ",method," for dimensional reduction, out dimension coerced to ",out_dim)
    }
    mapped <- mapped[ ,c(1:out_dim)]
    colnames(mapped) <- paste(method, c(1:out_dim), sep = "_")
    rownames(mapped) <- rnames
    cat("DONE\n")
    return(mapped)
} 
