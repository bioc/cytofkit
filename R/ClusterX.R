#' Fast clustering by automaticly search and find of density peaks 
#' 
#' This package implement the clustering algorithm described by Alex Rodriguez
#' and Alessandro Laio (2014) with improvements of automatic peak detection and 
#' parallel implementation
#' 
#' @param data A data matrix for clustering.
#' @param dimReduction Dimenionality reduciton method.
#' @param outDim Number of dimensions will be used.
#' @param dc Distance cutoff value.
#' @param gaussian If apply gaussian to esitmate the density.
#' @param alpha Signance level for peak detection.
#' @param detectHalos If detect the halos.
#' @param parallel If run the algorithm in parallel.
#' @param nCore Number of cores umployed for parallel compution.
#' 
#' @return a object of \code{ClusterX} class
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom pdist pdist
#' @importFrom plyr llply
#' @export
#' 
#' @author Chen Hao
#' 
#' @examples
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' ClusterXRes <- ClusterX(data)
ClusterX <- function(data, dimReduction = NULL, outDim=2, dc, gaussian=TRUE, alpha = 0.001, 
                     detectHalos = FALSE, parallel = FALSE, nCore = 4) {
    res <- NULL
    
    message("ClusterX currently is not avilable, but will come soon...!")
    
    return(res)
}


