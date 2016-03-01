#' Subset detection by clustering
#' 
#' Apply clustering algorithms to detect cell subsets. \code{densVM} and \code{densityClustX} clustering is 
#' based on the transformend ydata; \code{Rphenograph} is directly applied on the high dimemnional xdata. And 
#' \code{densVM} need the xdata to train the VM model.
#' 
#' @param ydata a matrix of the dimension reduced(transformed) data
#' @param xdata a matrix of the expression data
#' @param method cluster method including \code{densVM}, \code{densityClustX} and \code{Rphenograph}.
#' @return a vector of the clusters assigned for each row of the ydata
#' @export
#' @examples
#' d<-system.file('extdata', package='cytofkit')
#' fcsFile <- list.files(d, pattern='.fcs$', full=TRUE)
#' xdata <- cytof_exprsMerge(fcsFile, mergeMethod = 'fixed', fixedNum = 100)
#' ydata <- cytof_dimReduction(xdata, method = "tsne")
#' clusters <- cytof_cluster(ydata, xdata, method = "densVM")
cytof_cluster <- function(ydata = NULL, xdata = NULL, method = "densVM") {
    
    if(method == "densVM"){
        cat("  DensVM...")
        clusters <- densVM(ydata, xdata)$cluster$cluster
    } else if (method == "ClusterX"){
        cat("  ClusterX...")
        clusters <- ClusterX(ydata, gaussian=TRUE, alpha = 0.001, detectHalos = FALSE)$cluster
    } else if(method == "Rphenograph"){
        cat("  PhenoGraph...")
        clusters <- as.numeric(membership(Rphenograph(xdata, k=30)))
    } else if(is.null(method)){
        return(NULL)
    } else{
        stop("clusterMethod [",method,"] doesn't exist in cytofkit!")
    }
    
    if( length(clusters) != ifelse(is.null(ydata), nrow(xdata), nrow(ydata)) ){
        message("Cluster is not complete, cluster failed, try other cluster method!")
        return(NULL)
    }else{
        if(!is.null(xdata) && !is.null(row.names(xdata))){
            names(clusters) <- row.names(xdata)
        }else if(!is.null(ydata) && !is.null(row.names(ydata))){
            names(clusters) <- row.names(ydata)
        }
        cat(" DONE!\n")
        return(clusters)
    }
}

