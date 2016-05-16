#' cytofkit: an integrated analysis pipeline for mass cytometry data
#' 
#' This package is designed to facilitate the analysis workflow of mass cytometry data with 
#' automatic subset identification and mapping of cellular progression. Both command line and 
#' a GUI client are provided for executing the workflow easily.
#' 
#' This package integrates merging methods of multiple FCS files, dimension reduction methods (PCA, t-SNE and ISOMAP) 
#' and clustering methods (DensVM, densityClustX, and Rphenograph) for rapid subset detection. 
#' Cell subsets can be visualized in scatter plot and heat map. The method isomap is also provided to map the cellular progression. 
#' This workflow can be easily executed with the main function \code{\link{cytofkit}} or through the GUI client \code{\link{cytofkit_GUI}}.
#' 
#' Pre-processing
#' 
#' Using function \code{\link{cytof_exprsMerge}}, one or multiple FCS files will be loaded via the *read.FCS* 
#' function in the *flowCore* package. Then transformation was applied to the expression value 
#' of selected markers of each FCS file. Transformation methods include \code{auto_lgcl}, \code{fixed_lgcl}, 
#' \code{arcsin} and \code{biexp}, where \code{auto_lgcl} is the default.Then mutilple FCS files are 
#' merged using method \code{all}, \code{min}, \code{fixed} or \code{ceil}.
#' 
#' Dimensionality reduction
#' 
#' Using function \code{\link{cytof_dimReduction}}, t-Distributed Stochastic Neighbor Embedding (\code{tsne}) 
#' is suggested for dimensionality reduction although we also provide methods like \code{isomap} and \code{pca}.
#' 
#' Cluster 
#' 
#' Using function \code{\link{cytof_cluster}}, three cluster method are provided, \code{densVM}, \code{ClusterX}
#' and \code{Rphenograph}. \code{densVM}, \code{densityClustX} are performend on the dimension reducted data, while \code{Rphenograph}
#' is peroformed directed on the high dimensional expression data. 
#' 
#' Post-processing
#' 
#'  - Using function \code{\link{cytof_clusterPlot}} to visualize the cluster results in a catter plot, in which dots represent cells, colours 
#' indicate their assigned clusters and point shapes represent their belonging samples.
#' 
#'  - Using function \code{\link{cytof_heatmap}} to generate heat map to visualize the mean expression of every marker in every cluster. 
#'  This heat maps is useful to interrogate marker expression to identify each cluster's defining markers. 
#' 
#'  - Using function \code{\link{cytof_progressionPlot}} to visualize the expression patter of selected markers against the estimated 
#'  cellular progression order.
#'  
#'  - Using function \code{\link{cytof_addToFCS}} to add any dimension reduced data, cluster results, progression data into the original FCS files,
#'  new FCS files will be saved for easy checking with other softwares like FlowJo. 
#' 
#' All the above post processing can be automatically implemented and saved using one function \code{\link{cytof_writeResults}}.
#' 
#' 
#' @examples
#' 
#' ## Run on GUI
#' #cytofkit_GUI()  # remove the hash symbol to launch the GUI
#' 
#' ## Run on command
#' dir <- system.file('extdata',package='cytofkit')
#' file <- list.files(dir, pattern='.fcs$', full=TRUE)
#' parameters <- list.files(dir, pattern='.txt$', full=TRUE)
#' ## remove the hash symbol to run the following command
#' #cytofkit(fcsFile = file, markers = parameters, projectName = 'test')   
#' 
#' ## Checking the vignettes for more details 
#' if(interactive()) browseVignettes(package = 'cytofkit')
#' 
#' @seealso \code{\link{cytofkit}}, \code{\link{cytofkit_GUI}}
#' @references \url{http://signbioinfo.github.io/cytofkit/}
#' @docType package
#' @name cytofkit-package
#' 
NULL



#' cytofkit: an integrated analysis pipeline for mass cytometry data
#' 
#' A user friendly GUI is provided for easy usage of cytofkit, \code{\link{cytofkit_GUI}}.
#' 
#' \code{cytofkit} provides a workflow for one or multiple CyTOF data analysis, 
#' including data preprocess with merging methods of multiple fcs file, expression data transformation, 
#' dimension reduction with PCA, isomap or tsne(default), clustering methods(densVM, ClusterX, Rphenograph) 
#' for subpopulation detection, and estimation of cellular progression with isomap. The analysis results 
#' can be visualized with scatter plot, heatmap plot or progression plot. Moreover theses results can be saved back to 
#' FCS files. By default the results will be automatically saved for further annotation. An interactive web application is
#' provided for interactive exploration of the analysis results, \code{cytofkitShinyAPP}.
#' 
#' 
#' @param fcsFiles it can be either the name of the path where stores your FCS files or a vector of FCS file names. 
#' @param markers it can be either a text file that specifies the makers to be used for analysis or a vector of the marker names.
#' @param projectName a prefix that will be added to the names of result files.
#' @param mergeMethod when multiple fcs files are selected, cells can be combined using 
#' one of the four different methods including \code{ceil}, \code{all}, \code{min}, \code{fixed}. 
#' The default option is \code{ceil}, up to a fixed number (specified by \code{fixedNum}) of cells are sampled 
#' without replacement from each fcs file and combined for analysis.
#' \code{all}: all cells from each fcs file are combined for analysis. 
#' \code{min}: The minimum number of cells among all the selected fcs files are sampled from each fcs file and combined for analysis. 
#' \code{fixed}: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than 
#' fixedNum) from each fcs file and combined for analysis.
#' @param fixedNum up to fixedNum of cells from each fcs file are used for analysis.
#' @param ifCompensation Boolean value to decide if do compensation. This will be applied to flow cytometry data.
#' @param transformMethod dat transformation method, either \code{auto_lgcl}, \code{fixed_lgcl}, \code{arcsin} or \code{biexp}.
#' @param dimReductionMethod the method used for dimensionality reduction, including \code{tsne}, \code{pca} and \code{isomap}.
#' @param clusterMethods the clustering method(s) used for subpopulation detection, including \code{densVM}, \code{ClusterX} and \code{Rphenograph}. Multiple selection are accepted.
#' @param visualizationMethods the method(s) used for visualize the cluster data, including \code{tsne}, \code{pca} and \code{isomap}. Multiple selection are accepted.
#' @param progressionMethod use the first ordination score of \code{isomap} to estimated the preogression order of cells, choose \code{NULL} to ignore.
#' @param uniformClusterSize the uniform size of each cluster.
#' @param resultDir the directory where result files will be generated.
#' @param saveResults if save the results, and the post-processing results including scatter plot, heatmap, and statistical results.
#' @param saveObject save the resutls into RData objects for loading back to R for further analysis
#' @param saveToFCS save the results back to the FCS files, new FCS files will be generated.
#' @param scaleTo scale the expression values to the same scale after transformation, default is NULL, should be a vector of two numbers if scale.
#' @param q quantile of negative values removed for auto w estimation in logicle transformation, default is 0.05.
#' @param ... more arguments contral the logicle transformation
#' 
#' @return a list containing \code{expressionData}, \code{dimReductionMethod}, \code{visualizationMethods}, \code{dimReducedRes}, \code{clusterRes} and \code{progressionRes}. If choose 'saveResults = TRUE', results will be saved into files under \code{resultDir}
#' @author Chen Jinmiao, Chen Hao
#' @references \url{http://signbioinfo.github.io/cytofkit/}
#' @seealso \code{\link{cytofkit}}, \code{\link{cytofkit_GUI}}
#' @useDynLib cytofkit
#' @export
#' @examples
#' dir <- system.file('extdata',package='cytofkit')
#' file <- list.files(dir, pattern='.fcs$', full=TRUE)
#' parameters <- list.files(dir, pattern='.txt$', full=TRUE)
#' ## remove the hash symbol to run the following command
#' #cytofkit(fcsFile = file, markers = parameters, projectName = 'test')   
cytofkit <- function(fcsFiles = getwd(), markers = NULL, 
                     projectName = "cytofkit", 
                     mergeMethod = "ceil", 
                     fixedNum = 10000, 
                     ifCompensation = FALSE, 
                     transformMethod = "auto_lgcl", 
                     dimReductionMethod = "tsne", 
                     clusterMethods = "ClusterX", 
                     visualizationMethods = "tsne", 
                     progressionMethod = NULL, 
                     uniformClusterSize = 500,
                     resultDir = getwd(), 
                     saveResults = TRUE, saveObject = TRUE, saveToFCS = TRUE, 
                     scaleTo = NULL, q = 0.05, ...) {
    
    
    ## arguments checking
    if (is.null(fcsFiles)) {
        fcsFiles <- list.files(path = getwd(), pattern = ".fcs$", 
            full.names = TRUE)
        rawFCSdir <- getwd()
    } else if (length(fcsFiles) == 1 && file.info(fcsFiles)$isdir) {
        fcsFiles <- list.files(path = fcsFiles, pattern = ".fcs$", 
            full.names = TRUE)
        rawFCSdir <- fcsFiles
    } else{
        if(dirname(fcsFiles[1]) == "."){
            rawFCSdir <- getwd()
        }else{
            rawFCSdir <- dirname(fcsFiles[1])  
        }
    }
    
    if (is.null(fcsFiles) || length(fcsFiles) < 1) {
        stop("No FCS file found, please select your fcsFiles!")
    } else if (!all(file.exists(fcsFiles))) {
        stop("Can not find file(s):", fcsFiles[which(!file.exists(fcsFiles))])
    }
    
    if (is.null(markers)) 
        markers <- "parameter.txt"
    
    if (length(markers) == 1 && file.exists(markers)) {
        markers <- as.character(read.table(markers, sep = "\t", 
                                           header = TRUE)[, 1])
    }
        
    if (is.null(markers) || length(markers) < 1) 
        stop("no marker selected!")
    
    if (!is.null(mergeMethod) && !(mergeMethod %in% c("ceil", "all", "min", "fixed"))) 
        stop("mergeMethod doesn't exist in cytofkit!")
    
    if (!is.null(transformMethod) && !(transformMethod %in% c("auto_lgcl", "fixed_lgcl", "arcsin", 
        "biexp"))) 
        stop("transformMethod doesn't exist in cytofkit!")
    
    if (!is.null(dimReductionMethod) && !(dimReductionMethod %in% c("tsne", "pca", "isomap"))) 
        stop("dimReductionMethod doesn't exist in cytofkit!")
    
    ## force dimReductionMethod to tSNE
    if(length(dimReductionMethod) > 1 || is.null(dimReductionMethod)) 
        dimReductionMethod <- "tsne"
    
    if (!is.null(clusterMethods) && !all(clusterMethods %in% c("densVM", 
        "ClusterX", "Rphenograph"))) 
        stop("clusterMethods doesn't exist in cytofkit!")
    
    if (is.null(visualizationMethods)) {
        visualizationMethods <- dimReductionMethod
    } else if (!is.null(visualizationMethods) && !(all(visualizationMethods %in% 
        c("tsne", "pca", "isomap")))) 
        stop("visualizationMethods doesn't exist in cytofkit!")
    
    if (!is.null(progressionMethod) && !(all(progressionMethod %in% 
        c("tsne", "pca", "isomap")))) 
        stop("progressionMethod doesn't exist in cytofkit!")
    
    if (!is.null(uniformClusterSize) && !(is.numeric(uniformClusterSize))) 
        stop("uniformClusterSize must be a numeric number!")
    
    
    ## print arguments for user info
    message("Input arguments:")
    cat("* Input FCS files for analysis:\n ")
    cat(paste0("  -", basename(fcsFiles), "\n"))
    cat("* Makrers:\n ")
    cat(paste0("  -", markers, "\n"))
    cat("* File merging method: ")
    cat(mergeMethod, "\n")
    cat("* Data transformation method: ")
    cat(transformMethod, "\n")
    cat("* Dimensional reduction method: ")
    cat(dimReductionMethod, "\n")
    cat("* Data clustering method(s): ")
    cat(clusterMethods, "\n")
    cat("* Data visualization method(s): ")
    cat(visualizationMethods, "\n")
    cat("* Subset progression analysis method: ")
    cat(progressionMethod, "\n\n")
    
    
    ## get marker-filtered, transformed, combined exprs data
    message("Extract expression data...")
    ## match_markers to extract expression of all markers(saved for visualization in shiny web app)
    all_marker_names <- match_markers(fcsFiles[1])
    
    exprs_data_all <- tryCatch(
        {
            cytof_exprsMerge(fcsFiles, comp = ifCompensation, verbose = FALSE, 
                             markers = all_marker_names, transformMethod = transformMethod, 
                             scaleTo = scaleTo, q = q, mergeMethod = mergeMethod, 
                             fixedNum = fixedNum)
        }, error=function(cond) {
            message("jump")
            return(NULL)
        }
    )    
    
    if(is.null(exprs_data_all)){
        exprs_data_all <- cytof_exprsMerge(fcsFiles, comp = ifCompensation, verbose = FALSE, 
                                   markers = markers, transformMethod = transformMethod, 
                                   scaleTo = scaleTo, q = q, mergeMethod = mergeMethod, 
                                   fixedNum = fixedNum)
        exprs_data <- exprs_data_all
    }else{
        exprs_data <- cytof_exprsMerge(fcsFiles, comp = ifCompensation, verbose = FALSE, 
                                       markers = markers, transformMethod = transformMethod, 
                                       scaleTo = scaleTo, q = q, mergeMethod = mergeMethod, 
                                       fixedNum = fixedNum)
    }
    cat("  ", nrow(exprs_data), " x ", ncol(exprs_data), " data was extracted!\n")
    
    
    ## dimension reduced data, a list
    message("Dimension reduction...")
    alldimReductionMethods <- unique(c(visualizationMethods, dimReductionMethod))
    allDimReducedList <- lapply(alldimReductionMethods, 
                                cytof_dimReduction, data = exprs_data)
    names(allDimReducedList) <- alldimReductionMethods
    
    
    ## cluster results, a list
    message("Run clustering...")
    head(allDimReducedList[[dimReductionMethod]])
    head(exprs_data)
    cluster_res <- lapply(clusterMethods, cytof_cluster, 
                          ydata = allDimReducedList[[dimReductionMethod]], 
                          xdata = exprs_data)
    names(cluster_res) <- clusterMethods
    
    
    ## progression analysis results, a list  
    ## NOTE, currently only the first cluster method resutls 
    ## are used for preogression visualization(by default: cluster_res[[1]])
    message("Progression analysis...")   
    if(!(is.null(uniformClusterSize)) && !is.null(progressionMethod) && 
       progressionMethod %in% alldimReductionMethods){
           progression_res <- list(sampleData = exprs_data, 
                                   sampleCluster = cluster_res[[1]], 
                                   progressionData = allDimReducedList[[progressionMethod]])
    }else if(!is.null(progressionMethod)){
        progression_res <- cytof_progression(data = exprs_data, 
                                             cluster = cluster_res[[1]], 
                                             method = progressionMethod,
                                             uniformClusterSize = uniformClusterSize)
    }else{
        progression_res <- NULL
    }
    
    
    ## wrap the results
    analysis_results <- list(expressionData = exprs_data,
                             dimReductionMethod = dimReductionMethod,
                             visualizationMethods = visualizationMethods,
                             dimReducedRes = allDimReducedList,
                             clusterRes = cluster_res, 
                             progressionRes = progression_res,
                             allExpressionData = exprs_data_all)
     
    
    ## save the results
    message("Analysis DONE, saving the reuslts...") 
    if(saveObject){
        objFile <- paste0(resultDir, .Platform$file.sep, projectName, ".RData")
        save(analysis_results, file = objFile)
        cat("Analysis obejct is saved in ", objFile, "\n")
    }
    
    if (saveResults == TRUE) {
        cat("Writing results\n")
        cytof_writeResults(analysis_results = analysis_results, 
                           projectName=projectName, 
                           resultDir=resultDir,
                           saveToFCS = saveToFCS,
                           rawFCSdir=rawFCSdir)
    } else {
        return(analysis_results)
    }
}

# this is for save all markers for visualization in the shiny web app
# save all markers may cause error in transformation
match_markers <- function(fcsFile, markers=NULL){
    fcs <- suppressWarnings(read.FCS(fcsFile))
    pd <- fcs@parameters@data
    
    if (!(is.null(markers))) {
        right_marker <- markers %in% pd$desc || markers %in% pd$name
        if (!(right_marker)) {
            stop("\n Selected marker(s) is not in the input fcs files \n please check your selected markers! \n")
        } else {
            desc_id <- match(markers, pd$desc)
            name_id <- match(markers, pd$name)
            mids <- c(desc_id, name_id)
            marker_id <- unique(mids[!is.na(mids)])
        }
        return(marker_id)
    } else {
        pdname <- as.character(pd$name)
        tl_channel <- pdname %in% c("Time", "Event_length", "Length", "length", 
                                     "time", "event_length", "Event", "Event #") 
        return(pdname[!tl_channel])
    }
}





#' A Shiny app to interactively visualize the analysis results 
#' 
#' Load the RData object saved by cytofkit, explore the analysis results with interactive control
#'
#' @import shiny
#' @export
#' @examples 
#' if (interactive()) cytofkit::cytofkitShinyAPP()
cytofkitShinyAPP = function() {
    shiny::runApp(system.file('shiny', package = 'cytofkit'))
}



#' check the package update news
#' 
#' @export
cytofkitNews <- function() 
{
    newsfile <- file.path(system.file(package = "cytofkit"), 
                          "NEWS")
    file.show(newsfile)
}


