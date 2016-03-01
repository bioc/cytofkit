#' Merge the transformed expression data of FCS file(s) of selected markers
#' 
#' Apply transformation of selected markers of each FCS file, arcsin, biexponential, auto logicle
#' transformation and fixed logicle transformation are provided, then mutilple 
#' FCS files are merged using method \code{all}, \code{min}, \code{fixed} or \code{ceil}
#' 
#' @param fcsFiles the input fcsFiles (usually more than 1 file)
#' @param comp Boolean value tells if do compensation
#' @param verbose Boolean value detecides if print the massage details
#' @param markers Selected markers for analysis, either from names or from description
#' @param transformMethod transformation method, \code{auto_lgcl}, \code{fixed_lgcl}, \code{arcsin} or \code{biexp}
#' @param scaleTo scale the expression to same scale, default is NULL, should be a vector of two numbers if scale
#' @param mergeMethod merge method for mutiple FCS expression data. cells can be combined using 
#' one of the four different methods including \code{ceil}, \code{all}, \code{min}, \code{fixed}. 
#' The default option is \code{ceil}, up to a fixed number (specified by \code{fixedNum}) of cells are sampled 
#' without replacement from each fcs file and combined for analysis.
#' \code{all}: all cells from each fcs file are combined for analysis. 
#' \code{min}: The minimum number of cells among all the selected fcs files are sampled from each fcs file and combined for analysis. 
#' \code{fixed}: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than 
#' fixedNum) from each fcs file and combined for analysis.
#' @param fixedNum the fixed number of cells for merging multiple FCSs
#' @param w Linearization width in asymptotic decades
#' @param t Top of the scale data value
#' @param m Full width of the transformed display in asymptotic decades
#' @param a Additional negative range to be included in the display in asymptotic decades
#' @param q quantile of negative values removed for auto w estimation, default is 0.05
#' 
#' @return Merged FCS expression data matrix of selected markers after transformation
#' @export
#' @examples
#' d<-system.file('extdata',package='cytofkit')
#' fcsFiles <- list.files(d,pattern='.fcs$',full=TRUE)
#' merged <- cytof_exprsMerge(fcsFiles)
cytof_exprsMerge <- function(fcsFiles, comp = FALSE, verbose = FALSE, 
    markers = NULL, transformMethod = "auto_lgcl", scaleTo = NULL, 
    mergeMethod = "ceil", fixedNum = 10000, w = 0.1, t = 4000, 
    m = 4.5, a = 0, q = 0.05) {
    
    exprsL <- mapply(cytof_exprsExtract, fcsFiles, MoreArgs = list(comp = comp, 
        verbose = verbose, markers = markers, transformMethod = transformMethod, 
        scaleTo = scaleTo, w = w, t = t, m = m, a = a, q = q), 
        SIMPLIFY = FALSE)
    
    if (mergeMethod == "all") {
        merged <- do.call(rbind, exprsL)
    } else if (mergeMethod == "min") {
        minSize <- min(sapply(exprsL, nrow))
        mergeFunc <- function(x) {
            x[sample(nrow(x), size = minSize, replace = FALSE), ]
        }
        merged <- do.call(rbind, lapply(exprsL, mergeFunc))
    } else if (mergeMethod == "fixed") {
        mergeFunc <- function(x) {
            x[sample(nrow(x), size = fixedNum, replace = ifelse(nrow(x) < 
                fixedNum, TRUE, FALSE)), ]
        }
        merged <- do.call(rbind, lapply(exprsL, mergeFunc))
    } else if (mergeMethod == "ceil") {
        mergeFunc <- function(x) {
            if (nrow(x) < fixedNum) {
                x
            } else {
                x[sample(nrow(x), size = fixedNum, replace = FALSE), 
                  ]
            }
        }
        merged <- do.call(rbind, lapply(exprsL, mergeFunc))
    } else {
        stop("mergeMethod [", mergeMethod, "] doesn't exit for cytofkit!")
    }
    
    return(merged)
}


#' Extract the expression matrix of the FCS data
#' 
#' Extract the FCS expresssion data and apply the transformation 
#' 
#' @param fcsFile The name of the FCS file
#' @param comp Boolean value tells if do compensation
#' @param verbose Boolean value detecides if print the massage details
#' @param markers Selected markers for analysis, either from names or from description
#' @param transformMethod transformation method, \code{auto_lgcl}, \code{fixed_lgcl}, \code{arcsin} or \code{biexp}
#' @param scaleTo scale the expression to same scale, default is NULL, should be a vector of two numbers if scale
#' @param w Linearization width in asymptotic decades
#' @param t Top of the scale data value
#' @param m Full width of the transformed display in asymptotic decades
#' @param a Additional negative range to be included in the display in asymptotic decades
#' @param q quantile of negative values removed for auto w estimation, default is 0.05
#' 
#' @return The transformend expression data matrix with selected markers
#' @importFrom flowCore read.FCS compensate estimateLogicle logicleTransform parameters transformList arcsinhTransform biexponentialTransform
#' @importMethodsFrom flowCore transform
#' @importClassesFrom flowCore transformList
#' @export
#' @examples
#' d<-system.file('extdata',package='cytofkit')
#' fcsFile <- list.files(d,pattern='.fcs$',full=TRUE)
#' transformed <- cytof_exprsExtract(fcsFile)
cytof_exprsExtract <- function(fcsFile, comp = FALSE, verbose = FALSE, 
    markers = NULL, transformMethod = "auto_lgcl", scaleTo = NULL, 
    w = 0.1, t = 4000, m = 4.5, a = 0, q = 0.05) {
    
    ## load FCS files
    name <- sub(".fcs", "", basename(fcsFile))
    if (verbose) {
        fcs <- read.FCS(fcsFile)
    } else {
        fcs <- suppressWarnings(read.FCS(fcsFile))
    }
    
    ## compensation if provided in the data
    if (comp == TRUE) {
        if (comp && !is.null(fcs@description$SPILL)) {
            fcs <- applyComp(fcs, "SPILL")
        } else if (comp && !is.null(fcs@description$SPILLOVER)) {
            fcs <- applyComp(fcs, "SPILLOVER")
        } else if (comp && !is.null(fcs@description$COMP)) {
            fcs <- applyComp(fcs, "COMP")
        }
    }
    
    ## marker check, allow mix marker names
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
    } else {
        marker_id <- 1:ncol(fcs@exprs)
    }
    
    ## exprs transformation
    if (transformMethod == "auto_lgcl") {
        trans <- auto_lgcl(fcs, channels = colnames(fcs@exprs)[marker_id])
        transformed <- flowCore::transform(fcs, trans)
        exprs <- transformed@exprs[, marker_id]
    } else if (transformMethod == "fixed_lgcl") {
        trans <- flowCore::logicleTransform(w = w, t = t, m = m, a = a)
        exprs <- apply(fcs@exprs[, marker_id], 2, trans)
    } else if (transformMethod == "arcsin") {
        trans <- flowCore::arcsinhTransform(a=1, b=1, c=1)
        exprs <- apply(fcs@exprs[, marker_id], 2, trans)
    } else if(transformMethod == "biexp") {
        trans <- flowCore::biexponentialTransform(a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0)
        exprs <- apply(fcs@exprs[, marker_id], 2, trans)
    }else{
        stop("transformMethod [", transformMethod, "] doesn't exist for cytofkit!")
    }
    
    ## rescale data
    if (!is.null(scaleTo)) {
        exprs <- apply(exprs, 2, function(x) scaleData(x, scaleTo))
    }
     
    ## add rownames and colnames   
    col_names <- paste0(pd$name, "<", pd$desc,">")
    colnames(exprs) <- col_names[marker_id]
    row.names(exprs) <- paste(name, 1:nrow(exprs), sep = "_")
    
    return(exprs)
}


#' apply compensation on the FCS expression data
#' 
#' @param fcs FCS file.
#' @param keyword Keywords.
applyComp <- function(fcs, keyword) {
    comp_fcs <- compensate(fcs, fcs@description[[keyword]])
}

#' rescale the data
#' 
#' @param x data.
#' @param range The range of the data.
scaleData <- function(x, range = c(0, 4.5)) {
    (x - min(x))/(max(x) - min(x)) * (range[2] - range[1]) + range[1]
}

#' a modified version of "estimateLogicle" from flowCore
#' 
#' 
#' @param x Data.
#' @param channels Channel names.
#' @param m Para m.
#' @param q Para q.
auto_lgcl <- function(x, channels, m = 4.5, q = 0.05) {
    if (!is(x, "flowFrame")) 
        stop("x has to be an object of class \"flowFrame\"")
    if (missing(channels)) 
        stop("Please specify the channels to be logicle transformed")
    indx <- channels %in% colnames(x@exprs)
    if (!all(indx)) 
        stop(paste("Channels", channels[!indx], "were not found in the FCS file.\n ", 
            sep = " "))
    rng <- range(x)
    trans <- lapply(channels, function(p) {
        lgclParaEstimate(x@exprs[, p], m = m, q = q)
    })
    transformList(channels, trans)
}

lgclParaEstimate <- function(data, m = 4.5, q = 0.05, type = "instrument") {
    t <- max(data)
    ndata <- data[data < 0]
    w <- 0
    a <- 0
    
    if (missing(m)) {
        if (type == "instrument") 
            m <- 4.5
        else m <- log10(t) + 1
    }
    
    if (length(ndata)) {
        r <- .Machine$double.eps + quantile(ndata, q)
        if (10^m * abs(r) <= t) {
            w <- 0  ## this special check to avoid failure of negative w
        } else {
            w <- (m - log10(t/abs(r)))/2
            if(w>2) {
                w <- 0.5
                cat("w is coerced to 0.5\n")
            }
        }
    }
    # cat(paste('sign_autoLgcl parameters:', ' w=', w, ' t=',t,'
    # m=',m,' a=',a,'\n', sep = ''))
    logicleTransform(w = w, t = t, m = m, a = a)
} 
