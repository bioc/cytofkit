
library(ggplot2)
library(gplots)
library(reshape2)
library(reshape)
library(plyr)
library(VGAM)


visuaPlot <- function(obj, xlab, ylab, zlab, pointSize=1, 
                      addLabel=TRUE, labelSize=1, selectSamples,
                      removeOutlier = TRUE){
    
    data <- cbind(obj$allExpressionData, do.call(cbind, obj$dimReducedRes))
    data <- as.data.frame(data)
    clusterMethods <- names(obj$clusterRes)
    for(cname in clusterMethods){
        data[[cname]] <- as.factor(obj$clusterRes[[cname]])
    }
    row.names(data) <- row.names(obj$expressionData)
    samples <- sub("_[0-9]*$", "", row.names(obj$expressionData))
    
    data <- data[samples %in% selectSamples, ]
    nsamples <- samples[samples %in% selectSamples]
    data$sample <- nsamples
    sample_num <- length(unique(nsamples))
    if (sample_num >= 8) {
        shape_value <- LETTERS[1:sample_num]
    } else {
        shape_value <- c(1:sample_num) + 15
    }
    
    if(zlab %in% clusterMethods){
        cluster_num <- length(unique(data[[zlab]]))
        col_legend_row <- ceiling(cluster_num/15)
        size_legend_row <- ceiling(sample_num/4)
        shapeLab <- "sample"
        
        gp <- ggplot(data, aes_string(x=xlab, y=ylab, colour = zlab, shape = shapeLab)) + 
            geom_point(size = pointSize) + scale_shape_manual(values = shape_value) + 
            scale_colour_manual(values = rainbow(cluster_num)) + 
            xlab(xlab) + ylab(ylab) + #coord_fixed() + 
            guides(colour = guide_legend(nrow = col_legend_row, override.aes = list(size = 4)), 
                   shape = guide_legend(nrow = size_legend_row, override.aes = list(size = 4))) + 
            theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        
        if(addLabel){
            edata <- data[ ,c(xlab, ylab, zlab)]
            colnames(edata) <- c('x', "y", "z")
            center <- aggregate(cbind(x,y) ~ z, data = edata, median)
            gp <- gp + annotate("text", label = center[,1], x=center[,2], y = center[,3], 
                                size = labelSize, colour = "black")
        }
        
        gp <- gp + theme(legend.position = "bottom", axis.text=element_text(size=14),
                         axis.title=element_text(size=18,face="bold"))
    
    }else{
        ## ggplot aes cannot recognize marker names with symbol <>
        title <- zlab
        data <- data[,c(xlab, ylab, zlab)]
        if(removeOutlier)
            data[,zlab] <- remove_outliers(data[,zlab])
        zlab <- "Expression"
        colnames(data) <- c(xlab, ylab, zlab)
        gp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = zlab)) + 
            geom_point(size = pointSize) + theme_bw() + #coord_fixed() + 
            scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=median(data[[zlab]])) +
            theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
    }
    
    return(gp)
}


heatMap <- function(data, clusterMethod = "DensVM", type = "mean", selectSamples,
                    cex_row_label = 1, cex_col_label = 1, scaleMethod = "none") {
    exprs <- data$expressionData
    samples <- sub("_[0-9]*$", "", row.names(exprs))
    exprs <- exprs[samples %in% selectSamples, ]
    ifMultiFCS <- length(selectSamples) > 1
    dataj <- data$clusterRes[[clusterMethod]][samples %in% selectSamples]
    exprs_cluster <- data.frame(exprs, cluster = dataj)
    
    if(type == "mean"){
        cluster_stat <- aggregate(. ~ cluster, data = exprs_cluster, mean)
        rownames(cluster_stat) <- paste("cluster_", cluster_stat$cluster, sep = "")
        cluster_stat <- cluster_stat[ ,-which(colnames(cluster_stat) == "cluster")]
    }else if (type == "median"){
        cluster_stat <- aggregate(. ~ cluster, data = exprs_cluster, median)
        rownames(cluster_stat) <- paste("cluster_", cluster_stat$cluster, sep = "")
        cluster_stat <- cluster_stat[ ,-which(colnames(cluster_stat) == "cluster")]
    }else if(type == "percentage" && ifMultiFCS){
        sampleName <- sub("_[0-9]*$", "", row.names(exprs))
        clusterCounts <- as.data.frame(table(sampleName, dataj))
        colnames(clusterCounts) <- c("sample", "cluster", "cellCount")
        sampleCellCount <- as.data.frame(table(sampleName))
        colnames(sampleCellCount) <- c("sample", "totalCellCount")
        clust_cellCount <- merge(clusterCounts, sampleCellCount, by = "sample")
        clust_cellCount$percentage <- round(clust_cellCount$cellCount/clust_cellCount$totalCellCount * 100, 2)
        cluster_stat <- reshape::cast(clust_cellCount, sample ~ cluster, value = "percentage")
        percColNames <- cluster_stat$sample
        cluster_stat <- cluster_stat[, -which(colnames(cluster_stat) == "sample")]
        percRowNames <- paste("cluster_", colnames(cluster_stat), sep = "")
        cluster_stat <- t(as.matrix(cluster_stat))
        row.names(cluster_stat) <- percRowNames
        colnames(cluster_stat) <- percColNames
        
    }else{
        return(NULL)
    }

    cluster_stat <- as.matrix(cluster_stat)
    heatmap.2(cluster_stat, col = bluered, trace = "none", 
              symbreaks = FALSE, scale = scaleMethod, 
              margins = c(8, 8),
              cexRow = cex_row_label, 
              cexCol = cex_col_label, 
              srtCol = 30, symkey = FALSE, 
              keysize = 1, key.par=list(mgp=c(2, 1, 0),mar=c(4, 3, 4, 0)),
              main = paste(clusterMethod, type, "heatmap", sep = " "))
}



progressionPlot <- function(data, orderCol="isomap_1", clusterCol = "cluster", 
                            trend_formula="expression ~ sm.ns(Pseudotime, df=3)"){
    
    progressionData <- data$progressionRes
    if(!is.null(progressionData)){
        data <- do.call(cbind, progressionData) 
        markers <- colnames(progressionData[[1]]) 
        colnames(data) <- c(markers, "cluster", colnames(progressionData[[3]]))
       
        if(!is.data.frame(data)) data <- data.frame(data, check.names = FALSE)
        if(!all(markers %in% colnames(data))) stop("Unmatching markers found!")
        if(!(length(orderCol)==1 && orderCol %in% colnames(data)))
            stop("Can not find orderCol in data")
        if(!(length(clusterCol)==1 && clusterCol %in% colnames(data)))
            stop("Can not find clusterCol in data")
        
        orderValue <- data[[orderCol]]
        data <- data[order(orderValue), c(markers, clusterCol)]
        data$Pseudotime <- sort(orderValue)
        
        mdata <- melt(data, id.vars = c("Pseudotime", clusterCol))
                      # variable.name = "markers", variable.name= "expression")
        colnames(mdata) <- c("Pseudotime", clusterCol, "markers", "expression")
        mdata$markers <- factor(mdata$markers)
        mdata[[clusterCol]] <- factor(mdata[[clusterCol]])
        min_expr <- min(mdata$expression)
        
        ## tobit regression
        vgamPredict <- ddply(mdata, .(markers), function(x) { 
            fit_res <- tryCatch({
                vg <- suppressWarnings(vgam(formula = as.formula(trend_formula), 
                                            family = VGAM::tobit(Lower = min_expr, lmu = "identitylink"), 
                                            data = x, maxit=30, checkwz=FALSE))
                res <- VGAM::predict(vg, type="response")
                res[res < min_expr] <- min_expr
                res
            }
            ,error = function(e) {
                print("Error!")
                print(e)
                res <- rep(NA, nrow(x))
                res
            }
            )
            expectation = fit_res
            data.frame(Pseudotime=x$Pseudotime, expectation=expectation)
        })
        
        color_by <- clusterCol
        plot_cols <- round(sqrt(length(markers)))
        cell_size <- 1
        x_lab <- orderCol
        y_lab <- "Expression"
        legend_title <- clusterCol
        
        ## copied from monocle package
        monocle_theme_opts <- function(){
            theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
                theme(panel.border = element_blank(), axis.line = element_line()) +
                theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
                theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
                theme(panel.background = element_rect(fill='white')) +
                theme(legend.position = "right") +
                theme(axis.title = element_text(size = 15)) }
        
        q <- ggplot(aes(Pseudotime, expression), data=mdata) 
        q <- q + geom_point(aes_string(color=color_by), size=I(cell_size))
        q <- q + geom_line(aes(Pseudotime, expectation), data=vgamPredict)
        q <- q + facet_wrap(~markers, ncol=plot_cols, scales="free_y")
        q <- q + ylab(y_lab) + xlab(x_lab) + theme_bw()
        q <- q + guides(colour = guide_legend(title = legend_title, override.aes = list(size = cell_size*3)))
        q <- q + monocle_theme_opts() 
        
        return(q)
    }else{
        return(NULL)
    }
}


remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- qnt[1] - H
    y[x > (qnt[2] + H)] <- qnt[2] + H
    y
}
