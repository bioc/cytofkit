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
    if(!is.null(dimReduction))
        data <- cytof_dimReduction(data, method=dimReduction, out_dim = outDim)
    if(missing(dc))
        dc <- estimateDc(data, sampleSize = 10000)
    if(parallel){
        cl <- makeCluster(nCore)  
        registerDoParallel(cl)
    }
    rho <- localDensity(data, dc, gaussian=gaussian, ifParallel = parallel)
    deltaWid <- minDistToHigher(data, rho, ifParallel = parallel)
    delta <- deltaWid[[1]]
    higherID <- deltaWid[[2]]
    peakID <- peakDetect(rho, delta, alpha)
    cluster <- clusterAssign(peakID, higherID, rho)
    clusTable <- as.vector(table(cluster))
    if(sum(clusTable < length(rho)*0.0005) > 0){
        cat("Noise cluster removing, ")
        peakID <- peakID[clusTable >= length(rho)*alpha]
        cluster <- clusterAssign(peakID, higherID, rho)
    }
    
    if(detectHalos){
        halo <- haloDetect(data, rho, cluster, peakID, dc)
    }else{
        halo <- NULL
    }
    
    if(parallel){
        stopCluster(cl)
    }
    
    if(ncol(data) < 3){
        plotData <- data
    }else{
        plotData <- data[ ,c(1,2)]
    }
        
    res <- list(cluster = cluster, dc = dc, rho = rho, delta = delta, peakID = peakID, 
                higherID = higherID, halo = halo, plotData = plotData)
    
    #class(res) <- 'ClusterX'
    return(res)
}


#' Estimate the distance cutoff (density neighbourhood) from down-sampled data
#' 
#' This function estimate a distance cutoff value from the down-samples data,
#' wchich meet the criteria that the average neighbor rate (number of points 
#' within the distance cutoff value) fall between the provided range. 
#' 
#' @param data Numeric matrix of data or data frame.
#' @param sampleSize The size of the down-sampled data.
#' @param neighborRateLow The lower bound of the neighbor rate (default 0.01).
#' @param neighborRateHigh The upper bound of the neighbor rate (default 0.15).
#' 
#' @return A numeric value giving the estimated distance cutoff value.
estimateDc <- function(data, sampleSize = 10000, neighborRateLow=0.01, neighborRateHigh=0.02) {
    data <- as.matrix(data)
    dataLens <- nrow(data)
    if(dataLens > sampleSize){
        sample <- data[sample(1:dataLens, sampleSize, replace = FALSE), ]
    }else{ sample <- data }
    
    comb <- as.matrix(dist(sample, method = "euclidean"))
    size <- nrow(comb)
    dc <- min(comb)
    dcMod <- median(comb)*0.05
    
    while(TRUE) {
        neighborRate <- mean((apply(comb < dc, 1, sum)-1)/size)
        #neighborRate <- mean(apply((exp(-(comb/dc)^2)), 1, sum) - 1) / size
        if(neighborRate > neighborRateLow && neighborRate < neighborRateHigh) break  
        if(neighborRate >= neighborRateHigh) {
            dc <- dc - dcMod
            dcMod <- dcMod/2
        }else{
            dc <- dc + dcMod
        }
    }
    cat('Distance cutoff calculated to', dc, '\n')
    dc
}


#' Computes the local density of points in a data matrix
#' 
#' This function calculate the local density for each point in the matrix. 
#' With a rowise implementation of the pairwise distance calculation, makes 
#' the local density estimation faster and memory efficient. A big benifit 
#' is the aviliability for big data. Parallel computing is supported for 
#' fast calculation. The computation can either be done using a simple summation 
#' of the points with the distance cutoff for each observation, or by applying 
#' a gaussian kernel scaled by the distance cutoff (more robust for low-density data)
#' 
#' @param data Numeric matrix of data or data frame.
#' @param dc A numeric value specifying the distance cutoff.
#' @param gaussian Logical value decide if a gaussian kernel be used to estimate the density (defaults to TRUE).
#' @param ifParallel A boolean decides if run parallelly
#' 
#' @return A vector of local density values with index matching the row names of data.
localDensity <- function(data, dc, gaussian=FALSE, ifParallel = FALSE) {
    cat("calculate local density...")
    splitFactor <- splitFactorGenerator(nrow(data))
    dataFolds <- split.data.frame(data, splitFactor)

    rholist <- llply(dataFolds, function(datai, data, dc, gaussian) {
        
        suppressWarnings(idist <- as.matrix(pdist(datai, data)))
        
        if(gaussian){
            apply((exp(-(idist/dc)^2)), 1, sum) - 1
        }else{
            apply(idist < dc, 1, sum) - 1 } 
        }, data = data, dc = dc, gaussian = gaussian, .parallel = ifParallel)
    
    rho <- do.call(base::c, rholist)
    if(is.null(row.names(data))) {
        names(rho) <- NULL
    } else {
        names(rho) <- row.names(data)
    }
    cat("DONE!\n")
    rho
}


#' Calculate distance to closest observation of higher density
#' 
#' This function finds, for each observation, the minimum distance to an 
#' observation of higher local density. With a rowise implementation of 
#' the pairwise distance calculation, makes the local density estimation 
#' faster and memory efficient. A big benifit is the aviliability for big 
#' data. Parallel computing is supported for fast calculation.
#' 
#' @param data Numeric matrix of data or data frame.
#' @param rho A vector of local density values as outputted by \code{localDensity}.
#' @param ifParallel A boolean decides if run parallelly.
#' 
#' @return A list of distances to closest observation of higher density and the ID
minDistToHigher <- function(data, rho, ifParallel = FALSE) {
    cat("Search nearest neighobour with higher density...")
    splitFactor <- splitFactorGenerator(nrow(data))
    dataFolds <- split.data.frame(data, splitFactor)
    rhoFolds <- split(rho, splitFactor)
    
    dataRhoList <- mapply(function(datai, rhoi) {
        list(datai, rhoi)}, dataFolds, rhoFolds, SIMPLIFY = FALSE)
    
    deltaWidList <- llply(dataRhoList, function(x, data, rho){
        datai <- x[[1]]
        rhoi <- x[[2]]
        suppressWarnings(datai2dataDist <- as.matrix(pdist(datai, data)))
        
        rhoi2rhoComp <- do.call(rbind, lapply(rhoi, function(x) x < rho)) 
        sapply(seq_len(nrow(datai2dataDist)), function(i){
            distToAllHigherPoints <- datai2dataDist[i, rhoi2rhoComp[i, ]]
            if(length(distToAllHigherPoints) == 0) {
                c(max(datai2dataDist[i, ]), which.max(datai2dataDist[i, ]))
            } else {
                c(min(distToAllHigherPoints), 
                  which(rhoi2rhoComp[i, ] == TRUE)[which.min(distToAllHigherPoints)] )
            }} )
    }, data = data, rho = rho, .parallel = ifParallel)
    
    deltaWid <- do.call(cbind, deltaWidList)
    delta <- deltaWid[1, ]
    names(delta) <- names(rho)
    id <- deltaWid[2, ]
    names(id) <- names(rho)
    cat("DONE!\n")
    return(list(delta = delta, higherID = id))
}



#' Automatic peak detection
#' 
#' Automatic detect peaks by searching high denisty point with anomalous large distance to
#' higher denisty peaks. rho and delta are transformed to one index, and the anomalous peaks
#' are detected using generalized ESD method.
#' 
#' @param rho A vector of the local density, outout of \code{localDensity}
#' @param delta A vector of distance to closest observation of higher density
#' @param alpha The level of statistical significance for peak detection.
#' 
#' @return a vector containing the indexes of peaks
peakDetect <- function(rho, delta, alpha = 0.001){
    cat("Peak detection...")
    delta[is.infinite(delta)] <- max(delta[!(is.infinite(delta))])^2
    rdIndex <- scale01(rho) * delta   ## transform delta, important for big data
    #rdIndex <- log(rho + 1) * delta
    peakID1 <- detect_anoms_sd(rdIndex, direction = "pos", alpha = alpha)
    peakID2 <- detect_anoms_sd(delta, direction = "pos", alpha = alpha)
    peakID <- intersect(peakID1, peakID2)
    cat("DONE!\n")
    peakID
}

scale01 <- function(x){
    x <- (x-min(x))/(max(x)-min(x)) 
    x
}

truncTrans <- function(x, a, b=1){
    if(a < b)
        x[x<=a] <- b
    x
}

#' assign clusters to non-peak points
#' 
#' @param peakID A vector of the peak ID.
#' @param higherID A vector of IDs of the nearest neighbor with higer density. 
#' @param rho A vector of the density values.
#' 
#' @return the cluster ID
clusterAssign <- function(peakID, higherID, rho){
    cat("Cluster assigning...")
    runOrder <- order(rho, decreasing = TRUE)
    cluster <- rep(NA, length(rho))
    for(i in runOrder) {
        cluster[i] <- ifelse(i %in% peakID, match(i, peakID), cluster[higherID[i]] ) }
    cat("DONE!\n")
    cluster
}

#' differentiate halo form cores
#' 
#' @param data Input data.
#' @param rho Density values.
#' @param cluster The cluster IDs.
#' @param peakID The peak IDs.
#' @param dc The distance cutoff value.
#' 
#' @return the boolean value indicating if the point belongs to halo
haloDetect <- function(data, rho, cluster, peakID, dc){
    cat("diffenentiate halos from cores ...")
    clusterRhoThres <- sapply(1:length(peakID), function(clusteri){
        dataOfClusteri <- data[cluster == clusteri, ]
        otherData <- data[cluster != clusteri, ]
        rhoOfClusteri <- rho[cluster == clusteri]
        splitFactor <- splitFactorGenerator(nrow(dataOfClusteri), nrow(otherData))
        dataFolds <- split.data.frame(dataOfClusteri, splitFactor)
        haveColseNeighbourList <- lapply(dataFolds, function(x){
            suppressWarnings(dataFoldi2otherDataDist <- as.matrix(pdist(x, otherData)))
            apply(dataFoldi2otherDataDist < dc/2, 1, any)
        })
        checkRes <- do.call(base::c, haveColseNeighbourList)
        rhoThres <- ifelse(any(checkRes), max(rhoOfClusteri[checkRes]), min(rhoOfClusteri))
        rhoThres
    })
    halo <- rho < clusterRhoThres[cluster]
    cat("DONE!\n")
    halo
}



## generate split factors to split rowNum into folds, each of size around foldSize
splitFactorGenerator <- function(rowNum, colNum){
    if(missing(colNum)){
        colNum <- rowNum
    }
    foldSize <- round(32772663 / colNum)  ## each chunk with maxi 250Mb
    foldNum <- ceiling(rowNum / foldSize)
    lastfoldSize <- rowNum - (foldNum-1) * foldSize
    if(foldNum > 1){
        splitFactor <- c(rep(1:(foldNum-1), each = foldSize), rep(foldNum, lastfoldSize))
    }else{
        splitFactor <- rep(foldNum, lastfoldSize)
    }
    return(splitFactor)
}




#' Outlier detection
#' 
#' Using generialized ESD to detect outliers, iterate and remove point with ares higher than lamda 
#' in a univariate data set assumed to come from a normally distributed population.
#' 
#' @param data A vectors of boservations.
#' @param max_anoms Maximal percentile of anomalies.
#' @param alpha The level of statistical significance with which to accept or reject anomalies.
#' @param direction Directionality of the anomalies to be detected. Options are: 'pos' | 'neg' | 'both'.
#' 
#' @return A vector containing indexes of the anomalies (outliers).
#' 
detect_anoms_sd <- function(data, max_anoms=0.1, alpha = 0.01, direction='pos') {
    
    num_obs <- length(data)
    names(data) <- 1:num_obs
    data <- na.omit(data)
    max_outliers <- trunc(num_obs*max_anoms)
    
    anoms <- NULL
    for (i in 1L:max_outliers){
        ares <- switch(direction, 
                       pos = data - mean(data),
                       neg = mean(data) - data,
                       both = abs(data - mean(data)) )
        
        p <- switch(direction, 
                    pos = 1 - alpha/(num_obs-i+1),
                    neg = 1 - alpha/(num_obs-i+1),
                    both = 1 - alpha/(2*(num_obs-i+1)) )
        
        data_sigma <- sd(data)
        if(data_sigma == 0) break
        ares <- ares/data_sigma
        maxAres <- max(ares)
        
        t <- qt(p, (num_obs-i-1))
        lam <- t*(num_obs-i) / sqrt((num_obs-i-1+t**2)*(num_obs-i+1))
        
        if(maxAres > lam){
            maxAres_id <- which(ares == maxAres)
            anoms <- c(anoms, names(maxAres_id))
            data <- data[-maxAres_id ]
        }else
            break
    }
    
    return(as.numeric(anoms))
}

