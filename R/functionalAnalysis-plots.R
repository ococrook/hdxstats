## plotting for epitope mapping

manhattentplot <- function(params,
                           sequences,
                           difference = 0,
                           region = NULL, 
                           nrow = 1){
    
    butterfly <- matrix(0, ncol = length(unique(sequences)), nrow = 1)
    colnames(butterfly) <- unique(sequences)
    butterfly <- params$ebayesres$fdr[unique(sequences)]
    butterflydf <- as.data.frame(butterfly)
    butterflydf$Sequence <- rownames(butterflydf)
    butterflydf$protection <- 1*(difference > 0)
    butterflydf$protection[is.na(butterflydf$protection)] <- 0
    colnames(butterflydf)[1] <- "p_value"
    butterfly_long <- butterflydf
    butterfly_long$position <- rep(seq.int(nrow(butterfly_long)/1), each = 1)
    butterfly_long$region <- unique(region[, c("Start", "End")])
    
    plot.list <- list()
    r <- nrow
    for (i in 1:r) {
        
        xannot <- paste0("[", butterfly_long$region[c((i - 1)*nrow(butterfly_long)/r + 1):c(i*nrow(butterfly_long)/r ),1], ",",
                         butterfly_long$region[c((i - 1)*nrow(butterfly_long)/r + 1):c(i*nrow(butterfly_long)/r ),2],"]")
        
        plot.list[[i]] <- ggplot(butterfly_long, aes(x = position,
                                                     y = -log10(p_value),
                                                     color = factor(protection),
                                                     group = -log10(p_value))) + 
            geom_point(size = 3) + 
            scale_color_manual(values = brewer.pal(n = 3, name = "Set2")[1:2], labels = c("protected","deprotected")) + 
            theme_classic() + ylim(c(0, max(-log10(butterfly_long$p_value)))) + 
            ylab("-log10 pvalue") + xlab("Peptide") + xlim(c((i - 1)*nrow(butterfly_long)/r + 1, i*nrow(butterfly_long)/r )) + 
            scale_x_continuous(breaks = c((i - 1)*nrow(butterfly_long)/r + 1):c(i*nrow(butterfly_long)/r ), labels = xannot) + 
            ggtitle("Manhatten plot") + geom_hline(yintercept = 1.301, linetype = "dashed", colour = "red") + 
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 15)) + 
            labs(color = "Protection")
    }
    
    return(plot.list = plot.list)
    
}

plotEpitopeMap <- function(AAString, 
                         peptideSeqs,
                         numlines = 5,
                         maxmismatch = 0,
                         by = 5,
                         scores = NULL,
                         name = "-log10 p values", threshold = -log10(0.05)){
    
    # Test
    stopifnot("AAString must be an object of class AAString"= class(AAString) == "AAString")
    stopifnot("peptideSeqs must be a character vector"= is.character(peptideSeqs) == "TRUE")
    
    # Storage and global variables
    plot.list <- list()
    coverage <- matrix(0, ncol = length(AAString), nrow = 1)
    n <- ceiling(length(coverage)/numlines)
    colnames(coverage) <- strsplit(as.character(AAString), "")[[1]]
    
    # Compute AA stringset and match to dictionary
    peptideset <- AAStringSet(x = peptideSeqs)
    allPatterns <- matchPDict(pdict = peptideset,
                              subject = AAString,
                              max.mismatch = maxmismatch)  
    
    # Compute coverage numbers
    for (i in seq_along(allPatterns)) {
        
        begin <- allPatterns[[i]]@start[1]
        end <- allPatterns[[i]]@start[1] - 1 + allPatterns[[i]]@width[1]
        coverage[, seq.int(begin, end)] <- coverage[, seq.int(begin, end)] + 1
    }
    
    ncov <- max(coverage)
    
    start <- sapply(allPatterns, function(x) x@start) - 1
    end <- start + sapply(allPatterns, function(x) x@width)
    
    
    peptideMap <- matrix(0, ncol = length(AAString), nrow = ncov + 3)
    colnames(peptideMap) <- strsplit(as.character(AAString), "")[[1]]
    rownames(peptideMap) <- seq.int(1:nrow(peptideMap))
    
    if (is.null(scores) == TRUE){
        for (i in seq_along(start)){
            
            if (i == 1) {
                peptideMap[1, start[i]:end[i]] <- 1
            } else {
                j <- which.min(rowSums(peptideMap[, start[i]:end[i]] > 0))
                peptideMap[j, start[i]:end[i]] <- 1
            }
        }
    } else {
        for (i in seq_along(start)){
            
            if (i == 1) {
                peptideMap[1, start[i]:end[i]] <- scores[i]
            } else {
                j <- which.min(rowSums(peptideMap[, start[i]:end[i]] > 0))
                peptideMap[j, start[i]:end[i]] <- scores[i]
            }
        }
    }
    
    
    
    plot.list <- list()
    if (is.null(scores) == TRUE){
        for (i in 1:ceiling(length(coverage)/n)) {
            
            if (i < ceiling(length(coverage)/n)){    
                mygrid <- expand.grid(X = factor(1:n), Y = rownames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]))
                mygrid$Z <- as.factor(c(t(peptideMap[, n * (i - 1) + 1:n])))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + 
                    scale_fill_manual(breaks = c(0, 1), values = c("white", alpha("#1B7837", 0.7))) +
                    theme_classic() + scale_x_discrete(breaks = 1:n, labels = colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE])) +
                    theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            } else {
                mygrid <- expand.grid(X = factor(1:(ncol(peptideMap)%%n + 1)), Y = rownames(peptideMap[, (n * (i-1)):ncol(peptideMap), drop = FALSE]))
                mygrid$Z <- as.factor(c(t(peptideMap[, (n * (i-1)):ncol(peptideMap)])))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + 
                    scale_fill_manual(breaks = c(0,1), values = c("white", alpha("#1B7837", 0.7)))  +
                    theme_classic() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1), labels = colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)])) +
                    theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            }
        }
    } else{
        
        sc <- scale_fill_manual(name, 
                                values = c("white", brewer.pal(n = 2, name = "Set2")),
                                labels = c(" ", "Not Signifcant", "Significant"),
                                breaks = levels(factor(peptideMap))[seq(1, nlevels(factor(peptideMap)), by = by)], drop = FALSE)
        
        for (i in 1:ceiling(length(coverage)/n)) {
            
            
            if (i < ceiling(length(coverage)/n)){    
                mygrid <- expand.grid(X = factor(1:n), Y = rownames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]))
                mygrid$Z <- factor(c(t(peptideMap[, n * (i - 1) + 1:n])), levels = levels(factor(peptideMap)))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc +
                    theme_classic() + scale_x_discrete(breaks = 1:n, labels = colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE])) +
                    theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            } else {
                mygrid$Z <- factor(c(t(peptideMap[, (n * (i-1)):ncol(peptideMap)])), levels = levels(factor(peptideMap)))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc + 
                    theme_classic() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1), labels = colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)])) +
                    theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            }
        }
    }
    return(plot.list)
    
}

plotEpitopeMapFdr <- function(AAString, 
                              peptideSeqs,
                              numlines = 5,
                              maxmismatch = 0,
                              by = 5,
                              scores = NULL,
                              name = "-log10 p values",
                              threshold = -log10(0.05)){
    
    # Test
    stopifnot("AAString must be an object of class AAString"= class(AAString) == "AAString")
    stopifnot("peptideSeqs must be a character vector"= is.character(peptideSeqs) == "TRUE")
    
    # Storage and global variables
    plot.list <- list()
    coverage <- matrix(0, ncol = length(AAString), nrow = 1)
    n <- ceiling(length(coverage)/numlines)
    colnames(coverage) <- strsplit(as.character(AAString), "")[[1]]
    
    # Compute AA stringset and match to dictionary
    peptideset <- AAStringSet(x = peptideSeqs)
    allPatterns <- matchPDict(pdict = peptideset,
                              subject = AAString,
                              max.mismatch = maxmismatch)  
    
    # Compute coverage numbers
    for (i in seq_along(allPatterns)) {
        
        begin <- allPatterns[[i]]@start[1]
        end <- allPatterns[[i]]@start[1] - 1 + allPatterns[[i]]@width[1]
        coverage[, seq.int(begin, end)] <- coverage[, seq.int(begin, end)] + 1
    }
    
    ncov <- max(coverage)
    
    start <- sapply(allPatterns, function(x) x@start) - 1
    end <- start + sapply(allPatterns, function(x) x@width)
    
    
    peptideMap <- matrix(0, ncol = length(AAString), nrow = ncov + 3)
    colnames(peptideMap) <- strsplit(as.character(AAString), "")[[1]]
    rownames(peptideMap) <- seq.int(1:nrow(peptideMap))
    
    if (is.null(scores) == TRUE){
        for (i in seq_along(start)){
            
            if (i == 1) {
                peptideMap[1, start[i]:end[i]] <- 1
            } else {
                j <- which.min(rowSums(peptideMap[, start[i]:end[i]] > 0))
                peptideMap[j, start[i]:end[i]] <- 1
            }
        }
    } else {
        for (i in seq_along(start)){
            
            if (i == 1) {
                peptideMap[1, start[i]:end[i]] <- scores[i]
            } else {
                j <- which.min(rowSums(peptideMap[, start[i]:end[i]] > 0))
                peptideMap[j, start[i]:end[i]] <- scores[i]
            }
        }
    }
    
    
    
    plot.list <- list()
    if (is.null(scores) == TRUE){
        for (i in 1:ceiling(length(coverage)/n)) {
            
            if (i < ceiling(length(coverage)/n)){    
                mygrid <- expand.grid(X = factor(1:n), Y = rownames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]))
                mygrid$Z <- as.factor(c(t(peptideMap[, n * (i - 1) + 1:n])))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + 
                    scale_fill_manual(breaks = c(0, 1), values = c("white", alpha("#1B7837", 0.7))) +
                    theme_classic() + scale_x_discrete(breaks = 1:n, labels = colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE])) +
                    theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            } else {
                mygrid <- expand.grid(X = factor(1:(ncol(peptideMap)%%n + 1)), Y = rownames(peptideMap[, (n * (i-1)):ncol(peptideMap), drop = FALSE]))
                mygrid$Z <- as.factor(c(t(peptideMap[, (n * (i-1)):ncol(peptideMap)])))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + 
                    scale_fill_manual(breaks = c(0,1), values = c("white", alpha("#1B7837", 0.7)))  +
                    theme_classic() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1), labels = colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)])) +
                    theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            }
        }
    } else{
        
        sc <- scale_fill_manual(name, 
                                values = c("white", viridis(n = nlevels(factor(peptideMap))), alpha = 0.7),
                                labels = substring(levels(factor(peptideMap))[seq(1, nlevels(factor(peptideMap)), by = by)], 1, 4),
                                breaks = levels(factor(peptideMap))[seq(1, nlevels(factor(peptideMap)), by = by)], drop = FALSE)
        
        for (i in 1:ceiling(length(coverage)/n)) {
            
            
            if (i < ceiling(length(coverage)/n)){    
                mygrid <- expand.grid(X = factor(1:n), Y = rownames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]))
                mygrid$Z <- factor(c(t(peptideMap[, n * (i - 1) + 1:n])), levels = levels(factor(peptideMap)))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc +
                    theme_classic() + scale_x_discrete(breaks = 1:n, labels = colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE])) +
                    theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            } else {
                mygrid$Z <- factor(c(t(peptideMap[, (n * (i-1)):ncol(peptideMap)])), levels = levels(factor(peptideMap)))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc + 
                    theme_classic() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1), labels = colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)])) +
                    theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            }
        }
    }
    return(plot.list)
    
}

plotEpitopeMapResidue <- function(AAString, 
                              peptideSeqs,
                              numlines = 5,
                              maxmismatch = 0,
                              by = 5,
                              scores = NULL,
                              name = "-log10 p values",
                              threshold = -log10(0.05)){
    
    # Test
    stopifnot("AAString must be an object of class AAString"= class(AAString) == "AAString")
    stopifnot("peptideSeqs must be a character vector"= is.character(peptideSeqs) == "TRUE")
    
    # Storage and global variables
    plot.list <- list()
    coverage <- matrix(0, ncol = length(AAString), nrow = 1)
    n <- ceiling(length(coverage)/numlines)
    colnames(coverage) <- strsplit(as.character(AAString), "")[[1]]
    
    # Compute AA stringset and match to dictionary
    peptideset <- AAStringSet(x = peptideSeqs)
    allPatterns <- matchPDict(pdict = peptideset,
                              subject = AAString,
                              max.mismatch = maxmismatch)  
    
    # Compute coverage numbers
    for (i in seq_along(allPatterns)) {
        
        begin <- allPatterns[[i]]@start[1]
        end <- allPatterns[[i]]@start[1] - 1 + allPatterns[[i]]@width[1]
        coverage[, seq.int(begin, end)] <- coverage[, seq.int(begin, end)] + 1
    }
    
    ncov <- max(coverage)
    
    start <- sapply(allPatterns, function(x) x@start) - 1
    end <- start + sapply(allPatterns, function(x) x@width)
    
    
    peptideMap <- matrix(0, ncol = length(AAString), nrow = ncov + 3)
    colnames(peptideMap) <- strsplit(as.character(AAString), "")[[1]]
    rownames(peptideMap) <- seq.int(1:nrow(peptideMap))
    
    if (is.null(scores) == TRUE){
        for (i in seq_along(start)){
            
            if (i == 1) {
                peptideMap[1, start[i]:end[i]] <- 1
            } else {
                j <- which.min(rowSums(peptideMap[, start[i]:end[i]] > 0))
                peptideMap[j, start[i]:end[i]] <- 1
            }
        }
    } else {
        for (i in seq_along(start)){
            
            if (i == 1) {
                peptideMap[1, start[i]:end[i]] <- scores[i]
            } else {
                j <- which.min(rowSums(peptideMap[, start[i]:end[i]] > 0))
                peptideMap[j, start[i]:end[i]] <- scores[i]
            }
        }
    }
    peptideMap[peptideMap == 0] <- NA
    averageMap <- apply(peptideMap, 2, function(x) 1/mean(1/x, na.rm = TRUE))
    averageMap[is.nan(averageMap)] <- 1
    averageMap <- -log10(averageMap)
    averageMap <- t(as.matrix(averageMap))
    rownames(averageMap) <- paste0(seq.int(1:nrow(averageMap)))
    
    sc <- scale_fill_manual(name, 
                            values = c("white", viridis(n = nlevels(factor(averageMap))), alpha = 0.7),
                            labels = substring(levels(factor(averageMap))[seq(1, nlevels(factor(averageMap)), by = by)], 1, 4),
                            breaks = levels(factor(averageMap))[seq(1, nlevels(factor(averageMap)), by = by)], drop = FALSE)
    
    for (i in 1:ceiling(length(coverage)/n)) {
        
        
        if (i < ceiling(length(coverage)/n)){    
            mygrid <- expand.grid(X = factor(1:n), Y = rownames(averageMap[, n * (i - 1) + 1:n, drop = FALSE]))
            mygrid$Z <- factor(c(t(averageMap[, n * (i - 1) + 1:n])), levels = levels(factor(averageMap)))
            
            plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc +
                theme_classic() + scale_x_discrete(breaks = 1:n, labels = colnames(averageMap[, n * (i - 1) + 1:n, drop = FALSE])) +
                theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
        } else {
            mygrid$Z <- factor(c(t(averageMap[, (n * (i-1)):ncol(averageMap)])), levels = levels(factor(averageMap)))
            
            plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc + 
                theme_classic() + 
                scale_x_discrete(breaks = 1:(ncol(averageMap)%%n + 1), labels = colnames(averageMap[, (n * (i-1)):ncol(averageMap), drop = FALSE])) +
                theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
        }
    }
    
    return(plot.list)
}


ComputeAverageMap <- function(AAString, 
                              peptideSeqs,
                              numlines = 5,
                              maxmismatch = 0,
                              by = 5,
                              scores = NULL,
                              name = "-log10 p values",
                              threshold = -log10(0.05)){
    
    # Test
    stopifnot("AAString must be an object of class AAString"= class(AAString) == "AAString")
    #stopifnot("peptideSeqs must be a character vector"= is.character(peptideSeqs) == "TRUE")
    
    # Storage and global variables
    plot.list <- list()
    coverage <- matrix(0, ncol = length(AAString), nrow = 1)
    n <- ceiling(length(coverage)/numlines)
    colnames(coverage) <- strsplit(as.character(AAString), "")[[1]]
    
    # Compute AA stringset and match to dictionary
    peptideset <- AAStringSet(x = peptideSeqs)
    allPatterns <- matchPDict(pdict = peptideset,
                              subject = AAString,
                              max.mismatch = maxmismatch)  
    
    # Compute coverage numbers
    for (i in seq_along(allPatterns)) {
        
        begin <- allPatterns[[i]]@start[1]
        end <- allPatterns[[i]]@start[1] - 1 + allPatterns[[i]]@width[1]
        coverage[, seq.int(begin, end)] <- coverage[, seq.int(begin, end)] + 1
    }
    
    ncov <- max(coverage)
    
    start <- sapply(allPatterns, function(x) x@start) - 1
    end <- start + sapply(allPatterns, function(x) x@width)
    
    
    peptideMap <- matrix(0, ncol = length(AAString), nrow = ncov + 3)
    colnames(peptideMap) <- strsplit(as.character(AAString), "")[[1]]
    rownames(peptideMap) <- seq.int(1:nrow(peptideMap))
    
    if (is.null(scores) == TRUE){
        for (i in seq_along(start)){
            
            if (i == 1) {
                peptideMap[1, start[i]:end[i]] <- 1
            } else {
                j <- which.min(rowSums(peptideMap[, start[i]:end[i]] > 0))
                peptideMap[j, start[i]:end[i]] <- 1
            }
        }
    } else {
        for (i in seq_along(start)){
            
            if (i == 1) {
                peptideMap[1, start[i]:end[i]] <- scores[i]
            } else {
                j <- which.min(rowSums(peptideMap[, start[i]:end[i]] > 0))
                peptideMap[j, start[i]:end[i]] <- scores[i]
            }
        }
    }
    peptideMap[peptideMap == 0] <- NA
    averageMap <- apply(peptideMap, 2, function(x) 1/mean(1/x, na.rm = TRUE))
    averageMap[is.nan(averageMap)] <- 1
    averageMap <- -log10(averageMap)
    averageMap <- t(as.matrix(averageMap))
    rownames(averageMap) <- paste0(seq.int(1:nrow(averageMap)))
    
   
    return(averageMap = averageMap)
}    

plotAverageMaps <- function(averageMaps,
                            name = "-log10 p value",
                            numlines = 2,
                            by = 5){
    
    stopifnot(class(averageMaps) == "list")
    plot.list <- list()
    
    map <- do.call(rbind, averageMaps)
    coverage <- ncol(map)
    n <- ceiling(coverage/numlines)
    
    sc <- scale_fill_manual(name, 
                            values = c("white", viridis(n = nlevels(factor(map))), alpha = 0.7),
                            labels = substring(levels(factor(map))[seq(1, nlevels(factor(map)), by = by)], 1, 4),
                            breaks = levels(factor(map))[seq(1, nlevels(factor(map)), by = by)], drop = FALSE)

    for (i in 1:ceiling(coverage/n)) {
        
        
        if (i < ceiling(coverage/n)){    
            mygrid <- expand.grid(X = factor(1:n), Y = rownames(map[, n * (i - 1) + 1:n, drop = FALSE]))
            mygrid$Z <- factor(c(t(map[, n * (i - 1) + 1:n])), levels = levels(factor(map)))
            
            plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 1) + sc +
                theme_classic() + scale_x_discrete(breaks = 1:n, labels = colnames(map[, n * (i - 1) + 1:n, drop = FALSE])) +
                theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + scale_y_discrete(breaks = levels(mygrid$Y), labels = levels(mygrid$Y)) + 
                ylab("dAb") + xlab("AA sequence")
        } else {
            mygrid$Z <- factor(c(t(map[, (n * (i-1)):ncol(map)])), levels = levels(factor(map)))
            
            plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 1) + sc + 
                theme_classic() + theme(legend.position = "none") + 
                scale_x_discrete(breaks = 1:(ncol(map)%%n + 1), labels = colnames(map[, (n * (i-1)):ncol(map), drop = FALSE])) +
                scale_y_discrete(breaks = levels(mygrid$Y), labels = levels(mygrid$Y)) + 
                ylab("dAb") + xlab("AA sequence")
                
        }
    }
    
    return(plot.list)
}
