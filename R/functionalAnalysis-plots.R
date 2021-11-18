##' A manhatten plot for epitope mapping
##' @title manhatten plot
##' @param params An object of class `HdxStatRes`
##' @param sequences A character vector containing the measured peptide sequences
##' @param difference A numeric vector with deuterium differences for each peptide
##' @param region The start and end of the sequences provided.
##' @param nrow The number of rows to plot the manhatten plot over. Useful for larger
##' proteins
##' @return use for side effect which returns manhatten plot 
##' @md 
##'
##' @rdname functional-plots
manhattenplot <- function(params,
                          sequences,
                          difference = 0,
                          region = NULL, 
                          nrow = 1){
    
    butterfly <- matrix(0, ncol = length(unique(sequences)), nrow = 1)
    colnames(butterfly) <- unique(sequences)
    butterfly <- params@results$ebayes.fdr[unique(sequences)]
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

##' Generate an epitope map by plotting the signifcantly changed peptides
##' with respect to a differential HDX-MS experiment
##' @title Epitope Map
##' @param AAString An object of class `AAString` for the protein of interest
##' @param peptideSeqs A character vector of peptide sequences
##' @param numlines The number of lines to plot the protein over. Useful for larger
##'  proteins. Default is 5.
##' @param maxmismatch A numeric indicating if incorrect mapping is allowed. Number 
##'  indicated the number of mismatched amino acids. Default is 0.
##' @param by A value to indicate the legend breaks. Default is NULL.
##' @param scores A numeric vector indicating score to be used for plotting
##' @param name The name of the legend for the score plotting. 
##'  Default is "-log10 p values".
##' @param threshold The threshold used to determine significance. Default is
##' `-log10(0.05)`. Note the log scale.
##' @md
##' 
##' @rdname functional-plots
plotEpitopeMap <- function(AAString, 
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
##' Generate an epitope map by plotting the signifcantly changed peptides
##' with respect to a differential HDX-MS experiment. Plots FDR rather than
##' thresholding.
##' @param AAString An object of class `AAString` for the protein of interest
##' @param peptideSeqs A character vector of peptide sequences
##' @param numlines The number of lines to plot the protein over. Useful for larger
##'  proteins. Default is 5.
##' @param maxmismatch A numeric indicating if incorrect mapping is allowed. Number 
##'  indicated the number of mismatched amino acids. Default is 0.
##' @param by A value to indicate the legend breaks. Default is NULL.
##' @param scores A numeric vector indicating score to be used for plotting
##' @param name The name of the legend for the score plotting. 
##'  Default is "-log10 p values".
##' @param threshold The threshold used to determine significance. Default is
##' `-log10(0.05)`. Note the log scale.
##' @md
##' 
##' @rdname functional-plots
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
                                labels = substring(levels(factor(peptideMap)), 1, 4),
                                breaks = levels(factor(peptideMap)), drop = FALSE)
        
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
##' Generate an epitope barcode on the residue scale using a harmonic mean
##' averageing approach.
##' @param AAString An object of class `AAString` for the protein of interest
##' @param peptideSeqs A character vector of peptide sequences
##' @param numlines The number of lines to plot the protein over. Useful for larger
##'  proteins. Default is 5.
##' @param maxmismatch A numeric indicating if incorrect mapping is allowed. Number 
##'  indicated the number of mismatched amino acids. Default is 0.
##' @param by A value to indicate the legend breaks. Default is NULL.
##' @param scores A numeric vector indicating score to be used for plotting. Most 
##'  likely adusted p-values.
##' @param name The name of the legend for the score plotting. 
##'  Default is "-log10 p values".
##' @param threshold The threshold used to determine significance. Default is
##' `-log10(0.05)`. Note the log scale.
##' @md
##' 
##' @rdname functional-plots
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
                            labels = substring(levels(factor(averageMap)), 1, 4),
                            breaks = levels(factor(averageMap)), drop = FALSE)
    
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
##' Underlying computation for compute residue level mappings. Will performed 
##' harmonic mean averaging to obtain residue level results.
##' @param AAString An object of class `AAString` for the protein of interest
##' @param peptideSeqs A character vector of peptide sequences
##' @param numlines The number of lines to plot the protein over. Useful for larger
##'  proteins. Default is 5.
##' @param maxmismatch A numeric indicating if incorrect mapping is allowed. Number 
##'  indicated the number of mismatched amino acids. Default is 0.
##' @param by A value to indicate the legend breaks. Default is NULL.
##' @param scores A numeric vector indicating score to be used for plotting. Most 
##'  likely adusted p-values.
##' @param name The name of the legend for the score plotting. 
##'  Default is "-log10 p values".
##' @param threshold The threshold used to determine significance. Default is
##' `-log10(0.05)`. Note the log scale.
##' @md
##' 
##' @rdname functional-plots
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
##' Plotting for comparing average maps. This function will simultaneous plot
##'  several barcodes ontop of each other so the comparison is easier
##' @param averageMaps A list of average maps generated by the `computeAverageMaps`
##'  function
##' @param name The name of the legend for the score plotting. 
##'  Default is "-log10 p values".
##' @param numlines The number of lines to plot the protein over. Useful for larger
##'  proteins. Default is 2.
##' @param by A value to indicate the legend breaks. Default is NULL.
##' @md
##' 
##' @rdname functional-plots
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
                            values = c("white", viridis(n = nlevels(factor(map)), alpha = 0.7)), drop = FALSE,
                            labels = substring(levels(factor(map)), 1, 4),
                            breaks = levels(factor(map)))

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

##' Generate a forest plot
##' 
##' @param params An object of class `HdxStatRes`
##' @param condition If there are multiple conditions which ones to plot. Default
##' is `c(1,2)`
##' @return Side effect produces a forest plot
##' @md
##' 
##' @rdname functional-plots
forestPlot <- function(params, condition = c(1,2)) {
    
    stopifnot("params must be an HdxStatModel"=class(params)=="HdxStatModel")
    
    model1 <- params@alternative@nlsmodels[[1]]
    model2 <- params@alternative@nlsmodels[[2]]
    
    tmp <- summary(model1)$parameters
    tmp2 <- summary(model1)$parameters
    
    df <- data.frame(tmp)
    df$rownames <- rownames(df)
    df2 <- data.frame(tmp2)
    df2$rownames <- rownames(df2)
    df$condition <- condition[1]
    df2$condition <- condition[2]
    
    # now need predictions
    data <- params@vis$data
    
    # compute quantities of interest
    out <- data %>% group_by(timepoint) %>% do(tidy(t.test(value ~ condition, data = .)))
    df3 <- data.frame(cbind(out$estimate, out$conf.low, out$conf.high))
    df3$rownames <- paste0("Timepoint ", out$timepoint)
    df3$condition <- rep(c("Deuterium Difference"), nrow(df3))
    colnames(df3) <- c("Estimate", "confL", "confU", "rownames", "condition")
    
    #Compute confidence intervals
    df$confL <- df$Estimate - 1.96*df$Std..Error
    df$confU <- df$Estimate + 1.96*df$Std..Error
    df2$confL <- df2$Estimate - 1.96*df2$Std..Error
    df2$confU <- df2$Estimate + 1.96*df2$Std..Error
    
    # data frame or ggplot
    mydf <- rbind(df[,colnames(df3)], df2[,colnames(df3)], df3)
    
    gg <- ggplot(mydf, aes(x = Estimate, y = rownames, xmin = confL,
                           xmax = confU, col = factor(condition))) + 
        geom_point(size = 2, position = position_dodge(width = 0.5)) + 
        theme_classic() + labs(x = "Effect Size", y = "Effect", col = "condition", title = "Forest plot for Effect sizes") + 
        geom_errorbarh(height=0.3, position = position_dodge(width = 0.5), lwd = 1.5) + 
        scale_color_viridis(discrete = TRUE, begin = 0.2, direction = 1) + theme(text = element_text(size = 16)) + 
        geom_vline(xintercept = 0, linetype = "dashed")
    
    print(gg)
}
