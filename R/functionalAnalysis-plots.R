##' A manhatten plot for a multi-level testing prodcedure based on the 
##' haromonic mean p-value
##' @title hmpWindow (Currently in development and not used)
##' @param params An object of class HdxStatRes
##' @param sequences A character vector containing the measured peptide sequences
##' @param region The start and end of the sequences provided.
##' @param interval interval size to average over.
##' @return use for side effect which returns manhattan plot 
##' @md 
##'
##' @rdname functional-plots
hmpWindow <- function(params,
                      sequences,
                      region,
                      interval = 20){
    
    stopifnot("params must be a class of HdxStatRes"=is(params, "HdxStatRes"))
    
    butterfly <- matrix(0, ncol = length(unique(sequences)), nrow = 1)
    colnames(butterfly) <- unique(sequences)
    butterfly <- params@results$ebayes.pvals[unique(sequences)]
    butterflydf <- as.data.frame(butterfly)
    butterflydf$Sequence <- rownames(butterflydf)
    colnames(butterflydf)[1] <- "p_value"
    butterfly_long <- butterflydf
    butterfly_long$position <- rep(seq.int(nrow(butterfly_long)/1), each = 1)
    butterfly_long$region <- unique(region[, c("Start", "End")])
    
    # compute harmonic mean p-value
    phmp <- vector(mode = "numeric", nrow(butterfly_long) - interval)
    region <- matrix(NA, nrow = nrow(butterfly_long) - interval, ncol = 2)
    for (i in seq.int(nrow(butterfly_long) - interval)){
        phmp[i] <- p.hmp(butterfly_long[seq.int(i, i + interval), "p_value"], L = nrow(butterfly_long))
        region[i, ] <- c(min(butterfly_long[seq.int(i, i + interval), "region"]), 
                         max(butterfly_long[seq.int(i, i + interval), "region"]))
    }
    
    
    plot.list <- list()
    r <- nrow
    xannot <- paste0("[", region[,1], ",", region[,2], "]")
    df <- data.frame(p_value = phmp, region = xannot)
    df$position <- seq.int(nrow(df))
    for (i in seq.int(r)) {
        
        plot.list[[i]] <- ggplot(df, aes(x = position,
                                         y = -log10(p_value/((interval + 1)/nrow(butterfly_long))),
                                         group = -log10(p_value))) + 
            geom_point(size = 3, color = brewer.pal(n = 3, name = "Set2")[2]) + 
            theme_classic() + ylim(c(0, max(-log10(df$p_value)))) + 
            ylab("-log10 harmonic adjusted p-value") + xlab("region") + xlim(c(1, nrow(df))) + 
            scale_x_continuous(breaks = 1:nrow(df), labels = xannot) + 
            ggtitle("Manhatten hmp plot") + geom_hline(yintercept = 1.301, linetype = "dashed", colour = "red") + 
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 15))
    }
    
    return(plot.list = plot.list)
    
}
##' A manhattan plot for epitope mapping
##' @title manhattan plot
##' @param params An object of class `HdxStatRes`
##' @param sequences A character vector containing the measured peptides sequences,
##' if there are different charge states the vector must identify the charge state.
##' @param difference A numeric vector with deuterium differences for each peptide
##' @param region The start and end of the sequences provided. Columns must be 
##' called "Start" and "End".
##' @param nrow The number of rows to plot the manhatten plot over. Useful for larger
##' proteins
##' @return use for side effect which returns manhatten plot 
##' @md 
##'
##' @rdname functional-plots
manhattanplot <- function(params,
                          sequences,
                          difference = 0,
                          region = NULL, 
                          nrow = 1){
    
    stopifnot("Column names of region incorrect"=all(c("Start", "End") %in% colnames(region)))
    
    if(nrow(params@results) != length(unique(sequences))){
       warning("Sequences provided don't match the results from the statistical
             analysis. Double check sequences and subset if necessary.
             Results may not be coherent. Automatic subsetting has been applied.")
        # order important do not switch
        region <- region[which(sequences %in% rownames(params@results)), c("Start", "End")]
        sequences <- sequences[sequences %in% rownames(params@results)]

    }
    
    
    
    butterfly <- matrix(0, ncol = length(unique(sequences)), nrow = 1)
    colnames(butterfly) <- unique(sequences)
    butterfly <- params@results$ebayes.fdr[unique(sequences)]
    butterflydf <- as.data.frame(butterfly)
    butterflydf$Sequence <- rownames(butterflydf)
    butterflydf$protection <- 1*(difference[butterflydf$Sequence] > 0)
    butterflydf$protection[is.na(butterflydf$protection)] <- 0
    colnames(butterflydf)[1] <- "p_value"
    butterfly_long <- butterflydf
    butterfly_long$position <- rep(seq.int(nrow(butterfly_long)/1), each = 1)
    
    butterfly_long$region <-  unique(region[sequences %in% unique(sequences), c("Start", "End")])
    
    
    
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
    
    start <- sapply(allPatterns, function(x) x@start)
    end <- start + sapply(allPatterns, function(x) x@width) - 1
    
    
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
                    theme_bw() + 
                    scale_x_discrete(breaks = seq.int(n), 
                                     labels = paste0(colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]), "-",
                                                     seq.int(from = n* (i - 1) + 1, to = n *(i - 1)  + n, by = 1)
                                     )) +
                    theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                    ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            } else {
                mygrid <- expand.grid(X = factor(1:(ncol(peptideMap)%%n + 1)), Y = rownames(peptideMap[, (n * (i-1)):ncol(peptideMap), drop = FALSE]))
                mygrid$Z <- as.factor(c(t(peptideMap[, (n * (i-1)):ncol(peptideMap)])))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + 
                    scale_fill_manual(breaks = c(0,1), values = c("white", alpha("#1B7837", 0.7)))  +
                    theme_bw() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1),
                                                  labels = paste0(colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)]), "-",
                                                                  seq.int(n*(i-1), ncol(peptideMap)))) +
                    theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                    ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            }
        }
    } else{
        
        sc <- scale_fill_manual(name, 
                                values = c("white", brewer.pal(n = 3, name = "Set2")),
                                labels = c(" ", "Not Signifcant", "Significant"),
                                breaks = levels(factor(peptideMap))[seq(1, nlevels(factor(peptideMap)), by = by)], drop = FALSE)
        
        for (i in 1:ceiling(length(coverage)/n)) {
            
            
            if (i < ceiling(length(coverage)/n)){    
                mygrid <- expand.grid(X = factor(1:n), Y = rownames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]))
                mygrid$Z <- factor(c(t(peptideMap[, n * (i - 1) + 1:n])), levels = levels(factor(peptideMap)))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc +
                    theme_bw() + 
                    scale_x_discrete(breaks = seq.int(n), 
                                     labels = paste0(colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]), "-",
                                                                       seq.int(from = n* (i - 1) + 1, to = n *(i - 1)  + n, by = 1)
                        )) +
                    theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
                    ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y)))) + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))

            } else {
                mygrid$Z <- factor(c(t(peptideMap[, (n * (i-1)):ncol(peptideMap)])), levels = levels(factor(peptideMap)))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc + 
                    theme_bw() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1),
                                                  labels = paste0(colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)]), "-",
                                                                  seq.int(n*(i-1), ncol(peptideMap)))) +
                    theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                    ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
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
    
    start <- sapply(allPatterns, function(x) x@start)
    end <- start + sapply(allPatterns, function(x) x@width) - 1
    
    
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
                    theme_bw() + 
                    scale_x_discrete(breaks = seq.int(n), 
                                     labels = paste0(colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]), "-",
                                                     seq.int(from = n* (i - 1) + 1, to = n *(i - 1)  + n, by = 1)
                                     )) +
                    theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                    ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            } else {
                mygrid <- expand.grid(X = factor(1:(ncol(peptideMap)%%n + 1)), Y = rownames(peptideMap[, (n * (i-1)):ncol(peptideMap), drop = FALSE]))
                mygrid$Z <- as.factor(c(t(peptideMap[, (n * (i-1)):ncol(peptideMap)])))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + 
                    scale_fill_manual(breaks = c(0,1), values = c("white", alpha("#1B7837", 0.7)))  +
                    theme_bw() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1),
                                                  labels = paste0(colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)]), "-",
                                                                  seq.int(n*(i-1), ncol(peptideMap)))) +
                    theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
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
                    theme_bw() + 
                    scale_x_discrete(breaks = seq.int(n), 
                                     labels = paste0(colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]), "-",
                                                     seq.int(from = n* (i - 1) + 1, to = n *(i - 1)  + n, by = 1)
                                     )) +
                    theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                    ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                    scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
            } else {
                mygrid$Z <- factor(c(t(peptideMap[, (n * (i-1)):ncol(peptideMap)])), levels = levels(factor(peptideMap)))
                
                plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc + 
                    theme_bw() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1),
                                                  labels = paste0(colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)]), "-",
                                                                  seq.int(n*(i-1), ncol(peptideMap)))) +
                    theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
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
                            values = c("white", viridis(n = nlevels(factor(averageMap)), alpha = 0.7)),
                            labels = substring(levels(factor(averageMap)), 1, 4),
                            breaks = levels(factor(averageMap)), drop = FALSE)
    
    for (i in 1:ceiling(length(coverage)/n)) {
        
        
        if (i < ceiling(length(coverage)/n)){    
            mygrid <- expand.grid(X = factor(1:n), Y = rownames(averageMap[, n * (i - 1) + 1:n, drop = FALSE]))
            mygrid$Z <- factor(c(t(averageMap[, n * (i - 1) + 1:n])), levels = levels(factor(averageMap)))
            
            plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc +
                theme_bw() + 
                scale_x_discrete(breaks = seq.int(n), 
                                 labels = paste0(colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]), "-",
                                                 seq.int(from = n* (i - 1) + 1, to = n *(i - 1)  + n, by = 1)
                                 )) +
                theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
                scale_y_discrete(breaks = 1:length(levels(mygrid$Y)), labels = rep("", length(levels(mygrid$Y))))
        } else {
            mygrid$Z <- factor(c(t(averageMap[, (n * (i-1)):ncol(averageMap)])), levels = levels(factor(averageMap)))
            
            plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 0.7) + sc + 
                theme_bw() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1),
                                              labels = paste0(colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)]), "-",
                                                              seq.int(n*(i-1), ncol(peptideMap)))) +
                theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                ylim(c(0, max(coverage) + 1)) + ylab("") + xlab("AA sequence") + 
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
                theme_bw() + 
                scale_x_discrete(breaks = seq.int(n), 
                                 labels = paste0(colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]), "-",
                                                 seq.int(from = n* (i - 1) + 1, to = n *(i - 1)  + n, by = 1)
                                 )) +
                theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                ylim(c(0, max(coverage) + 1)) + scale_y_discrete(breaks = levels(mygrid$Y), labels = levels(mygrid$Y)) + 
                ylab("dAb") + xlab("AA sequence")
        } else {
            mygrid$Z <- factor(c(t(map[, (n * (i-1)):ncol(map)])), levels = levels(factor(map)))
            
            plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 1) + sc + 
                theme_bw() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1),
                                              labels = paste0(colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)]), "-",
                                                              seq.int(n*(i-1), ncol(peptideMap)))) +
                theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
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
    tmp2 <- summary(model2)$parameters
    
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
    
    return(gg)
}
##' Underlying computation for compute residue level differences plots. Uses
##' average differences at each residue. 
##' @param object An object of class `QFeatures` contains the hdx data.
##' @param AAString An object of class `AAString` for the protein of interest
##' @param peptideSeqs A character vector of peptide sequences
##' @param numlines The number of lines to plot the protein over. Useful for larger
##'  proteins. Default is 5.
##' @param maxmismatch A numeric indicating if incorrect mapping is allowed. Number 
##'  indicated the number of mismatched amino acids. Default is 0.
##' @param by A value to indicate the legend breaks. Default is NULL.
##' @param scores A numeric vector indicating score to be used for plotting. Most 
##'  likely adjusted p-values.
##' @param cols Columns for which to compute the difference. The difference is
##' relative to the first entry.
##' @param name The name of the legend for the score plotting. 
##'  Default is "-log10 p values".
##' @param threshold The threshold used to determine significance. Default is
##' `-log10(0.05)`. Note the log scale.
##' @md
##' 
##' @rdname functional-plots
hdxdifference <- function(object, 
                          AAString, 
                          peptideSeqs,
                          numlines = 5,
                          maxmismatch = 0,
                          by = 5,
                          cols = c(1,4),
                          scores = NULL,
                          name = "-log10 p value"){
    
    # Test
    stopifnot("AAString must be an object of class AAString"= class(AAString) == "AAString")
    #stopifnot("peptideSeqs must be a character vector"= is.character(peptideSeqs) == "TRUE")
    stopifnot("cols must be length 2"=length(cols) == 2)
    
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
    
    ## compute differences
    diff <- assay(object)[, cols[2]] - assay(object)[, cols[1]]
    
    for (i in seq_along(start)){
            
            if (i == 1) {
                peptideMap[1, start[i]:end[i]] <- diff[i]
            } else {
                j <- which.min(rowSums(peptideMap[, start[i]:end[i]] > 0))
                peptideMap[j, start[i]:end[i]] <- diff[i]
            }
    }
    
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
    
    
    
    peptideMap[peptideMap == 0] <- NA
    diffMap <- apply(peptideMap, 2, function(x) mean(x, na.rm = TRUE))
    diffMap[is.nan(diffMap)] <- 0
    diffMap <- t(as.matrix(diffMap))
    rownames(diffMap) <- paste0(seq.int(1:nrow(diffMap)))
    
    
    return(list(averageMap = averageMap, diffMap = diffMap))
}    
##' Plotting for comparing difference maps. This function will simultaneous plot
##'  several difference barcodes on top of each other so the comparison is easier
##' @param averageMaps A list of average maps generated by the `computeAverageMaps`
##'  function
##' @param diffMaps A list of difference maps generated by the `hdxdifference`
##' function 
##' @param name The name of the legend for the score plotting. 
##'  Default is "-log10 p values (singed)". Indicated the significance and direction
##' @param numlines The number of lines to plot the protein over. Useful for larger
##'  proteins. Default is 2.
##' @param by A value to indicate the legend breaks. Default is NULL.
##' @md
##' 
##' @rdname functional-plots
hdxheatmap <- function(averageMaps,
                       diffMaps,
                       name = "-log10 p value (signed)",
                       numlines = 2,
                       by = 5){

    stopifnot(class(averageMaps) == "list")
    stopifnot(class(diffMaps) == "list")
    plot.list <- list()
    
    map <- do.call(rbind, averageMaps)
    diff <- do.call(rbind, diffMaps)
    coverage <- ncol(map)
    n <- ceiling(coverage/numlines)
    
    map <- map * sign(diff)
    
    values <-  c(colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(11, 7)], alpha = 0.9)(sum(unique(c(map)) < 0)),
                 "white",
                 colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(6, 1)], alpha = 0.9)(sum(unique(c(map)) > 0)))
                 
    names(values) <- factor(sort(unique(c(map))))
    
    sc <- scale_fill_manual(name, 
                            values = values,
                            drop = FALSE,
                            labels = round(sort(unique(c(map))),3),
                            breaks = levels(factor(map)))
    
    for (i in 1:ceiling(coverage/n)) {
        
        
        if (i < ceiling(coverage/n)){    
            mygrid <- expand.grid(X = factor(1:n), Y = rownames(map[, n * (i - 1) + 1:n, drop = FALSE]))
            mygrid$Z <- factor(c(t(map[, n * (i - 1) + 1:n])), levels = levels(factor(map)))
            
            plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 1) + sc +
                theme_bw() + 
                scale_x_discrete(breaks = seq.int(n), 
                                 labels = paste0(colnames(peptideMap[, n * (i - 1) + 1:n, drop = FALSE]), "-",
                                                 seq.int(from = n* (i - 1) + 1, to = n *(i - 1)  + n, by = 1)
                                 )) +
                theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                ylim(c(0, max(coverage) + 1)) +
                theme(legend.position = "none") + ylim(c(0, max(coverage) + 1)) + scale_y_discrete(breaks = levels(mygrid$Y), labels = levels(mygrid$Y)) + 
                ylab("dAb") + xlab("AA sequence")
        } else {
            mygrid$Z <- factor(c(t(map[, (n * (i-1)):ncol(map)])), levels = levels(factor(map)))
            
            plot.list[[i]] <- ggplot(mygrid, aes(x = X, y = Y, fill = Z), show.legend = FALSE) + geom_tile(height = 1) + sc + 
                theme_bw() + scale_x_discrete(breaks = 1:(ncol(peptideMap)%%n + 1),
                                              labels = paste0(colnames(peptideMap[, (n * (i-1)):ncol(peptideMap)]), "-",
                                                              seq.int(n*(i-1), ncol(peptideMap)))) +
                theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                scale_y_discrete(breaks = levels(mygrid$Y), labels = levels(mygrid$Y)) + 
                ylab("dAb") + xlab("AA sequence")
            
        }
    }
    
    return(plot.list)
}

