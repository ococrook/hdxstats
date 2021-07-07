## uptake kinetics 
differentialUptakeKinetics <- function(object,
                                       feature = NULL,
                                       design = NULL,
                                       formula = NULL,
                                       start = list(a = NULL, b = 0.001,  d = NULL, p = 1)){
    
    .out <- NULL
    
    if (!is.null(design)){
        colnames(object)[[1]] <- design
    }
    
    ## move object to Long Format
    ldf <- longFormat(object)
    
    ## add time points using regex
    ldf$timepoint <- as.numeric(str_match(ldf$colname, "X\\s*(.*?)\\s*rep")[, 2])
    
    ## add replicate number using regex
    ldf$replicates <- as.numeric(str_match(ldf$colname, "rep\\s*(.*)\\s*cond")[, 2])
    
    #add charge number using regex
    ldf$chargestate <- as.numeric(str_match(ldf$rowname, "_\\s*(.*)")[, 2])
    
    #add conditions using regex
    ldf$condition <- as.factor(str_match(ldf$colname, "cond\\s*(.*)")[, 2])
    
    if (length(feature) ==1){
        if (!is.null(feature)){
            if (is.null(formula)){
                formula <- value ~ a * (1 - exp(-b*(timepoint)^p)) + d
            }
            
            data <- data.frame(ldf[ldf$rowname == feature, ])
            data <- data[,!colSums(is.na(data)) == nrow(data)]
            data <- na.omit(data)
            if (nrow(data) < 3) {
                print("too few data points to fit model")
                return()
            }
            
            if (is.null(start$a)){
                start$a <- max(data$value)
            }
            
            if (is.null(start$d) & length(start) > 2){
                start$d <- min(data$value)
            }

            nonlin_mod <- try(nlsLM(data = data, 
                                  formula = formula, 
                                  start = start,
                                  control = nls.lm.control(maxiter = 500, ftol = 10^{-8}),
                                  trace = FALSE, 
                                  lower = rep(0, length(start)), algorithm = "LM", na.action = na.exclude))
            
            if(inherits(nonlin_mod, "try-error")){
                print("model fit failed, likely exessive missing values")
                return(nonlin_mod)
            }
            
            
            myPredict1 <- expand.grid(timepoint = seq(min(data$timepoint, 0), max(data$timepoint), 0.5))  
            #expand.grid here in case your model has more than one variable
            #Caution, extrapolating well beyond the data
            myPredict1$fit <- predict(nonlin_mod, newdata = myPredict1) 
            
            datalist <- group_split(data, condition)
            
            nlmod <- lapply(datalist,  function(x){nonlin_mod <- 
                try(nlsLM(data = x, 
                        formula = formula, 
                        start = start,
                        control = nls.lm.control(maxiter = 500, ftol = 10^{-8}),
                        trace = FALSE, 
                        lower = rep(0, length(start)), algorithm = "LM", na.action = na.exclude))})
            
            if (any(sapply(nlmod, function(x) inherits(x, "try-error")))){
                print("Could not fit model, likely exessive missing values")
                return(nlmod)
            }
            
            myPredict <- lapply(1:length(nlmod), function(j)
                expand.grid(timepoint = seq(min(data$timepoint, 0), max(data$timepoint), 0.5)))  
            #expand.grid here in case your model has more than one variable
            #Caution, extrapolating well beyond the data
            myPredict <- lapply(1:length(nlmod), function(j){
                myPredict[[j]] <- predict(nlmod[[j]], newdata = myPredict[[j]])
                return(myPredict[[j]])}) 
            
            myPredictjoin <- do.call(cbind, myPredict)
            colnames(myPredictjoin) <- sapply(datalist, function(x) x$condition[1])
            myPredictjoin <- data.frame(myPredictjoin)
            myPredictjoin$timepoint <- unlist(expand.grid(timepoint = seq(min(data$timepoint, 0), max(data$timepoint), 0.5)))
            myPredict_long <- myPredictjoin %>% pivot_longer(cols = 1:length(datalist))
            colnames(myPredict_long) <- c("timepoint", "condition", "fit")
            df <- data.frame(myPredict_long)
            
            df$condition <- gsub("X", "", df$condition)
            
            gg <- ggplot(data, aes(y = value,
                                   x = timepoint,
                                   color = condition)) + geom_point(size = 4) + theme_classic() + 
                geom_line(data = myPredict1, aes(x = timepoint, y = fit),
                          inherit.aes = FALSE, size = 2, col = mycolours[2]) + 
                geom_line(data = df, aes(x = timepoint, y = fit, color = condition), inherit.aes = FALSE,  size = 2) + 
                labs(x = "Exposure", y = "Deuterium Incoperation", color = "condition",
                     title = paste0(data$rowname[1])) + 
                scale_color_viridis(alpha = 0.7, discrete = TRUE)
            
            
            
            
            .out <- list(model_null = nonlin_mod, model_alt = nlmod, gg = gg)
            
        }
        
    }
    
    return(.out)
}

computeRSS <- function(nlmod_null,
                       nlmod_alt){
    
    # compute residual sum of squares
    RSS0 <- sum(resid(nlmod_null)^2)
    RSS1 <- sum(unlist(sapply(nlmod_alt, function(x) resid(x)))^2)
    
    # compute number of parameters
    p_null <- length(summary(nlmod_null)$parameters[,1])
    p_alt <- 2 * p_null
    
    # compute first degree of freedom
    d1 <- p_alt - p_null
    
    # compute second degrees of freedom
    d2 <- stats::nobs(nlmod_null) - p_alt
    
    .res <- list(RSS0 = RSS0, RSS1 = RSS1, d1 = d1, d2 = d2)
    
    return(.res)
}


computeFstat <- function(RSS0,
                         RSS1,
                         d1,
                         d2){
    
    numerator <- (RSS0 - RSS1)/d1
    denomenator <- RSS1/d2
    Fstat <- numerator/denomenator
 
    return(list(Fstat = Fstat, numerator = numerator, denomenator = denomenator))   
}


computePval <- function(Fstat,
                       d1,
                       d2){
    
    res <- pf(q = Fstat, df1 = d1, df2 = d2, lower.tail = FALSE)
    
    return(res)
}

hdxEbayes <- function(RSS0,
                      RSS1,
                      d1,
                      d2){
    
    numerator <- (RSS0 - RSS1)/d1
    denomenator <- RSS1/d2
    
    # squeeze variances
    newvar <- limma::squeezeVar(var = denomenator, df = d2)
    
    # compute moderated F statistic
    modFstat <- numerator/newvar$var.post
    
    # Estimate new degrees of freedom empirically

    # compute pvalues
    pvalues <- pf(q = modFstat, df1 = d1, df2 = d2, lower.tail = FALSE)
    
    # compute FDR
    fdr <- p.adjust(pvalues, method = "BH")
    
    return(list(pvalues = pvalues, fdr = fdr))
}

ttestUptakeKinetics <- function(object,
                                feature = NULL,
                                design = NULL,
                                formula = NULL){
    
    .out <- NULL
    
    if (!is.null(design)){
        colnames(object)[[1]] <- design
    }
    
    ## move object to Long Format
    ldf <- longFormat(object)
    
    ## add time points using regex
    ldf$timepoint <- as.numeric(str_match(ldf$colname, "X\\s*(.*?)\\s*rep")[, 2])
    
    ## add replicate number using regex
    ldf$replicates <- as.numeric(str_match(ldf$colname, "rep\\s*(.*)\\s*cond")[, 2])
    
    #add charge number using regex
    ldf$chargestate <- as.numeric(str_match(ldf$rowname, "_\\s*(.*)")[, 2])
    
    #add conditions using regex
    ldf$condition <- as.factor(str_match(ldf$colname, "cond\\s*(.*)")[, 2])

    datalist <- group_split(data.frame(ldf), c(rowname))
    tres <- list()
    for (j in seq_along(datalist)){
        
        data <- na.omit(datalist[[j]]) %>% group_by(timepoint) %>% filter(n() > 3)
        if (nrow(data) > 0){
            tres[[j]] <- na.omit(datalist[[j]]) %>% group_by(timepoint) %>% filter(n() > 3) %>% do(tidy(t.test(value ~ condition, data = .)))
            tres[[j]]$rownames <- datalist[[j]]$rowname[1]  
        }
    }
    
    res <- do.call(rbind, tres)
    res$fdr <- p.adjust(res$p.value, method = "BH")
    
    return(res)
}

lmUptakeKinetics <- function(object,
                             feature = NULL,
                             design = NULL,
                             formula = NULL){
    
    .out <- NULL
    
    if (!is.null(design)){
        colnames(object)[[1]] <- design
    }
    
    ## move object to Long Format
    ldf <- longFormat(object)
    
    ## add time points using regex
    ldf$timepoint <- as.numeric(str_match(ldf$colname, "X\\s*(.*?)\\s*rep")[, 2])
    
    ## add replicate number using regex
    ldf$replicates <- as.numeric(str_match(ldf$colname, "rep\\s*(.*)\\s*cond")[, 2])
    
    #add charge number using regex
    ldf$chargestate <- as.numeric(str_match(ldf$rowname, "_\\s*(.*)")[, 2])
    
    #add conditions using regex
    ldf$condition <- as.factor(str_match(ldf$colname, "cond\\s*(.*)")[, 2])
    
    ## convert timepoitn to factor
    ldf$timepoint <- factor(ldf$timepoint)
    
    datalist <- group_split(data.frame(ldf), c(rowname))
    
    lmres <- list()
    for (j in seq_along(datalist)){
        lmres[[j]] <- data.frame(summary(datalist[[j]] %>% group_by(timepoint) %>% 
                                lmerTest::lmer(value ~ condition + timepoint + condition:timepoint + (1|replicates), data = .))$coefficients)
        colnames(lmres[[j]])[5] <- "p.value"
        lmres[[j]]$rownames <- datalist[[j]]$rowname[1]
    }
    
    pvals <- lapply(lmres, function(x) x[,"p.value"])
    pvals <- sapply(pvals, '[', seq(max(sapply(pvals, length))))
    
    id <- min(which(sapply(lmres, function(x) nrow(x) == nrow(pvals))))
    rownames(pvals) <- rownames(lmres[[id]])
    colnames(pvals) <- sapply(lmres, function(x) x$rownames[1])
    
    fdr <- matrix(p.adjust(c(pvals), method = "BH"), nrow = nrow(pvals))
    rownames(fdr) <- rownames(pvals)
    colnames(fdr) <- colnames(pvals)
    
    return(fdr)
}


processFunctional <- function(object,
                              qfeature){
    

    # object compute from running non-linear models models
    res <- lapply(1:length(object),
                  function(z) try(computeRSS(nlmod_null = object[[z]]$model_null,
                                             nlmod_alt = object[[z]]$model_alt)))
    
    # get rownames
    rw <- rownames(qfeature)[[1]][which(!sapply(res, function(x) class(x)) == "try-error")]

    # remove data with error
    res_filtered <- res[which(!sapply(res, function(x) class(x)) == "try-error")]
    
    # Compute F Statistics
    Fstats <- lapply(1:length(res_filtered),
                     function(z) computeFstat(RSS0 = res_filtered[[z]]$RSS0,
                                              RSS1 = res_filtered[[z]]$RSS1,
                                              d1 = res_filtered[[z]]$d1,
                                              d2 = res_filtered[[z]]$d2))
    
    # Compute p-values
    pvals <- lapply(1:length(Fstats), function(z) computePval(Fstat = Fstats[[z]]$Fstat, 
                                                              d1 = res_filtered[[z]]$d1, 
                                                              d2 = res_filtered[[z]]$d2))

    # fdr correction
    fdr <- p.adjust(unlist(pvals), method = "BH")

    # empirical Bayes analysis
    ebayesres <- hdxEbayes(RSS0 = sapply(res_filtered, function(x) x$RSS0),
                           RSS1 = sapply(res_filtered, function(x) x$RSS1),
                           d1 = sapply(res_filtered, function(x) x$d1),
                           d2 = sapply(res_filtered, function(x) x$d2))

    # get names for emprical bayes    
    names(ebayesres$fdr) <- names(ebayesres$pvalues) <- rw

    Fstatres <- lapply(Fstats, unlist)
    names(Fstatres) <- rw
    names(pvals) <- rw
    names(ebayesres$pvalues) <- rw
    
    .out <- list(Fstats = lapply(Fstats, unlist),
                 pvals = pvals,
                 fdr = fdr,
                 ebayesres = ebayesres,
                 fitcomplete = which(!sapply(res, function(x) class(x)) == "try-error"))
    return(.out)
    
}

forestPlot <- function(params, condition = c(1,2)) {
    
    
    tmp <- summary(params$model_alt[[1]])$parameters
    tmp2 <- summary(params$model_alt[[2]])$parameters
    
    df <- data.frame(tmp)
    df$rownames <- rownames(df)
    df2 <- data.frame(tmp2)
    df2$rownames <- rownames(df2)
    df$condition <- condition[1]
    df2$condition <- condition[2]
    
    # now need predictions
    data <- params$gg$data
    
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






