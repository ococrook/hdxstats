##' Generics for hdxstats package
##' 
##' @exportMethod vcov
##' @param object An instance of `HdxStatModel`
##' @return The variance-covariance matrix of the parameters
##' 
##' @rdname hdxstat-methods
setMethod("vcov", "HdxStatModel", 
          function(object){
            .nullvcov <- vcov(object@nullmodel)
            .altvcov <- lapply(object@alternative@nlsmodels, vcov)
            .out <- list(nullvcov = .nullvcov, altvcov = .altvcov)
            .out
            })
##' @exportMethod anova
##' @param object An instance of `HdxStatModel`
##' @return An analysis of variance
##' 
##' @rdname hdxstat-methods
setMethod("anova", "HdxStatModel",
          function(object){
              do.call("anova", c(list(object@nullmodel), object@alternative@nlsmodels))  
          })
##' @exportMethod logLik
##' @param object An instance of `HdxStatModel`
##' @return The log likelihod of the fitted model assuming normally distributed residuals.
##' 
##' @rdname hdxstat-methods
setMethod("logLik", "HdxStatModel",
          function(object){
              .out <- sapply(c(list(object@nullmodel), object@alternative@nlsmodels), stats::logLik)
              names(.out) <- c("null", paste0("alt", seq.int(length(object))))
              .out
          })
##' @exportMethod likRatio
##' @param object An instance of `HdxStatModel`
##' @return The log likelihood ratio between the alternative and null models
##' 
##' @rdname hdxstat-methods
setMethod("likRatio", "HdxStatModel",
          function(object){
              .loglik <- hdxstats::logLik(object)
              .lr <- 2 * (sum(.loglik[-1]) - .loglik[1])
              names(.lr) <- "logLR"
              .lr
          })
##' @exportMethod wilk
##' @param object An instance of `HdxStatModel`
##' @return Applies Wilk's theorem to generate a p-value based on the likelihood ratio test
##' 
##' @rdname hdxstat-methods
setMethod("wilk", "HdxStatModel",
          function(object){
            .lr <- likRatio(object)
            palt <- sum(sapply(object@alternative@nlsmodels, function(x) summary(x)$df[1]))
            pnull <- summary(object@alternative@nlsmodels[[1]])$df[1]
            .pval <- pchisq(.lr, df = palt - pnull, lower.tail = FALSE)
            names(.pval) <- "p-value"
            .pval
          })
##' @exportMethod coef
##' @param object An instance of `HdxStatModel`
##' @return Returns the coefficients of the null and alternative models
##' 
##' @rdname hdxstat-methods
setMethod("coef", "HdxStatModel",
          function(object) {
              .out <- t(sapply(c(list(object@nullmodel), object@alternative@nlsmodels), coef))
              rownames(.out) <- c("null", paste0("alt", seq.int(length(object))))
              .out
          })
##' @exportMethod deviance
##' @param object An instance of `HdxStatModel`
##' @return Returns the deviance of the fitted models
##' 
##' @rdname hdxstat-methods
setMethod("deviance", "HdxStatModel",
          function(object) {
              .out <- sapply(c(list(object@nullmodel), object@alternative@nlsmodels), deviance)
              names(.out) <- c("null", paste0("alt", seq.int(length(object))))
              .out
          })
##' @exportMethod residuals
##' @param object An instance of `HdxStatModel`
##' @return The residuals from the fitted models
##' 
##' @rdname hdxstat-methods
setMethod("residuals", "HdxStatModel",
          function(object) {
              .out <- lapply(c(list(object@nullmodel), object@alternative@nlsmodels), residuals)
              names(.out) <- c("null", paste0("alt", seq.int(length(object))))
              .out
          })
##' @exportMethod summary
##' @param object An instance of `HdxStatmodel`
##' @return Returns a summary of the fitted models.
##' 
##' @rdname hdxstat-methods
setMethod("summary", "HdxStatModel",
          function(object) {
              .out <- lapply(c(list(object@nullmodel), object@alternative@nlsmodels), summary)
              names(.out) <- c("null", paste0("alt", seq.int(length(object))))
              .out
          })
##' @exportMethod fitUptakeKinetics
##' @param  object An instance of class `QFeatures`
##' @param feature The rowname of feature to be modelled
##' @param design The design defining conditions and replicates. Default is NULL.
##' @param formula The formula for the non-linear fit.
##' @param start The initial guess for the parameters. Parameters must match formula
##' @return Returns an instance of `HdxStatModels`
##' 
##' @rdname hdxstat-methods
setMethod("fitUptakeKinetics", "QFeatures",
          function(object, 
                   feature = NULL,
                   design = NULL,
                   formula = NULL,
                   start = list(a = NULL, b = 0.001,  d = NULL, p = 1),
                   maxAttempts = 5){
              .res <- lapply(feature,
                             function(x) differentialUptakeKinetics(object = object,
                                                                    feature = x,
                                                                    start = start,
                                                                    formula = formula,
                                                                    design = design, 
                                                                    maxAttempts = maxAttempts))
              .res <- .res[which(!sapply(.res, function(x) class(x)) == "try-error")]
              .res <- .res[which(sapply(.res, function(x) class(x) == "HdxStatModel"))]
              .res2 <- .hdxstatmodels(statmodels = .res)
              .res2
          })

