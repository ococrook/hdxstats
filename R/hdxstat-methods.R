

setMethod("vcov", "HdxStatModel", 
          function(object){
            .nullvcov <- vcov(object@nullmodel)
            .altvcov <- lapply(object@alternative@nlsmodels, vcov)
            .out <- list(nullvcov = .nullvcov, altvcov = .altvcov)
            .out
            })

setMethod("anova", "HdxStatModel",
          function(object){
              do.call("anova", c(list(object@nullmodel), object@alternative@nlsmodels))  
          })

setMethod("logLik", "HdxStatModel",
          function(object){
              .out <- sapply(c(list(object@nullmodel), object@alternative@nlsmodels), logLik)
              names(.out) <- c("null", paste0("alt", seq.int(length(object))))
              .out
          })

setMethod("likRatio", "HdxStatModel",
          function(object){
              .loglik <- logLik(object)
              .lr <- 2 * (sum(.loglik[-1]) - .loglik[1])
              names(.lr) <- "logLR"
              .lr
          })

setMethod("wilk", "HdxStatModel",
          function(object){
            .lr <- likRatio(object)
            palt <- sum(sapply(.out@alternative@nlsmodels, function(x) summary(x)$df[1]))
            pnull <- summary(.out@alternative@nlsmodels[[1]])$df[1]
            .pval <- pchisq(.lr, df = palt - pnull, lower.tail = FALSE)
            names(.pval) <- "p-value"
            .pval
          })

setMethod("coef", "HdxStatModel",
          function(object) {
              .out <- t(sapply(c(list(object@nullmodel), object@alternative@nlsmodels), coef))
              rownames(.out) <- c("null", paste0("alt", seq.int(length(object))))
              .out
          })

setMethod("deviance", "HdxStatModel",
          function(object) {
              .out <- sapply(c(list(object@nullmodel), object@alternative@nlsmodels), deviance)
              names(.out) <- c("null", paste0("alt", seq.int(length(object))))
              .out
          })

setMethod("residuals", "HdxStatModel",
          function(object) {
              .out <- lapply(c(list(.out@nullmodel), .out@alternative@nlsmodels), residuals)
              names(.out) <- c("null", paste0("alt", seq.int(length(object))))
              .out
          })

setMethod("summary", "HdxStatModel",
          function(object) {
              .out <- lapply(c(list(.out@nullmodel), .out@alternative@nlsmodels), summary)
              names(.out) <- c("null", paste0("alt", seq.int(length(object))))
              .out
          })
