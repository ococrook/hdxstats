##' @name nls The nls class
##' @md
##' @rdname hdxstat-class
setOldClass("nls")
##' @name gg The ggplot class
##' @md
##' @rdname hdxstat-class
setOldClass("gg")
##' Classes for hdxstats package
##' 
##' @slot nlsmodels A `list()` containing nls models results in `nls` instances.
##' Each element must be a valid `nls` instance.
##' @param object An instance of class `nlsmodel`
##' @md
##' @rdname hdxstat-class
.nlsList <- setClass("nlsList",
                     slots = c(nlsmodels = "list"),
                     validity = function(object){
                        msg <- validMsg(NULL, NULL)
                        if (!is(class(object[[1]]), "nls"))
                            mgs <- validMsg(msg, "Not an nls object")
                        if (is.null(msg)) TRUE
                        else msg
                     })
##' @slot nullmodel A `nls` object describing the fitted null model
##' @slot alternative An `nlsList` object described the fitted alternative models
##' @slot vis A `ggplot` object visualising the fitted objects
##' @slot method A `character` string indicating the model fitted
##' @slot formula An instance of `formula` class
##' @md
##' @rdname hdxstat-class
.hdxstatmodel <- setClass("HdxStatModel",
                          slots = c(nullmodel = "nls",
                                    alternative = "nlsList",
                                    vis = "gg",
                                    method = "character",
                                    formula = "formula"),
                          validity = function(object){
                              msg <- validMsg(NULL, NULL)
                              sl <- sapply(object@alternative@nlsmodels, function(x) inherits(x, "nls"))
                              if (!all(sl))
                                  msg <- validMsg(msg, "Not all alternative models are nls")
                              if (is.null(msg)) TRUE
                              else msg
                          })
##' @slot statmodels A `list` of models. Each instance must a valid `HdxStatModel`.
##' @md
##' @rdname hdxstat-class
.hdxstatmodels <- setClass("HdxStatModels", 
                           slots = c(statmodels = "list"),
                           validity = function(object){
                               msg <- validMsg(NULL, NULL)
                               ms <- sapply(object@statmodels, function(x) inherits(x, "HdxStatModel"))
                               if (!all(ms))
                                   msg <- validMsg(msg, "Not all elements are valid HdxStatModels")
                               if (is.null(ms)) TRUE
                               else msg
                           })

##' @slot results A `DataFrame` contains summarised results from an hdx experiment.
##' @slot method The statistical method applied method applied. 
##' @md
##' @rdname hdxstat-class
.hdxstatres <- setClass("HdxStatRes", 
                         slots = c(results = "DataFrame",
                                   method = "character"),
                         validity = function(object){
                             msg <- validMsg(NULL, NULL)
                             ms <- inherits(object@results, "DataFrame")
                             if (!all(ms))
                                 msg <- validMsg(msg, "Results must be provided as a DataFrame")
                             if (is.null(ms)) TRUE
                             else msg
                           })
##' @md
##' @rdname hdxstat-class
setMethod("show", "HdxStatModel",
          function(object){
              cat("Object of class \"", class(object), "\"\n", sep = "")
              cat("Method:", object@method, "\n")
              cat("Fitted", length(object@alternative), "\n")
              invisible(NULL)
          })

##' @md
##' @rdname hdxstat-class
setMethod("show", "HdxStatModels",
          function(object){
              cat("Object of class \"", class(object), "\"\n", sep = "")
              cat("Number of models", length(object@statmodels), "\n")
              invisible(NULL)
          })

##' @md
##' @rdname hdxstat-class
setMethod("show", "HdxStatRes",
          function(object){
              cat("Object of class \"", class(object), "\"\n", sep = "")
              cat("Analysed using", object@method, "\n")
              invisible(NULL)
          })


##' @rdname hdxstat-class
setMethod("length", "nlsList", 
          function(x) length(x@nlsmodels))

##' @rdname hdxstat-class
setMethod("length", "HdxStatModel",
          function(x) length(x@alternative))

##' @rdname hdxstat-class
setMethod("length", "HdxStatModels",
          function(x) length(x@statmodels))

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname hdxstat-class
setMethod("[[", "nlsList",
          function(x, i, j = "missing", drop = "missing") x@nlsmodels[[i]])

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname hdxstat-class
setMethod("[", "nlsList",
          function(x, i, j = "missing", drop = "missing") {
              if (any(i > length(x)))
                  stop("Index out of bounds. Only ", length(x), " param(s) available.")
              x@nlsmodels <- x@nlsmodels[i]
              x
          })
