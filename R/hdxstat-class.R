##' @slot nlsmodels A `list()` containing nls models results in `nls` instances.
##' Each element must be a valid `nls` instance.
##' @md
##' @rdname hdxstat-class
.nlsList <- setClass("nlsList",
                     slots = c(nlsmodels = "list"),
                     validity = function(object){
                        msg <- validMsg(NULL, NULL)
                        if (class(object[[1]]) != "nls")
                            mgs <- validMsg(msg, "Not an nls object")
                        if (is.null(msg)) TRUE
                        else msg
                     })
##' @slot nullmodel A `nls` object describing the fitted null model
##' @slot alternative An `nlsList` object described the fitted alternative models
##' @slot vis A `ggplot` object visualising the fitted objects
##' @slot method A `character` string indicating the model fitted
##' @md
##' @rdname hdxstat-class
.hdxstatmodel <- setClass("HdxStatModel",
                          slots = c(nullmodel = "nls",
                                    alternative = "nlsList",
                                    vis = "ggplot",
                                    method = "character"),
                          validity = function(object){
                              msg <- validMsg(NULL, NULL)
                              sl <- sapply(object@nlsList, function(x) inherits(x, "nls"))
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


##' @md
##' @rdname hdxstat-class
setMethod("show", "HdxStatModel",
          function(object){
              cat("Object of class \"", class(object), "\"\n", sep = "")
              cat("Method:", object@method, "\n")
              cat("Fitted", length(objective@alternative), "\n")
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

