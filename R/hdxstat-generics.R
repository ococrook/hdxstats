##' Generics for hdxstats package
##' @param object An object of class `HdxStatModel`
##' @param ... Additional arguments to `likRatio`
##' @rdname hidden_aliases
setGeneric("likRatio", function(object, ...)
    standardGeneric("likRatio"))
##' @rdname hidden_aliases
setGeneric("wilk", function(object, ...)
    standardGeneric("wilk"))
##' @rdname  hidden_aliases
setGeneric("fitUptakeKinetics", function(object, ...)
    standardGeneric("fitUptakeKinetics"))
##'@rdname hidden_aliases
setGeneric("logLik", function(object, ...)
    standardGeneric("logLik"))
setGeneric("residuals", function(object, ...)
    standardGeneric("residuals"))
setGeneric("summary", function(object, ...)
    standardGeneric("residuals"))    
setGeneric("vcov", function(object, ...)
    standardGeneric("residuals"))        