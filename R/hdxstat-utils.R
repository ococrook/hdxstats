##' Computes the number of exchangeable amides based on the sequnece
##' @title Compute exchangeable amides.
##' @param sequence The sequence of the peptide
##' @return Returns a numeric indicating the number of exchangeable amides
##' @md
##' 
##' @rdname hdxstats-utils
exchangeableAmides <- function(sequence) {
    
    n <- length(sequence)
    x <- vector(mode = "numeric", length = n)
    
    for(i in 1:n) {
        seq_vector <- strsplit(as.character(sequence[i]), split = "")[[1]]
        x[i] <- length(na.omit(sub("P", NA, seq_vector))) - 2
    }
    
    return(x)	
}
##' Normalise hdx data by putting on percentage incoperation scale or
##' using back-exchange correction
##' @title Normalise hdx data
##' @param object An object of class QFeatures
##' @param sequences Sequences corresponding to the hdxdata rows
##' @param method Either "pc" or "bc". Percentgae incoporation or 
##' back-exchange correction. Default is "pc"
##' @param correction A numeric vector indicated complete incoperation values 
##' for correction. Must be provided if method is "bc"
##' @return Returns a normalised 'QFeatures' object.
##' @md
##' 
##' @rdname hdxstats-utils
normalisehdx <- function(object,
                         sequences = NULL,
                         method = "pc",
                         correction = NULL){
    
    # checks
    stopifnot("Object is not an instance of QFeatures"=class(object) == "QFeatures")
    stopifnot("method is not one of pc or bc"=method %in% c("pc", "bc"))
    
    if (method == "bc"){
        stopifnot("correction must be numeric"=class(correction) == "numeric")
        stopifnot("correction must have compatible dimensions"=length(correction) == nrow(object))
    }
    
    if (method == "pc"){
        num_exch_sites <- exchangeableAmides(sequences)
        x <- t(vapply(1:nrow(assay(object)),
                            function(n) assay(object)[n,]/max(num_exch_sites[n], 1),
                            FUN.VALUE = numeric(ncol(assay(object)))))
    } else if (method == "bc"){
        
        x <- t(vapply(1:nrow(assay(object)),
                                function(n) assay(object)[n,]/correction[n],
                                FUN.VALUE = numeric(ncol(assay(object)))))
    }
    
    # parse as qFeatures object
    x <- DataFrame(x)
    x$rownames <- rownames(object)[[1]]
    qFeat <- readQFeatures(data.frame(x), ecol = 1:ncol(assay(object)), name = names(object), fnames = "rownames")
    
    return(qFeat)
}
