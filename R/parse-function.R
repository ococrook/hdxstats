##' Parse deuterium data into QFeature objects
##' @param object An object of class DataFrame
##' @param assayname The name of assay. The Default is incoperation.
##' @param rownames The rownames of the features. Default is NULL, in which case
##' the rownames are extract from the `sequence` and `charge` in the object
##' @param quantcol The columns wiht the quantitative data.
##' @param filter Currently unused filter
##' @param filterScore Currently usused filter score
##' @param sequence The name of the column where the peptide sequence is stored.
##' Default is "pep_sequence"
##' @param charge The name of the column where the peptide charge is stored.
##' Default is "pep_chrage"
##' @param  design The design which will become the column names
##' @return An instance of class `QFeatures` storing the quantitative mass-spectrometry
##' data
##' @md
##' 
##' @rdname parse-function
parseDeutData <- function(object,
                          assayname = "incoperation",
                          rownames = NULL,
                          quantcol = NULL,
                          filter = NULL,
                          filterScore = 0.9,
                          sequence = "pep_sequence",
                          charge = "pep_charge",
                          design = NULL){
    
    stopifnot("object must be of class DataFrame"= class(object) == "DFrame")
    stopifnot("You must provide a design"= is.null(design) == FALSE)

    # if rownames are null use sequence and charge 
    if(is.null(rownames)){
        row.names <- paste0(object[, sequence], "_", object[, charge])
        rownames(object) <- row.names
    }else {
        row.names <- rownames
        rownames(object) <- row.names
    }
    
    
    # generate storage
    qDF <- DataFrame(rownames = row.names, matrix(NA, ncol = length(design), nrow = nrow(object)))
    colnames(qDF)[-1] <- design
    rownames(qDF) <- qDF$rownames
    
    ## add quantitative data
    qDF[rownames(qDF), -1] <- object[rownames(qDF), quantcol]
    
    # parse as qFeatures object
    qFeat <- readQFeatures(data.frame(qDF), ecol = 2:ncol(qDF), name = assayname, fnames = "rownames")
    
    return(qFeat)   
}
