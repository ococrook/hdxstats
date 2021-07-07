# parse deuterium data into QFeature objects
parseDeutData <- function(object,
                          assayname = "incoperation",
                          rownames = NULL,
                          quantcol = NULL,
                          filter = NULL,
                          filterScore = 0.9,
                          sequence = "pep_sequence",
                          charge = "pep_charge",
                          design = NULL,
                          ecol = NULL){
    
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
