##' Color function
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param dataset The dataset for which numerical values will be colour enconded
##' @param scale_limits You can force a range of numerical values to be mapped only
##' @return Returns a color function
##' @md
##' 
##' @rdname pdb-visualisation
define_color_function <- function(dataset, 
                                  scale_limits = NULL,
                                  cmap_name = NULL){
  
    if (is.null(scale_limits)) {
      # By default work out the numerical domain from the data
      scale_limits <- c(min(dataset[!is.na(dataset)]), 
                        max(dataset[!is.na(dataset)]))
    }
  
    message(paste("INFO: Your scale limits are ", round(scale_limits, 2)))
  
    if (!is.null(cmap_name) && cmap_name == "ProtDeprot") {
        message("INFO: Negative values will be coloured in Blue, and positive ones on Red")
        
        colormap <- colorRamp(brewer.pal(8, "Blues"))
        domain <- c(0.0, abs(scale_limits[1]))
        color_function_deprotected <- col_bin(colormap, domain, na.color="#808080")
        
        colormap <- colorRamp(brewer.pal(8, "Reds"))
        domain <- c(0.0, abs(scale_limits[2]))
        color_function_protected <- col_bin(colormap, domain, na.color="#808080")
        
        output_function <- function(x){
          if (!is.na(x) && x >= 0.0){
            return(color_function_protected(abs(x)))
          }else{
            return(color_function_deprotected(abs(x)))
          }
        }
        
        warning("WARNING: NA values will be coloured in grey")
        
    }else if (is.null(cmap_name) || cmap_name == "viridis") {
        message("INFO: Your values will be coloured using Viridis")
      
        n_values <- length(unique(sort(dataset)))
        col_pal = c("white", viridis(n_values))
        output_function <- col_bin(col_pal, scale_limits, na.color="#808080")
        
        warning("WARNING: NA values will be coloured in grey")
    }
    
    return(output_function)
}

##' Map estimated (de)protection values to Blue-White-Red colormap values
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param dataset A list of estimated protection (positive) and deprotection (negative) values per residue
##' @param pdb_filepath The path to the PDB file for mapping of (de)protection values
##' @param scale_limits You can force a range of numerical values to be mapped only
##' @param cmap_name Specifies the name of the color map. Current options: "ProtDeprot" and "viridis"
##' @return Returns a list of colours and residue numbers to be inputted into NGLVieweR, 
##' with protection values in Red and deprotection values in Blue, and NA values in Grey.
##' @md
##' 
##' @rdname pdb-visualisation
hdx_to_pdb_colours <- function(dataset,
                               pdb_filepath,
                               scale_limits = NULL,
                               cmap_name = NULL){
  
    if (!file.exists(pdb_filepath)){
      stop(paste("ERROR: PDB filepath does not exist:", pdb_filepath))
    }
  
    # Extract residue numbers from available residues in PDB coordinates ---
    
    pdb_content <- read.pdb(pdb_filepath)
    sequence_from_pdb <- pdbseq(pdb_content)
    sequence_residue_numbers_pdb <- strtoi(row.names(data.frame(sequence_from_pdb)))
    
    # Report available data ---
    
    n_residues_from_pdb <- nchar(paste(sequence_from_pdb, collapse=""))
    message(paste("INFO: Your HDX input dataset has", length(colnames(dataset)), "entries"))
    message(paste("INFO: And excluding NA data you only have", length(dataset[!is.na(dataset)]), "entries"))
    message(paste("INFO: However, your input PDB has only", n_residues_from_pdb, "residues in total"))
    
    # Work out residues and (de)protection values that can be mapped onto PDB ---
    # Note: some residues are likely to be missing in the PDB
    
    sequence_residue_numbers_hdx <- seq(length(dataset)) # NEED TO TWEAK THIS!!!
    residue_numbers_hdx_pdb <- intersect(sequence_residue_numbers_hdx, sequence_residue_numbers_pdb)
    
    dataset_for_pdb_mapping <- dataset[residue_numbers_hdx_pdb]
    
    aa_sequence_from_pdb <- sequence_from_pdb[sequence_residue_numbers_pdb %in% residue_numbers_hdx_pdb]
    aa_sequence_from_pdb <- as.vector(aa_sequence_from_pdb)
    
    pdb_viewer_data <- list("values"= dataset_for_pdb_mapping,
                            "residues"= residue_numbers_hdx_pdb,
                            "aa"= aa_sequence_from_pdb)
    
    # Define colormap according to scale_limits--
    
    color_function <- define_color_function(dataset, scale_limits, cmap_name)
    
    residue_selections <- pdb_viewer_data$residues
    df <- data.frame(x=unlist(lapply(pdb_viewer_data$values, color_function)), y=residue_selections)
    color_parameters <- to_list(for(i in seq_len(length(residue_selections))) c(df$x[i], df$y[i]))

    return(color_parameters)
}
