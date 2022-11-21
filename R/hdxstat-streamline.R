##' Accessor for vis slot.
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param object of class "HdxStatModel"
##' @return "ggplot2" object
##' @md
##'
get_graphics_object <- function(object){ return(object@vis) }

##' Accessor for statmodels slot. Use: get_statmodels(results$fitted_models)
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param object of class "HdxStatModel"
##' @return "list" of objects of class "HdxStatModel"
##' @md
##'
get_statmodels <- function(object){ return(object@statmodels) }

##' Accessor for results slot. Use: get_ebayes_results(results$functional_analysis)
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param object of class "HdxStatModel"
##' @return "DataFrame" of results from Empirical Bayes functional analysis
##' @md
##'
get_ebayes_results <- function(object){ object@results }

##' Fit and apply Empirical Bayes (functional analysis) to an input QFeatures data object
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param data Input QFeatures data object or a selection 
##' @param peptide_selection Peptide string chain, e.g., 'DKELKAKGKSAL_4'
##' @param method Options: 
##'               "fit"  for a regular fitting
##'               "dfit" for differential fitting
##' @param starting_parameters To optimise for kinetic model based on 'formula'
##' @param formula Functional form or the fitting curve to optimise its parameters given 'starting_parameters'
##'                e.g., `formula = value ~ a * (1 - exp(-0.07*(timepoint))))`
##' @param maxAttempts The number of restarts attempts to get better model fits. Default is 5.               
##' @return Results from Empirical Bayes on fitted kinetic functionals or curves of Deuterium uptake
##' @md
##'  
analyse_kinetics <- function(data, 
                             peptide_selection = NULL,
                             method = NULL,
                             starting_parameters = list(a = NULL, b = 0.001,  d = NULL, p = 1),
                             formula = NULL,
                             maxAttempts = 5){
  
  if (method == "fit"){
    message("INFO: Performing fitting of Deuterium uptake kinetics. Method: 'hdxstats::fitUptakeKinetics' ")
    
    fitting_method <- fitUptakeKinetics
    tryCatch({
      fitted_models <- fitting_method(object = data,
                                      feature = peptide_selection,
                                      start = starting_parameters,
                                      formula = formula,
                                      maxAttempts = maxAttempts)
      }, error = function(e){ message(paste("ERROR:", conditionMessage(e), ". 'hdxstats::fitUptakeKinetics' returned some errors","\n")) })
    
    tryCatch({
      functional_analysis <- processFunctional(object = data,
                                               params = fitted_models)
      }, error = function(e){ message(paste("ERROR:", conditionMessage(e), ". 'hdxstats::processFunctional' returned some errors","\n")) })
    
    results <- list("fitted_models" = fitted_models, "functional_analysis" = functional_analysis, "method" = "fitUptakeKinetics")
  }
  
  if (method == "dfit"){
    message("INFO: Performing differential fitting of Deuterium uptake kinetics. Method: 'hdxstats::differentialUptakeKinetics' ")
    
    fitting_method <- differentialUptakeKinetics
    
    if (is.null(formula)){
      message("INFO: You did not specify a 'formula' for your fitting model.")
      message("INFO: Fitting will be performed for (default): 'formula <- value ~ a * (1 - exp(-b*(timepoint)^p)) + d' ")
      
      fitted_models <- fitting_method(object = data,
                                      feature = peptide_selection,
                                      start = starting_parameters)
      results <- list("fitted_models" = fitted_models, "method" = "differentialUptakeKinetics")
    }
    
    else{
      message("INFO: You specified your own 'formula' for your fitting model.")
      
      fitted_models <- fitting_method(object = data,
                                      feature = peptide_selection,
                                      start = starting_parameters,
                                      formula = formula)
      results <- list("fitted_models" = fitted_models, "method" = "differentialUptakeKinetics")
    }
    
  }
  
  return(results)
}

##' Extract data provided a CSV filepath
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param results List of functional analysis results 
##' @param data_selection QFeatures object selection
##' @param type Visualisation type for functional analysis results
##'             Options:
##'             'kinetics' Visualise fittings to Deu uptake kinetics data
##'             'forest' Visualise quality of fitting parameters for model
##'             'manhattan' Visualise ...
##'             'epitope' Visualise epitope map to see data overlap
##'             'protection' Visualise protection/deprotection plots
##' @param level Options:
##'             'peptide' applicable to 'type = epitope' 
##'             'residue' applicable to 'type = epitope' and 'type = protection'
##' @param fasta Path to FASTA file containing the sequence of protein of interest
##' @param pdb Path to PDB coordinate file of protein of interest
##' @return ggplot object(s) or NGLVieweR object (if pdb not NULL)
##' @md
##' 
visualise_hdx_data <- function(results,
                               data_selection = NULL,
                               type = NULL,
                               level = NULL,
                               fasta = NULL,
                               pdb = NULL){
  
  if (type == "kinetics" & results$method == "fitUptakeKinetics"){
    n_models <- length(get_statmodels(results$fitted_models))
    mssg <- paste("INFO: I found ", n_models, " models in your results data", sep = " ")
    message(mssg)
    message("INFO: You selected 'kinetics' to visualise from your results")
    
    graphics = list()
    for (i in seq_len(n_models)){graphics[[i]] = get_graphics_object(get_statmodels(results$fitted_models)[[i]])}
    return(graphics)
    message("INFO: I appended all ggplot output objects to a list")
  }
  
  else if (type == "kinetics" & results$method == "differentialUptakeKinetics"){
    n_models <- length(results$fitted_models)
    mssg <- paste("INFO: I found ", n_models, " models in your results data.", sep = " ")
    message(mssg)
    
    graphics = get_graphics_object(results$fitted_models)
    return(graphics)
    message("INFO: I appended all ggplot output objects to a list")
  }
  
  else if (type == "forest"){
    n_models <- length(get_statmodels(results$fitted_models))
    mssg <- paste("INFO: I found ", n_models, " models in your results data.", sep = " ")
    message(mssg)
    message("INFO: You selected 'forest' to visualise from your results")
    
    graphics = list()
    for (i in seq_len(n_models)){
      tryCatch({
        graphics[[i]] <- forestPlot(params = get_statmodels(results$fitted_models)[[i]])
      }, error = function(e){cat("ERROR:",conditionMessage(e), ". Failed for model number: ", i, "\n")})
    }
    return(graphics)
    message("INFO: I appended all 'forestPlot' output objects to a list")
  }
  
  else if (type == "manhattan"){
    n_cols <- length(as.vector(colnames(data_selection))$incoperation)
    
    mssg <- paste("INFO: I found",n_cols,"columns in your data selection. I will split your data selection into two and take their difference.")
    message(mssg)
    message("INFO: First half")
    message(paste(colnames(assay(data_selection)[,seq_len((n_cols/2))]), sep="\n"))
    message("INFO: Second half")
    message(paste(colnames(assay(data_selection)[,seq.int((1+n_cols/2), n_cols)]), sep="\n"))
    
    data_diff <- assay(data_selection)[,(1+n_cols/2):n_cols] - assay(data_selection)[,1:(n_cols/2)]
    successful_results<- rownames(data_selection)[["incoperation"]] %in% rownames(get_ebayes_results(results$functional_analysis))
    region_def <- as.data.frame(unique(rowData(data_selection)[["incoperation"]][, c("Start", "End")][successful_results,]))
    
    graphics <- list()
    for (i in seq_len((n_cols/2))){
      tryCatch({
        graphics[i] <- manhattanplot(params = results$functional_analysis,
                                     sequences = rownames(get_ebayes_results(results$functional_analysis)), 
                                     region = region_def,
                                     difference = data_diff[,i],
                                     nrow = 1)
      }, error = function(e){ message(paste("ERROR:",conditionMessage(e), ". Failed for model number: ", i)) })
    }
    
    mssg <- paste("INFO: You have ", length(graphics), "Manhattan plots")
    message(mssg)
    message("INFO: You selected 'manhattan' to visualise from your results")
    
    return(graphics)
    message("INFO: I appended all ggplot output objects to a list")
  }
  else if (type == "epitope"){
    
    message("INFO: You selected 'epitope' to visualise from your results")
    
    scores <- get_ebayes_results(results$functional_analysis)$ebayes.fdr
    peptide_charge_names <- rownames(get_ebayes_results(results$functional_analysis))
    peptide_sequences <- unlist(lapply(strsplit(peptide_charge_names, split="_"), function(x) head(x,n=1)))
    # NOTE: What about the charged states? Do they get usually ignored?
    
    if (!file.exists(fasta)){
      stop(paste("ERROR: FASTA file does not exist: ", fasta))
    }
    
    if (level == "peptide" & !is.null(fasta)){
      
      fasta_data <- readAAStringSet(filepath = fasta, "fasta")
      mssg <- paste("INFO: You input FASTA file contains", length(fasta_data), " sequences. I will take the first entry by default.")
      print(mssg)
      
      graphics <- plotEpitopeMap(AAString = fasta_data[[1]],
                                 peptideSeqs = peptide_sequences,
                                 numlines = 2,
                                 maxmismatch = 2,
                                 by = 1, # NOTE: What's the role of this arg that's never called in function?
                                 scores = 1 * (-log10(scores[unique(peptide_charge_names)])  > -log10(0.05)) + 0.0001,
                                 name = "significant")
      
      mssg <- paste("INFO: You have ", length(graphics), "parts for your Epitope map")
      message(mssg)
      return(graphics)
      message("INFO: I appended all ggplot output objects to a list")
      
    }else if (level == "residue" & !is.null(fasta) & is.null(pdb)){
      
      fasta_data <- readAAStringSet(filepath = fasta, "fasta")
      message <- paste("INFO: You input FASTA file contains", length(fasta_data), ". I will take the first entry by default.")
      print(message)
      
      graphics <- plotEpitopeMapResidue(AAString = fasta_data[[1]],
                                        peptideSeqs = peptide_sequences,
                                        numlines = 2,
                                        maxmismatch = 1,
                                        by = 5,
                                        scores = scores[unique(peptide_charge_names)],
                                        name = "-log10 p value")
      
      mssg <- paste("INFO: You have ", length(graphics), "parts for your Epitope map")
      message(mssg)
      return(graphics)
      message("INFO: I appended all ggplot output objects to a list")
      
    }else if (level == "residue" & !is.null(fasta) & !is.null(pdb)){
      
      fasta_data <- readAAStringSet(filepath = fasta, "fasta")
      message <- paste("INFO: You input FASTA file contains", length(fasta_data), ". I will take the first entry by default.")
      print(message)
      
      graphics_data <- ComputeAverageMap(AAString = fasta_data[[1]],
                                         peptideSeqs = unique(peptide_sequences),
                                         numlines = 2, maxmismatch = 1,
                                         by = 10, scores = scores[unique(peptide_charge_names)],
                                         name = "-log10 p value")
      
      mycolor_parameters <- hdx_to_pdb_colours(graphics_data, pdb)
      graphics <- NGLVieweR(pdb_filepath) %>%
        stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
        addRepresentation("cartoon") %>%
        addRepresentation("cartoon", param = list(color=mycolor_parameters, backgroundColor="white"))
      
      return(graphics)
    }
  }
  else if (type == "protection"){
    
    message("INFO: You selected 'protection' to visualise from your results")
    
    scores <- get_ebayes_results(results$functional_analysis)$ebayes.fdr
    peptide_charge_names <- rownames(get_ebayes_results(results$functional_analysis))
    peptide_sequences <- unlist(lapply(strsplit(peptide_charge_names, split="_"), function(x) head(x,n=1)))
    
    if (level == "residue" & !is.null(fasta) & is.null(pdb)){
      
      fasta_data <- readAAStringSet(filepath = fasta, "fasta")
      n_cols <- length(as.vector(colnames(data_selection))$incoperation)
      
      mssg <- paste("INFO: I found",n_cols,"columns in your data selection. I will split your data selection into two and take their difference.")
      message(mssg)
      print(colnames(assay(data_selection)[, seq_len(n_cols/2)]))
      print(colnames(assay(data_selection)[, seq.int((1+n_cols/2), n_cols)]))
      
      hdx_average <- ComputeAverageMap(AAString = fasta_data[[1]],
                                       peptideSeqs = unique(peptide_sequences),
                                       numlines = 2, 
                                       maxmismatch = 1,
                                       by = 10, 
                                       scores = scores[unique(peptide_charge_names)],
                                       name = "-log10 p value")
      hdx_diff <- list()
      graphics <- list()
      qDF <- data_selection
      for (i in 1:(n_cols/2)){
        hdx_diff[[i]] <- hdxdifference(object = data_selection,
                                       AAString = fasta_data[[1]],
                                       peptideSeqs = unique(peptide_sequences),
                                       numlines = 2,
                                       maxmismatch = 1,
                                       by = 10,
                                       scores = scores[unique(peptide_charge_names)],
                                       cols = c(i,(n_cols/2 + i)),
                                       name = "-log10 p value (signed)")
        
        graphics[[i]] <- hdxheatmap(averageMaps = list(hdx_average), diffMaps = list(hdx_diff[[i]]$diffMap)) 
      }
      
      mssg <- paste("INFO: You have ", length(graphics), "Protection/Deprotection heatmaps")
      message(mssg)
      return(graphics)
      message("INFO: I appended all output heatmaps to a list")
    }
    if (level == "residue" & !is.null(fasta) & !is.null(pdb)){
      
      fasta_data <- readAAStringSet(filepath = fasta, "fasta")
      n_cols <- length(as.vector(colnames(data_selection))$incoperation)
      
      message <- paste("INFO: I found",n_cols,"columns in your data selection. I will split your data selection into two and take their difference.")
      print(message)
      print(colnames(assay(data_selection)[, seq_len((n_cols/2))]))
      print(colnames(assay(data_selection)[, seq.int((1+n_cols/2),n_cols)]))
      
      hdx_average <- ComputeAverageMap(AAString = fasta_data[[1]],
                                       peptideSeqs = unique(peptide_sequences),
                                       numlines = 2, 
                                       maxmismatch = 1,
                                       by = 10, 
                                       scores = scores[unique(peptide_charge_names)],
                                       name = "-log10 p value")
      hdx_diff <- list()
      graphics <- list()
      qDF <- data_selection
      for (i in seq.int((n_cols/2))){
        hdx_diff[[i]] <- hdxdifference(object = data_selection,
                                       AAString = fasta_data[[1]],
                                       peptideSeqs = unique(peptide_sequences),
                                       numlines = 2,
                                       maxmismatch = 1,
                                       by = 10,
                                       scores = scores[unique(peptide_charge_names)],
                                       cols = c(i, (n_cols/2 + i)),
                                       name = "-log10 p value (signed)")
      }
      
      graphics <- list()
      for (i in seq.int((n_cols/2))){
        graphics_data <- hdx_average + sign(hdx_diff[[i]]$diffMap)
        mycolor_parameters <- hdx_to_pdb_colours(graphics_data, pdb, cmap_name="ProtDeprot")
        graphics[[i]] <- NGLVieweR(pdb_filepath) %>%
          stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
          addRepresentation("cartoon") %>%
          addRepresentation("cartoon", param = list(color=mycolor_parameters, backgroundColor="white"))
      }
      
      return(graphics)
    }
    
  }else {
    return(NULL)
    stop("FATAL: Specify what 'type' of visualiation you want. Available options: 'kinetics', 'forest', 'manhattan', 'epitope', 'protection' ")
  }
}
