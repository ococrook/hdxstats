
make_parameter_file <- function(data, 
                                save = FALSE) {
  
  #Print column names
  print("INFO: I found these columns in your input CSV file")
  data_columns <- colnames(data)
  message <- paste(colnames(data))
  print(message)

  print("INFO: Specify the column name indicating the starting peptide residue numbers...")
  Start <- readline(prompt = "Start (residue number) = ")
  while (is.null(data[[Start]])){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    Start <- readline(prompt = "Start (residue number) = ")
  }
  
  print("INFO: Specify the column name indicating the ending peptide residue numbers...")
  End <- readline(prompt = "End (residue number) = ")
  while (is.null(data[[End]])){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    End <- readline(prompt = "End (residue number) = ")
  }
  
  print("INFO: Specify the column name indicating the peptide sequences...")
  Sequence <- readline(prompt = "Sequence (peptide) = " )
  while (is.null(data[[Sequence]])){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    Sequence <- readline(prompt = "Sequence (peptide) = " )
  }
  
  print("INFO: Specify the column name indicating the peptide charge state...")
  Charge <- readline(prompt = "Charge = ")
  while (is.null(data[[Charge]])){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    Charge <- readline(prompt = "Charge = ")
  }
  
  print("INFO: Specify the column name indicating the Deuterium uptake values ...")
  Deu_Uptake <- readline(prompt = "Deu_Uptake = ")
  while (is.null(data[[Deu_Uptake]])){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    Deu_Uptake <- readline(prompt = "Deu_Uptake = ")
  }
  
  print("INFO: Specify the column name indicating the Deuterium exposure timepoints...")
  Exposure_Time <- readline(prompt = "Exposure_Time = ")
  while (is.null(data[[Exposure_Time]])){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    Exposure_Time <- readline(prompt = "Exposure_Time = ")
  }
  
  print("INFO: Specify column names indicating relevant experimental conditions ...")
  print("INFO: IMPORTANT. You can provide more than one column name separared by commas (,) - I will merge them into a single label though.")
  Conditions <- readline(prompt = "Conditions = ")
  column_in_set <- unlist(strsplit(toString(gsub(" ", "", Conditions, fixed = TRUE)), split=",")) %in% data_columns
  while (!all(column_in_set)){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    Conditions <- readline(prompt = "Conditions = ")
  }
  
  print("INFO: Specify the column name indicating experimental replicates ...")
  Replicate <- readline(prompt = "Replicate = ")
  while (is.null(data[[Replicate]])){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    Replicate <- readline(prompt = "Replicate = ")
  }
  
  print("INFO: OPTIONAL. Specify the columns you want to ignore. Otherwise, leave blank.")
  print("INFO: IMPORTANT. You can provide more than one column name separared by commas (,)")
  Ignore <- readline(prompt = "Ignore = ")
  
  print("INFO: OPTIONAL. Specify other column names you want to tag along - I will merge these into a single string chain. Otherwise, leave blank.")
  Other <- readline(prompt = "Other = ")
  
  parameters <- list("Start" = Start, 
                     "End" = End, 
                     "Sequence" = Sequence,
                     "Charge" = Charge,
                     "Deu_Uptake" = Deu_Uptake,
                     "Exposure_Time" = Exposure_Time,
                     "Conditions" = Conditions,
                     "Replicate" = Replicate,
                     "Ignore" = Ignore,
                     "Other" = Other)
  
  if (file.exists(dirname(save))){
    saveRDS(parameters, file = save)
    message = paste("INFO: Saved your parameters in ", save)
    print(message)
    
  } else{
    print("ERROR: Your parameter_file was not saved. Provide a valid output path with 'save = outfile_path'")
    print(paste("The output file path that you provided was: ", save))
    }
  
  return(parameters)
}

##' Pre-process data and output a QFeatures instance
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param data
##' @param normalise
##' @param save
##' @return preprocessed data
##' @md
##' 
preprocess_data <- function(data, 
                            normalise = FALSE,
                            save = NULL,
                            parameter_file = NULL,
                            interactive = FALSE) {
  # Print column names
  print("INFO: I found these columns in your input CSV file...")
  data_columns <- colnames(data)
  message <- paste(colnames(data))
  print(message)

  if (interactive == TRUE){
    print("INFO: You chose 'interactive' mode to parse the columns from your CSV content and define parameters to format your output QFeatures data object.")
    parameters <- make_parameter_file(data)
    if (!is.null(parameter_file)){
      print("WARNING: You enabled 'interactive' as TRUE. This will override any 'parameter_file' you provided.")
    }
  }
  else if (!is.null(parameter_file) & file_test("-f", parameter_file)){
    print("INFO: You provided a 'parameter_file' form which I extract parameters to format your output QFeatures data object.")
    parameters <- readRDS(parameter_file)
    
  }
  else{
    print("ERROR: You either provided a invalid 'parameter_file'. I will quit pre-processing.")
    quit()
  }
  
  print("INFO: Reformatting your data to a wide format...")
  
  # Set default delimiters: X, rep, cond.
  delimiter.Exposure_Time <- "X"
  delimiter.Replicate <- "rep"
  delimiter.Conditions <- "cond"
  # Set column selections for pivot_wider 
  columns_names <- c(parameters$Exposure_Time, parameters$Replicate, parameters$Conditions)
  columns_fixed <- c(parameters$Sequence, parameters$Charge) 
  columns_values <- parameters$Deu_Uptake
  # Add delimiters to column entries
  data[[parameters$Replicate]] <- paste0(delimiter.Replicate, data[[parameters$Replicate]])
  data[[parameters$Exposure_Time]] <- paste0(delimiter.Exposure_Time, data[[parameters$Exposure_Time]])
  
  data_wide <- pivot_wider(data.frame(data),
                           values_from = columns_values,
                           names_from = all_of(columns_names),
                           id_cols = columns_fixed,
                           names_sep = "_")
  
  # Remove NA values
  data_wide <- data_wide[, colSums(is.na(data_wide)) != nrow(data_wide)]
  # Take all column names except 'columns_fixed'
  columns_to_remove <- 1:length(columns_fixed) # Remove columns_fixed
  old_columns_names <- colnames(data_wide)[-columns_to_remove]
  # Replace "_" with default delimiters and remove trailing strings
  new_object.colnames <- gsub(paste("_", delimiter.Replicate, sep=""), delimiter.Replicate, old_columns_names)
  new_object.colnames <- gsub("_", delimiter.Conditions, new_object.colnames)
  new_object.colnames <- gsub(" .*", "", new_object.colnames)
  
  # Parse data for selected columns
  print("INFO: Parsing your data as a qDF object class instance. Method: parseDeutData")
  
  initial_column <- length(columns_fixed)+1 # Fixed value by default
  last_column <- length(columns_fixed)+length(new_object.colnames) # Change to length value
  data_qDF <- parseDeutData(object = DataFrame(data_wide),
                            design = new_object.colnames,
                            quantcol = initial_column:last_column)
  
  # Normalise data 
  if (normalise) {
    print("INFO: Normalising data ... Method: normalisehdx")
    
    data_qDF <- normalisehdx(data_qDF,
                            sequence = unique(data[[parameters$Sequence]]),
                            method = "pc")
    
  }
  else{
    print("WARNING: Your output data is not normalised.")
  }
  
  # Save data
  if (!is.null(save)) {
    if (file.exists(dirname(save))) {
      
    saveRDS(data_qDF, file = save)
    print(paste("INFO: Saved output data in ", save))
    
    }
  } else {
    print("WARNING: Your output data was not saved.")
    print(paste("You provided the path ", save))
  }
  
  return(data_qDF)
}

##' Extract data provided a filepath
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param data_path
##' @return extract data provided an input filepath
##' @md
##' 
extract_hdx_data <- function(data_path,
                             normalise = FALSE,
                             save = NULL,
                             parameter_file = NULL,
                             interactive = FALSE) {

  if (file.exists(data_path)){
    if (file_ext(data_path) == "csv") {
      print("INFO: You gave me a CSV file of your HDX-MSM data")
      
      data <- read_csv(data_path)
      data.type <- "csv"
      }
    else if (file_ext(data_path) == "rsd"){
      print("INFO: You gave me a RSD file for your HDX-MSM data")
      print("INFO: I will assume your input data has alreayd been pre-processed")
      
      data <- readRDS(data_path)
      return(data)
    }
    else{
      print("ERROR: You provided an input file format that I cannot recognise")
      print("ERROR: Provide a valid input. I will provide a NULL output")
      return(NULL)
    }
  }
  else{
    print("ERROR: This is not a valid path. Try again.")
    return(NULL)
  }
    
  # 2. Pre-process data
  if (data.type == "csv"){
    print("INFO: I will pre-process your data parse it using QFeatures ...")
    data <- preprocess_data(data, 
                            normalise = normalise, 
                            save = save, 
                            parameter_file = parameter_file, 
                            interactive = interactive)
    
    print("INFO: I pre-processed you input CSV data content and now it's available as a QFeatures instance")
    
    return(data)
  }
  
}

##' Apply functional analysis to an input QFeatures data object
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param data Input QFeatures data object or a selection 
##' @param peptide_selection Peptide string chain, e.g., 'DKELKAKGKSAL_4'
##' @param method Options: 
##'               "fit"  for a regular fitting
##'               "dfit" for differential fitting
##' @param starting_parameters To optimise for kinetic model based on 'formula'
##' @param formula Functional form or the fitting curve to optimise its parameters given 'starting_parameters'
##'                e.g., 'formula = value ~ a * (1 - exp(-0.07*(timepoint))))'
##' @return Results from Empirical Bayes on fitted kinetic functionals or curves of Deuterium uptake
##' @md
##'  
analyse_kinetics <- function(data, 
                             peptide_selection = NULL,
                             method = NULL,
                             starting_parameters = list(a = NULL, b = 0.001,  d = NULL, p = 1),
                             formula = NULL){
  
  if (method == "fit"){
    print("INFO: Performing fitting of Deuterium uptake kinetics. Method: 'hdxstats::fitUptakeKinetics' ")
    
    fitting_method <- fitUptakeKinetics
    fitted_models <- fitting_method(object = data,
                                    feature = peptide_selection,
                                    start = starting_parameters)
    functional_analysis <- processFunctional(object = data,
                                             params = fitted_models)
    
    results <- list("fitted_models" = fitted_models, "functional_analysis" = functional_analysis, "method" = "fitUptakeKinetics")
  }

  else if (method == "dfit"){
    print("INFO: Performing differential fitting of Deuterium uptake kinetics. Method: 'hdxstats::differentialUptakeKinetics' ")
    
    fitting_method <- differentialUptakeKinetics
    
    if (is.null(formula)){
      print("INFO: You did not specify a 'formula' for your fitting model.")
      print("INFO: Fitting will be performed for (default): 'formula <- value ~ a * (1 - exp(-b*(timepoint)^p)) + d' ")
      
      fitted_models <- fitting_method(object = data,
                                      feature = peptide_selection,
                                      start = starting_parameters)
      results <- list("fitted_models" = fitted_models, "method" = "differentialUptakeKinetics")
    }
    
    else{
      print("INFO: You specified your own 'formula' for your fitting model.")
      
      fitted_models <- fitting_method(object = data,
                                      feature = peptide_selection,
                                      start = starting_parameters,
                                      formula = formula)
      results <- list("fitted_models" = fitted_models, "method" = "differentialUptakeKinetics")
    }
    
  }
  
  return(results)
}

visualise_hdx_data <- function(results,
                               type = NULL,
                               level = NULL,
                               pdb = NULL){
  
  if (type == "kinetics" & results$method == "fitUptakeKinetics"){
    n_models <- length(results$fitted_models@statmodels)
    message <- paste("INFO: I found ", n_models, " models in your results data.", sep = " ")
    print(message)
    
    graphics = list()
    for (i in 1:n_models){graphics[[i]] = results$fitted_models@statmodels[[i]]@vis}
    return(graphics)
    print("INFO: I appended all ggplot output objects to a list")
  }
  
  else if (type == "kinetics" & results$method == "differentialUptakeKinetics"){
    n_models <- length(results$fitted_models)
    message <- paste("INFO: I found ", n_models, " models in your results data.", sep = " ")
    print(message)
    
    graphics = results$fitted_models@vis
    return(graphics)
    print("INFO: I appended all ggplot output objects to a list")
  }
  
  else if (type == "forest"){
    n_models <- length(results$fitted_models@statmodels)
    message <- paste("INFO: I found ", n_models, " models in your results data.", sep = " ")
    print(message)
    
    graphics = list()
    for (i in 1:n_models){graphics[[i]] = forestPlot(params = results$fitted_models@statmodels[[i]])}
    return(graphics)
    print("INFO: I appended all 'forestPlot' output objects to a list")
  }
}
  # 
  # else if (type == "manhatten"){
  #   return(NULL)
  # }
  # 
  # else if (type == "peptide"){
  #   return(NULL)
  # }
  # 
  # else if (type == "epitope"){
  #   return(NULL)
  # }
  # 
  # else if (type == "protection"){
  #   return(NULL)
  # }
  # 
  # else {
  #   print("FATAL: Specify what 'type' of visualiation you want. Available options: 'kinetics', 'forest', 'manhatten', 'peptide', 'epitope', 'protection' ")
  #   return(NULL)
  # }
}

graphics <- visualise_hdx_data(results, type="kinetics") # READY
graphics <- visualise_hdx_data(results, type="forest") # READY 
graphics <- visualise_hdx_data(results, type="manhatten")# <<<<--- NEXT
graphics <- visualise_hdx_data(results, type="peptide")
graphics <- visualise_hdx_data(results, type="epitope", level="peptide")
graphics <- visualise_hdx_data(results, type="epitope", level="residue") # Return heatmap
graphics <- visualise_hdx_data(results, type="protection", level="peptide") # Return heatmap
graphics <- visualise_hdx_data(results, type="protection", level="residue") # Return heatmap
graphics <- visualise_hdx_data(results, type="protection", level="residue", pdb="my_pdb_path")


########################################################################
# TUTORIAL 1 & 2
#data_filepath <- "/homes/sanjuan/R/x86_64-pc-linux-gnu-library/4.1/hdxstats/extdata/MBP.csv"
#data_filepath <- system.file("extdata", csv_filename, package = "hdxstats")
#data_filepath <- "vignettes/data/MBPqDF.rsd" # or, CSV file

# TUTORIAL 3 & 4 (HOIP assays)
#data_filepath <- system.file("extdata", "N64184_1a2_state.csv", package = "hdxstats")

# TUTORIAL secA-casestudy
#data_path <- "inst/extdata/Project_2_SecA_Cluster_Data.csv"

# TUTORIAL flexible-fits.Rmd
#MBPpath <- system.file("extdata", "MBP.csv", package = "hdxstats")

#######################################################
data_filepath <- "vignettes/data/MBPqDF.rsd"
hdx_data <- extract_hdx_data(data_filepath) # DONE


# TEST 1
# INPUT
data_selection <- hdx_data[,1:24]
all_peptides <- rownames(data_selection)[[1]]
starting_parameters <- list(a = NULL, b = 0.001,  d = NULL, p = 1)
# OUTPUT
results <- analyse_kinetics(data = data_selection, 
                            method = "fit", 
                            peptide_selection = all_peptides, 
                            start = starting_parameters)
# TEST 2
# INPUT
data_selection <- hdx_data[,61:100]
all_peptides <- rownames(data_selection)[[1]]
starting_parameters <- list(a = NULL, b = 0.001,  d = NULL, p = 1)
# OUTPUT
results <- analyse_kinetics(data = data_selection, 
                            method = "fit", 
                            peptide_selection = all_peptides, 
                            start = starting_parameters)
# TEST 3
# INPUT
data_selection <- hdx_data[,1:100]
all_peptides <- rownames(data_selection)[[1]] # get all peptides
starting_parameters <- list(a = NULL, b = 0.0001,  d = NULL, p = 1)
# OUTPUT
results <- analyse_kinetics(data = data_selection, 
                            method = "dfit", 
                            peptide_selection = all_peptides[37], 
                            start = starting_parameters)

#######################################################
graphics <- visualise_hdx_data(results, type="kinetics") # This applies to both regular and differential fits
graphics <- visualise_hdx_data(results, type="forest") # This only applies to regular fits
graphics <- visualise_hdx_data(results, type="manhatten") # This only applies to differential fits
graphics <- visualise_hdx_data(results, type="epitope", level="peptide")
graphics <- visualise_hdx_data(results, type="epitope", level="residue") # Return heatmap
graphics <- visualise_hdx_data(results, type="protection", level="peptide") # Return heatmap
graphics <- visualise_hdx_data(results, type="protection", level="residue") # Return heatmap
graphics <- visualise_hdx_data(results, type="protection", level="residue", pdb="my_pdb_path")

# GUI: Make GUI by assembling these building blocks
