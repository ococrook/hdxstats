
make_parameter_file <- function(data, 
                                save = FALSE) {
  
  #Print column names
  print("INFO: I found these columns in your input CSV file")
  data_columns <- colnames(data)
  column_message <- paste(colnames(data))
  print(column_message)


  print("INFO: Specify the column name indicating the starting peptide residue numbers... OR, enter NA")
  Start <- readline(prompt = "Start (residue number) = ")
  while (is.null(data[[Start]]) & Start != "NA"){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    print(column_message)
    
    Start <- readline(prompt = "Start (residue number) = ")
  }
  
  print("INFO: Specify the column name indicating the ending peptide residue numbers... OR, enter NA")
  End <- readline(prompt = "End (residue number) = ")
  while (is.null(data[[End]]) & End != "NA"){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    print(column_message)
    
    End <- readline(prompt = "End (residue number) = ")
  }
  
  print("INFO: Specify the column name indicating the peptide sequences... OR, enter NA")
  Sequence <- readline(prompt = "Sequence (peptide) = " )
  while (is.null(data[[Sequence]]) & Sequence != "NA"){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    print(column_message)
    
    Sequence <- readline(prompt = "Sequence (peptide) = " )
  }
  
  print("INFO: Specify the column name indicating the peptide charge state... OR, enter NA")
  Charge <- readline(prompt = "Charge = ")
  while (is.null(data[[Charge]]) & Charge != "NA"){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    print(column_message)
    
    Charge <- readline(prompt = "Charge = ")
  }
  
  print("INFO: Specify the column name indicating the Deuterium uptake values ... OR, enter NA")
  Deu_Uptake <- readline(prompt = "Deu_Uptake = ")
  while (is.null(data[[Deu_Uptake]]) & Deu_Uptake != "NA"){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    print(column_message)
    
    Deu_Uptake <- readline(prompt = "Deu_Uptake = ")
  }
  
  print("INFO: Specify the column name indicating the Deuterium exposure timepoints... OR, enter NA")
  Exposure_Time <- readline(prompt = "Exposure_Time = ")
  while (is.null(data[[Exposure_Time]]) & Exposure_Time != "NA"){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    print(column_message)
    
    Exposure_Time <- readline(prompt = "Exposure_Time = ")
  }
  
  print("INFO: Specify column names indicating relevant experimental conditions ... OR, enter NA")
  print("INFO: IMPORTANT. You can provide more than one column name separared by commas (,) - I will merge them into a single label though.")
  Conditions <- readline(prompt = "Conditions = ")
  column_in_set <- unlist(strsplit(toString(gsub(" ", "", Conditions, fixed = TRUE)), split=",")) %in% data_columns
  while (!all(column_in_set) & Conditions != "NA"){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    print(column_message)
    
    Conditions <- readline(prompt = "Conditions = ")
  }
  
  print("INFO: Specify the column name indicating experimental replicates ... OR, enter NA")
  Replicate <- readline(prompt = "Replicate = ")
  while (is.null(data[[Replicate]]) & Replicate != "NA"){
    print("ERROR: Not a valid column name in your input CSV data. Try again.")
    print(column_message)
    
    Replicate <- readline(prompt = "Replicate = ")
  }
  
  print("INFO: OPTIONAL. Specify the columns you want to ignore. Otherwise, leave blank.")
  print("INFO: IMPORTANT. You can provide more than one column name separared by commas (,)")
  Ignore <- readline(prompt = "Ignore = ")
  
  print("INFO: OPTIONAL. Specify other column names you want to tag along - I will merge these into a single string chain. Otherwise, leave blank.")
  Other <- readline(prompt = "Other = ")
  
  print("INFO: Indicate whether I should convert your 'Exposure_Time' values. Options: TRUE or FALSE")
  convert_time <- readline(prompt = "convert_time = ")
  if (convert_time) {
    print("INFO: what are the original time units of your data? Available units: h (Hours), m (Minutes), s (Seconds).")
    
    original_time_units <- readline(prompt = "original_time_units = ")
    while (!original_time_units %in% c("s", "m", "h")) {
      print("ERROR: Not a valid time unit. Available units: h (Hours), m (Minutes), s (Seconds).")
      original_time_units <- readline(prompt = "original_time_units = ")
    }
  }
  else{
    convert_time <- FALSE
    original_time_units <- "NA"
  }
  
  parameters <- list("Start" = Start, 
                     "End" = End, 
                     "Sequence" = Sequence,
                     "Charge" = Charge,
                     "Deu_Uptake" = Deu_Uptake,
                     "Exposure_Time" = Exposure_Time,
                     "Conditions" = Conditions,
                     "Replicate" = Replicate,
                     "Ignore" = Ignore,
                     "Other" = Other,
                     "convert_time" = convert_time,
                     "original_time_units" = original_time_units)
  
  if (save != FALSE) {
    if (file.exists(dirname(save))) {
      saveRDS(parameters, file = save)
      
      message = paste("INFO: Saved your parameters in ", save)
      print(message)
    }
    
  } else{
    print("WARNING: Your parameters were not saved. Provide a valid output path with 'save = outfile_path'")
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
                            parameters = NULL,
                            parameter_file = NULL,
                            interactive = FALSE) {

  if (interactive == TRUE){
    print("INFO: You chose 'interactive' mode to parse the columns from your CSV content and define parameters to format your output QFeatures data object.")
    parameters <- make_parameter_file(data, save = FALSE)
    if (!is.null(parameter_file)){
      print("WARNING: You enabled 'interactive' as TRUE. This will override any 'parameter_file' or 'parameters' you provided.")
    }
  }
  
  if (!is.null(parameter_file)){
    print("INFO: You provided a 'parameter_file', I will extract parameters from this to format your output QFeatures data object.")
    if (file_test("-f", parameter_file)){
      parameters <- readRDS(parameter_file)
    }
  }
  
  if (!is.null(parameters)){
    if (is.list(parameters)){
      print("INFO: You provided a list of 'parameters', I will extract parameters from this to format your output QFeatures data object.")
    }
  }
  
  if (is.null(parameter_file) & is.null(parameters)) {
    print("ERROR: You either provided a invalid 'parameter_file' or list of 'parameters'. I will quit pre-processing.")
  }
  
  print("INFO: Stripped your 'Exposure_Time' values from non-numeric characters.")
  data[[parameters$Exposure_Time]] <- as.numeric(gsub("[^0-9.-]", "", data[[parameters$Exposure_Time]]))
  
  if (parameters$convert_time){
    if (parameters$original_time_units == 'h') {
      print("INFO: Your original_time_units == 'h'. I will convert your 'Exposure_Time' values to seconds (s).")
      data[[parameters$Exposure_Time]] <- 3600*data[[parameters$Exposure_Time]]
    }
    else if (parameters$original_time_units == 'm') {
      print("INFO: Your original_time_units == 'm'. I will convert your 'Exposure_Time' values to seconds (s).")
      data[[parameters$Exposure_Time]] <- 60*data[[parameters$Exposure_Time]]
    }
    else if (parameters$original_time_units == 's') {
      print("INFO: Your original_time_units == 's'. I will not convert your 'Exposure_Time' values.")
    }
  }
  
  if (parameters$Replicate == "NA"){
    print("INFO: Your 'Replicate' column appears to be NA. I will add this column with 1 values just to label your data.")
    data$Replicate <- 1
    parameters$Replicate <- "Replicate"
  }
  
  if (parameters$Charge == "NA"){
    print("INFO: Your 'Charge' column appears to be NA. I will add this column with 0 values just to label your data.")
    data$Charge <- 0
    parameters$Charge <- "Charge"
  }
  
  print("INFO: Reformatting your data to a wide format...")
  # Set default delimiters: X, rep, cond.
  delimiter.Exposure_Time <- "X" # <T>
  delimiter.Replicate <- "rep" # <R>
  delimiter.Conditions <- "cond" # <C>
  
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
                           names_sep = "<>")
  
  # Remove NA values
  data_wide <- data_wide[, colSums(is.na(data_wide)) != nrow(data_wide)]
  # Take all column names except 'columns_fixed'
  columns_to_remove <- 1:length(columns_fixed) # Remove columns_fixed
  old_columns_names <- colnames(data_wide)[-columns_to_remove]
  # Replace "_" with default delimiters and remove trailing strings
  new_object.colnames <- gsub(paste("<>", delimiter.Replicate, sep=""), delimiter.Replicate, old_columns_names)
  new_object.colnames <- gsub("<>", delimiter.Conditions, new_object.colnames)
  new_object.colnames <- gsub(" .*", "", new_object.colnames) # ??????
  
  # Parse data for selected columns
  print("INFO: Parsing your data as a qDF object class instance. Method: parseDeutData")
  
  initial_column <- length(columns_fixed)+1 # Fixed value by default
  last_column <- length(columns_fixed)+length(new_object.colnames) # Change to length value
  data_qDF <- parseDeutData(object = DataFrame(data_wide),
                            design = new_object.colnames,
                            quantcol = initial_column:last_column,
                            rownames = data_wide[[parameters$Sequence]],
                            sequence = parameters$Sequence,
                            charge = parameters$Charge)
  
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
    print("WARNING: Your output data was not saved. You can provide an output path with 'save = my_path'")
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
                             parameters = NULL,
                             parameter_file = NULL,
                             interactive = FALSE) {

  if (file.exists(data_path)){
    if (file_ext(data_path) == "csv") {
      print("INFO: You gave me a CSV file of your HDX-MSM data")
      
      data <- read_csv(data_path, show_col_types = FALSE)
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
                            parameters = parameters,
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
                               reference = NULL,
                               level = NULL,
                               fasta = NULL,
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

