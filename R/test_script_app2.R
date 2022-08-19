##' Pre-process data and output a QFeatures instance
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param data
##' @param normalise
##' @param save
##' @return preprocessed data
##' @md
##' 
preprocess_data <- function(data, 
                            normalise = TRUE,
                            save = TRUE) {
  # Print column names
  print("INFO: I found these columns in your input CSV file")
  data_columns <- colnames(data)
  message <- paste(colnames(data))
  print(message)
  
  # Transform data
  print("INFO: I reformatted your data to a wide format given your selected columns")
  
  columns_names <- c("hx_time", "replicate_cnt", "hx_sample") # to mash up
  columns_fixed <- c("pep_sequence", "pep_charge") # to keep fixed
  columns_values <- c("d")
  
  data_wide <- pivot_wider(data.frame(data),
                           values_from = columns_values,
                           names_from = columns_names,
                           id_cols = columns_fixed)
  
  # Remove NA's
  print("INFO: I removed NA values")
  
  data_wide <- data_wide[, colSums(is.na(data_wide)) != nrow(data_wide)]
  
  # Subtract columns from youir data
  print("INFO: I remove your selected columns from your output data")
  
  columns_to_remove <- c(1,2) # subtract names from first two columns
  old_columns_names <- colnames(data_wide)[-columns_to_remove]

  # add X and rep markers to new column names
  print("INFO: I reformatted your column labels")
  
  new_object.colnames <- paste0("X", old_columns_names)
  new_object.colnames <- gsub("0_", "0rep", new_object.colnames)
  new_object.colnames <- gsub("_", "cond", new_object.colnames)
  new_object.colnames <- gsub("%", "", new_object.colnames) 
  new_object.colnames <- gsub(" .*", "", new_object.colnames)
  
  # Parse data for selected columns
  print("INFO: I parsed your data as a qDF object class instance")
  
  initial_column <- 3
  last_column <- 102
  data_qDF <- parseDeutData(object = DataFrame(data_wide),
                            design = new_object.colnames,
                            quantcol = initial_column:last_column)
  
  # Normalise data 
  if (normalise) {
    print("INFO: Normalised data")
    
    data_qDF_normalised <- normalisehdx(data_qDF,
                                        sequence = unique(data$pep_sequence),
                                        method = "pc")
    data_qDF <- data_qDF_normalised
  }
  else{
    print("WARNING: Your output data is not normalised")
  }
  
  # Save data
  if (save) {
    print("INFO: Saved output data in")
    
    saveRDS(data_qDF, file='data/MBPqDF.rsd')
  } else{
    print("WARNING: Your output data was not saved")
  }
  
  return(data_qDF)
}

##' Extract data provided a filepath
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param data_path
##' @return extract data provided an input filepath
##' @md
##' 
extract_hdx_data <- function(data_path) {

  if (file_test("-f", data_path)){
    if (file_ext(data_path) == "csv") {
      print("You gave me a CSV file for your HDX-MSM data")
      
      data <- read_csv(data_path)
      data.type <- "csv"
      }
    else if (file_ext(data_path) == "rsd"){
      print("You gave me a RSD file for your HDX-MSM data")
      print("I will assume your input data has alreayd been pre-processed")
      
      data <- readRDS(data_path)
      return(data)
    }
    else{
      print("You provided an input file format that I cannot recognise")
      print("Provide a valid input. I will not continue execution.")
      return(NULL)
    }
  }
  else{
    print("This is not a valid path. Try again.")
    return(NULL)
  }
    
  # 2. Pre-process data
  if (data.type == "csv"){
    print("INFO: I will pre-process your data parse it using QFeatures ...")
    data <- preprocess_data(data, normalise = FALSE, save = FALSE)
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
      
      functional_analysis <- fitting_method(object = data,
                                feature = peptide_selection,
                                start = starting_parameters)
      results <- list("functional_analysis" = functional_analysis, "method" = "differentialUptakeKinetics")
    }
    
    else{
      print("INFO: You specified your own 'formula' for your fitting model.")
      
      functional_analysis <- fitting_method(object = data,
                                feature = peptide_selection,
                                start = starting_parameters,
                                formula = formula)
      results <- list("functional_analysis" = functional_analysis, "method" = "differentialUptakeKinetics")
    }
    
  }
  
  return(results)
}

visualise_hdx_data <- function(results,
                               type = NULL,
                               level = NULL,
                               pdb = NULL){
  
  if (type == "kinetics"){
    n_models <- length(results$fitted_models@statmodels)
    message <- paste("INFO: I found ", n_models, " models in your results data.", sep = " ")
    print(message)
    
    graphics = list()
    for (i in 1:n_models){graphics[[i]] = results$fitted_models@statmodels[[i]]@vis}
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

graphics <- visualise_hdx_data(results, type="kinetics") # <<<<--- NEXT
graphics <- visualise_hdx_data(results, type="forest")
graphics <- visualise_hdx_data(results, type="manhatten")
graphics <- visualise_hdx_data(results, type="epitope", level="peptide")
graphics <- visualise_hdx_data(results, type="epitope", level="residue") # Return heatmap
graphics <- visualise_hdx_data(results, type="protection", level="peptide") # Return heatmap
graphics <- visualise_hdx_data(results, type="protection", level="residue") # Return heatmap
graphics <- visualise_hdx_data(results, type="protection", level="residue", pdb="my_pdb_path")

# GUI: Make GUI by assembling these building blocks


# This condenses Tutorial 1
# Note that differentialUptake can only be applied to all data in relation to a single selected peptide
# If feature = all_peptides, this will procude a NULL output
# THi ALREADY PROCESS functionals by default
#all_peptides <- rownames(data_qDF)[[1]] # get all peptides
#peptide_selection <- all_peptides[37]