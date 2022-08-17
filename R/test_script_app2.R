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

  if (file_test("-f", data_path)) {
    if (file_ext(data_path) == "csv") {
      print("You gave me a CSV file for your HDX-MSM data")
      
      data <- read_csv(data_path)
      data.type <- "csv"
      }
    else if (file_ext(data_path) == "rsd") {
      print("You gave me a RSD file for your HDX-MSM data")
      print("I will assume your input data has alreayd been pre-processed")
      
      data <- readRDS(data_path)
      return(data)
    }
    else {
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
  if (data.type == "csv") {
    print("INFO: I will pre-process your data parse it using QFeatures ...")
    data <- preprocess_data(data, normalise = FALSE, save = FALSE)
    print("INFO: I pre-processed you input CSV data content and now it's available as a QFeatures instance")
    
    return(data)
  }
  
}

##' Visualise data


##' Apply functional analysis
analyse_data <- function(a, b){
  return()
}

##'  


# 1. Extract the data, and pre-process if necessary
#data_path <- "inst/extdata/Project_2_SecA_Cluster_Data.csv"
data_path <- "/homes/sanjuan/R/x86_64-pc-linux-gnu-library/4.1/hdxstats/extdata/MBP.csv"

data_filepath <- "vignettes/data/MBPqDF.rsd"
hdx_data <- extract_hdx_data(data_filepath)
results <- analyse_hdx_data(hdx_data)
graphics <- visualise_hdx_data(results, type="kinetics")
graphics <- visualise_hdx_data(results, type="manhatten")
graphics <- visualise_hdx_data(results, type="epitope", level="peptide")
graphics <- visualise_hdx_data(results, type="epitope", level="residue")
graphics <- visualise_hdx_data(results, type="protection", level="peptide")
graphics <- visualise_hdx_data(results, type="protection", level="residue")
graphics <- visualise_hdx_data(results, type="protection", level="residue", pdb="my_pdb_path")

# 2. Visualise the data

# Create heatmaps to visualise content and a few data points

# 3. Analyse data (statistic modelling)
