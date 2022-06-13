# Load essential libraries
library("NGLVieweR")
library("tidyverse")
library("RColorBrewer")
library("scales")
library("comprehenr")
library("dplyr")
library("bio3d")
library("psych") # to compute harmonic.mean
###################################################################
# MODIFY THIS FOR A DIFFERENT Ab
###################################################################
# Load input CSV files
averaged_protection <- read.csv("vignettes/data/averaged_protection_both.csv")
averaged_deprotection <- read.csv("vignettes/data/averaged_deprotection_both.csv")

# Select column name from input CSV datafiles
#column_name <- "dAb41_3"
#column_name <- "dAb25_1"
column_name <- "dAb6_1"
#column_name <- "dAb27_1"
###################################################################

# List to store data for visualisation
pdb_viewer_data <- list()

# Extract residue numbers from available residues in PDB coordinates
pdb_filepath <- "vignettes/data/5edv_chainA_clean_renumbered.pdb"
pdb_content <- read.pdb(pdb_filepath)

sequence_from_pdb <- pdbseq(pdb_content)
sequence_residue_numbers_pdb <- strtoi(row.names(data.frame(sequence_from_pdb)))
n_residues_from_pdb <- nchar(paste(sequence_from_pdb, collapse=""))

# Workout mean likelihood for protection values and available residue numbers in HDX data
dataset <- averaged_protection
X <- dataset$X # residue ranges
Y <- dataset[[column_name]] # protection likelihood

output_list = list()
for (i in 1:length(X)){
  # generate residue numbers sequence from range
  residue_number_range <- eval(parse(text=paste("seq", X[i], sep = "")))
  for (resn in residue_number_range){
    key <- paste("residue-", as.character(resn),sep = "")
    # gather all protection likelihood values per residue number
    output_list[[key]] <- append(output_list[[key]], Y[i])
  }
}

#geometric_mean_values <- lapply(output_list, function(x) exp(mean(log(x))))
geometric_mean_values <- lapply(output_list, function(x) harmonic.mean(x))
geometric_mean_values <- unname(unlist(geometric_mean_values)) # extract numerical values only

n_residues_hdx <- length(geometric_mean_values)
sequence_residue_numbers_hdx <- strtoi(sub("residue-", "", names(output_list)))

# Work out residues and protection values to map onto PDB
residue_numbers_hdx_pdb <- intersect(sequence_residue_numbers_hdx, sequence_residue_numbers_pdb)
geometric_mean_values_pdb <- geometric_mean_values[sequence_residue_numbers_hdx %in% sequence_residue_numbers_pdb]

# Save data into list
pdb_viewer_data[["protection"]] <- list("values"= geometric_mean_values_pdb,"residues"= residue_numbers_hdx_pdb)

# Workout mean likelihood for deprotection values and available residue numbers in HDX data
dataset <- averaged_deprotection
X <- dataset$X # residue ranges
Y <- dataset[[column_name]] # deprotection likelihood

output_list = list()
for (i in 1:length(X)){
  # generate residue numbers sequence from range
  residue_number_range <- eval(parse(text=paste("seq", X[i], sep = "")))
  for (resn in residue_number_range){
    key <- paste("residue-", as.character(resn),sep = "")
    # gather all protection likelihood values per residue number
    output_list[[key]] <- append(output_list[[key]], Y[i])
  }
}

geometric_mean_values <- lapply(output_list, function(x) exp(mean(log(x))))
geometric_mean_values <- unname(unlist(geometric_mean_values)) # extract numerical values only

n_residues_hdx <- length(geometric_mean_values)
sequence_residue_numbers_hdx <- strtoi(sub("residue-", "", names(output_list)))

# Work out residues and deprotection values to map onto PDB
residue_numbers_hdx_pdb <- intersect(sequence_residue_numbers_hdx, sequence_residue_numbers_pdb)
geometric_mean_values_pdb <- geometric_mean_values[sequence_residue_numbers_hdx %in% sequence_residue_numbers_pdb]

# Save data into list
pdb_viewer_data[["deprotection"]] <- list("values"= geometric_mean_values_pdb,"residues"= residue_numbers_hdx_pdb)

# Extract HDX values and assign protection/deprotection labels per residue
extracted_values <- c()
assigned_labels <- c()
for (i in 1:length(pdb_viewer_data$protection$residues)) {
  protection_value <- pdb_viewer_data$protection$values[i]
  deprotection_value <- pdb_viewer_data$deprotection$values[i]
  if (protection_value > deprotection_value) {
    assigned_labels <- append(assigned_labels, "protected")
    extracted_values <- append(extracted_values, protection_value)
  } else {
    assigned_labels <- append(assigned_labels, "deprotected")
    extracted_values <- append(extracted_values, deprotection_value)
  }
}

# Save data in dataframe and write as CSV
output_data <- data.frame("resnum"=pdb_viewer_data$protection$residues,"labels"=assigned_labels,"hdx_values"=extracted_values)
output_file <- paste("vignettes/data/hdx_", column_name, "_labelled_residues_pdb_5edv_chainA.csv", sep="")
write_csv(output_data, output_file)

# Work out colour parameters 
# Protection color function
residue_selections <- sprintf("%s", pdb_viewer_data$protection$residues)
values_per_residue <- pdb_viewer_data$protection$values
colormap <- colorRamp(brewer.pal(8, "Reds"))
domain <- c(min(values_per_residue), max(values_per_residue))
color_function_protected <- col_bin(colormap, domain = domain)

# Deprotection color function
residue_selections <- sprintf("%s", pdb_viewer_data$deprotection$residues)
values_per_residue <- pdb_viewer_data$deprotection$values
colormap <- colorRamp(brewer.pal(8, "Blues"))
domain <- c(min(values_per_residue), max(values_per_residue))
color_function_deprotected <- col_bin(colormap, domain = domain)

# Define general colour function
color_function <- function(extracted_values, assigned_labels) {
  colours <- c()
  for (i in 1:length(extracted_values)) {
    if (assigned_labels[i] == "protected") {
      colours <- append(colours, color_function_protected(extracted_values[i]))
    }
    else if (assigned_labels[i] == "deprotected") {
      colours <- append(colours, color_function_deprotected(extracted_values[i]))
    }
  }
  return(colours)
}

# Define selection according to label
residue_selections <- pdb_viewer_data$protection$residues
df <- data.frame(x=color_function(extracted_values, assigned_labels), y=residue_selections)
color_parameters <- to_list(for(i in 1:length(residue_selections)) c(df$x[i], df$y[i]))

# View mapped colour parameters onto PDB
NGLVieweR(pdb_filepath) %>%
  stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
  addRepresentation("cartoon") %>%
  addRepresentation("cartoon", param = list(color=color_parameters, backgroundColor="white"))
