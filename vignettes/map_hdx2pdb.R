#!/usr/bin/env Rscript

# Load essential libraries
library("shiny")
library("NGLVieweR")
library("tidyverse")
library("RColorBrewer")
library("scales")
library("comprehenr")
library("dplyr")
library("bio3d")
library("psych") # to compute harmonic.mean

#' Return harmonic mean protection/deprotection probability per residue
#' 
get_mean_values_per_residue <- function(X, Y) {
    output_list = list()
    for (i in 1:length(X)){
        # generate residue numbers sequence from range
        residue_number_range <- eval(parse(text=paste("seq", X[i], sep = "")))
        for (resn in residue_number_range){
            key <- paste("residue-", as.character(resn),sep = "")
            # gather all likelihood values per residue number
            output_list[[key]] <- append(output_list[[key]], Y[i])
        }
    }
    
    sequence_residue_numbers_hdx <- strtoi(sub("residue-", "", names(output_list)))
    
    #mean_values <- lapply(output_list, function(x) exp(mean(log(x))))
    mean_values <- lapply(output_list, function(x) harmonic.mean(x))
    mean_values <- unname(unlist(mean_values)) # extract numerical values only
    
    results <- list("sequence_residue_numbers_hdx"= sequence_residue_numbers_hdx, 
                    "mean_values"= mean_values
                    )
    return(results)
}

#' Label residues according to deuterium protection/deprotection probability
#' 
#' assigns a protection or deprotection label according to max probability
#' outputs a dataframe 
label_residues <- function(pdb_viewer_data) {
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
    
    output_labelled_data <- data.frame("resnum"=pdb_viewer_data$protection$residues,
                              "labels"=assigned_labels,
                              "hdx_values"=extracted_values
                              )
    return(output_labelled_data)
}

#' Colour residues according to protection/deprotection probability
#' 
#' given labelled mean probability values
#' protected values are coloured in red
#' deprotected values are coloured in blue
color_function <- function(output_labelled_data) {
    extracted_values <- output_labelled_data$hdx_values
    assigned_labels <- output_labelled_data$labels
    
    colormap <- colorRamp(brewer.pal(8, "Reds"))
    domain <- c(0.03, 1) # min and max values
    color_function_protected <- col_bin(colormap, domain)
    
    colormap <- colorRamp(brewer.pal(8, "Blues"))
    domain <- c(0.03, 1) # min and max values
    color_function_deprotected <- col_bin(colormap, domain)
    
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

# Load input CSV files
averaged_protection <- read.csv("data/averaged_protection_both.csv")
averaged_deprotection <- read.csv("data/averaged_deprotection_both.csv")

# Select column name from input CSV datafiles
antibody <- "dAb41_3"
#antibody <- "dAb25_1"
#antibody <- "dAb6_1"
#antibody <- "dAb27_1"

#pdb_viewer_data <- list() # List to store data for visualisation
pdb_filepath <- "data/5edv_chainA_clean_renumbered.pdb"
pdb_content <- read.pdb(pdb_filepath)

# Extract residue numbers from available residues in PDB coordinates
sequence_from_pdb <- pdbseq(pdb_content)
sequence_residue_numbers_pdb <- strtoi(row.names(data.frame(sequence_from_pdb)))
n_residues_from_pdb <- nchar(paste(sequence_from_pdb, collapse=""))

# Iterate over raw protection/deprotection data to generated labelled data
dataset <- list("protection"=averaged_protection, "deprotection"=averaged_deprotection)

pdb_viewer_data <- list() # List to store data for visualisation
for (x in names(dataset)){
    X <- dataset[[x]]$X # residue ranges
    Y <- dataset[[x]][[antibody]] # deprotection likelihood
    
    data <- get_mean_values_per_residue(X, Y)
    sequence_residue_numbers_hdx <- data$sequence_residue_numbers_hdx
    mean_values <- data$mean_values
    
    # Work out residues and deprotection values to map onto PDB
    residue_numbers_hdx_pdb <- intersect(sequence_residue_numbers_hdx, sequence_residue_numbers_pdb)
    mean_values_pdb <- mean_values[sequence_residue_numbers_hdx %in% sequence_residue_numbers_pdb]
    
    pdb_viewer_data[[x]] <- list("values"= mean_values_pdb,"residues"= residue_numbers_hdx_pdb)
}

# Save labelled data per residue in CSV file
output_labelled_data <- label_residues(pdb_viewer_data)
output_file <- paste("data/hdx_", antibody, "_labelled_residues_pdb_5edv_chainA.csv", sep="")
write_csv(output_labelled_data, output_file)

# Define selection according to label
residue_selections <- pdb_viewer_data$protection$residues
df <- data.frame(x=color_function(output_labelled_data), y=residue_selections)
color_parameters <- to_list(for(i in 1:length(residue_selections)) c(df$x[i], df$y[i]))

# View mapped colour parameters onto PDB
ui <- fluidPage(NGLVieweROutput("structure"))
server <- function(input, output) {
    output$structure <- renderNGLVieweR({
    
    NGLVieweR(pdb_filepath) %>%
    stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
    addRepresentation("cartoon") %>%
    addRepresentation("cartoon", param = list(color=color_parameters, backgroundColor="white")) %>%
    setQuality("high") %>%
    setFocus(0) #%>%
    #setSpin(TRUE)
    })
}
shinyApp(ui, server)