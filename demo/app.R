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
library("data.table") # to transpose table
library("shinydashboard")
library("DT")

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
  
  #resn_min <- min(sequence_residue_numbers_hdx)
  #resn_max <- max(sequence_residue_numbers_hdx)
  #sequence_residue_numbers_hdx_full <- seq(resn_min, resn_max)
  
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
                                     "aa"=pdb_viewer_data$protection$aa,
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
color_function <- function(output_labelled_data, scale_limits = c(0.03, 1, 0.03, 1)) {
  extracted_values <- output_labelled_data$hdx_values
  assigned_labels <- output_labelled_data$labels
  
  colormap <- colorRamp(brewer.pal(8, "Reds"))
  domain <- c(scale_limits[1], scale_limits[2]) # min and max values
  color_function_protected <- col_bin(colormap, domain)
  
  colormap <- colorRamp(brewer.pal(8, "Blues"))
  domain <- c(scale_limits[3], scale_limits[4]) # min and max values
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

map_hdx2pdb <- function(dataset, antibody, pdb_filepath, scale_limits = c(0.03, 1, 0.03, 1)) {
  pdb_content <- read.pdb(pdb_filepath)
  
  # Extract residue numbers from available residues in PDB coordinates
  sequence_from_pdb <- pdbseq(pdb_content)
  sequence_residue_numbers_pdb <- strtoi(row.names(data.frame(sequence_from_pdb)))
  n_residues_from_pdb <- nchar(paste(sequence_from_pdb, collapse=""))
  
  # Iterate over raw protection/deprotection data to generated labelled data
  #dataset <- list("protection"=averaged_protection, "deprotection"=averaged_deprotection)
  
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
    
    aa_sequence_from_pdb <- sequence_from_pdb[sequence_residue_numbers_pdb %in% residue_numbers_hdx_pdb]
    aa_sequence_from_pdb <- as.vector(aa_sequence_from_pdb)
    #residues_pdb_hdx_diff <- setdiff(sequence_residue_numbers_pdb, sequence_residue_numbers_hdx)
    
    # NOTE: We only visualise those residues for which HDX data and PDB coords are available
    pdb_viewer_data[[x]] <- list("values"= mean_values_pdb,"residues"= residue_numbers_hdx_pdb, "aa"= aa_sequence_from_pdb)
  }
  
  # Save labelled data per residue in CSV file
  output_labelled_data <- label_residues(pdb_viewer_data)
  #output_file <- paste("data/hdx_", antibody, "_labelled_residues_pdb_5edv_chainA.csv", sep="")
  #write_csv(output_labelled_data, output_file)
  
  # Define selection according to label
  residue_selections <- pdb_viewer_data$protection$residues
  df <- data.frame(x=color_function(output_labelled_data, scale_limits = scale_limits), y=residue_selections)
  color_parameters <- to_list(for(i in 1:length(residue_selections)) c(df$x[i], df$y[i]))
  
  return(color_parameters)
}
#################################################
# Load input CSV files
averaged_protection <- read.csv("data/averaged_protection_both.csv")
averaged_deprotection <- read.csv("data/averaged_deprotection_both.csv")
dataset <- list("protection"=averaged_protection, "deprotection"=averaged_deprotection)

# List of Antibodies present in HDX data
antibody_selection <- colnames(averaged_protection)[c(2:length(averaged_protection))]

# Antigen PDB coordinate file
pdb_filepath <- "data/5edv_chainA_clean_renumbered.pdb"
#################################################

# Define shinny app
ui <- fluidPage(
  titlePanel("HDX protection/deprotection viewer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("antibody", "Antibody", antibody_selection),
      selectInput("representation", "Representation", c("cartoon", "backbone", "licorice", "ball+stick", "spacefill", "surface")),
      sliderInput("prot_range", "Protection colour-scale limits",min = 0, max = 100, value = c(0,100)),
      sliderInput("deprot_range", "Deprotection colour-scale limits",min = 0, max = 100, value = c(0,100)),
    ),
    
    mainPanel(NGLVieweROutput("structure"))
  ),

  dashboardBody(
    fluidRow(
      column(width = 12,
             box(title = "Sequence Viewer",
                 width = 12,
                 status = "primary",
                 div(style = 'overflow-x: scroll',
                     #tableOutput('trace_table')
                     DT::dataTableOutput('trace_table')
                     )
                 )
             )
      )
  ),
)


server <- function(input, output) {
  
  output$structure <- renderNGLVieweR({
    
  scale_limits <- c(input$prot_range[1], input$prot_range[2], input$deprot_range[1], input$deprot_range[2])/100
  mycolor_parameters <- map_hdx2pdb(dataset, input$antibody, pdb_filepath, scale_limits)
  
  #colour antigen PDB according to protection/deprotection probability
  NGLVieweR(pdb_filepath) %>%
      stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
      addRepresentation(input$representation) %>%
      addRepresentation(input$representation, param = list(name="sele1", color=mycolor_parameters, backgroundColor="yellow")) #%>%
      #setQuality("high")
    
    })
  
  #render protection/deprotection probability for chosen antibody
  observeEvent(input$antibody, {
    # load precomputed protection/deprotection probabilities
    filename_csv <- paste("data/hdx_", input$antibody, "_labelled_residues_pdb_5edv_chainA.csv", sep="")
    df <- data.frame(read_csv(filename_csv))
    df[,'hdx_values'] = round(df[,'hdx_values'],2) # round probability values to 2 significantd digits
    cut_vals <- df$resnum
    
    df <- transpose(subset(df, select = -c(labels))) # ignore labels column and transpose dataframe
    scale_limits <- c(input$prot_range[1], input$prot_range[2], input$deprot_range[1], input$deprot_range[2])/100
    mycolor_parameters <- map_hdx2pdb(dataset, input$antibody, pdb_filepath, scale_limits)
    
    output$trace_table <- DT::renderDataTable({
      datatable(df,
                #fontWeight =  c(5, 'normal', 'bold'),
                options = list(dom='t', ordering=F), # no search, no sorting
                colnames = rep("", ncol(df)), # no column names
                rownames = c('<b>resnumber</b>','<b>aa</b>','<b>probability</b>'), # rename rows
                escape = FALSE
                ) %>% formatStyle(columns=names(df),
                                  backgroundColor=styleEqual(cut_vals, sapply(mycolor_parameters,"[[",1))
        ) 
    })
  })
}

shinyApp(ui, server)