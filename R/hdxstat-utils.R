##' Computes the number of exchangeable amides based on the sequnece
##' @title Compute exchangeable amides.
##' @param sequence The sequence of the peptide
##' @return Returns a numeric indicating the number of exchangeable amides
##' @md
##' 
##' @rdname hdxstats-utils
exchangeableAmides <- function(sequence) {
    
    n <- length(sequence)
    x <- vector(mode = "numeric", length = n)
    
    for(i in seq.int(n)) {
        seq_vector <- strsplit(as.character(sequence[i]), split = "")[[1]]
        x[i] <- length(na.omit(sub("P", NA, seq_vector[-seq.int(2)])))
    }
    
    return(x)	
}
##' Normalise hdx data by putting on percentage incoperation scale or
##' using back-exchange correction
##' @title Normalise hdx data
##' @param object An object of class QFeatures
##' @param sequences Sequences corresponding to the hdxdata rows
##' @param method Either "pc" or "bc". Percentgae incoporation or 
##' back-exchange correction. Default is "pc"
##' @param correction A numeric vector indicated complete incoperation values 
##' for correction. Must be provided if method is "bc"
##' @return Returns a normalised 'QFeatures' object.
##' @md
##' 
##' @rdname hdxstats-utils
normalisehdx <- function(object,
                         sequences = NULL,
                         method = "pc",
                         correction = NULL){
    
    # checks
    stopifnot("Object is not an instance of QFeatures"=class(object) == "QFeatures")
    stopifnot("method is not one of pc or bc"=method %in% c("pc", "bc", "undeuterated", "intercept"))
    
    if (method == "bc"){
        stopifnot("correction must be numeric"=class(correction) == "numeric")
        stopifnot("correction must have compatible dimensions"=length(correction) == nrow(object))
    }
    
    if (method == "pc"){
        num_exch_sites <- exchangeableAmides(sequences)
        x <- t(vapply(seq.int(nrow(assay(object))),
                      function(n) assay(object)[n,]/max(num_exch_sites[n], 1),
                      FUN.VALUE = numeric(ncol(assay(object)))))
        
        # parse as qFeatures object
        x <- DataFrame(x)
        x$rownames <- rownames(object)[[1]]
        qFeat <- readQFeatures(data.frame(x), ecol = 1:ncol(assay(object)), name = names(object), fnames = "rownames")
        rowData(qFeat)[["incoperation"]] <- rowData(object)[["incoperation"]]

        return(qFeat)
    }
    
    if (method == "undeuterated"){
        # Subtract row minimum value, regardless of replicate and condition
        x <- DataFrame(data.frame(assay(object) - apply(assay(object), 1, function(x) min(x, na.rm = TRUE))))
        
        # parse as qFeatures object
        x <- DataFrame(x)
        x$rownames <- rownames(object)[[1]]
        qFeat <- readQFeatures(data.frame(x), ecol = 1:ncol(assay(object)), name = names(object), fnames = "rownames")
        rowData(qFeat)[["incoperation"]] <- rowData(object)[["incoperation"]]

        return(qFeat)
        
    }
    
    if (method == "intercept"){
        mssg <-paste(" You have", length(rownames(assay(object))), "peptide-charge paired values")
        rlog::log_info(mssg)
        
        ldf_new <- rbind()
        for (peptide_charge in rownames(assay(object))){
            peptide_charge_data <- as.data.frame(assay(object))[peptide_charge, ]
            peptide_charge_data <- longFormat(peptide_charge_data)
            peptide_charge_data$condition <- as.factor(str_match(peptide_charge_data$colname, "cond\\s*(.*)")[, 2])
            
            Deu_min_global <- apply(assay(object), 1, function(x) min(x, na.rm = TRUE))[[peptide_charge]]
            
            peptide_charge_conditions <- unique(peptide_charge_data$condition)
            mssg <- paste(" For ", peptide_charge, ", you have", length(peptide_charge_conditions), "conditions")
            rlog::log_info(mssg)
            
            #ldf_new <- rbind()
            for (state in peptide_charge_conditions){
                ldf <- peptide_charge_data %>% subset(condition == state)
                ldf$timepoint <- as.numeric(str_match(ldf$colname, "X\\s*(.*?)\\s*rep")[, 2])
                ldf$replicates <- as.factor(str_match(ldf$colname, "rep\\s*(.*)\\s*cond")[, 2])
                
                ldf$replicates <- unlist(lapply(strsplit(as.vector(as.factor(str_match(ldf$colname, "rep\\s*(.*)\\s*cond")[, 2])), split="_"), function(x) tail(x, n=1)))
                
                # Subtract Deu uptake value at 0 timepoint
                mssg <- paste(" You have", length(unique(ldf$replicates)), "replicates, for", state)
                rlog::log_info(mssg)
                
                for (n_replicate in unique(ldf$replicates)){
                    single_replicate_data <- ldf %>% subset(replicates == n_replicate)
                    x <- single_replicate_data %>% subset(timepoint == 0)
                    
                    if (all(is.na(x$value))){
                        mssg <- " All Deu uptake values for the zero timepoint are NA. I will take the minimum across all conditions."
                        rlog::log_info(mssg)
                        single_replicate_data$value <- single_replicate_data$value - Deu_min_global
                        ldf_new <- rbind(ldf_new, single_replicate_data)
                        
                    }else{
                        mssg <- " At least one Deu uptake values for the zero timepoint is NA. I will take the minimum of all zero timepoints"
                        rlog::log_info(mssg)
                        Deu_min <- min(x$value, na.rm = TRUE)
                        single_replicate_data$value <- single_replicate_data$value - Deu_min
                        ldf_new <- rbind(ldf_new, single_replicate_data)
                        
                    }
                }
            }
        }
        
        x <- DataFrame(ldf_new)
        x_wide <- pivot_wider(data.frame(x), values_from = value, id_cols = rowname, names_from = colname)
        x_wide_df <- DataFrame(x_wide)
        x_wide_df$rownames <- x_wide$rowname
        qFeat <- readQFeatures(data.frame(x_wide_df), ecol = 2:ncol(assay(object)), name = names(object), fnames = "rownames")
        rowData(qFeat)[["incoperation"]] <- rowData(object)[["incoperation"]]

        return(qFeat)
        
    } else if (method == "bc"){
        
        x <- t(vapply(seq.int(nrow(assay(object))),
                      function(n) assay(object)[n, ]/correction[n],
                      FUN.VALUE = numeric(ncol(assay(object)))))
        
        # parse as qFeatures object
        x <- DataFrame(x)
        x$rownames <- rownames(object)[[1]]
        qFeat <- readQFeatures(data.frame(x), ecol = 1:ncol(assay(object)), name = names(object), fnames = "rownames")
        rowData(qFeat)[["incoperation"]] <- rowData(object)[["incoperation"]]

        return(qFeat)
    }
}

##' Function to define parameters corresponding to CSV column names to feed into data curation pipeline
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param data Data content from read from CSV file
##' @param save_parameters Output path for list with parameters for data curation pipeline
##' @return list of parameters to CSV column names
##' @md
##'  
##' @rdname hdxstat-utils
make_parameter_file <- function(data, 
                                save_parameters = FALSE) {
    
    #Print column names
    rlog::log_info(" I found these columns in your input CSV file")
    data_columns <- colnames(data)
    column_message <- paste(colnames(data))
    rlog::log_info(column_message)
    
    
    print(" Specify the column name indicating the starting peptide residue numbers... OR, enter NA")
    Start <- readline(prompt = "Start (residue number) = ")
    while (is.null(data[[Start]]) & Start != "NA"){
        rlog::log_warn("Not a valid column name in your input CSV data. Try again.")
        print(column_message)
        
        Start <- readline(prompt = "Start (residue number) = ")
    }
    
    rlog::log_info(" Specify the column name indicating the ending peptide residue numbers... OR, enter NA")
    End <- readline(prompt = "End (residue number) = ")
    while (is.null(data[[End]]) & End != "NA"){
        rlog::log_error("Not a valid column name in your input CSV data. Try again.")
        print(column_message)
        
        End <- readline(prompt = "End (residue number) = ")
    }
    
    print(" Specify the column name indicating the peptide sequences... OR, enter NA")
    Sequence <- readline(prompt = "Sequence (peptide) = " )
    while (is.null(data[[Sequence]]) & Sequence != "NA"){
        rlog::log_error("Not a valid column name in your input CSV data. Try again.")
        print(column_message)
        
        Sequence <- readline(prompt = "Sequence (peptide) = " )
    }
    
    print(" Specify the column name indicating the peptide charge state... OR, enter NA")
    Charge <- readline(prompt = "Charge = ")
    while (is.null(data[[Charge]]) & Charge != "NA"){
        rlog::log_error("Not a valid column name in your input CSV data. Try again.")
        print(column_message)
        
        Charge <- readline(prompt = "Charge = ")
    }
    
    print(" Specify the column name indicating the Deuterium uptake values ... OR, enter NA")
    Deu_Uptake <- readline(prompt = "Deu_Uptake = ")
    while (is.null(data[[Deu_Uptake]]) & Deu_Uptake != "NA"){
        rlog::log_error("Not a valid column name in your input CSV data. Try again.")
        print(column_message)
        
        Deu_Uptake <- readline(prompt = "Deu_Uptake = ")
    }
    
    print(" Specify the column name indicating the Deuterium exposure timepoints... OR, enter NA")
    Exposure_Time <- readline(prompt = "Exposure_Time = ")
    while (is.null(data[[Exposure_Time]]) & Exposure_Time != "NA"){
        rlog::log_error("Not a valid column name in your input CSV data. Try again.")
        print(column_message)
        
        Exposure_Time <- readline(prompt = "Exposure_Time = ")
    }
    
    print(" Specify column names indicating relevant experimental conditions ... OR, enter NA")
    print(" IMPORTANT. You can provide more than one column name separared by commas (,) - I will merge them into a single label though.")
    Conditions <- readline(prompt = "Conditions = ")
    column_in_set <- unlist(strsplit(toString(gsub(" ", "", Conditions, fixed = TRUE)), split=",")) %in% data_columns
    while (!all(column_in_set) & Conditions != "NA"){
        rlog::log_error("Not a valid column name in your input CSV data. Try again.")
        print(column_message)
        
        Conditions <- readline(prompt = "Conditions = ")
    }
    
    print(" Specify the column name indicating experimental replicates ... OR, enter NA")
    Replicate <- readline(prompt = "Replicate = ")
    while (is.null(data[[Replicate]]) & Replicate != "NA"){
        rlog::log_error("Not a valid column name in your input CSV data. Try again.")
        print(column_message)
        
        Replicate <- readline(prompt = "Replicate = ")
    }
    
    print(" OPTIONAL. Specify the columns you want to ignore. Otherwise, leave blank.")
    print(" IMPORTANT. You can provide more than one column name separared by commas (,)")
    Ignore <- readline(prompt = "Ignore = ")
    
    print(" OPTIONAL. Specify other column names you want to tag along - I will merge these into a single string chain. Otherwise, leave blank.")
    Other <- readline(prompt = "Other = ")
    
    print(" Indicate whether I should convert your 'Exposure_Time' values. Options: TRUE or FALSE")
    convert_time <- readline(prompt = "convert_time = ")
    if (convert_time) {
        print(" what are the original time units of your data? Available units: h (Hours), m (Minutes), s (Seconds).")
        
        original_time_units <- readline(prompt = "original_time_units = ")
        while (!original_time_units %in% c("s", "m", "h")) {
            rlog::log_error("Not a valid time unit. Available units: h (Hours), m (Minutes), s (Seconds).")
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
    
    if (save_parameters != FALSE) {
        if (file.exists(dirname(save_parameters))) {
            saveRDS(parameters, file = save_parameters)
            
            mssg = paste(" Saved your parameters in ", save_parameters)
            rlog::log_info(mssg)
        }
        
    } else{
        rlog::log_warn("Your parameters were not saved. Provide a valid output path with 'save_parameters = outfile_path'")
    }
    
    return(parameters)
}

##' Pre-process data and output a QFeatures instance with curated data from CSV file
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param data Data content from read from CSV file
##' @param normalise Normalisation method to call `hdxstats::normalisehdx method`
##' @param save_qDF Output path for QFeatures object containing curated data
##' @param parameters List with parameters for data curation pipeline (OPTIONAL)
##' @param parameter_file Path for list with parameters for data curation pipeline (OPTIONAL)
##' @param save_parameters Output path for list with parameters for data curation pipeline
##' @param interactive If neither `parameters` nor `parameter_file` provided, `interactive` mode will call `make_parameters_file`
##' @return QFeatures object
##' @md
##' 
##' @rdname hdxstat-utils
preprocess_data <- function(data, 
                            normalise = FALSE,
                            save_qDF = NULL,
                            parameters = NULL,
                            parameter_file = NULL,
                            save_parameters = FALSE,
                            interactive = FALSE) {
    
    if (interactive == TRUE){
        print(" You chose 'interactive' mode to parse the columns from your CSV content and define parameters to format your output QFeatures data object.")
        parameters <- make_parameter_file(data, save_parameters = save_parameters)
        if (!is.null(parameter_file)){
            rlog::log_warn("You enabled 'interactive' as TRUE. This will override any 'parameter_file' or 'parameters' you provided.")
        }
    }
    
    if (!is.null(parameter_file)){
        rlog::log_info(" You provided a 'parameter_file', I will extract parameters from this to format your output QFeatures data object.")
        if (file_test("-f", parameter_file)){
            parameters <- readRDS(parameter_file)
        }
    }
    
    if (!is.null(parameters)){
        if (is.list(parameters)){
            rlog::log_info(" You provided a list of 'parameters', I will extract parameters from this to format your output QFeatures data object.")
        }
    }
    
    if (is.null(parameter_file) & is.null(parameters)) {
        stop("You either provided a invalid 'parameter_file' or list of 'parameters'. I will quit pre-processing.")
    }
    
    rlog::log_info(" Stripped your 'Exposure_Time' values from non-numeric characters.")
    data[[parameters$Exposure_Time]] <- as.numeric(gsub("[^0-9.-]", "", data[[parameters$Exposure_Time]]))
    
    if (parameters$convert_time){
        if (parameters$original_time_units == 'h') {
            rlog::log_info(" Your original_time_units == 'h'. I will convert your 'Exposure_Time' values to seconds (s).")
            data[[parameters$Exposure_Time]] <- 3600*data[[parameters$Exposure_Time]]
        }
        else if (parameters$original_time_units == 'm') {
            rlog::log_info(" Your original_time_units == 'm'. I will convert your 'Exposure_Time' values to seconds (s).")
            data[[parameters$Exposure_Time]] <- 60*data[[parameters$Exposure_Time]]
        }
        else if (parameters$original_time_units == 's') {
            rlog::log_info(" Your original_time_units == 's'. I will not convert your 'Exposure_Time' values.")
        }
    }
    
    if (parameters$Replicate == "NA"){
        rlog::log_info(" Your 'Replicate' column appears to be NA. I will add this column with 1 values just to label your data.")
        data$Replicate <- 1
        parameters$Replicate <- "Replicate"
    }
    
    if (parameters$Charge == "NA"){
        rlog::log_info(" Your 'Charge' column appears to be NA. I will add this column with 0 values just to label your data.")
        data$Charge <- 0
        parameters$Charge <- "Charge"
    }
    
    rlog::log_info(" Reformatting your data to a wide format...")
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
    rlog::log_info(" Removing NA values from your data")
    data_wide <- data_wide[, colSums(is.na(data_wide)) != nrow(data_wide)]
    # Take all column names except 'columns_fixed'
    columns_to_remove <- 1:length(columns_fixed) # Remove columns_fixed
    old_columns_names <- colnames(data_wide)[-columns_to_remove]
    # Replace "_" with default delimiters and remove trailing strings
    new_object.colnames <- gsub(paste("<>", delimiter.Replicate, sep=""), delimiter.Replicate, old_columns_names)
    new_object.colnames <- gsub("<>", delimiter.Conditions, new_object.colnames)
    new_object.colnames <- gsub(" .*", "", new_object.colnames) # ?
    
    # Parse data for selected columns
    rlog::log_info(" Parsing your data as a qDF object class instance. Method: parseDeutData")
    
    initial_column <- length(columns_fixed)+1 # Fixed value by default
    last_column <- length(columns_fixed)+length(new_object.colnames) # Change to length value
    data_qDF <- parseDeutData(object = DataFrame(data_wide),
                              design = new_object.colnames,
                              quantcol = initial_column:last_column,
                              #rownames = data_wide[[parameters$Sequence]],
                              sequence = parameters$Sequence,
                              charge = parameters$Charge)
    
    rlog::log_info(" Saving a list of 'Start' and 'End' residue numbers as part of qDF object rowData")
    peptide_names_original <- paste0(data[[parameters$Sequence]], "_", data[[parameters$Charge]])
    peptide_names_qDF <- rownames(assay(data_qDF))
    first_matches <- match(unique(peptide_names_qDF), peptide_names_original)
    rowData(data_qDF)[["incoperation"]][["Start"]] <- data[[parameters$Start]][first_matches]
    rowData(data_qDF)[["incoperation"]][["End"]] <- data[[parameters$End]][first_matches]

    # Normalise data 
    if (normalise) {
        rlog::log_info(" Normalising data ... Method: normalisehdx")
        
        data_qDF <- normalisehdx(data_qDF,
                                 sequences = unique(data[[parameters$Sequence]]),
                                 method = "pc")
    }
    else{
        rlog::log_warn("Your output data is not normalised.")
    }
    
    # Save data
    if (!is.null(save_qDF)) {
        if (file.exists(dirname(save_qDF))) {
            
            saveRDS(data_qDF, file = save_qDF)
            rlog::log_info(paste(" Saved output data in ", save_qDF))
            
        }
    } else {
        rlog::log_warn("Your output data was not saved. You can provide an output path with 'save = my_path'")
    }
    
    return(data_qDF)
}

##' Extract data provided a CSV filepath
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' @param data_path CSV data filepath
##' @param normalise Normalisation method to call `hdxstats::normalisehdx method`
##' @param save_qDF Output path for QFeatures object containing curated data
##' @param parameters List with parameters for data curation pipeline (OPTIONAL)
##' @param parameter_file Path for list with parameters for data curation pipeline (OPTIONAL)
##' @param save_parameters Output path for list with parameters for data curation pipeline
##' @param interactive If neither `parameters` nor `parameter_file` provided, `interactive` mode will call `make_parameters_file`
##' @return QFeatures object
##' @md
##' 
##' @rdname hdxstat-utils
extract_hdx_data <- function(data_path,
                             normalise = FALSE,
                             save_qDF = NULL,
                             parameters = NULL,
                             parameter_file = NULL,
                             save_parameters = FALSE,
                             interactive = FALSE) {
    
    if (file.exists(data_path)){
        if (xfun::file_ext(data_path) == "csv") {
            rlog::log_info(" You gave me a CSV file of your HDX-MSM data")
            
            data <- read_csv(data_path, show_col_types = FALSE)
            data.type <- "csv"
        }
        else if (xfun::file_ext(data_path) == "rsd"){
            rlog::log_info(" You gave me a RSD file for your HDX-MSM data")
            rlog::log_info(" I will assume your input data has alreayd been pre-processed")
            
            data <- readRDS(data_path)
            return(data)
        }
        else{
            return(NULL)
            stop("You provided an input file format that I cannot recognise")
            
        }
    }
    else{
        return(NULL)
        stop("This is not a valid path. Try again.")
        
    }
    
    # Pre-process data
    if (data.type == "csv"){
        rlog::log_info(" I will pre-process your data parse it using QFeatures ...")
        data <- preprocess_data(data, 
                                normalise = normalise, 
                                save_qDF = save_qDF, 
                                parameters = parameters,
                                parameter_file = parameter_file,
                                save_parameters = save_parameters,
                                interactive = interactive)
        
        rlog::log_info(" I pre-processed you input CSV data content and now it's available as a QFeatures instance")
        
        return(data)
    }
}