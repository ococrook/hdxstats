context("data parsing consistency")


test_that("data parse", {
    
    
    data(hdx_data) 
    
    # INPUT
    csv_filepath <- csv_filepath <- system.file("extdata", "MBP.csv", package = "hdxstats") # input csv with HDX-MS datas
    parameter_file <- "vignettes/tests/MBP_qDF.hdxp"
    
    # OUTPUT
    output_file = "vignettes/tests/MBP_qDF.rsd" # QFeatures object output path
    
    # Run curation pipeline provided 'parameter_file'
    hdx_data2 <- extract_hdx_data(csv_filepath, save_qDF = output_file, parameter_file = parameter_file)
    
    expect_equal(hdx_data, hdx_data2)

})