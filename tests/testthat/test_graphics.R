context("test visualisations")


test_that("visualisations", {
    
    
    data(hdx_data)
    m <- 5
    # INPUT
    data_selection <- hdx_data[, seq.int(24)]
    all_peptides <- rownames(data_selection)[[1]][seq.int(m)]
    starting_parameters <- list(a = NULL, b = 0.001,  d = NULL, p = 1)
    
    # OUTPUT
    results <- analyse_kinetics(data = data_selection,
                                method = "fit",
                                peptide_selection = all_peptides,
                                start = starting_parameters)
    
    
    graphics_kinetics <- visualise_hdx_data(results, type="kinetics")
    graphics_forest <- visualise_hdx_data(results, type="forest")
    
    
    expect_equal(length(graphics_kinetics), m)
    expect_equal(length(graphics_forest), m)
    expect_equal(vapply(graphics_kinetics, function(x) x$data$rowname[1], "character"),
                 all_peptides)
    expect_equal(vapply(graphics_forest, function(x) x$data$rowname[1], "character"),
                 rep("a", m))
    
    
    
})
