context("Functional fitting tests")

test_that("functional fits", {
    
    data(hdx_data) 
    
    # INPUT
    data_selection <- hdx_data[,1:100]
    all_peptides <- rownames(data_selection)[[1]][1]
    starting_parameters <- list(a = NULL, b = 0.001,  d = NULL, p = 1)
    
    # OUTPUT
    results <- analyse_kinetics(data = data_selection,
                                method = "dfit",
                                peptide_selection = all_peptides,
                                start = starting_parameters)
    
    data_selection <- hdx_data[,1:100]
    all_peptides <- rownames(data_selection)[[1]]
    starting_parameters <- list(a = NULL, b = 0.001,  d = NULL, p = 1)
    
    # OUTPUT
    results2 <- analyse_kinetics(data = data_selection,
                                method = "dfit",
                                peptide_selection = all_peptides[1],
                                start = starting_parameters)
    expect_equal(results, results2)
    expect_equal(class(results$fitted_models)[1], "HdxStatModel")
    
    # correct number of models fitted
    num_models <- length(unique(sub(".*cond", "",colnames(data_selection)[[1]])))
    expect_equal(num_models, length(results$fitted_models@alternative))
    expect_equal("Functional Model", results$fitted_models@method)
    
    # did plotting work
    expect_equal("gg", class(results$fitted_models@vis)[1])
    
})
