context("Test statistics")

test_that("statistics", {
    
    data(hdx_data)
    
    # INPUT
    n <- 5
    data_selection <- hdx_data[,seq.int(24)]
    all_peptides <- rownames(data_selection)[[1]][seq.int(n)]
    starting_parameters <- list(a = NULL, b = 0.001,  d = NULL, p = 1)
    
    # OUTPUT
    results <- analyse_kinetics(data = data_selection,
                                method = "fit",
                                peptide_selection = all_peptides,
                                start = starting_parameters)
    
    expect_equal(nrow(results$functional_analysis@results), n)
    expect_true(all(results$functional_analysis@results$ebayes.fdr <= 1))
    expect_true(all(results$functional_analysis@results$ebayes.fdr >= 0))
    expect_equal(rownames(results$functional_analysis@results), all_peptides)
    
    results2 <- analyse_kinetics(data = data_selection,
                                method = "dfit",
                                peptide_selection = all_peptides[1],
                                start = starting_parameters)
    
    
    
    expect_equal(length(anova(results2$fitted_models)$Res.Df), 3) 
    expect_equal(names(logLik(results2$fitted_models)), c("null", "alt1", "alt2"))
    expect_equal(length(residuals(results2$fitted_models)), 3)
    expect_equal(length(vcov(results2$fitted_models)), 2)
    expect_true(likRatio(results2$fitted_models) > 0 )
    expect_true(wilk(results2$fitted_models) <= 1)
    expect_equal(ncol(coef(results2$fitted_models)), length(starting_parameters))
    expect_equal(deviance(results2$fitted_models), 3)
    
    # check default formula
    formula <- value ~ a * (1 - exp(-b * (timepoint)^p)) + d
    out <- summary(results2$fitted_models)
    expect_equal(out$null$formula, formula)
    
    
})