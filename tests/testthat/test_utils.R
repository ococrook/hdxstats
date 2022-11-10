context("test utils")

test_that("utils", {
    
    out <- exchangeableAmides("ABPPPPP")
    out2 <- exchangeableAmides("ABCDEFG")
    out3 <- exchangeableAmides("APBCDEF")
    out4 <- exchangeableAmides("PPPPPP")
    
    expect_equal(out, 0)
    expect_equal(out2, 5)
    expect_equal(out3, 5)
    expect_equal(out4, 0)
    
    
})