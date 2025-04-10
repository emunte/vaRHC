test_that("wrong_input", {
  expect_error(vaR())
  expect_error(vaR(gene="sidsd", variant = "c.1A>G"))
  expect_error(vaR(gene = "MLH1", variant= "c.3A>G"))
  #expect_error(vaR(gene = "CDKN2A", variant= "c.1A>G", spliceai.program = TRUE))
  })


test_that("right_input",{
  expect_no_error(vaR(gene="BRCA1", variant = "c.211A>G"))
})

test_that("no_spliceai_precalculated",{
  expect_error(vaR(gene = "RAD51C", variant = "c.1A>G", spliceai.program = FALSE))
})

#test_that("dels_warning",{
#  expect_warning(vaR(gene = "MLH1", variant = "c.1232_1233del", spliceai.program = FALSE, remote = FALSE))
#})
