
data("ex_vaRbatch")
ex_vaRbatch[] <- lapply(ex_vaRbatch, as.character)#convert to character

test_that("wrong_input", {
  expect_error(vaRbatch())
  expect_error(vaRbatch(all.variants="sidsd"))
  expect_error(vaRbatch(all.variants = ex_vaRbatch, remote=TRUE, browser= "FIraxofa" ))
})


test_that("no_output",{
  expect_error(all <- vaRbatch( all.variants = ex_vaRbatch, spliceai.program = FALSE))
})

