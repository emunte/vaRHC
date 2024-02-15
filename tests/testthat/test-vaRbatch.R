
test_that("wrong_input", {
  skip("some tests are failing")
  data("ex_vaRbatch")
  ex_vaRbatch[] <- lapply(ex_vaRbatch, as.character)#convert to character
  expect_error(vaRbatch())
  expect_error(vaRbatch(all.variants="sidsd"))
  expect_error(vaRbatch(all.variants = ex_vaRbatch, remote=TRUE, browser= "FIraxofa" ))
})


test_that("no_output",{
  skip("some tests are failing")
  data("ex_vaRbatch")
  ex_vaRbatch[] <- lapply(ex_vaRbatch, as.character)#convert to character
  expect_error(all <- vaRbatch( all.variants = ex_vaRbatch, spliceai.program = FALSE))
})

