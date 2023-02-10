vcf.file <- system.file("extdata", "vcf_test.vcf", package="vaRHC")
test_that("example",{
  expect_no_error(variants <- vcfToCodingDNA(file.vcf = vcf.file, nrows = 6000, assembly = "hg19"))
})
