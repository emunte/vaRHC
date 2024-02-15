test_that("example",{
  skip("file does not exist")
  # aixo tenia nrows = 6000
  vcf.file <- system.file("extdata", "vcf_test.vcf", package="vaRHC")
  expect_no_error(variants <- vcfToCodingDNA(file.vcf = vcf.file, assembly = "hg19"))
})
