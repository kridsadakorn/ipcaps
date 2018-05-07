context("test cal.eigen.fit")

test_that("test cal.eigen.fit",{

  data(example_SNP)

  PCs <- cal.pc.linear(simsnp$snp, no.pc = 3)

  expect_length(PCs, 2)
  expect_equal(fst.pairwise[2], -0.00195062)
  expect_type(ls(snp), "character")

  save.file <- file.path(tempdir(),"new_SNP")
  file.remove(file1)
})
