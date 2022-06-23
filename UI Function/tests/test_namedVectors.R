library(testthat)
library(bain)

# set this to local path
source("C:/Users/e_lib/OneDrive/Documents/GitHub/bayesynth/UI Function/pbf.R")


# with only estimates to work with, the sample information must first be converted to bain objects
# pbf() does not work with raw estimates
sample1 <- list(
  ests = c(x = 0.2, z = -0.3),
  covmat = diag(c(x = 0.04^2, z = 0.03^2)),
  N = 500
)

sample2 <- list(
  ests = c(x = 0.1, z = 0),
  covmat = diag(c(x = 0.06^2, z = 0.02^2)),
  N = 100
)

alldata <- list(sample1, sample2)


test_that("pbf works with named vector",{
  bains <- lapply(alldata, function(sample){
    bain(x = sample[['ests']], Sigma = sample[['covmat']], n = sample[['N']], hypothesis = "z=0")
  })
  expect_error(pbf(bains), NA)
  expect_equivalent(pbf(bains)[1,1], 0)
})

