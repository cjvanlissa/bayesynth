library(testthat)
library(bain)

# set this to local path
source("C:/Users/e_lib/OneDrive/Documents/GitHub/bayesynth/UI Function/pbf.R")


set.seed(100)
ttests <- lapply(1:4, function(i){
  tt = as.data.frame(cbind(y = rnorm(1000,0,1), x = rnorm(1000, 0.2,1)))
  t_test(tt$y, tt$x)
})

test_that("pbf works for t_tests, single hypothesis", {
  expect_error({pbf(ttests, "x=y")}, NA)
  res <- pbf(ttests, "x=y")
  expect_equivalent(res[1,1], 1.644e-14)
})

test_that("pbf works and replicates for t_tests, multiple hypotheses", {
  expect_error({pbf(ttests, "x=y;x>0;y=0")}, NA)
  for(i in c(1:5)){
    res <- pbf(ttests, "x=y;x>0;y=0")
    expect_equivalent(res[1,1], 1.644e-14) 
    expect_equivalent(res[2,1], 0.2393395)
    expect_equivalent(res[3,1], 2.014e-34)
  }
})

