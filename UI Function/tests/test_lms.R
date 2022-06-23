library(testthat)
library(bain)

# set this to local path
source("C:/Users/e_lib/OneDrive/Documents/GitHub/bayesynth/UI Function/pbf.R")


set.seed(100)
lms <- lapply(1:4, function(i){
  tt = as.data.frame(cbind(y = rnorm(1000,0,1), x = rnorm(1000, 0.2,1), z = rnorm(1000,0.1,1)))
  lm(tt$y ~ tt$x + tt$z)
})


test_that("pbf works for lms, single hypothesis", {
  expect_error({pbf(lms, "x>z")}, NA)
  res <- pbf(lms, "x>z")
  expect_equivalent(round(res[1,1],5), round(2.022794,5))
})

# NOTE: testing hypotheses with multiple inequality signs does NOT exactly replicate, shown in fourth test.
test_that("pbf works for lms, multiple hypotheses", {
  expect_error({pbf(lms, "x>z;0<z<1")}, NA)
  res <- pbf(lms, "x>z;0<z<1")
  expect_equivalent(round(res[1,1],4), round(2.022794,4))
  expect_equivalent(round(res[2,1],4), round(0.09675034,4))
})

test_that("pbf does replicate for lms, single inequality or equality sign", {
  expect_error({pbf(lms, "x>z;z=0.1")}, NA)
  res_ini <- pbf(lms, "x>z;z=0.1")
  for(i in c(1:5)){
    res <- pbf(lms, "x>z;z=0.1")
    expect_equivalent(res[1,1], res_ini[1,1])
    expect_equivalent(res[2,1], res_ini[2,1])
  }
})

# The replication issue seems to be in lapplying bain to all lms objects
# replicating a single (in)equality constraint is no issue.
# so pbf() works fine, but bain seems to not exactly replicate findings if multiple constraints are used in same hypothesis
test_that("pbf does NOT replicate for lms, multiple (in)equality signs", {
  expect_error({pbf(lms, "0<z<x")}, NA)
  res_ini <- pbf(lms, "0<z<x")
  for(i in c(1:4)){
    res <- pbf(lms, "0<z<x")
    expect_equivalent(round(res[1,1],6), round(res_ini[1,1],6))
  }
})

