library(testthat)
library(bain)

# set this to local path
source("C:/Users/e_lib/OneDrive/Documents/GitHub/bayesynth/UI Function/pbf.R")


set.seed(100)
lms <- lapply(1:4, function(i){
  tt = as.data.frame(cbind(y = rnorm(1000,0,1), x = rnorm(1000, 0.2,1), z = rnorm(1000,0.1,1)))
  lm(tt$y ~ tt$x + tt$z)
})

test_that("pbf works and replicates on bain objects, single same hypothesis", {
  bains <- lapply(lms, bain, "x<z")
  expect_error(pbf(bains), NA)
  res_ini <- pbf(bains)
  expect_equivalent(res_ini[1,1], 0.4943658)
  for(i in c(1:4)){
    res <- pbf(bains)
    expect_equivalent(res[1,1], res_ini[1,1])
  }
})

test_that("pbf does replicate for bains, multiple (in)equality signs", {
  bains <- lapply(lms, bain, "0<z<x")
  expect_error({pbf(bains)}, NA)
  res_ini <- pbf(bains, "0<z<x")
  for(i in c(1:4)){
    res <- pbf(bains, "0<z<x")
    expect_equivalent(round(res[1,1],6), round(res_ini[1,1],6))
  }
})

# Might be an idea to calculate pbf for all hypotheses that occur in bain objects anyway
# but log which objects contained which hypotheses
# so in this case H: x<z will have its pbf calculated for sample 1 and 2
# and x>z for sample 3 and 4
test_that("pbf throws error if hypotheses differ within list of bain objects", {
  bains_hypdiffs <- c(
    lapply(lms[1:2], bain, "x<z"),
    lapply(lms[3:4], bain, "x>z")
  )
  expect_error(pbf(bains_hypdiffs))
})

test_that("pbf works with difference in hypothesis ordering within bain objects", {
  bains_same_order <- lapply(lms, bain, "x=0;x<z")
  expect_error(pbf(bains_same_order),NA)
  
  bains_diff_order <- c(
    lapply(lms[1:2], bain, "x<z;x=0"),
    lapply(lms[3:4], bain, "x=0;x<z")
    )
  expect_error(pbf(bains_diff_order),NA)
  
  expect_equivalent(pbf(bains_same_order), pbf(bains_diff_order))
})

test_that("pbf only retains hypotheses that are available in all objects", {
  lms2 <- lapply(1:4, function(i){
    tt = as.data.frame(cbind(y = rnorm(1000,0,1), x = rnorm(1000, 0.2,1), z = rnorm(1000,0.1,1),
                             a = rnorm(1000,0.3,1)))
    lm(tt$y ~ tt$x + tt$z + tt$a)
  })
  
  bains_single_common_hyp <- c(
    list(bain(lms2[[1]], "x<z;x=0")),
    list(bain(lms2[[2]], "x<z;x=0;a=3")),
    lapply(lms2[3:4], bain, "x=0;z>0.1")
  )
  
  expect_error(pbf(bains_single_common_hyp),NA)
  res <- pbf(bains_single_common_hyp)
  expect_equivalent(rownames(res), "H1: x=0")
})
