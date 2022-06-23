library(testthat)
library(bain)
library(lavaan)


# set this to local path
source("C:/Users/e_lib/OneDrive/Documents/GitHub/bayesynth/UI Function/pbf.R")


model1 <- 'A =~ Ab + Al + Af + An + Ar + Ac 
           B =~ Bb + Bl + Bf + Bn + Br + Bc'
hypotheses1 <-
  "A=~Ab > .6 & A=~Al > .6 & A=~Af > .6 & A=~An > .6 & A=~Ar > .6 & A=~Ac >.6& 
   B=~Bb > .6 & B=~Bl > .6 & B=~Bf > .6 & B=~Bn > .6 & B=~Br > .6 & B=~Bc >.6"
fit1 <- sem(model1, data = sesamesim[1:80,], std.lv = TRUE)
fit2 <- sem(model1, data = sesamesim[81:160,], std.lv = TRUE)
fit3 <- sem(model1, data = sesamesim[161:nrow(sesamesim),], std.lv = TRUE)

lavaans <- c(fit1,fit2,fit3)


# again, does not replicate exactly when pbf() is run multiple times due to multiple constraints in the hypotheses
test_that("pbf works with lavaan objects", {
  expect_error(pbf(lavaans, hypotheses1, standardize=T), NA)
  res <- pbf(lavaans, hypotheses1, standardize = T)
  expect_equal(res[1,1], 1883568.366)
})
