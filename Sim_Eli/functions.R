


BFs <- function(es, n, hyp_val, k, errorsd){
  
  # create dataframe for each group
  dfs <- lapply(1:k, function(i){
    df <- rmvnorm(n, sigma = matrix(c(1, es, es, 1), nrow = 2))
    df + rnorm(2*n, sd = errorsd)
  })
  
  # obtain estimates for correlations and their standard errors for every dataset
  res <- sapply(dfs, function(x){
    est <- cor(x)[2,1]                        #estimate for correlation between x and y
    se_est <- sqrt((1 - est^2)/(n - 2)) #estimate for standard error = sqrt((1-r^2)/df) with df = N - 2
    c(est, se_est)               #return estimate, se and sample size of particular set
  })
  
  # Classic approach
  classic_allsig <- all(apply(res, 2, function(x){
    x[1]-1.96*x[2] > hyp_val
  }))
  
  # necessary naming for bain and further preparing
  colnames(res) <- paste0('r', 1:k)
  sig <- lapply(res[2,], matrix)    # make list of covariance matrices for the datasets, as shown in the Bain vignette
  #ngroup <- rep(n, k)       # obtain sample size per group
  
  #run bf_individual to extract product bf and geometric product bf
  bf_individual <- lapply(paste0(colnames(res), ">", hyp_val), # for every group, we hypothesize that r_k > hyp_val
                          bain,                 # call bain
                          x = res[1,],          # estimates
                          Sigma = sig,          # all rho's are assumed to be independent
                          n = rep(n, k),           # pass the named vector of sample sizes per group to bain
                          group_parameters = 1, # every group k has 1 parameter which is rho_xy
                          joint_parameters = 0) # they do not share parameters 

  # extract BF_ic and BF_iu for the parameter of every group
  BFs <- t(sapply(bf_individual, function(x){
    c(x$fit$BF.c[1], x$fit$BF.u[1])
  }))

  # obtain geometric product and regular product
  gp_and_prod <- apply(BFs, 2, function(x){
    prod_bf <- prod(x)         #obtain product bf
    c(prod_bf^(1/k), prod_bf)  #concatenate geometric product and regular product
  })
  
  # create bf_together object to obtain the BFs for the (group1, group2, group3) > hyp_val
  bf_together <- bain(res[1,], 
                      hypothesis = gsub("(r1)", "r1", paste0("(", paste0(colnames(res), collapse = ", "), ") > ", hyp_val), fixed = TRUE), 
                      n = sum(rep(n, k)), # n = sum of sample sizes over all groups.
                      Sigma = diag(res[2,], ncol = ncol(res)))      # assume independence between groups
  
  # returns in order: gpbf_ic, gpbf_iu, prodbf_ic, prodbf_iu, tbf_ic, tbf_iu
  return(c(
    c(0,10)[classic_allsig+1],
    gp_and_prod[1,],
    gp_and_prod[2,], 
    c(bf_together$fit$BF.c[1], bf_together$fit$BF.u[1])))
}




#product BF
# gPBF <- function(BFs){
#   N  <- ifelse(is.null(nrow(BFs)), length(BFs), nrow(BFs)) #to distinguish between vectors and matrices
#   
#   res <- apply(BFs, 2, function(x){
#     GP <- prod(x) ^ (1 / N) #Geometric mean
#     ER <- abs((GP < 1) - sum(x > 1)/N) #Evidence Rate
#     SR <- ifelse(GP < 1, #Stability Rate
#                  sum(x < GP) / N, 
#                  sum(x > GP) / N)
#     c(GP, ER, SR)
#   })
#   
#   rownames(res) <- c("Geometric Product", "Evidence Rate", "Stability Rate")
#   out <- list("GPBF" = res, "BFs" = BFs, "N" = N)
#   class(out) <- "gPBF"
#   return(out)
# }
# 
# 
# prop_correct <- function(x, BF_threshold, var_n){
#   if(length(BF_threshold) == 1 && BF_threshold == 3){
#     sum(x >= BF_threshold) / (var_n)
#   } else if(length(BF_threshold) == 1 && BF_threshold == 0.33){
#     sum(x <= 0.33) / (var_n)
#   }else{
#     sum(x > min(BF_threshold) & x < max(BF_threshold)) / (var_n)
#   }
# }
