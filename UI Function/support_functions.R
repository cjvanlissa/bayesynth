# Support functions

##### Prepares data for bain() -------------------------------------------------
preprocess_data <- function(data, hypothesis, grouping){
  # obtain vector of unique groups
  unique_groups <- unique(data[grouping])[,1]    
  k <- length(unique_groups)
  if(!k > 1){stop("At least two groups are required to aggregate Bayes Factors")}
  
  # obtain the hypothesized value
  hyp_val <- as.numeric(unlist(regmatches(hypothesis, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",hypothesis))))
  constraint <- unlist(regmatches(hypothesis, gregexpr("(>|<|=)", hypothesis)))
  
  # obtain per group the data.frame to calculate BF on.
  df_groups <- lapply(unique_groups, function(group){
    df_group <- subset(data[data[grouping] == group, ], select = -eval(parse(text = grouping)))
    sapply(df_group, as.numeric)
  })
  
  # get the correlation coefficient and standard error per group 
  res <- sapply(df_groups, function(df_group){
    est <- cor(df_group)[2,1]                        #estimate for correlation between x and y
    se_est <- sqrt((1 - est^2)/(nrow(df_group) - 2)) #estimate for standard error = sqrt((1-r^2)/df) with df = N - 2
    c(est, se_est)               #return estimate, se and sample size of particular set
  })
  colnames(res) <- paste0('r', 1:k)
  
  return(list(hyp_val = hyp_val, res = res, k = k, constraint = constraint))
}

###### Calculates allsig -------------------------------------------------------
# Do not know how this function precisely works.... Copied it from the simulation script.
allsig <- function(res, hyp_val){
  classic_allsig <- all(apply(res, 2, function(x){
    (x[1]-hyp_val)/x[2] > 1.644854
  }))
  c(0,10)[classic_allsig+1]
}


###### Calculates tbf_iu and tbf_ic --------------------------------------------
tbf <- function(res, hyp_val, n, k, constraint){
  bf_together <- bain(res[1,], 
                      hypothesis = gsub("(r1)", "r1", paste0("(", paste0(colnames(res), collapse = ", "), ")", constraint, hyp_val), fixed = TRUE), 
                      n = sum(rep(n, k)),
                      Sigma = diag(res[2,], ncol = ncol(res)))
  return(c(tbf_ic = bf_together$fit$BF.c[1], tbf_iu = bf_together$fit$BF.u[1]))
}

###### Calculates prodbf_iu, prodbf_ic, gpbf_iu and gpbf_ic --------------------
prod_and_gpbf <- function(res, hyp_val, n, k, constraint){
  sig <- lapply(res[2,], matrix)
  bf_individual <- lapply(paste0(colnames(res), constraint, hyp_val), 
                          bain,       
                          x = res[1,],
                          Sigma = sig,        
                          n = rep(n, k),           
                          group_parameters = 1, 
                          joint_parameters = 0) 
  
  BFs <- t(sapply(bf_individual, function(x){
    c(x$fit$BF.c[1], x$fit$BF.u[1])
  }))
  
  gp_and_prod <- apply(BFs, 2, function(x){
    prod_bf <- prod(x)         
    c(prod_bf^(1/k), prod_bf)  
  })
  
  return(c(gpbf_ic = gp_and_prod[1,][1], gpbf_iu = gp_and_prod[1,][2],
           prodbf_ic = gp_and_prod[2,][1], prodbf_iu = gp_and_prod[2,][2]))
}

##### To generate an example data.frame ----------------------------------------
example_df <- function(){
  set.seed(6164900)
  y <- rnorm(1000)
  x <- y + rnorm(1000,0,3)
  z <- rnorm(1000)
  M <- as.data.frame(cbind(y = y,        # generate outcome variable from std normal dist
             x = x,                      # generate predictor from normal dist
             z = z,
             k = rep(c("A", "B", "C", "D"), each = 250)))   # make 4 groups
  M[,c(1:3)] <- sapply(M[,c(1:3)], as.numeric)
  return(M)
}
 