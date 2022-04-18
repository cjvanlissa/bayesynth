# Support functions

allsig <- function(res, hyp_val){
  classic_allsig <- all(apply(res, 2, function(x){
    (x[1]-hyp_val)/x[2] > 1.644854
  }))
  c(0,10)[classic_allsig+1]
}


tbf <- function(res, hyp_val, n, k){
  bf_together <- bain(res[1,], 
                      hypothesis = gsub("(r1)", "r1", paste0("(", paste0(colnames(res), collapse = ", "), ") > ", hyp_val), fixed = TRUE), 
                      n = sum(rep(n, k)),
                      Sigma = diag(res[2,], ncol = ncol(res)))
  return(c(tbf_ic = bf_together$fit$BF.c[1], tbf_iu = bf_together$fit$BF.u[1]))
}

prod_and_gpbf <- function(res, hyp_val, n, k){
  sig <- lapply(res[2,], matrix)
  bf_individual <- lapply(paste0(colnames(res), ">", hyp_val), 
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