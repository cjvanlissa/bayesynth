library(bain)

# pbf moet ook werken voor meerdere hypotheses en de posterior model probabilities
# moeten gereturned worden.
# zodra functie compleet is, fork bain en implementeer
pbf <- function(x, hypothesis, ...){
  
  if(!all(sapply(x, inherits, what = "bain"))){ 
    cl <- match.call()
    cl[[1L]] <- quote(bain)  
    for(i in (1:length(x))){
      cl[['x']] <- x[[i]]
      x[[i]] <- eval.parent(cl)
    }
    # cl[['x']] <- lapply(x, eval.parent, cl)
    cl[['x']] <- x
    cl[[1L]] <- quote(pbf)
    eval.parent(cl)
  }
  
  hypotheses <- sapply(x, function(y){y$hypotheses}) # check if hypotheses are equal
  if(length(unique(hypotheses)) != 1){
    stop("hypotheses of bain objects are not all equal")}
  BFs <- sapply(x, function(y){y$fit$BF.c[1]})
  res <- list(BFs = BFs, pbf = prod(BFs)) # obtain pbf ic, might need to change dependent on alternative hyp
  return(res)
}

# try out with non-bain objects
ttests <- lapply(1:4, function(i){
  tt = as.data.frame(cbind(y = rnorm(1000,0,1), x = rnorm(1000, 0.2,1)))
  t_test(tt$y, tt$x)
})
pbf(ttests, "x=y")

# and with bain-objects
bains <- lapply(ttests,bain, hypothesis = "x=y")
pbf(bains, "x=y")
