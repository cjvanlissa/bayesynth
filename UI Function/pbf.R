library(bain)

# pbf moet ook werken voor meerdere hypotheses en de posterior model probabilities
# moeten gereturned worden.
# zodra functie compleet is, fork bain en implementeer
pbf <- function(x, ...){
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
  # Merge the hypotheses from list item 1 and 2 into object merged
  if(length(x) > 1){
    hyps <- x[[1]]$hypotheses
    for(i in length(x)-1){
      browser()
      hyps <- c(hyps, x[[i+1]]$hypotheses)
      # Drop all non-duplicated hypotheses from merged
      hyps <- hyps[duplicated(hyps)]
      # If merged now has length 0, throw error
      if(length(hyps) == 0){
        stop("The objects passed to pbf() have no hypotheses in common.")
      }
      # Else, go back to step 1, but now merge merged with list item 3
    }
  }
  BFs <- sapply(x, function(y){y$fit$BF.c[which(y$hypotheses %in% hyps)]})
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
