# Maybe I am thinking to hard about this, but it seems that the function does not yet deal with
# hypotheses in bain objects potentially having different locations.
# I'll show an example

library(bain)

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
  BFs <- do.call(cbind, lapply(x, function(y){y$fit$BF.c[which(y$hypotheses %in% hyps)]}))
  res <- list(BFs = BFs, pbf = apply(BFs, 1, prod)) # obtain pbf ic, might need to change dependent on alternative hyp
  return(res)
}


# we create four lm objects
set.seed(6164)
lms <- lapply(1:4, function(i){
  tt = as.data.frame(cbind(y = rnorm(1000,0,1), x = rnorm(1000, 0.2,1), z=rnorm(1000, 0.2,1)))
  lm(tt$y ~ tt$z + tt$x)
})
pbf(lms, "x=z;z<0") #pbfs of 32957 and 529




# now I call bain on the lm objects, but the change the order of the two hypotheses for the first and last two lm objects
bains_order1 <- lapply(lms[1:2], bain, hypothesis = 'x=z;z<0')
bains_order2 <- lapply(lms[3:4], bain, hypothesis = 'z<0;x=z')
bains <- c(bains_order1,bains_order2)




# now running pbf(bains) should give the same output as when the pbf was called on the lm objects, but it does not
pbf(bains)  # pbfs of 324 and 53815
# this is because the values in column 3 and 4, belonging to objects 3 and 4 are now swapped


## suggestion for improvement: Line 93 - 98
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
  
  BFs <- do.call(cbind, lapply(x, function(y){
    hyps_in <- which(y$hypotheses %in% hyps)           # find which hypotheses exist in y
    hyps_indices <- match(hyps[hyps_in], y$hypotheses) # find index of those hypotheses
    y$fit$BF.c[hyps_indices]                           # get hypotheses in correct order
  }))
  rownames(BFs) <- paste0(sprintf('H%d: ', 1:length(hyps)),hyps) # give names
  
  res <- list(BFs = BFs, pbf = apply(BFs, 1, prod)) # obtain pbf ic, might need to change dependent on alternative hyp
  return(res)
}

pbf(lms, "x=z;z<0") # pbfs of 32957 and 529
pbf(bains)          # pbfs of 529 and 32957, only the order is different but the names give away which hypothesis is evaluated

