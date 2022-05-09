renamefactors <- function(x){
  out <- gsub("moderators", "M", x, fixed = T) #this changes 'moderators' to M
  out <- gsub("tau2", "$\\tau^2$", out, fixed = T) #nothing
  out <-gsub("es", "$\\beta$", out, fixed = T) #changes 'es'to '$beta$ but I do not know why
  out<-gsub("mean_n", "$n$", out, fixed = T) #nothing
  out<-gsub("model", "Model ", out, fixed = T) #capitalizes 'model'
  out<-gsub("k_train", "$k$ ", out, fixed = T)
  out<-gsub("alpha_mod", "$\\omega$ ", out, fixed = T)
  out
}

interpret <- function(x){
  thelng <- (length(x)-1)
  if(is.na(x[1])) return(NA)
  x <- sum(sign(diff(x)))
  if(x == thelng){
    return("positive")
  }
  if(-1*x == thelng){
    return("negative")
  }
  return("other")
}



#creates traceplots for metrics
Traceplot <- function(test, train, param, alg){
  plot(x = 1:length(test), y = test, type = 'n', main = paste('Traceplot Test and Train', param , ' for', alg), xlab = 'iterations', ylab = paste('values for ', param)) #empty plot with correct dimensions for the traced values
  lines(x = 1:length(test),  y = test, col = 'red') #superimpose traced values for testing r2
  lines(x= 1: length(train), y = train, col = 'blue') #superimpose traced values for training r2
  legend('bottomleft', c('Train R2', 'Test R2'), lty = 1, col = c('blue', 'red'), bg = 'white', cex = 0.75) #add informative legend
}


#creates densities for metric
plotdens <- function(df){
  df <- as.data.frame(df)
  for(column in (lc+1):ncol(df)){ #plot densities for all algorithms
    plot(density(df[,column]),
         main = paste0('histogram of ', colnames(df)[column])
         #xlim = c(-30,1)
    )
  }
}

#Check for significant moderators using eta squared
EtaSq<-function (x)
{
  anovaResults <- summary.aov(x)[[1]]
  anovaResultsNames <- rownames(anovaResults)
  SS <- anovaResults[,2] #SS effects and residuals
  k <- length(SS) - 1  # Number of factors
  ssResid <- SS[k + 1]  # Sum of Squares Residual
  ssTot <- sum(SS)  # Sum of Squares Total
  SS <- SS[1:k] # takes only the effect SS
  anovaResultsNames <- anovaResultsNames[1:k]
  etaSquared <- SS/ssTot # Should be the same as R^2 values
  partialEtaSquared <- SS/(SS + ssResid)
  res <- cbind(etaSquared, partialEtaSquared)
  colnames(res) <- c("Eta^2", "Partial Eta^2")
  rownames(res) <- anovaResultsNames
  return(res)
}


