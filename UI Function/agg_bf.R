# this functions main input must be a formula, for example y ~ x as well as a hypothesis, for example x > 0.2.
# the user must also specify the variable that represents the groups
# this function than tests the hypothesis for multiple groups and calculates an aggregated BF.
# the user can specify which algorithm to use (for example, "allsig", "gpbf_ic", "prodbf_ic", "prodbf_iu" and "tbf_ic").


# for example, an input matrix can look as follows
set.seed(6164900)
y <- rnorm(1000)
x <- y + rnorm(1000,0,3)
cor(x,y)
M <- cbind(y = y,            # generate outcome variable from std normal dist
           x = x,        # generate predictor from normal dist
           k = rep(c("A", "B", "C", "D"), each = 250))   # make 4 groups
# View(M)

# suppose the hypothesis is that the correlation coefficient of x and y > 0.2
# then for every k, the BF is calculated for the formula y ~ x 
# finally it is aggregated using some algorithm specified and this aggregated BF is returned

# this function below is an example of what it would look like.
# it relies on three other functions that occur in the simulation
source(file.path("UI Function", "support_functions.R"))

agg_bf <- function(formula = F, grouping, data, algorithm, hypothesis){
  
  # see if bain is can be loaded
  tryCatch(library(bain), error = function(e){stop("installation of package 'bain' is required")})
  
  # see if chosen algorithm is supported
  supported_algorithms <- c("allsig", "gpbf_ic", "gpbf_iu", "prodbf_ic", "prodbf_iu", "tbf_ic", "tbf_iu")
  if(!algorithm %in% supported_algorithms){
    stop(paste("not a valid algorithm, please use one of the supported algorithms: ", 
               paste(supported_algorithms, collapse = ", ")))
  }
  
  # convert to data frame for easy subsetting and obtain vector unique groups
  data <- as.data.frame(data)                    
  unique_groups <- unique(data[grouping])[,1]    
  k <- length(unique_groups)
  if(!k > 1){stop("At least two groups are required to aggregate Bayes Factors")}
  
  # obtain the hypothesized value
  hyp_val <- as.numeric(unlist(regmatches(hypothesis, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",hypothesis))))
  
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
  
  # return bf based on chosen algorithm
  out <- switch(algorithm,
      allsig ={allsig(res, hyp_val)},
      tbf_ic = {tbf(res, hyp_val, floor(nrow(data)/k), k)['tbf_ic']},
      tbf_iu = {tbf(res, hyp_val, floor(nrow(data)/k), k)['tbf_iu']},
      prodbf_ic ={prod_and_gpbf(res, hyp_val, floor(nrow(data)/k), k)['prodbf_ic']},
      prodbf_iu ={prod_and_gpbf(res, hyp_val, floor(nrow(data)/k), k)['prodbf_iu']},
      gpbf_ic ={prod_and_gpbf(res, hyp_val, floor(nrow(data)/k), k)['gpbf_ic']},
      gpbf_iu ={prod_and_gpbf(res, hyp_val, floor(nrow(data)/k), k)['gpbf_iu']})
  
  return(out)
}

agg_bf(F, grouping = 'k', data = M, algorithm = 'prodbf_ic', hypothesis = 'x>0.2')
agg_bf(F, grouping = 'k', data = M, algorithm = 'tbf_ic', hypothesis = 'x>0.2')

