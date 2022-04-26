# agg_bf() is a function that can aggregated Bayes Factors for hypotheses.
# For now, only correlation coefficients are supported, but this could be expanded to formulas and multiple hypotheses.
# Suppose the hypothesis is that the correlation coefficient of y and x > 0.2.
# Then for every k, the BF is calculated for the formula y ~ x based on the algorithm the user selects (for example, "allsig", "gpbf_ic", "prodbf_ic", "prodbf_iu" and "tbf_ic").
# Finally it is aggregated using some algorithm specified and this aggregated BF is returned.
# This means the user must also specify the variable that represents the groups.


# The function below is an example of what it would look like.
# It relies on four other functions.
# 1.) preprocess_data(): function to make data ready to be used by bain
# 2.) allsig(): function to calculate allsig BF (do not know this algorithm, however and might be incorrect)
# 3.) tbf(): calculates tbf_ic and tbf_iu
# 4.) prod_and_gpbf: calculated prodbf_ic, prodbf_iu, gpbf_ic, gpbf_iu

# for example, an input data frame can look as follows, where z remains unused for now
source(file.path("UI Function", "support_functions.R"))
M <- example_df()

# the function
agg_bf <- function(formula, grouping, data, algorithm, hypothesis){
  
  # see if bain is can be loaded
  tryCatch(library(bain), error = function(e){stop("installation of package 'bain' is required")})
  
  # see if chosen algorithm is supported
  supported_algorithms <- c("allsig", "gpbf_ic", "gpbf_iu", "prodbf_ic", "prodbf_iu", "tbf_ic", "tbf_iu")
  if(!algorithm %in% supported_algorithms){
    stop(paste("not a valid algorithm, please use one of the supported algorithms: ", 
               paste(supported_algorithms, collapse = ", ")))
  }
  
  # convert data to data.frame and subset for only variables used in analysis
  data <- tryCatch(as.data.frame(data), error = function(e){stop('Make sure data can be converted to data.frame')})
  data <- cbind(stats::model.frame(formula, data), data[grouping])
  if(ncol(data) < 3){stop("data must have at least 3 columns")}
  
  # prepare data for bain and obtain parameters to be passed to tbf() and prod_and_gpbf() functions
  prep <- preprocess_data(data, hypothesis, grouping)
  
  # return bf based on chosen algorithm
  out <- switch(algorithm,
      allsig    = {allsig(prep$res, prep$hyp_val)},
      tbf_ic    = {tbf(prep$res, prep$hyp_val, floor(nrow(data)/prep$k), prep$k, prep$constraint)['tbf_ic']},
      tbf_iu    = {tbf(prep$res, prep$hyp_val, floor(nrow(data)/prep$k), prep$k, prep$constraint)['tbf_iu']},
      prodbf_ic = {prod_and_gpbf(prep$res, prep$hyp_val, floor(nrow(data)/prep$k), prep$k, prep$constraint)['prodbf_ic']},
      prodbf_iu = {prod_and_gpbf(prep$res, prep$hyp_val, floor(nrow(data)/prep$k), prep$k, prep$constrain)['prodbf_iu']},
      gpbf_ic   = {prod_and_gpbf(prep$res, prep$hyp_val, floor(nrow(data)/prep$k), prep$k, prep$constrain)['gpbf_ic']},
      gpbf_iu   = {prod_and_gpbf(prep$res, prep$hyp_val, floor(nrow(data)/prep$k), prep$k, prep$constrain)['gpbf_iu']})
  
  return(out)
}

# test for five algorithms and three different constraints
(c(true_effect = cor(M$x, M$y)))
agg_bf(y~x, grouping = 'k', data = M, algorithm = 'prodbf_ic', hypothesis = '>0.2')
agg_bf(y~x, grouping = 'k', data = M, algorithm = 'tbf_ic', hypothesis = '>0.2')
agg_bf(y~x, grouping = 'k', data = M, algorithm = 'gpbf_ic', hypothesis = '<0.2')

#the 'equals' constraint seems to malfunction
agg_bf(y~x, grouping = 'k', data = M, algorithm = 'gpbf_iu', hypothesis = '= 0') 

#also does not seem to work yet for negative hypothesized values
agg_bf(y~x, grouping = 'k', data = M, algorithm = 'prodbf_iu', hypothesis = '<-0.5')


