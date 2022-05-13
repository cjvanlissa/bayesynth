#clear console and env
rm(list=ls(all.names = T))
cat("\014")

# Delete folder where previous results are in and create new 'Results' folder
unlink("results", recursive = T)
dir.create("results")

#load necessary packages
dependencies <- c('MASS', 'bain', 'metafor', 'lme4', 'mvtnorm')
lapply(dependencies, function(x){library(x, character.only = T)})

# Check package versions
versions <- c(
  compareVersion(as.character(packageVersion("bain")), "0.2.8"),
  compareVersion(as.character(packageVersion("MASS")), "7.3.55"))
if(!all(versions == 0)) stop("Using the incorrect version of one or more packages.")

# Load simulation functions from source -----------------------------------
# source('sim/functions.R')

# set conditions for simulation
hyper_parameters<-list(
  ndataset = 1:1000,            # number of replications per condition
  es = c(0.1, .2),        # true effect size = true correlation with outcome
  errorsd = c(0.81, .5, 0),   # corresponds to reliability of .6, .8, and 1
  n = c(20, 80, 200, 500),             # mean sample size per group
  k = c(2, 3, 10),               # number of groups
  hyp_val = c(0.1)    # thresholds for informative hypotheses
)

# Create hypergrid with simulation parameters and save it as .RData file extension
summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)
set.seed(6164900)
summarydata$seed <- sample(1:.Machine$integer.max, nrow(summarydata))
saveRDS(summarydata, file = "./sim/summarydata.RData")
# summarydata<-readRDS("./sim/summarydata.RData")

# prepare parallel processing
library(doSNOW)
nclust <- parallel::detectCores() 
cl <- makeCluster(nclust) 
registerDoSNOW(cl) 

# add progression bar
pb <- txtProgressBar(min = 0, max = nrow(summarydata), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

# run simulation
tab <- foreach(rownum = 1:nrow(summarydata), .options.snow = opts, .packages = c("bain", "mvtnorm", "lme4", "metafor"), .combine = rbind) %dopar% {
  # Set seed
  attach(summarydata[rownum, ])
  set.seed(seed)
  
  # create dataframe for each group
  dfs <- lapply(1:k, function(i){
    df <- rmvnorm(n, sigma = matrix(c(1, es, es, 1), nrow = 2))
    df + rnorm(2*n, sd = errorsd)
  })
  
  # IPD meta
  df_join <- scale(do.call(rbind, dfs))
  df_join <- data.frame(cbind(df_join, rep(1:k, each = n)))
  res_ipd <- lmer(formula = X1 ~ X2 + (1|X3),
                             data    = df_join) #to run the model
  ci_ipd <- confint(res_ipd, level = .90)
  
  # obtain estimates for correlations and their standard errors for every dataset
  res <- sapply(dfs, function(x){
    est <- cor(x)[2,1]                        #estimate for correlation between x and y
    se_est <- sqrt((1 - est^2)/(n - 2)) #estimate for standard error = sqrt((1-r^2)/df) with df = N - 2
    c(est, se_est)               #return estimate, se and sample size of particular set
  })
  
  df_meta <- data.frame(t(res))
  
  res_met <- rma(yi = df_meta[,1], sei = df_meta[,2], level = 90)

  # Classic approach
  classic_votecount <- round(mean(apply(res, 2, function(x){
    (x[1]-hyp_val)/x[2] > 1.644854
  })))
  
  # necessary naming for bain and further preparing
  colnames(res) <- paste0('r', 1:k)
  sig <- lapply(res[2,], function(x){matrix(x^2)})    # make list of covariance matrices for the datasets, make sure to square the standard errors
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
                      Sigma = diag(res[2,]^2, ncol = ncol(res))) # assume independence between groups, square standard errors
  
  # returns in order: gpbf_ic, gpbf_iu, prodbf_ic, prodbf_iu, tbf_ic, tbf_iu
  fits <- c(rownum,
    c(0,10)[(ci_ipd[4,1] > hyp_val) +1], #return 10 if IPD is significantly greater than hyp, 0 if false
    c(0,10)[(res_met$ci.lb > hyp_val) +1], #return 10 if RMA is significantly greater than hyp, 0 if false
    c(0,10)[classic_votecount+1], #return 10 if classic_votecount = T, 0 if false
    gp_and_prod[1,],
    gp_and_prod[2,], 
    c(bf_together$fit$BF.c[1], bf_together$fit$BF.u[1]))
  write.table(x = t(fits), file = sprintf("results/results_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
  NULL
}

#Close cluster
stopCluster(cl)


# End of simulation -------------------------------------------------------
stop("End of simulation")

# Merge files -------------------------------------------------------------
library(data.table)

# read in the simulation conditions 
res <- readRDS("./sim/summarydata.RData")
setDT(res)

f <- list.files("results/", full.names = TRUE)

tab <- rbindlist(lapply(f, fread))
setorderv(tab, cols = "V1", order=1L, na.last=FALSE)

# algorithms (geometric product Bayes Factor, product Bayes Factor, together Bayes Factor)
algorithms <- c("gpbf", "prodbf", "tbf")
hyps <- c("_ic", "_iu")  # using both complementary and unconstrained
alg_names <- c("ipd", "rma", "votecount", paste0(rep(algorithms, each = length(hyps)), hyps))

conditions <- colnames(res)

# make sure results are same length as conditions
if(!(tab[1,1] == 1 & tab[nrow(tab), 1] == nrow(res) & length(unique(tab$V1)) == nrow(res))){
  stop("Results not the same length as number of simulation iterations")
}

# give appropriate names to the simulation results and omit identification variable 'V1'
names(tab) <- c("V1", alg_names)
tab[, "V1" := NULL]

# cbind conditions and results
res<- cbind(res, tab)
rm(tab)

# write results to .RData and .csv extension and delete .txt files in the results folder.
fwrite(res, file.path("sim", paste0("sim_results_", Sys.Date(), ".csv")))
saveRDS(res, file.path("sim", paste0("sim_results_", Sys.Date(), ".RData")))

# END OF FILE
tabres <- res[, lapply(.SD, function(x){mean(x > 3)}), .SDcols = alg_names, by = c("es", "errorsd", "n", "k")]
write.csv(tabres, "tabres.csv", row.names = FALSE)
saveRDS(tabres, "tabres.RData")
#colMeans(res[, .SD > 3, .SDcols = alg_names, .gro])