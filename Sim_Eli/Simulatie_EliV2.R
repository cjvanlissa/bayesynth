#clear console and env
rm(list=ls(all.names = T))
cat("\014")

#load necessary packages
dependencies <- c('MASS', 'bain')
lapply(dependencies, function(x){library(x, character.only = T)})

# Check package versions
versions <- c(
  compareVersion(as.character(packageVersion("bain")), "0.2.8"),
  compareVersion(as.character(packageVersion("MASS")), "7.3.55"))
if(!all(versions == 0)) stop("Using the incorrect version of one or more packages.")

# Load simulation functions from source -----------------------------------
source('Sim_Eli/functions.R')

# set conditions for simulation
hyper_parameters<-list(
  ndataset = 1:2,            # number of replications per condition
  es = c(0, 0.1),        # true effect size = true correlation with outcome
  tau2 = c(0, .1),            # reliability of outcome (standard deviation of outcome) (NOTE)
  n = c(30, 300),             # mean sample size per group
  k = c(3, 10),               # number of groups
  hyp_val = c(0, 0.1)    # thresholds for informative hypotheses
)

# Create hypergrid with simulation parameters and save it as .RData file extension
summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)
saveRDS(summarydata, file = "./Sim_Eli/summarydata.RData")
# summarydata<-readRDS("./Sim_Eli/summarydata.RData")


# Delete folder where previous results are in and create new 'Results' folder
unlink("./sim_Eli/Results", recursive = T)
dir.create("./Sim_Eli/Results")

# number of algorithms used
nalg <- 6

# prepare parallel processing
library(doSNOW)
nclust <- parallel::detectCores() 
cl <- makeCluster(nclust) 
registerDoSNOW(cl) 

# add progression bar
pb <- txtProgressBar(min = 0, max = nrow(summarydata), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

# export nalg and vector of seeds to all clusters
set.seed(6164900)
seeds <- sample(1:.Machine$integer.max, nrow(summarydata))
clusterExport(cl, c('nalg', 'seeds'))

# run simulation
foreach(rownum = 1:nrow(summarydata), .options.snow = opts, .packages = c("bain", "MASS")) %dopar% {
  
  set.seed(seeds[rownum])
  res <-  tryCatch(
    c(rownum, 
      BFs(es = summarydata[rownum, "es"],
             n = summarydata[rownum, "n"],
             hyp_val = summarydata[rownum, "hyp_val"],
             k =summarydata[rownum, "k"],
             tau2 = summarydata[rownum, "tau2"])),
          error = function(e){c(rownum, rep(NA, nalg))}
    )
  
  # every cluster writes its own .txt file where the results for every iteration are stored
  write.table(x = t(res), file = sprintf("./Sim_Eli/Results/results%d.txt" , Sys.getpid()), 
              sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}

#Close cluster
stopCluster(cl)


# End of simulation -------------------------------------------------------
stop("End of simulation")

# Merge files -------------------------------------------------------------
library(data.table)

# algorithms (geometric product Bayes Factor, product Bayes Factor, together Bayes Factor)
algorithms <- c("gpbf", "prodbf", "tbf")
hyps <- c("_ic", "_iu")  # using both complementary and unconstrained
alg_names <- paste0(rep(algorithms, each = length(hyps)), hyps)

# read in the simulation conditions 
res <- readRDS("./Sim_Eli/summarydata.RData")
conditions <- colnames(res)

# read in the simulation results
f <- list.files(path = "./Sim_Eli/Results", pattern = "result", full.names = T)
tab <- setorderv(
  do.call(rbind, lapply(f, function(file){
    fread(file)
  })),
  cols = "V1", order = 1L, na.last = F)

# make sure results are same length as conditions
if(!(tab$V1[1] == 1 & tail(tab$V1, 1) == nrow(res) & length(unique(tab$V1)) == nrow(res))){
  stop("Results not the same length as number of simulation iterations")
}

# give appropriate names to the simulation results and omit identification variable 'V1'
names(tab) <- c("V1", alg_names)
tab[, "V1" := NULL]

# cbind conditions and results
res<- cbind(res, tab)
rm(tab)

# write results to .RData and .csv extension and delete .txt files in the results folder.
fwrite(res, file.path("Sim_Eli", paste0("sim_results_", Sys.Date(), ".csv")))
saveRDS(res, file.path("Sim_Eli", paste0("sim_results_", Sys.Date(), ".RData")))
f <- list.files("./Sim_Eli/Results", full.names = TRUE)
file.remove(f)

# END OF FILE