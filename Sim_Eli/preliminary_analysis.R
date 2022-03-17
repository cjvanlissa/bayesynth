########################
# Preliminary Analyses #
########################
library(data.table)
rm(list=ls())

# algorithms we used
algorithms <- c("gpbf", "prodbf", "tbf")
hyps <- c("_ic", "_iu")
alg_names <- paste0(rep(algorithms, each = length(hyps)), hyps)

# some preparation
summarydata <- readRDS("./Sim_Eli/summarydata.RData")
conditions <- c("ndataset", "es", "tau2", "n", "k", "hyp_val")
res <- readRDS(file.path("Sim_Eli", "sim_results_2022-03-14.RData"))
lc <- length(conditions)

# obtain unique conditions
unique_conds <- unique(summarydata[,2:lc])
unique_conds <- t(apply(unique_conds, 1, function(x){
  paste(paste0(colnames(unique_conds), " == ", x), collapse = " & ")
}))

# convert res to datatable and make outcome values numeric
res <- as.data.table(res)
res[, c(lc:ncol(res)):=lapply(.SD, as.numeric), .SDcols=c(lc:ncol(res))]

# we will sample 16 conditions to plot
N <- 16
colors <- c('#c15cf7', '#dda7fa', '#5c8df7', '#acc5fa', '#f5bd3b', '#ffe7b0')

# Sample 16 integers from 1 to length(unique_conds)
set.seed(6164900)
sample_cond <- sample(1:length(unique_conds), 16)


# function to plot histogram of BFs of a condition which is called in the next function
# while subsetting the res datatable.
plotcond <- function(outcome, name, condition, col){
  condition <- unlist(strsplit(condition, "&"))
  condition <- gsub("==", "=", condition)
  condition <- paste(condition, collapse = " ")
  
  propBF3 <- sum(outcome > 3) / length(outcome)
  hist(outcome, main = condition, col = col,
       xlab = paste0("median: ", round(median(outcome),2),
                     ";    sd: ", round(sd(outcome),2), 
                     ";    prop BF > 3: ", round(propBF3,2)),
       cex.lab = 1.3)
}


# Plot histogram s of the calculated BayesFactors from 16 random conditions in a single grid
plot16 <- function(name, line = -1.5, colr){
  par(mfrow = c(4,4))
  for(i in seq_along(sample_cond)){
    res[eval(parse(text=unique_conds[sample_cond[i]])),
        .(plotcond(eval(parse(text = name)), name, unique_conds[sample_cond[i]], col = colr))]
  }
  mtext(paste0("Density of ", name, " for 16 randomly sampled conditions"), side = 3, line = line, outer = T, cex = 1)
}

# Save image to pdf
pdf(file = paste0("Sim_Eli/first_results.pdf"), width = 15, height = 10)
for(i in 1:length(alg_names)){plot16(alg_names[i], colr = colors[i])}
dev.off()
