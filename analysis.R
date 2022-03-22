#dat <- readRDS("tabres.RData")

library(data.table) # load in library
dat <- readRDS(file.path("Sim_Eli", "sim_results_2022-03-19.RData")) # read in results

# some preparation
varsout <- c("allsig", "gpbf_ic", "gpbf_iu", "prodbf_ic", "prodbf_iu", "tbf_ic", "tbf_iu")
varspred <- c("es", "errorsd", "n", "k", "hyp_val")

#########################
# MARGINAL DESCRIPTIVES #
#########################

# First we obtain some descriptives over the marginal conditions.
descriptives <- function(datf){ 
  descs <- as.matrix(
    datf[, lapply(.SD, function(var){
    medians <- median(var) # obtain overall medians per algorithm over medians per condition
    mads <- mad(var)       # obtain Median absolute deviations per conditions
    sums <- sum(var>3)     # obtain sums where the median is > 3
    maxvals <- max(var)    # obtain the highest median
    c(medians, mads, sums, maxvals)
  }), .SDcols = varsout])
  rownames(descs) <- c("median", "mad", "sum>3", "max")
  descs
}

(marginal_descs <- descriptives(dat[ ,lapply(.SD, median), .SDcols = varsout, by = varspred]))
# In theory, the algorithms should have 108/3 = 36 conditions where the median BF > 3. This is because BF > 3 only
# in the conditions where true_es > hyp_val, which is one third of the time.
# The marginal median for all algorithms is < 1, suggesting that the BF in most conditions is indeed < 3.
# The Median Absolute Deviations are also not too severe, giving more proof for reliable estimates.
# Interestingly, the prodbf_ic does have a notable higher mad.
# The finding that the algorithms are rather conservative is further supported by the number of times
# the median BF does exceed 3. Both gpbf algorithms count 0 conditions where median BF > 3, while for allsig this
# number is 1. For the prodbf_iu and tbf algorithms 5 or 6 times is counted and only the prodbf seems to have a 
# more notable number of times the median BF > 3, namely 16.

(es_gt_hyp_val <- descriptives(dat[es == .2 ,lapply(.SD, median), .SDcols = varsout, by = varspred])) # when true BF > 3
all(es_gt_hyp_val['sum>3',] == marginal_descs['sum>3',]) 
# even when we subset for when es > hyp_val no median gets above 3. Also, all cases where the median BF > 3
# was indeed when es > hyp_val.

descriptives(dat[es < .2 ,lapply(.SD, median), .SDcols = varsout, by = varspred]) # when true BF < 3
# also interesting is that the maximum median BF in the cases where es <= hyp_val all are ~1. 
# this provides more support for conservative algorithms.



###############################
# CONFUSION MATRIX STATISTICS #
###############################

# obtain TP, FP, TN and FN for every algorithm in wide format
lapply_at <- function(var, truees) {
  results <- sapply(var, function(var) {
    table(ordered(truees > .1, levels = c("FALSE", "TRUE")), ordered(var > 3, levels = c("FALSE", "TRUE")))
  })
  names(results) <- vapply(names(var), paste, c("TN", "FN", "FP", "TP"), sep = "_", 
                           FUN.VALUE = character(4),
                           USE.NAMES = FALSE)
  as.list(results)
}

# in long format
metrics <- function(datf){
  res <- datf[!es < .1, sapply(.SD, function(var) {
    table(ordered(es > .1, levels = c("FALSE", "TRUE")), ordered(var > 3, levels = c("FALSE", "TRUE")))
  }), .SDcols = varsout]
  rownames(res) <- c("TN", "FN", "FP", "TP")
  return(res)
}

# Example: per condition
metrics(dat[n == 20,])
metrics(dat[n == 80,])
by_n <- dat[!es < .1, lapply_at(.SD, truees = es), .SDcols = varsout, by = c("n")]


# obtain metrics per algorithm()
get_stats <- function(res){ 
  apply(res, 2, function(x){
    attach(as.list(x))
    alpha = FP / (FP + TN) #False positive rate 
    beta = FN / (TP + FN) # False negative rate
    sensitivity = 1-beta
    specificity = 1-alpha
    c(False_positive_rate_alpha = alpha,
      False_negative_rate_beta = beta,
      sensitivity = sensitivity,
      specificity = specificity, 
      pos_lr = sensitivity / (1 - specificity),
      neg_lr = (1 - sensitivity) / specificity)
  })
}


library(ggplot2)
plot_df <-function(stats){ 
  data.frame(Method = colnames(stats),
        Sensitivity = stats["sensitivity", ],
        Specificity = stats["specificity", ],
        FPR = stats["False_positive_rate_alpha", ],
        FNR = stats["False_negative_rate_beta", ],
        Pos_lr = stats["pos_lr", ],
        Neg_lr = stats["neg_lr", ])
}

#create plot df for subset of data
create_df_plot <- function(datf){
  metrics_long <- metrics(datf)
  stats <- get_stats(metrics_long)
  df_plot <- plot_df(stats)
  return(df_plot)
}

# marginal
df_plot <- create_df_plot(dat)
        
ggplot(df_plot, aes(x = Sensitivity, y = Specificity, colour = Method)) +
  geom_point() +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1)) +
  theme_bw()

# conditional on n
ndf <- rbind(create_df_plot(dat[n==20]),
           create_df_plot(dat[n==80]),
           create_df_plot(dat[n==200]),
           create_df_plot(dat[n==500]))
ndf$cond <- c(rep("n = 20", 7), rep("n = 80", 7), rep("n = 200", 7), rep("n = 500", 7))

ggplot(ndf, aes(x = Sensitivity, y = Specificity, colour = Method, shape = cond)) +
  geom_point() +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0.75,1)) +
  theme_bw()


# conditional on errorsd
errorsddf <- rbind(create_df_plot(dat[errorsd == 0]),
             create_df_plot(dat[errorsd == 0.5]),
             create_df_plot(dat[errorsd == 0.81]))
errorsddf$cond <- c(rep("errorsd = 0", 7), rep("errorsd = 0.5", 7), rep("errorsd = 0.81", 7))

ggplot(errorsddf, aes(x = Sensitivity, y = Specificity, colour = Method, shape = cond)) +
  geom_point() +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1)) +
  theme_bw()
