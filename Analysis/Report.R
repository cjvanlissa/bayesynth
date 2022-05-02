## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
library(data.table) # load in libraries
library(DT)
library(ggplot2)
dat <- readRDS(file.path("Sim_Eli", "sim_results_2022-03-19.RData")) #
dat <- dat[!es < .1]


## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
varsout <- c("allsig", "gpbf_ic", "gpbf_iu", "prodbf_ic", "prodbf_iu", "tbf_ic", "tbf_iu")
varspred <- c("es", "errorsd", "n", "k", "hyp_val")


## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
descriptives <- function(datf){ 
  descs <- as.matrix(
    datf[, lapply(.SD, function(var){
    quants <- quantile(var, c(0.25, 0.50, 0.75, 0.90, 0.95))# obtain overall quantiles
    mads <- mad(var)       # obtain Median absolute deviations per conditions
    sums <- sum(var>3)     # obtain sums where the median is > 3
    maxvals <- max(var)    # obtain the highest median
    minvals <- min(var)    # obtain lowest median
    c(quants[1], quants[2], quants[3], quants[4], quants[5], mads, sums, maxvals, minvals)
  }), .SDcols = varsout])
  rownames(descs) <- c('25 %', '50 %',  '75 %','90 %', '95 %', "mad", "sum>3", "max", "min")
  descs
}
marginal_descs <- descriptives(dat[ ,lapply(.SD, median), .SDcols = varsout, by = varspred])
datatable(round(marginal_descs,5), options = list(dom = 't'))


## --------------------------------------------------------------------------------------------------------------------------------------------------
all.equal(dat$tbf_iu, dat$prodbf_iu)
cor(dat$tbf_iu, dat$prodbf_iu)


## --------------------------------------------------------------------------------------------------------------------------------------------------
dat <- dat[,tbf_iu:=NULL]
varsout <- c("allsig", "gpbf_ic", "gpbf_iu", "prodbf_ic", "prodbf_iu", "tbf_ic")


## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
es_gt_hyp_val <- descriptives(dat[es == .2 ,lapply(.SD, median), .SDcols = varsout, by = varspred])
datatable(round(es_gt_hyp_val,5),options = list(dom = "t"))


## --------------------------------------------------------------------------------------------------------------------------------------------------
all(es_gt_hyp_val['sum>3',varsout] == marginal_descs['sum>3',varsout]) 


## ---- echo =F--------------------------------------------------------------------------------------------------------------------------------------
es_st_hyp_val <- descriptives(dat[es == .1 ,lapply(.SD, median), .SDcols = varsout, by = varspred])
datatable(round(es_st_hyp_val,5), options = list(dom = 't'))


## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
prodbf_gt_3 <- dat[,lapply(.SD, median), .SDcols = varsout, by = varspred][prodbf_ic > 3][order(k, prodbf_ic, decreasing = T)]
datatable(round(prodbf_gt_3,3), options = list(pageLength = 5))


## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
# obtain confusion matrix metrics with a dimension of 1 by alg*metrics (so 1x28 in this datatable)
lapply_at <- function(var, truees) {
  results <- sapply(var, function(var) {
    table(ordered(truees > .1, levels = c("FALSE", "TRUE")), ordered(var > 3, levels = c("FALSE", "TRUE")))
  })
  names(results) <- vapply(names(var), paste, c("TN", "FN", "FP", "TP"), sep = "_", 
                           FUN.VALUE = character(4),
                           USE.NAMES = FALSE)
  as.list(results)
}

# in long format gives a 4x7 matrix with 4 metrics and 7 algorithms
metrics <- function(datf){
  res <- datf[, sapply(.SD, function(var) {
    table(ordered(es > .1, levels = c("FALSE", "TRUE")), ordered(var > 3, levels = c("FALSE", "TRUE")))
  }), .SDcols = varsout]
  rownames(res) <- c("TN", "FN", "FP", "TP")
  return(res)
}

# obtain metrics per algorithm based on FP, TN, FN and TP() and collect in a matrix
get_stats <- function(res){ 
  suppressMessages(
  apply(res, 2, function(x){
    attach(as.list(x))
    alpha = FP / (FP + TN) #False positive rate 
    beta = FN / (TP + FN) # False negative rate
    sensitivity = 1-beta
    specificity = 1-alpha
    correct = (TN + TP) /(FP + FN + TN + TP)
    pos_lr = sensitivity / (1 - specificity)
    neg_lr = (1 - sensitivity) / specificity
    c(correct = correct,
      False_positive_rate_alpha = alpha,
      False_negative_rate_beta = beta,
      sensitivity = sensitivity,
      specificity = specificity, 
      pos_lr = ifelse(is.na(pos_lr), 0, pos_lr),
      neg_lr = ifelse(is.na(neg_lr), 0, neg_lr))
  })
  )
}

# create df from stats matrix for plotting purposes
plot_df <-function(stats){ 
  data.frame(Method = colnames(stats),
        correct = stats["correct", ],
        Sensitivity = stats["sensitivity", ],
        Specificity = stats["specificity", ],
        FPR = stats["False_positive_rate_alpha", ],
        FNR = stats["False_negative_rate_beta", ],
        Pos_lr = stats["pos_lr", ],
        Neg_lr = stats["neg_lr", ])
}

# all in one function. Create plot df for a subset of the data table.
create_df_plot <- function(datf){
  metrics_long <- metrics(datf)
  stats <- get_stats(metrics_long)
  df_plot <- plot_df(stats)
  return(df_plot)
}

#set colours to use for different algorithms
colours <- c("black", "orange" , "red", "#79aaf7", "blue", "#00a62f" )

# obtain numeric values each conditions take on
cond_num <- function(cond){return(as.numeric(names(table(dat[,..cond]))))}

#obtain condition names. For k it returns c("k = 1", "k = 3", "k = 10")
cond_names <- function(cond){return(paste0(cond, " = ",names(table(dat[,..cond]))))}


## ---- echo = F, message=F--------------------------------------------------------------------------------------------------------------------------
df_plot <- create_df_plot(dat)
marginal_metrics <- rbind(df_plot[,2:ncol(df_plot)], 
      mean = apply(df_plot[,2:ncol(df_plot)], 2 ,mean))
datatable(round(marginal_metrics, 5), options = list(dom = 't'))


## --------------------------------------------------------------------------------------------------------------------------------------------------
all(dat$gpbf_iu < 3)


## --------------------------------------------------------------------------------------------------------------------------------------------------
dat <- dat[,gpbf_iu:=NULL]
varsout <- c("allsig", "gpbf_ic", "prodbf_ic", "prodbf_iu", "tbf_ic")


## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
sd_df <- do.call(rbind, lapply(cond_num('errorsd'), function(x){create_df_plot(dat[errorsd == x])})) 
sd_df$condition <- rep(cond_names('errorsd'), each = length(varsout))

ggplot(sd_df, aes(x = Sensitivity, y = Specificity, colour = Method, shape = condition)) +
  geom_point(size = 2.5) +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0.85,1)) +
  scale_color_manual(values = colours) +
  theme_bw()


## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
n_df <- do.call(rbind, lapply(cond_num('n'), function(x){create_df_plot(dat[n == x])})) 
n_df$condition <- rep(cond_names('n'), each = length(varsout))

ggplot(n_df, aes(x = Sensitivity, y = Specificity, colour = Method, shape = condition)) +
  geom_point(size = 2.5) +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0.80,1)) +
  scale_color_manual(values = colours) +
  theme_bw()


## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
k_df <- do.call(rbind, lapply(cond_num('k'), function(x){create_df_plot(dat[k == x])})) 
k_df$condition <- rep(cond_names('k'), each = length(varsout))

ggplot(k_df, aes(x = Sensitivity, y = Specificity, colour = Method, shape = condition)) +
  geom_point(size = 2.5) +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0.85,1)) +
  scale_color_manual(values = colours) +
  theme_bw()

# THIS MUST RUN AGAIN WITH NEW DATA!
## ---- eval = F, echo = F---------------------------------------------------------------------------------------------------------------------------
## kn_df <- do.call(rbind, lapply(cond_num('k'), function(x){
##   do.call(rbind, lapply(cond_num('n'), function(y){
##     temp <- create_df_plot(dat[k == x & n == y])
##     temp$condition1 <- paste0('k = ', x)
##     temp$condition2 <- paste0('n = ', y)
##     temp
##     }))
##   }))
## 
## ke_df <- do.call(rbind, lapply(cond_num('k'), function(x){
##   do.call(rbind, lapply(cond_num('errorsd'), function(y){
##     temp <- create_df_plot(dat[k == x & errorsd == y])
##     temp$condition1 <- paste0('k = ', x)
##     temp$condition2 <- paste0('errorsd = ', y)
##     temp
##     }))
##   }))
## 
## en_df <- do.call(rbind, lapply(cond_num('errorsd'), function(x){
##   do.call(rbind, lapply(cond_num('n'), function(y){
##     temp <- create_df_plot(dat[errorsd == x & n == y])
##     temp$condition1 <- paste0('errorsd = ', x)
##     temp$condition2 <- paste0('n = ', y)
##     temp
##     }))
##   }))
## 
## saveRDS(kn_df, file = "./Analysis/plot_kn.RData")
## saveRDS(ke_df, file = "./Analysis/plot_ke.RData")
## saveRDS(en_df, file = "./Analysis/plot_en.RData")
## 


## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
kn_df <- readRDS("./Analysis/plot_kn.RData")
ke_df <- readRDS("./Analysis/plot_ke.RData")
en_df <- readRDS("./Analysis/plot_en.RData")



## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
ggplot(kn_df, aes(x = Sensitivity, y = Specificity, colour = Method, shape = condition1)) +
  geom_point(size = 2.5) +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0.70,1)) +
  scale_color_manual(values = colours) +
  facet_wrap(~condition2) +
  theme_bw()




## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
ggplot(ke_df, aes(x = Sensitivity, y = Specificity, colour = Method, shape = condition1)) +
  geom_point(size = 2.5) +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0.75,1)) +
  scale_color_manual(values = colours) +
  facet_wrap(~condition2) +
  theme_bw()




## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
ggplot(en_df, aes(x = Sensitivity, y = Specificity, colour = Method, shape = condition1)) +
  geom_point(size = 2.5) +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0.6,1)) +
  scale_color_manual(values = colours) +
  facet_wrap(~condition2) +
  theme_bw()




## ---- echo = F-------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
gc(full = T, reset = T)
cat("\014")

