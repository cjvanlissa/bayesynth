#dat <- readRDS("tabres.RData")
dat <- readRDS(file.path("Sim_Eli", "sim_results_2022-03-19.RData"))
lapply_at <- function(var, truees) {
  results <- sapply(var, function(var) {
    table(ordered(truees > .1, levels = c("FALSE", "TRUE")), ordered(var > 3, levels = c("FALSE", "TRUE")))
  })
  names(results) <- vapply(names(var), paste, c("TN", "FN", "FP", "TP"), sep = "_", 
                           FUN.VALUE = character(4),
                           USE.NAMES = FALSE)
  as.list(results)
}
varsout <- c("allsig", "gpbf_ic", "gpbf_iu", "prodbf_ic", "prodbf_iu", "tbf_ic", "tbf_iu")
varspred <- c("es", "errorsd", "n", "k", "hyp_val")

# Example: per condition
dat[!es < .1, lapply_at(.SD, truees = es), .SDcols = varsout, by = c("n", "errorsd")]

# Overall
#res <- dat[!es < .1, lapply_at(.SD, truees = es), .SDcols = varsout]
res <- dat[!es < .1, sapply(.SD, function(var) {
  table(ordered(es > .1, levels = c("FALSE", "TRUE")), ordered(var > 3, levels = c("FALSE", "TRUE")))
}), .SDcols = varsout]
rownames(res) <- c("TN", "FN", "FP", "TP")

stats <- apply(res, 2, function(x){
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

library(ggplot2)
df_plot <- data.frame(Method = colnames(stats),
        Sensitivity = stats["sensitivity", ],
        Specificity = stats["specificity", ])
        
ggplot(df_plot, aes(x = Sensitivity, y = Specificity, colour = Method, shape = Method)) +
  geom_point() +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1)) +
  theme_bw()