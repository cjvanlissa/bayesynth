#dat <- readRDS("tabres.RData")
library(data.table)
dat <- readRDS(file.path("Sim", "sim_results_2022-05-08.RData"))
lapply_at <- function(var, truees) {
  results <- sapply(var, function(var) {
    table(ordered(truees > .1, levels = c("FALSE", "TRUE")), ordered(var > 3, levels = c("FALSE", "TRUE")))
  })
  names(results) <- vapply(names(var), paste, c("TN", "FN", "FP", "TP"), sep = "_", 
                           FUN.VALUE = character(4),
                           USE.NAMES = FALSE)
  as.list(results)
}
varsout <- c("ipd", "rma", "allsig", "gpbf_ic", "gpbf_iu", "prodbf_ic", "prodbf_iu", "tbf_ic", "tbf_iu")
varspred <- c("errorsd", "n", "k", "hyp_val")

# Example: per condition
dat[!es < .1, lapply_at(.SD, truees = es), .SDcols = varsout, by = c("n", "errorsd")]

# Overall
#res <- dat[!es < .1, lapply_at(.SD, truees = es), .SDcols = varsout]
res <- dat[!es < .1, sapply(.SD, function(var) {
  table(ordered(es > .1, levels = c("FALSE", "TRUE")), ordered(var > 3, levels = c("FALSE", "TRUE")))
}), .SDcols = varsout]
rownames(res) <- c("TN", "FN", "FP", "TP")

tabres <- dat[!es < .1,  lapply(.SD, function(var) {
  as.vector(table(ordered(es > .1, levels = c("FALSE", "TRUE")), ordered(var > 3, levels = c("FALSE", "TRUE"))))
}), .SDcols = varsout, by = varspred]
tabres[, metric := rep(c("TN", "FN", "FP", "TP"), nrow(tabres)/4)]

tabres <- melt(tabres, measure.vars = varsout,
             variable.name = "alg", value.name = "val")
tabres <- dcast(tabres, errorsd + n + k + hyp_val + alg ~ metric, value.var = "val")
tabres[, alpha := FP / (FP + TN)]
tabres[, beta := FN / (TP + FN)]
tabres[, sensitivity := 1-beta]
tabres[, specificity := 1-alpha]
tabres[, accuracy := (TP + TN) / (TP+TN+FP+FN)]

saveRDS(tabres, "confusion_by_cond.RData")


# ANOVAs ------------------------------------------------------------------
source("analysis_functions.R")
df_acc <- dcast(tabres, errorsd + n + k + hyp_val ~ alg, value.var = "accuracy")
#the dependent variables (the test r2)
yvars <- varsout[!endsWith(varsout, "iu")]
conditions <- c("n", "k", "errorsd")
#creates a list for every Anova with the results for the anova and the effect sizes for the conditions on the algorithms
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, '~', paste(conditions, collapse = "+"))
  # form<-paste(yvar, '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2") #the ^2 signifies that we want all possible effects up until interactions
  thisaov<-aov(as.formula(form), data=df_acc) #change data according to dataframe you are using
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})

# Anova for the difference ------------------------------------------------

comps <- expand.grid(yvars, yvars)
comps <- comps[!comps$Var1 == comps$Var2, ]
comps <- t(apply(comps, 1, sort))
comps <- comps[!duplicated(comps), ]
#creates a list for every Anova with the results for the anova and the effect sizes for the conditions on the algorithms
diffanovas <- sapply(1:nrow(comps), function(i){
  form<-as.formula(paste("r2", '~ algo * (', paste(conditions, collapse = "+"), ")")) #the ^2 signifies that we want all possible effects up until interactions
  # form<-as.formula(paste("r2", '~ algo * ((', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2)")) #the ^2 signifies that we want all possible effects up until interactions
  tmp <- df_acc
  tmp <- tmp[, .SD, .SDcols = c(names(tmp)[1:3], comps[i, , drop = TRUE])]
  names(tmp) <- gsub("^(.+?)_r2", "r2_\\1", names(tmp))
  tmp = melt(tmp, id.vars = names(tmp)[1:3],
             measure.vars = names(tmp)[4:5],
             variable.name = "algo",
             value.name = "r2")
  thisaov<-aov(form, data=tmp) #change data according to dataframe you are using
  thisetasq<-EtaSq(thisaov)[ , 2]
  thisetasq <- thisetasq[startsWith(names(thisetasq), "algo")]
  thisetasq
})
colnames(diffanovas) <- paste0(comps[,1], " vs. ", comps[,2])
diffanovas <- data.frame(condition = gsub("algo:", "", rownames(diffanovas), fixed = T), diffanovas)
diffanovas$condition <- trimws(diffanovas$condition)
out <- list(difference = diffanovas[1, ])
#creates dataframe with effect sizes for all conditions for all algorithms
etasqs<-data.frame(sapply(anovas, `[[`, 2))
colnames(etasqs) <- yvars
etasqs$condition <- trimws(rownames(etasqs))
etasqs <- merge(etasqs, diffanovas, by = "condition", all.x = TRUE)

write.csv(etasqs, "effect_of_conditions.csv", row.names = FALSE)


stats <- apply(res, 2, function(x){
  attach(as.list(x))
  alpha = FP / (FP + TN) #False positive rate 
  beta = FN / (TP + FN) # False negative rate
  sensitivity = 1-beta
  specificity = 1-alpha
  accuracy = (TP + TN) / sum(x)
  c(False_positive_rate_alpha = alpha,
    False_negative_rate_beta = beta,
    sensitivity = sensitivity,
    specificity = specificity, 
    accuracy = accuracy,
    pos_lr = sensitivity / (1 - specificity),
    neg_lr = (1 - sensitivity) / specificity)
})
stats <- stats[, order(stats["accuracy",], decreasing = T)]
write.csv(stats, "confusion.csv")

library(ggplot2)
df_plot <- data.frame(Method = colnames(stats),
                      Sensitivity = stats["sensitivity", ],
                      Specificity = stats["specificity", ])
df_plot$Type <- "Frequentist"
df_plot$Type[endsWith(df_plot$Method, "iu")] <- "BF unconstrained"
df_plot$Type[endsWith(df_plot$Method, "ic")] <- "BF complement"
df_plot <- df_plot[!endsWith(df_plot$Method, "iu"), ]
ggplot(df_plot, aes(x = Sensitivity, y = Specificity, colour = Method, fill = Method, shape = Type)) +
  geom_point() +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1)) +
  geom_vline(xintercept = .8, linetype = 2) +
  geom_hline(yintercept = .95, linetype = 2) +
  theme_bw()