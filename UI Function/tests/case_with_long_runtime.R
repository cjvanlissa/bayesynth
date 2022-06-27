# I found a weird case where bain() takes a long time to run for certain hypothese, while it should not.
# weirdly enough, it seems to have to do with the column index in which certain variables exist in a df


# ttest
set.seed(100)
tt = as.data.frame(cbind(y = rnorm(1000,0,1), x = rnorm(1000, 0.2,1)))
ttest <- t_test(tt$y, tt$x)

# Bain seems run a long time for directional hypotheses for y in tt of the form a < y < b, but this works fine with x...

# it works fine with x
bain(ttest, hypothesis = "-2<x<0") # almost instant


# also for y, when only 1 constraint is imposed
bain(ttest, hypothesis = "y<0")


# but takes very long time for y if 2 constraints are imposed
bain(ttest, hypothesis = "-2<y<0")


# weirdly enough, it has to do with in what order x and y exist in the dataframe  tt'
# if we change their order it works fine for y, but not for x
set.seed(100)
tt = as.data.frame(cbind(x = rnorm(1000,0,1), y = rnorm(1000, 0.2,1))) # Note I changed order of x and y
ttest <- t_test(tt$y, tt$x)

# now x takes long time
bain(ttest, hypothesis = "-2<x<0")


#while y does not
bain(ttest, hypothesis = "-2<y<0")


