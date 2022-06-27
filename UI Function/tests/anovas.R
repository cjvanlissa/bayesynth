


sesamesim$site <- as.factor(sesamesim$site)

anov <- lm(postnumb~site-1,sesamesim[1:75,])
results <- bain(anov, "site2=site1=site3=site4=site5; site2>site1=site3=site4=site5; site2>site5>site1>site3>site4")
anov <- lm(postnumb~site-1,sesamesim[76:150,])
results2 <- bain(anov, "site2=site1=site3=site4=site5; site2=site1=site3>site4=site5; site2>site5>site1>site3>site4")
anov <- lm(postnumb~site-1,sesamesim[151:nrow(sesamesim),])
results3 <- bain(anov, "site2=site1=site3=site4=site5; site2=site1=site3=site4>site5; site2>site5>site1>site3>site4")

res_base <- pbf(list(results, results2, results3))

anov <- lm(postnumb~site-1,sesamesim[1:75,])
results <- bain(anov, "site2=site1>site3=site4=site5; site2>site1=site3=site4=site5; site2>site5>site1>site3>site4")
anov <- lm(postnumb~site-1,sesamesim[76:150,])
results2 <- bain(anov, "site2=site1=site3=site4=site5; site2=site1=site3>site4=site5; site2>site5>site1>site3>site4")
anov <- lm(postnumb~site-1,sesamesim[151:nrow(sesamesim),])
results3 <- bain(anov, "site2=site1=site3=site4=site5; site2=site1=site3=site4>site5; site2>site5>site1>site3>site4")

pbf(list(results, results2, results3))

anov <- lm(postnumb~site-1,sesamesim[1:75,])
results <- bain(anov, "site2=site1=site3=site4=site5; site2>site1=site3=site4=site5; site2>site5>site1>site3>site4")
anov <- lm(postnumb~site-1,sesamesim[76:150,])
results2 <- bain(anov, "site2=site1=site3>site4=site5; site2=site1=site3=site4=site5; site2>site5>site1>site3>site4")
anov <- lm(postnumb~site-1,sesamesim[151:nrow(sesamesim),])
results3 <- bain(anov, "site2=site1=site3=site4>site5; site2>site5>site1>site3>site4; site2=site1=site3=site4=site5")

res_changeorder <- pbf(list(results, results2, results3))

res_base
res_changeorder