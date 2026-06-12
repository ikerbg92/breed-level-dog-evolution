library(orthopolynom)
#Package version 2.29 of MCMCglmm is required to be able to include both the 
#estimate of common ancestry and haplotype-sharing in the analyses (following Garamszegi et al. 2020)
library(MCMCglmm)
library(PerformanceAnalytics)
library(dplyr) 
library(jtools)

#autocorrelation function to visualize model results
plot.acfs <- function(x) {
  n <- dim(x)[2]
  par(mfrow=c(1,1))
  for (i in 1:n) {
    acf(x[,i], 
        lag.max=100, 
        main=colnames(x)[i])
    grid()}}

#prior for models, based on Hadfield (2010)
prior1<-list(
  G=list(
    G1=list(V=1,nu=1),
    G2=list(V=1,nu=1)),
  R=list(V=1,nu=0.002)) 

#TABLE 1 OUTPUT:

#Data: life history traits for 428 dog breeds
data<-read.csv("data.csv", header=TRUE)

#Shared ancestry and haplotype sharing data for all breeds for which we have data for lifespan (dfl)
snpsdfl<-read.csv("snpsdfl.csv", row.names=1, header=TRUE)
hapdfl<-read.csv("hapdfl.csv", row.names=1, header=TRUE)

#Lifespan 
dfl<-subset(data, select=c(breed, abrevGaramszegi2019, lifespan))
dfl<- na.omit(dfl)

#To include genomic data as a random factor in the MCMCglmm model:
snpsdfl<-data.frame(sapply(snpsdfl, function(x) as.numeric(as.character(x))))
Snpsvd1 = svd(as.matrix(snpsdfl))
Snpsvd1 = Snpsvd1$v%*%(t(Snpsvd1$u)*sqrt(Snpsvd1$d))
dfl$phylo = Snpsvd1
hapdfl<-data.frame(sapply(hapdfl, function(x) as.numeric(as.character(x))))
Hapsvd2 = svd(as.matrix(hapdfl))
Hapsvd2 = Hapsvd2$v%*%(t(Hapsvd2$u)*sqrt(Hapsvd2$d))
dfl$haplo = Hapsvd2
row.names(dfl)<-dfl$abrevGaramszegi2019
dfl$phylo = Snpsvd1
dfl$haplo = Hapsvd2

#Analysis of the influence of shared ancestry and recent admixture on the among-breed variation in lifespan
lifespan1=MCMCglmm(lifespan ~ 1, 
                   random=~idv(phylo)+idv(haplo), 
                   data=dfl, 
                   family="gaussian", 
                   prior=prior1,
                   thin=50,
                   burnin=70000,
                   nitt=400000,
                   verbose=T)
plot(lifespan1)
summary(lifespan1)
autocorr.diag(lifespan1$VCV)
autocorr.diag(lifespan1$Sol)
plot.acfs(lifespan1$VCV)

#Shared ancestry and haplotype sharing data for all breeds for which we have data for adult weight (dfw)
snpsdfw<-read.csv("snpsdfw.csv", row.names=1, header=TRUE)
hapdfw<-read.csv("hapdfw.csv", row.names=1, header=TRUE)

#Adult weight 
dfw<-subset(data, select=c(breed, abrevGaramszegi2019, weight))
dfw<- na.omit(dfw)

#To include genomic data as a random factor in the MCMCglmm model:
snpsdfw<-data.frame(sapply(snpsdfw, function(x) as.numeric(as.character(x))))
Snpsvd1 = svd(as.matrix(snpsdfw))
Snpsvd1 = Snpsvd1$v%*%(t(Snpsvd1$u)*sqrt(Snpsvd1$d))
dfw$phylo = Snpsvd1
hapdfw<-data.frame(sapply(hapdfw, function(x) as.numeric(as.character(x))))
Hapsvd2 = svd(as.matrix(hapdfw))
Hapsvd2 = Hapsvd2$v%*%(t(Hapsvd2$u)*sqrt(Hapsvd2$d))
dfw$haplo = Hapsvd2
row.names(dfw)<-dfw$abrevGaramszegi2019
dfw$phylo = Snpsvd1
dfw$haplo = Hapsvd2

#Analysis of the influence of shared ancestry and recent admixture on the among-breed variation in adult weight
weight1=MCMCglmm(sqrt(weight) ~ 1, 
                 random=~idv(phylo)+idv(haplo), 
                 data=dfw, 
                 family="gaussian", 
                 prior=prior1,
                 thin=50,
                 burnin=70000,
                 nitt=400000,
                 verbose=T)
plot(weight1)
summary(weight1)
autocorr.diag(weight1$VCV)
autocorr.diag(weight1$Sol)
plot.acfs(weight1$VCV)

#Shared ancestry and haplotype sharing data for all breeds for which we have data for litter size (dfls)
snpsdfls<-read.csv("snpsdfls.csv", row.names=1, header=TRUE)
hapdfls<-read.csv("hapdfls.csv", row.names=1, header=TRUE)

#Litter size
dfls<-subset(data, select=c(breed, abrevGaramszegi2019, littersize))
dfls<- na.omit(dfls)

#To include genomic data as a random factor in the MCMCglmm model:
snpsdfls<-data.frame(sapply(snpsdfls, function(x) as.numeric(as.character(x))))
Snpsvd1 = svd(as.matrix(snpsdfls))
Snpsvd1 = Snpsvd1$v%*%(t(Snpsvd1$u)*sqrt(Snpsvd1$d))
dfls$phylo = Snpsvd1
hapdfls<-data.frame(sapply(hapdfls, function(x) as.numeric(as.character(x))))
Hapsvd2 = svd(as.matrix(hapdfls))
Hapsvd2 = Hapsvd2$v%*%(t(Hapsvd2$u)*sqrt(Hapsvd2$d))
dfls$haplo = Hapsvd2
row.names(dfls)<-dfls$abrevGaramszegi2019
dfls$phylo = Snpsvd1
dfls$haplo = Hapsvd2

#Analysis of the influence of shared ancestry and recent admixture on the among-breed variation in litter size
littersize1=MCMCglmm(littersize ~ 1, 
                     random=~idv(phylo)+idv(haplo), 
                     data=dfls, 
                     family="gaussian", 
                     prior=prior1,
                     thin=50,
                     burnin=70000,
                     nitt=400000,
                     verbose=T)
plot(littersize1)
summary(littersize1)
autocorr.diag(littersize1$VCV)
autocorr.diag(littersize1$Sol)
plot.acfs(littersize1$VCV)

#Shared ancestry and haplotype sharing data for all breeds for which we have data for birth weight (dfbw)
snpsdfbw<-read.csv("snpsdfbw.csv", row.names=1, header=TRUE)
hapdfbw<-read.csv("hapdfbw.csv", row.names=1, header=TRUE)

#Brith weight
dfbw<-subset(data, select=c(breed, abrevGaramszegi2019, bwkg))
dfbw<- na.omit(dfbw)

#To include genomic data as a random factor in the MCMCglmm model:
snpsdfbw<-data.frame(sapply(snpsdfbw, function(x) as.numeric(as.character(x))))
Snpsvd1 = svd(as.matrix(snpsdfbw))
Snpsvd1 = Snpsvd1$v%*%(t(Snpsvd1$u)*sqrt(Snpsvd1$d))
dfbw$phylo = Snpsvd1
hapdfbw<-data.frame(sapply(hapdfbw, function(x) as.numeric(as.character(x))))
Hapsvd2 = svd(as.matrix(hapdfbw))
Hapsvd2 = Hapsvd2$v%*%(t(Hapsvd2$u)*sqrt(Hapsvd2$d))
dfbw$haplo = Hapsvd2
row.names(dfbw)<-dfbw$abrevGaramszegi2019
dfbw$phylo = Snpsvd1
dfbw$haplo = Hapsvd2

#Analysis of the influence of shared ancestry and recent admixture on the among-breed variation in birth weight
birthweight1=MCMCglmm(sqrt(bwkg) ~ 1, 
                      random=~idv(phylo)+idv(haplo), 
                      data=dfbw, 
                      family="gaussian", 
                      prior=prior1,
                      thin=50,
                      burnin=70000,
                      nitt=400000,
                      verbose=T)
plot(birthweight1)
summary(birthweight1)
autocorr.diag(birthweight1$VCV)
autocorr.diag(birthweight1$Sol)
plot.acfs(birthweight1$VCV)

#TABLES 2 & 3 OUTPUT

#snps & hap: Shared ancestry and haplotype sharing among dog breeds based on genomic data from Parker et al (2017), following Garamszegi et al (2020)
snps<-read.csv("snps.csv", row.names=1, header=TRUE)
hap<-read.csv("hap.csv", row.names=1, header=TRUE)

#Subset lifehistory database to include only breeds for which we have data for all traits (n = 92 breeds) 
df92<-subset(data, select=c(breed,
                            abrevGaramszegi2019,
                            wAKC,
                            lifespan,
                            littersize,
                            bwkg,
                            bwkgxls,
                            growth))
df92<- na.omit(df92)
df92$abrevGaramszegi2019<-as.character(df92$abrevGaramszegi2019)
row.names(df92)<-df92$abrevGaramszegi2019

#Data visualization

##Weight
hist(sqrt(df92$wAKC))
##Lifespan
hist(df92$lifespan)
##Litter size
hist(df92$littersize)
##Birth weight
hist(sqrt(df92$bwkg))
##bwkgxls
hist(sqrt(df92$bwkgxls))
##(Adult weight - birth weight) / birth weight
hist(df92$growth)

#Data transformation
a1<-sqrt(df92$wAKC)
a2<-sqrt(df92$bwkg)
a3<-sqrt(df92$bwkgxls)
df92["weight"]<-a1
df92["birthweight"]<-a2
df92["reproductiveinvestment"]<-a3

#Dataframe with transformed data.
df<-subset(df92, select=c(breed,
                          abrevGaramszegi2019,
                          weight,
                          lifespan,
                          littersize,
                          birthweight,
                          reproductiveinvestment,
                          growth))

mcor1<-subset(df, select=c(weight,
                           lifespan,
                           littersize,
                           birthweight,
                           reproductiveinvestment,
                           growth))
res <- cor(mcor1)
round(res)
chart.Correlation(mcor1, 
                  histogram=TRUE, 
                  pch=20)

#To include genomic data as a random factor in the MCMCglmm model:
snps<-data.frame(sapply(snps, function(x) as.numeric(as.character(x))))
Snpsvd1 = svd(as.matrix(snps))
Snpsvd1 = Snpsvd1$v%*%(t(Snpsvd1$u)*sqrt(Snpsvd1$d))
df$phylo = Snpsvd1
hap<-data.frame(sapply(hap, function(x) as.numeric(as.character(x))))
Hapsvd2 = svd(as.matrix(hap))
Hapsvd2 = Hapsvd2$v%*%(t(Hapsvd2$u)*sqrt(Hapsvd2$d))
df$haplo = Hapsvd2
row.names(df)<-df$abrevGaramszegi2019
df$phylo = Snpsvd1
df$haplo = Hapsvd2

#MCMCglmm model: lifespan ~ weight * reproductive investment (Table 2)
Table2=MCMCglmm(scale(lifespan)~scale(reproductiveinvestment)*scale(weight),
                   random=~idv(phylo)+idv(haplo), 
                   data=df, 
                   family="gaussian",
                   prior=prior1, 
                   thin=50,
                   burnin=70000,
                   nitt=400000,
                   verbose=F)
plot(Table2)
summary(Table2)
autocorr.diag(Table2$VCV)
autocorr.diag(Table2$Sol)
plot.acfs(Table2$VCV)


#MCMCglmm model: lifespan ~ growth * reproductive investment (Table 3)
Table3=MCMCglmm(scale(lifespan)~scale(growth)*scale(reproductiveinvestment), 
                 random=~idv(phylo)+idv(haplo), 
                 data=df, 
                 family="gaussian",
                 prior=prior1, 
                 thin=50,
                 burnin=70000,
                 nitt=400000,
                 verbose=T)
plot(Table3)
summary(Table3)
autocorr.diag(Table3$VCV)
autocorr.diag(Table3$Sol)
plot.acfs(Table3$VCV)

### Supplementary material

#Confirmation of the relationship between lifespan and adult weight 
#We confirmed the previously reported negative association between lifespan and adult body size (Speakman et al. 2003; Fleming et al. 2011; Greer et al. 2011; Selman et al. 2013)

#MCMCglmm model: lifespan ~ weight (table S5)
S5=MCMCglmm(scale(lifespan) ~ scale(weight), 
                random=~idv(phylo)+idv(haplo), 
                data=df, 
                family="gaussian",
                prior=prior1, 
                thin=50,
                burnin=70000,
                nitt=400000,
                verbose=T)
plot(S5)
summary(S5)
autocorr.diag(S5$VCV)
autocorr.diag(S5$Sol)
plot.acfs(S5$VCV)

#MCMCglmm model: lifespan ~ reproductive investment (S6)
S6=MCMCglmm(scale(lifespan)~scale(reproductiveinvestment), 
                random=~idv(phylo)+idv(haplo), 
                data=df, 
                family="gaussian",
                prior=prior1, 
                thin=50,
                burnin=70000,
                nitt=400000,
                verbose=T)
plot(S6)
summary(S6)
autocorr.diag(S6$VCV)
autocorr.diag(S6$Sol)
plot.acfs(S6$VCV)


## #MCMCglmm model: lifespan ~ growth rate * reproductive investment (Table S7)

#Extra data attending reviewer comments on growth and growth rate:
datagr<-read.csv("datagr.csv", header=TRUE)

#Shared ancestry and haplotype sharing data for all breeds for which we have data for growth rates (dfgr)
snpsdfgr<-read.csv("snpsdfgr.csv", row.names=1, header=TRUE)
hapdfgr<-read.csv("hapdfgr.csv", row.names=1, header=TRUE)

#Growth rates subset for model (S8)
dfgr<-subset(datagr, select=c(breed, abrevGaramszegi2019, lifespan, bwkgxls,growthrate))
dfgr<- na.omit(dfgr)

#To include genomic data as a random factor in the MCMCglmm model:
snpsdfgr<-data.frame(sapply(snpsdfgr, function(x) as.numeric(as.character(x))))
Snpsvd1 = svd(as.matrix(snpsdfgr))
Snpsvd1 = Snpsvd1$v%*%(t(Snpsvd1$u)*sqrt(Snpsvd1$d))
dfgr$phylo = Snpsvd1
hapdfgr<-data.frame(sapply(hapdfgr, function(x) as.numeric(as.character(x))))
Hapsvd2 = svd(as.matrix(hapdfgr))
Hapsvd2 = Hapsvd2$v%*%(t(Hapsvd2$u)*sqrt(Hapsvd2$d))
dfgr$haplo = Hapsvd2
row.names(dfgr)<-dfgr$abrevGaramszegi2019
dfgr$phylo = Snpsvd1
dfgr$haplo = Hapsvd2

#Data visualization

##Lifespan
hist(dfgr$lifespan)
##bwkgxls
hist(sqrt(dfgr$bwkgxls))
##growth slopes
hist(sqrt(dfgr$growthrate))

#Data transformation
b1<-sqrt(dfgr$bwkgxls)
b2<-sqrt(dfgr$growthrate)

dfgr["rep_invest"]<-b1
dfgr["growthrates"]<-b2


#Dataframe with transformed data.
df34<-subset(dfgr, select=c(breed,
                          abrevGaramszegi2019,
                          lifespan,
                          rep_invest,
                          growthrates))

#MCMCglmm model: lifespan ~ growth rate * reproductive investment (Table S7)
S7=MCMCglmm(scale(lifespan)~scale(growthrates)*scale(rep_invest), 
                random=~idv(phylo)+idv(haplo), 
                data=dfgr, 
                family="gaussian",
                prior=prior1, 
                thin=50,
                burnin=70000,
                nitt=400000,
                verbose=T)
plot(S7)
summary(S7)
autocorr.diag(S7$VCV)
autocorr.diag(S7$Sol)
plot.acfs(S7$VCV)