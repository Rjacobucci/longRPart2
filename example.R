#wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/rjacobuc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")

#wisc <- read.table(file.choose()) # use directory to find wisc4vpe.dat
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")

# install longRPart2 from source



wisc.verb <- wisc[,c(1:4,9)]

# create subset for plotting
ntot <- nrow(wisc.verb)    # total number of observations
wisc.verb.sel <- wisc.verb[sample(ntot, 30), ]

wisc.long <- reshape(wisc.verb, varying = c("V1", "V2", "V4", "V6"), v.names = "verbal",
                     times = c(1, 2, 4, 6), direction = "long")

wisc.long.sel <- reshape(wisc.verb.sel, varying = c("V1", "V2", "V4", "V6"),
                         v.names = "verbal", times = c(1, 2, 4, 6),
                         direction = "long")
head(wisc.long,3)
names(wisc.long)[2] <- "grade"
names(wisc.long.sel)[2] <- "grade"


## longRPart2


lcart.mod1 <- longRPart2(method="lme",
                         fixedFormula=verbal ~ grade,
                         rPartFormula = ~ Moeducat,
                         randomFormula=~1|id,
                         data=wisc.long)
lcart.mod1$rpart_out


summary(lcart.mod1)
plot(lcart.mod1);text(lcart.mod1)

# mean slope for mother's education = 0 is 4.1
# mean slope for mother's education = 1 is 5


library(semtree);library(OpenMx)

resVars <- mxPath( from=c("V1", "V2", "V4", "V6"), arrows=2,
                   free=TRUE, values = c(1,1,1,1),
                   labels=c("residual","residual","residual","residual") )
latVars<- mxPath( from=c("intercept","slope"), arrows=2, connect="unique.pairs",
                  free=TRUE, values=c(1,.4,1), labels=c("vari1","cov1","vars1") )
intLoads<- mxPath( from="intercept", to=c("V1", "V2", "V4", "V6"), arrows=1,
                   free=FALSE, values=c(1,1,1,1) )
sloLoads<- mxPath( from="slope", to=c("V1", "V2", "V4", "V6"), arrows=1,
                   free=FALSE, values=c(1,2,4,6) )
manMeans<- mxPath( from="one", to=c("V1", "V2", "V4", "V6"), arrows=1,
                   free=FALSE, values=c(0,0,0,0) )
latMeans<- mxPath( from="one", to=c("intercept","slope"), arrows=1,
                   free=TRUE,  values=c(0,1), labels=c("meani1","means1") )
dataRaw<- mxData( observed=wisc[,c(1:4,9)], type="raw" )
lgm.mod<- mxModel("LGM", type="RAM",
                  manifestVars=c("V1", "V2", "V4", "V6"),
                  latentVars=c("intercept","slope"),dataRaw,
                  resVars, latVars, intLoads, sloLoads, manMeans, latMeans)
mod.run <- mxRun(lgm.mod);coef(mod.run)


mytree <- semtree(mod.run,wisc[,c(1:4,9)]) # only moeducat as covariate
plot(mytree)

# exact same results as longRPart2
