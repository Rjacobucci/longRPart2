wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/rjacobuc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")




wisc.verb <- wisc[,c(1:4,9)]

wisc.verb$noise <- round(rnorm(nrow(wisc.verb)),2)

# create subset for plotting
ntot <- nrow(wisc.verb)    # total number of observations
wisc.verb.sel <- wisc.verb[sample(ntot, 30), ]

wisc.long <- reshape(wisc.verb, varying = c("V1", "V2", "V4", "V6"), v.names = "verbal",
                     times = c(1, 2, 4, 6), direction = "long")

wisc.long.sel <- reshape(wisc.verb.sel, varying = c("V1", "V2", "V4", "V6"),
                         v.names = "verbal", times = c(1, 2, 4, 6),
                         direction = "long")
head(wisc.long,3)
names(wisc.long)[3] <- "grade"
names(wisc.long.sel)[3] <- "grade"


fm1 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
            data = Loblolly,
            fixed = Asym + R0 + lrc ~ 1,
            random = Asym ~ 1,
            start = c(Asym = 103, R0 = -8.5, lrc = -3.3))



mix1 <- nlme(verbal ~ SSasymp(), random = ~ grade | id,start=c(1,1),
             fixed = b0i + b1i | 1,
            data = wisc.long, method="ML" )
summary(mix1) # get same estimates as in LGM, notice SD not VAR

## longRPart2


lme1 <- lme(verbal ~ grade,random=~1|id,wisc.long)



lcart.mod1 <- longRPart2(method="lme",
                         fixedFormula=verbal ~ grade,
                         rPartFormula = ~ Moeducat + noise,
                         randomFormula=~1|id,
                         data=wisc.long,
                         control=rpart.control(cp=.0025))

lrpCV(lcart.mod1)

summary(lcart.mod1)


printcp(lcart.mod1$rpart_out)
plot(lcart.mod1$rpart_out)

lrp2Plot(lcart.mod1)
lrpTreePlot(lcart.mod1,use.n=F)
