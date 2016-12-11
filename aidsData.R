library(longRPart)
data(aidsData)
aidsData$Treatment = as.factor(aidsData$Treatment)
aidsData$Gender = as.factor(aidsData$Gender)
aidsData$W = as.factor(aidsData$W)
aidsData$Subject = as.factor(aidsData$Subject)


tree = longRPart(cd4~0+W,~Treatment,~1|Subject,R=corAR1(form=~Day),data=aidsData,control = rpart.control(maxdepth=3,cp=1e-7))
tree2 = longRPart(cd4~0+W,~Treatment+Gender,~1|Subject,R=corAR1(form=~Day),data=aidsData,control = rpart.control(maxdepth=3,cp=1e-7))

pdf("aidsPlots.pdf")
lrpPlot(tree)
lrpTreePlot(tree)
dev.off()

CI2Test = lrpCI(tree2,3)
pVal2Test = lrpPVal(tree2,3)

#CI2 = lrpCI(tree2,1000)
#pVal2 = lrpPVal(tree2,1000)
