library(longRPart)
data(pbkphData)
pbkphData$Time=as.factor(pbkphData$Time)
model = longRPart(pbkph~0+Time,~age+gender,~1|Subject,pbkphData,R=corExp(form=~time),control=rpart.control(minbucket=80,xval=10))
modelA = longRPart(pbkph~Time,~gender+age,~1|Subject,pbkphData,R=corExp(form=~time),control=rpart.control(minbucket=80,xval=10))

pdf("pbkphPlots.pdf")
lrpPlot(model,color=rainbow(4))
lrpTreePlot(model,color=rainbow(4))
dev.off()

CI = lrpCI(model,1000)
pVal = lrpPVal(model,1000)

#CV = lrpCV(model)
save.image(file='pbkphAnalysis.RData')
### CROSS-VALIDATION
#trying to get non-convergence error in a CV
#devs = list()
#for(i in 1:100){
#  set.seed(i)
#  devs[[i]] = lrpCV(modelA)
#}
