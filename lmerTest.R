# Equivalent functions between lme and lmer
m1 <- lmer(pbkph~time+(1|Subject),data=pbkphData)
m1a <- lme(pbkph~time,random=~1|Subject,data=pbkphData,na.action=na.omit)

m2 <- lmer(pbkph~as.factor(time)+(1|Subject),data=pbkphData)
m2a <-lme(pbkph~as.factor(time),random=~1|Subject,data=pbkphData,na.action=na.omit)

m3 <- lmer(pbkph~as.factor(time)+(1|Subject),data=pbkphData)
m3a <-lme(pbkph~as.factor(time),random=~1|Subject,data=pbkphData,na.action=na.omit,correlation=corExp(form=~time))

Rprof()
for(i in 1:1000){
  m1 <- lmer(pbkph~time+(1|Subject),data=pbkphData)
  m1a <- lmer2(pbkph~time+(1|Subject),data=pbkphData)
  m1b <- lme(pbkph~time,random=~1|Subject,data=pbkphData,na.action=na.omit)

  m2 <- lmer(pbkph~as.factor(time)+(1|Subject),data=pbkphData)
  m2a <- lmer2(pbkph~as.factor(time)+(1|Subject),data=pbkphData)
  m2b <-lme(pbkph~as.factor(time),random=~1|Subject,data=pbkphData,na.action=na.omit)
}
Rprof(NULL)
head(summaryRprof()$by.total,20)
