library(nlme)
# performing an lme timing test: Need to run this at home
# using the cd4 data
dat = read.table("aidsData.csv",sep=',')names(dat) = c("Subject","Treatment","Age","Gender","Week","cd4")# Subject should be a factordat$Subject = as.factor(dat$Subject)# it looks like the week variable was converted from days, so# I'm going to convert back for nowdat$Day = round(dat$Week*7,0)dat$D = dat$Daydat$Day = as.factor(dat$D)# That's note QUITE small enough (seems like too many variables),# so I'm just going to round week to nearest intdat$Week = round(dat$Week,0)dat$W = dat$Weekdat$Week = as.factor(dat$W)### one problem in dataset: patient 911 came on days 112 and 114, which both### fall in week 16.  for now, will change second observation to week 17dat$W[3484] = 17dat$Week[3484] = 17# Treatment seems pretty well dist'd (expected from RCT)# table(dat$Treatment)dat$Treatment = as.factor(dat$Treatment)# Age ranges from 15-74, and seems pretty normal, slightly skewed # left (as long as skewed left means the tail is longer on the# right, I forget how that works)# hist(dat$Age)# rounding Age to the nearest year reduces splits from 1200 to 56dat$Age = round(dat$Age)# WAYYYYY more males than females (4475 to 561), so that variable's # probably out# table(dat$Gender)dat$Gender = as.factor(dat$Gender)

#5036 observations, going to use samples of 
#50, 100, 200, 400, 800, 1600, 3200, 5036 to 
#test the computation time for lme
x = c(100,200,400,800,1600,3200,5036)
ind = 1:(dim(dat)[1])
times = list()
for(j in 1:10){
  order = sample(x,length(x),replace=F)
  times[[j]] = matrix(nrow=8,ncol=2)
  for(i in order){	
    d = dat[sample(ind,i,replace=F),]
    t = proc.time()[3]
    m = try(lme(cd4~Week,dat=d,random=~1|Subject,correlation=corExp(form=~W)),silent=T)
    times[[j]][x==i,1] = i
    if(class(m)=='lme'){
      times[[j]][x==i,2] = proc.time()[3]-t
    }
    else{
      times[[j]][x==i,2] = NA
    }
  }
}