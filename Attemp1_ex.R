# Data Generation

z = runif(200, 0, 10)

b0i = matrix(NA,200,1)
b1i = matrix(NA,200,1)
y = matrix(NA,200,6)
id = matrix(NA,200,1)
beta_0 = matrix(NA,200,1)
beta_1 = matrix(NA,200,1)
psi2_0 = matrix(NA,200,1)
psi2_1 = matrix(NA,200,1)
sigma2_u = matrix(NA,200,1)


for(i in 1:200){
  if (z[i] > 5) {beta_0[i]   = 10}  else {beta_0[i]   = 20}
  if (z[i] > 5) {beta_1[i]   = 10}  else {beta_1[i]   = 5}
  if (z[i] > 5) {psi2_0[i]   = 100} else {psi2_0[i]   = 64}
  if (z[i] > 5) {psi2_1[i]   = 25}  else {psi2_1[i]   = 16}
  if (z[i] > 5) {sigma2_u[i] = 4}   else {sigma2_u[i] = 4}

  b0i[i] = beta_0[i] + rnorm(1,0,sqrt(psi2_0[i]))
  b1i[i] = beta_1[i] + rnorm(1,0,sqrt(psi2_1[i]))

  id[i] = i

	for(t in 1:6){
        y[i,t] = b0i[i] + b1i[i] * (t - 1)/5 + rnorm(1,0,sqrt(sigma2_u))
      }
  }


ex.data.1 = as.data.frame(cbind(id,y,z))
names(ex.data.1)= c('id','y1','y2','y3','y4','y5','y6','z')
head(ex.data.1)

# Convert to Long Form
ex.data.2 = reshape(ex.data.1, idvar='id',
        varying=c('y1','y2','y3','y4','y5','y6'),
        times=c(0,1,2,3,4,5),
        v.names='y', direction='long')

ex.data.2 = ex.data.2[order(ex.data.2$id, ex.data.2$time),]

head(ex.data.2)

### Beginning of what package would be

## Things to Specify

#1 longitudinal model (e.g., y~b0i+b1i*(time/5))
#2 fixed-effect parameters (e.g., fixed=b0i+b1i~1)
#3 random-effect parameters (e.g., random=b0i+b1i~1)
#4 data
#5 Participant ID variable (e.g., group=~id)
#6 Starting values for fixed-effects parameters (e.g., start=c(10,5))
#7 Covariates used for splitting
#8 Model Type ('Fixed', 'Fixed+Ran', 'Fixed+Ran+Res')
#9 Minimum number of participants per node (Nmin)
#10 Criterion (e.g., Change in -2LL, p-value, CV)

# Run NLME
require(nlme)

model.0 = nlme(y~b0i+b1i*time,
               data=ex.data.2,
               fixed=b0i+b1i~1,
               random=b0i+b1i~1,
               group=~id,
               start=c(10,5))

summary(model.0)

dev.root = -2 * model.0$logLik



lcart.mod1 <- longRPart2(method="nlme",
                         nlme.model=y~b0i+b1i*time,
                         fixedFormula=b0i+b1i~1,
                         rPartFormula = ~ z,
                         group= ~ id,
                         randomFormula=b0i+b1i~1,
                         data=ex.data.2,
                         start=c(10,5))

