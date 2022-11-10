setwd('I:/Hohenheim/haiyin/Elk')
Dat<-read.table('AlcesCounts.txt')
head(Dat)

## 1. Formulate likelihood function and use MLE to estimate (r,K,S)

# ricker.fn - stochastic Ricker model
# negloglik.ricker - Likelihood function
# n.prev - Nt,population in year t
# n.next - Nt+1, population in year t+1
# r - parameter rt
# k - parameter K 
# s - parameter size in NegBinom
### all parameters are in lowercase

# 1.1 Ricker Model
ricker.fn<-function(n.prev,r,k){
  n.pred <- n.prev * exp(r*(1-n.prev/k))
  return(n.pred)
}

# 1.2 Likelihood function
negloglik.ricker<-function(par,n.prev,n.next){
  # Elements in par: 'r', 'k','s'
  r<-par['r']
  k<-par['k']
  s<-par['s']
  # n.pred - predicted population of year t+1, different from n.next
  n.pred<-ricker.fn(n.prev=n.prev,r=r,k=k)
  nLL<- -sum(dnbinom(x=n.next,mu=n.pred,size=s,log=TRUE))
  return(nLL)
}

# 1.3 Data Preparation
# rearrange the data in pairs of observed population size for each transition
n.years<-nrow(Dat)
Dat.trans<-data.frame(n.prev=Dat$Count[-n.years], n.next=Dat$Count[-1])
# start parameters
start.par<-c(r=0.5,k=500,s=10)

# 1.4 Optimization
MLE<-optim(par=start.par,fn=negloglik.ricker, n.prev=Dat.trans$n.prev, n.next=Dat.trans$n.next)
MLE

par.ricker<-MLE$par
par.ricker

AIC.ricker<- 2*MLE$value+2*length(MLE$par)
AIC.ricker
  
## 2. Formulate alternative model without density-dependence
#     Compare the models

# 2.1 exponential growth model
expo.fn<-function(n.prev,r){
  n.next <- n.prev * exp(r)
  return(n.next)
}

# 2.2 likelihood function
negloglik.expo<-function(par,n.prev,n.next){
  # Elements in par: 'r' 's'
  r<-par['r']
  s<-par['s']
  # n.pred - predicted population of year t+1, different from n.next
  n.pred<-expo.fn(n.prev=n.prev,r=r)
  nLL<- -sum(dnbinom(x=n.next,mu=n.pred,size=s,log=TRUE))
  return(nLL)
}


# 2.3 data preparation
# rearrange the data in pairs of observed population size for each transition
n.years<-nrow(Dat)
Dat.trans<-data.frame(n.prev=Dat$Count[-n.years], n.next=Dat$Count[-1])
start.par2<-c(r=0.5, s=10)

# 2.4 optimization
MLE2<-optim(par=start.par2,fn=negloglik.expo,n.prev=Dat.trans$n.prev,n.next=Dat.trans$n.next)
MLE2


par.expo<-MLE2$par
par.expo

AIC.expo<- 2*MLE2$value+2*length(MLE2$par)
AIC.expo

#3. Simulate the dynamics 1994-2013
# 3.1 Model 1: Ricker
Ricker.Sim <- function(N0, r, k, s, t.max){
  # Arguments
  # N0 - initial population size
  # r - intrinsic population growth rate
  # K - carrying capacity
  # s - size parameter of the NegBinom distribution 
  #     for demographic and environmental stochasticity
  # t.max - number of simulated time steps
  N.vec <- numeric(t.max)
  N <- N0
  for(t in 1:t.max){
    mu.N <- ricker.fn(n.prev = N, r = r, k = k)
    N <- rnbinom(n = 1, mu = mu.N, size = s)
    N.vec[t] <- N
  }
  return(N.vec)

  }

# Simulate a single realisation for the years 1994-2013
n.years.pred <- 2013 - 1993
N.ricker.sim <- Ricker.Sim(N0 = Dat$Count[Dat$Year == 1993], 
                    r = par.ricker['r'], k = par.ricker['k'], 
                    s =par.ricker['s'], t.max = n.years.pred)
plot(N.ricker.sim, type = 'l', 
     xlab = 'Time step', ylab = 'Population size N')

# replicate 1000 times
n.rep<-1000
Cols<-rainbow(n.rep)
# matrix to store predictions from all replicates
Pred.Mat<-matrix(0,nrow=n.rep,ncol=n.years.pred)
# loop
for(i in 1:n.rep){
  N.sim <- Ricker.Sim(N0 = Dat$Count[Dat$Year == 1993], 
                    r = par.ricker['r'], k = par.ricker['k'], 
                    s =par.ricker['s'], t.max = n.years.pred)
  lines(N.sim, col=Cols[i])
  Pred.Mat[i,]<-N.sim
}

# mean
Pred.mean.ricker<-colMeans(Pred.Mat)
plot(Pred.mean.ricker,type='l', ylim=c(100,400),xlab='Time Step from 1994 to 2013', 
     ylab='Population Size',lty=19, col='black')

# 25% and 75% quantiles
Pred.Q25<-apply(X=Pred.Mat, MARGIN=2, FUN=quantile, prob=0.25)
lines(Pred.Q25, lty=2,col='red')
Pred.Q75<-apply(X=Pred.Mat, MARGIN=2, FUN=quantile, prob=0.75)
lines(Pred.Q75, lty=2,col='blue')
legend('topright',lty=c(19,2,2),col=c('black','red','blue'),legend=c('mean','25% quantile','75% quantile'))

# Standard Deviation
Pred.sd<-apply(X=Pred.Mat, MARGIN=2, FUN=sd)
plot(Pred.sd, type='l')

# Variance
Pred.var<-apply(X=Pred.Mat,MARGIN=2,FUN=var)
plot(Pred.var,type='l',xlab='Time step for 1994 to 2013', 
     ylab='variance')
max(Pred.var)


# 3.2 Model 2: Exponential
Expo.Sim <- function(N0, r, s, t.max){
  # Arguments
  # N0 - initial population size
  # r - intrinsic population growth rate
  # K - carrying capacity
  # s - size parameter of the NegBinom distribution 
  #     for demographic and environmental stochasticity
  # t.max - number of simulated time steps
  N.vec <- numeric(t.max)
  N <- N0
  for(t in 1:t.max){
    mu.N <- expo.fn(n.prev = N, r = r)
    N <- rnbinom(n = 1, mu = mu.N, size = s)
    N.vec[t] <- N
  }
  return(N.vec)
}

# Simulate a single realisation for the years 1994-2013
n.years.pred <- 2013 - 1993
N.expo.sim <- Expo.Sim(N0 = Dat$Count[Dat$Year == 1993], 
                    r = par.expo['r'], 
                    s =par.expo['s'], t.max = n.years.pred)

plot(N.expo.sim, type = 'l', 
     xlab = 'Time step', ylab = 'Population size N')


# replicate 1000 times
n.rep<-1000
Cols<-rainbow(n.rep)
# matrix to store predictions from all replicates
Pred.Mat.exp<-matrix(0,nrow=n.rep,ncol=n.years.pred)
# loop
for(i in 1:n.rep){
  N.sim <- Expo.Sim(N0 = Dat$Count[Dat$Year == 1993], 
                    r = par.expo['r'], 
                    s =par.expo['s'], t.max = n.years.pred)
  lines(N.sim, col=Cols[i])
  Pred.Mat.exp[i,]<-N.sim
}

# mean
Pred.mean.exp<-colMeans(Pred.Mat.exp)
plot(Pred.mean.exp,type='l', ylim=c(100,2000),xlab='Time Step from 1994 to 2013', 
     ylab='Population Size',lty=19, col='black')
# 25% and 75% quantiles
Pred.Q25.exp<-apply(X=Pred.Mat.exp, MARGIN=2, FUN=quantile, prob=0.25)
lines(Pred.Q25.exp, lty=2,col='red')
Pred.Q75.exp<-apply(X=Pred.Mat.exp, MARGIN=2, FUN=quantile, prob=0.75)
lines(Pred.Q75.exp, lty=2,col='blue')
legend('topleft',lty=c(19,2,2),col=c('black','red','blue'),legend=c('mean','25% quantile','75% quantile'))

# Standard Deviation
Pred.sd.exp<-apply(X=Pred.Mat.exp, MARGIN=2, FUN=sd)
plot(Pred.sd.exp)

# Variance
Pred.var.exp<-apply(X=Pred.Mat.exp,MARGIN=2,FUN=var)
plot(Pred.var.exp,type='l',xlab='Time step for 1994 to 2013', ylab='variance')

# 4. Bootstrapping

start.par.bs<-par.ricker
n.bs <- 1000
Par.Mat <- matrix(0,nrow=n.bs, ncol=length(start.par.bs))
Dat.trans.bs<-data.frame(n.prev=Dat$Count[-n.years], n.next=Dat$Count[-1])
for(i.bs in 1:n.bs){
  #generating resampled data set
  bs.ind <- sample.int(n=nrow(Dat.trans.bs),repl=TRUE)
  Dat.bs <- Dat.trans.bs[bs.ind,]  
 
  # Estimate parameters from the resampled data set
  MLE.bs <- optim(par=start.par.bs,fn=negloglik.ricker, 
                  n.prev=Dat.bs$n.prev, n.next=Dat.bs$n.next)
  Par.Mat[i.bs,] <- MLE.bs$par
}

# Plot the distribution of parameters
colnames(Par.Mat)<-c('r','k','s')

# r
hist(Par.Mat[,'r'],col='green', xlab='Estimation Value of r',main="")
# Mean of parameter estimates for r
r.mean<-mean(Par.Mat[,'r'])
abline(v=r.mean, lwd=2)
# Standard Error= standard deviation of the parameter sample
r.se<-sd(Par.Mat[,'r'])
#95% Confidence Interval
r.ci<-quantile(Par.Mat[,'r'],prob=c(0.025,0.975))
abline(v=r.ci,lty=2)

# k
hist(Par.Mat[,'k'],col='green',main="", xlab='Estimation Value of k',xlim=c(0,700),breaks=1000000)
# Mean of parameter estimates for k
k.mean<-mean(Par.Mat[,'k'])
abline(v=k.mean, lwd=2)
# Standard Error= standard deviation of the parameter sample
k.se<-sd(Par.Mat[,'k'])
#95% Confidence Interval
k.ci<-quantile(Par.Mat[,'k'],prob=c(0.025,0.975))
abline(v=k.ci,lty=2)

# s
hist(Par.Mat[,'s'],col='green', xlab='Estimation Value of s',main="")
# Mean of parameter estimates for s
s.mean<-mean(Par.Mat[,'s'])
abline(v=s.mean, lwd=2)
# Standard Error= standard deviation of the parameter sample
s.se<-sd(Par.Mat[,'s'])
#95% Confidence Interval
s.ci<-quantile(Par.Mat[,'s'],prob=c(0.025,0.975))
abline(v=s.ci,lty=2)

# make the graphs into one
# 2 up 1 down
layout(mat = matrix(c(1,1,2,2,
                      0,3,3,0), nrow = 2, byrow = TRUE))
layout.show(n = 3)

# r
hist(Par.Mat[,'r'],col='slateblue4', xlab='Estimation Value of r',main="")
# Mean of parameter estimates for r
r.mean<-mean(Par.Mat[,'r'])
abline(v=r.mean, lwd=2)
#95% Confidence Interval
r.ci<-quantile(Par.Mat[,'r'],prob=c(0.025,0.975))
abline(v=r.ci,lty=2)
legend('topright',lty=c(19,2),lwd=c(2,1),bty='n',legend=c('mean','95%CI'))

# k
hist(Par.Mat[,'k'],col='royalblue1',main="", xlab='Estimation Value of k',xlim=c(0,700),breaks=1000000)
# Mean of parameter estimates for k
k.mean<-mean(Par.Mat[,'k'])
abline(v=k.mean, lwd=2)
#95% Confidence Interval
k.ci<-quantile(Par.Mat[,'k'],prob=c(0.025,0.975))
abline(v=k.ci,lty=2)
legend('topright',lty=2,bty='n',legend='95%CI')


# s
hist(Par.Mat[,'s'],col='skyblue1', xlab='Estimation Value of s',main="")
# Mean of parameter estimates for s
s.mean<-mean(Par.Mat[,'s'])
abline(v=s.mean, lwd=2)
#95% Confidence Interval
s.ci<-quantile(Par.Mat[,'s'],prob=c(0.025,0.975))
abline(v=s.ci,lty=2)
legend('topright',lty=c(19,2),lwd=c(2,1),bty='n',legend=c('mean','95%CI'))

#5. Simulate with bootstrapped parameters
# 5.1 Model 1: Ricker
Ricker.Sim2 <- function(N0, r, k, t.max){
  # Arguments
  # N0 - initial population size
  # r - intrinsic population growth rate
  # K - carrying capacity
  # t.max - number of simulated time steps
  N.vec <- numeric(t.max)
  N <- N0
  for(t in 1:t.max){
    mu.N <- ricker.fn(n.prev = N, r = r, k = k)
    N<-mu.N
    N.vec[t] <- N
  }
  return(N.vec)
  
}

# Simulate a single realisation for the years 1994-2013
n.years.pred <- 2013 - 1993
N.mat<-matrix(nrow=n.bs,ncol=n.years.pred)

for(i in 1:n.bs){
  r=Par.Mat[i,'r']
  k=Par.Mat[i,'k']

  N.mat[i,] <- Ricker.Sim2(N0 = Dat$Count[Dat$Year == 1993], 
                             r = r, k = k,  t.max = n.years.pred)
}

# mean
Pred.mean.ricker2<-colMeans(N.mat)
plot(Pred.mean.ricker2,type='l',ylim=c(200,300),xlab='Time Step from 1994 to 2013', 
     ylab='Population Size',lty=19,col='black')
# median
Pred.median.ricker2<-apply(X=N.mat,MARGIN=2,FUN=median)
lines(Pred.median.ricker2, lty=2,col='black')

# 25% and 75% quantiles
Pred.Q25<-apply(X=N.mat, MARGIN=2, FUN=quantile, prob=0.25)
lines(Pred.Q25, lty=2,col='red')
Pred.Q75<-apply(X=N.mat, MARGIN=2, FUN=quantile, prob=0.75)
lines(Pred.Q75, lty=2,col='blue')
legend('topleft',lty=c(19,2,2,2),col=c('black','black','red','blue'),legend=c('mean','median','25% quantile','75% quantile'))

# Standard Deviation
Pred.sd<-apply(X=N.mat, MARGIN=2, FUN=sd)
plot(Pred.sd)
# Variance
Pred.var<-apply(X=N.mat,MARGIN=2,FUN=var)
plot(Pred.var,type='l',xlab='Time step for 1994 to 2013', ylab='variance')
max(Pred.var)



### testing stochastic with BS
Ricker.Sim3 <- function(N0, r, k, s, t.max){
  # Arguments
  # N0 - initial population size
  # r - intrinsic population growth rate
  # K - carrying capacity
  # s - size parameter of the NegBinom distribution 
  #     for demographic and environmental stochasticity
  # t.max - number of simulated time steps
  N.vec <- numeric(t.max)
  N <- N0
  for(t in 1:t.max){
    mu.N <- ricker.fn(n.prev = N, r = r, k = k)
    N<-rnbinom(n = 1, mu = mu.N, size = s)
    N.vec[t] <- N
  }
  return(N.vec)
  
}

n.years.pred <- 2013 - 1993
N.mat<-matrix(nrow=n.bs,ncol=n.years.pred)

for(i in 1:n.bs){
  r=Par.Mat[i,'r']
  k=Par.Mat[i,'k']
  s=Par.Mat[i,'s']
  N.mat[i,] <- Ricker.Sim3(N0 = Dat$Count[Dat$Year == 1993], 
                           r = r, k = k, s=s, t.max = n.years.pred)
}

# mean
Pred.mean.ricker3<-colMeans(N.mat)
plot(Pred.mean.ricker2,type='l',ylim=c(100,300),xlab='Time Step', ylab='Population Size',lty=19)
Pred.median.ricker3<-apply(X=N.mat,MARGIN=2,FUN=median)
lines(Pred.median.ricker3, lty=2,col='black')
# 25% and 75% quantiles
Pred.Q25<-apply(X=N.mat, MARGIN=2, FUN=quantile, prob=0.25)
lines(Pred.Q25, lty=2,col='blue')
Pred.Q75<-apply(X=N.mat, MARGIN=2, FUN=quantile, prob=0.75)
lines(Pred.Q75, lty=2,col='red')
legend('topleft',lty=c(19,2,2,2),col=c('black','black','blue','red'),bty='n',legend=c('mean','median','25%','75%'))
