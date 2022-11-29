library(tsDyn) # to call VAR.sim
library(lavaan) # to call clpmModel
library(restriktor) # to call GORICA
library(bain) # Bayes factors

library(devtools) # Make sure you have Rtools (and a version which is compatable with your R version).
install_github("rebeccakuiper/ICweights")
library(ICweights)


# save estimater & GORICA values & GORICA weights
set.seed(123)
nsim<-1000
estsim1<-array(NA, dim=c(nsim, 2)) #beta, gamma
estsimw1<-array(NA, dim=c(nsim, 4)) #alpha1,beta1,delta1,gamma1
goricasim<-array(NA, dim=c(nsim, 2)) # H1 and H2
colnames(goricasim) <- c("H1", "H2")
weightssim<-array(NA, dim=c(nsim, 2))
colnames(weightssim) <- c("H1", "H2")
AICsim<-array(NA, dim=c(nsim, 2))
colnames(AICsim) <- c("H0", "Hu")
AICweightssim<-array(NA, dim=c(nsim, 2))
colnames(AICweightssim) <- c("H0", "Hu")
Penalty1 <- array(NA, dim=c(nsim, 2))
colnames(Penalty1) <- c("H1", "H2")
PT_weights1 <- array(NA, dim=c(nsim,2))
colnames(PT_weights1) <- c("H1", "H2")
Penalty2 <- array(NA, dim=c(nsim, 2))
colnames(Penalty2) <- c("H0", "Hu")
PT_weights2 <- array(NA, dim=c(nsim,2))
colnames(PT_weights2) <- c("H0", "Hu")
bainH1H2 <- array(NA, dim=c(nsim,2))
colnames(bainH1H2) <- c("H1", "H2")
bainH0H1H2 <- array(NA, dim=c(nsim,3))
colnames(bainH0H1H2) <- c("H0", "H1", "H2")


# specify value for generate dataset
p <- 300  # number of persons and thus number of rows in data
q <- 2   # number of variables
M <- 3   # number of measurement occasions
n <- q*M # number of columns in data
N<-100 # Number of iterations from VAR(1) model per person - and use only last M of those N iterations
B2<-matrix(c(0.22, 0.1, 0.2, 0.28),byrow = T, 2) # alpha, beta, gamma, delta
varxy <- array(NA, c(p, n)) # create array with nrow=p and ncol=n in order to keep simulate data
varxy2 <- array(NA, c(p, n))
colnames(varxy)<-c("x1", "x2", "x3", "y1", "y2", "y3")

# hypotheses
H0 <- "beta == gamma"
H1 <- "beta < gamma"
H2 <- "beta > gamma"

#H0 <- "abs(beta) == abs(gamma)"
#H1 <- "abs(beta) < abs(gamma)"
#H2 <- "abs(beta) > abs(gamma)"

simteller<- 1

for (simteller in 1:nsim) {
  
  print(paste0("Iteration", simteller))
  
  for(i in 1:p){
    var1 <- VAR.sim(B=B2, n=N, include="none")
    varxy[i,]<-var1[(N-M+1):N,]
    varxy2 <- scale(varxy)
    
  }
  
  clpmModel<- # <- ' for specify clpmModel and eding with '
    '
  kappa =~ 1*x1 + 1*x2 + 1*x3
  omega =~ 1*y1 + 1*y2 + 1*y3
  
  x1 ~ mu1*1 #intercepts
  x2 ~ mu2*1
  x3 ~ mu3*1
  y1 ~ pi1*1
  y2 ~ pi2*1
  y3 ~ pi3*1
  
  # RICLPM become CLPM because of this part
  kappa ~~ 0*kappa #variance
  omega ~~ 0*omega #variance
  kappa ~~ 0*omega #covariance
  
  #latent vars for AR and cross-lagged effects
  p1 =~ 1*x1 #each factor loading set to 1
  p2 =~ 1*x2
  p3 =~ 1*x3
  q1 =~ 1*y1
  q2 =~ 1*y2
  q3 =~ 1*y3
  
  #constrain autoregression and cross-lagged effects to be the same across both lags.
  p3 ~ alpha*p2 + beta*q2
  p2 ~ alpha*p1 + beta*q1
  
  q3 ~ delta*q2 + gamma*p2
  q2 ~ delta*q1 + gamma*p1
  
  p1 ~~ p1 #variance
  p2 ~~ u*p2
  p3 ~~ u*p3
  q1 ~~ q1 #variance
  q2 ~~ v*q2
  q3 ~~ v*q3
  
  p1 ~~ q1 #p1 and q1 covariance
  p2 ~~ uv*q2 #p2 and q2 covariance should also be constrained to be the same
  p3 ~~ uv*q3 #p2 and q2 covariance'
  
  clpmConstrainedsim <- lavaan(clpmModel, data = varxy2,
                  missing = 'ML', 
                  int.ov.free = F,
                  int.lv.free = F,
                  auto.fix.first = F,
                  auto.fix.single = F,
                  auto.cov.lv.x = F,
                  auto.cov.y = F,
                  auto.var = F)
  summary(clpmConstrainedsim, standardized = T)
  stdClpmSim<-standardizedsolution(clpmConstrainedsim, type = "std.all", se = TRUE, zstat = TRUE, 
                                          pvalue = TRUE, ci = TRUE, level = 0.95, cov.std = TRUE, 
                                          remove.eq = TRUE, remove.ineq = TRUE, remove.def = FALSE, 
                                          partable = NULL, GLIST = NULL, est = NULL)
  
  # subtract for alpha, beta, delta, gamma
  # GORICA values and weights
  est<-coef(clpmConstrainedsim)[c(10,14)]
  names(est) <- c("beta", "gamma")
  estw1<-coef(clpmConstrainedsim)[c(9,10,13,14)]
  vcov<-lavInspect(clpmConstrainedsim, "vcov")[c(10,14), c(10,14)]
  
  
  goricaResult1 <- goric(est, VCOV = vcov, H1, H2, comparison = "none", type = "gorica")
  goricaResult2 <- goric(est, VCOV = vcov, H0, comparison = "unconstrained", type = "gorica")
  
  bainResults1 <- bain(clpmConstrainedsim, "beta < gamma; beta > gamma")
  bainResults2 <- bain(clpmConstrainedsim, "beta = gamma; beta < gamma; beta > gamma")
  
  #parameter estimates
  estsim1[simteller,]<-est
  meanEst<-apply(estsim1, 2, mean, na.rm=TRUE)
  names(meanEst) <- c("beta", "gamma")
  estsimw1[simteller,]<-estw1
  meanEstsimw1<-apply(estsimw1, 2, mean, na.rm=TRUE) # na.rm = Excluding Missing Values(NA) from Analyses
  names(meanEstsimw1) <- c("alpha1", "beta1", "delta1", "gamma1")
  
  #Hypothesis evaluation
  
  goricasim[simteller,]<-goricaResult1$result[,4] 
  weightssim[simteller,]<-goricaResult1$result[,5]
  AICsim[simteller,]<-goricaResult2$result[,4] 
  AICweightssim[simteller,]<-goricaResult2$result[,5]
  Penalty1[simteller,] <- 2*goricaResult1$result[,3]
  PT_weights1[simteller,]<-IC.weights(Penalty1[simteller,])$IC.weights
  Penalty2[simteller,] <- 2*goricaResult2$result[,3]
  PT_weights2[simteller,]<-IC.weights(Penalty2[simteller,])$IC.weights
  bainH1H2[simteller,] <- bainResults1$fit[1:2,8]
  bainH0H1H2[simteller,] <- bainResults2$fit[1:3,8]
  
  simteller<- simteller+1

}

#parameter estimates
estsim
meanEst
estsimw1
meanEstsimw1


#causal dominance hypo
goricasim #H1:"beta < gamma" vs H2:"beta > gamma" 
weightssim
sum(weightssim[,1]>weightssim[,2], na.rm=TRUE) #support H1 
sum(weightssim[,2]>weightssim[,1], na.rm=TRUE)
AICsim
AICweightssim #AIC weights based on GORICA
sum(AICweightssim[,1]>AICweightssim[,2], na.rm=TRUE)
sum(AICweightssim[,2]>AICweightssim[,1], na.rm=TRUE)

#Penalty
Penalty1
PT_weights1
Penalty2
PT_weights2

#Bayes factors (PMPa only)
bainH1H2
bainH0H1H2

write.csv(goricasim, "G:/My Drive/PhD/Sukpa001/My Documents/My project/Project 1/Final_project1/Revision/Simulations results/CLPM/Bivariate/nonabsolute values/2v3w/(0.1, 0.2)/p300.goriga.H1H2.csv", row.names = FALSE)
write.csv(weightssim, "G:/My Drive/PhD/Sukpa001/My Documents/My project/Project 1/Final_project1/Revision/Simulations results/CLPM/Bivariate/nonabsolute values/2v3w/(0.1, 0.2)/p300.goriga.weights.H1H2.csv", row.names = FALSE)

write.csv(AICsim, "G:/My Drive/PhD/Sukpa001/My Documents/My project/Project 1/Final_project1/Revision/Simulations results/CLPM/Bivariate/nonabsolute values/2v3w/(0.1, 0.2)/p300.AIC.H0Hu.csv", row.names = FALSE)
write.csv(AICweightssim, "G:/My Drive/PhD/Sukpa001/My Documents/My project/Project 1/Final_project1/Revision/Simulations results/CLPM/Bivariate/nonabsolute values/2v3w/(0.1, 0.2)/p300.AIC.weights.H0Hu.csv", row.names = FALSE)

write.csv(Penalty1, "G:/My Drive/PhD/Sukpa001/My Documents/My project/Project 1/Final_project1/Revision/Simulations results/CLPM/Bivariate/nonabsolute values/2v3w/(0.1, 0.2)/p300.Penalty.goriga.H1H2.csv", row.names = FALSE)
write.csv(PT_weights1, "G:/My Drive/PhD/Sukpa001/My Documents/My project/Project 1/Final_project1/Revision/Simulations results/CLPM/Bivariate/nonabsolute values/2v3w/(0.1, 0.2)/p300.Penalty.weights.goriga.H1H2.csv", row.names = FALSE)

write.csv(Penalty2, "G:/My Drive/PhD/Sukpa001/My Documents/My project/Project 1/Final_project1/Revision/Simulations results/CLPM/Bivariate/nonabsolute values/2v3w/(0.1, 0.2)/p300.Penalty.AIC.H0Hu.csv", row.names = FALSE)
write.csv(PT_weights2, "G:/My Drive/PhD/Sukpa001/My Documents/My project/Project 1/Final_project1/Revision/Simulations results/CLPM/Bivariate/nonabsolute values/2v3w/(0.1, 0.2)/p300.Penalty.weights.AIC.H0Hu.csv", row.names = FALSE)

write.csv(bainH1H2, "G:/My Drive/PhD/Sukpa001/My Documents/My project/Project 1/Final_project1/Revision/Simulations results/CLPM/Bivariate/nonabsolute values/2v3w/(0.1, 0.2)/p300.PMPa.H1H2.csv", row.names = FALSE)
write.csv(bainH0H1H2, "G:/My Drive/PhD/Sukpa001/My Documents/My project/Project 1/Final_project1/Revision/Simulations results/CLPM/Bivariate/nonabsolute values/2v3w/(0.1, 0.2)/p300.PMPa.H0H1H2.csv", row.names = FALSE)



