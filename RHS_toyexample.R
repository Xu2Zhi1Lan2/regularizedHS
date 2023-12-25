library(MASS)
library(glmnet)
library(horseshoe)
source("RHS.R")

## generate data from conditional autoregression 
set.seed(2023)
Time = 50
J = 8*8

JL = 8
A = matrix(0,J,J) 
for (i in 1:J) {
  if(i<=(J-JL)){
    if(i%%JL!=0){A[i,i+1]=1}
    A[i,i+JL]=1
  }
  else{if(i<J){A[i,i+1]=1}}
}
A = A + t(A) # adjacency matrix
D = diag(as.vector(A%*%rep(1,J))) # diagonal matrix of row sum of A

x = mvrnorm(Time,rep(0,J),solve(D-0.99*A))
beta.true = rep(0,J)
S = c(1,2,1+JL,2+JL)
beta.true[S] = 0.5
matrix(beta.true,JL,JL)
y = x%*%beta.true + rnorm(Time,0,1)

## ridge regression
alpha = 0 # use L2 penalty only
set.seed(2)
enet_cv <- cv.glmnet(x, y, family = "gaussian", nfolds = 5, intercept = F,
                     alpha = alpha, type.measure = "mse")

beta_ridge = coef(enet_cv,s=enet_cv$lambda.min)[2:(J+1)]
matrix(beta_ridge,JL,JL)

## elastic net
alpha_seq <- seq(0.9, 0.1, -0.1)
cvm_min <- lamb <- rep(0, length(alpha_seq))
for (i in 1: length(alpha_seq)){
  set.seed(2)
  enet_cv <- cv.glmnet(x, y, family = "gaussian", nfolds = 5, intercept = F,
                       alpha = alpha_seq[i], type.measure = "mse")
  lamb[i] <- enet_cv$lambda.1se
  cvm_min[i] <- enet_cv$cvm[which.min(enet_cv$lambda)]
}
par_index <- which.min(cvm_min)
alpha_enet <- alpha_seq[par_index]
alpha_enet # optimal alpha = 0.1

set.seed(2)
cvfit <- cv.glmnet(x, y, family="gaussian", nfolds=5, intercept = F,
                   type.measure="mse", alpha = alpha_enet)
beta_elasticnet = coef(cvfit,s=cvfit$lambda.min)[2:(J+1)]
matrix(beta_elasticnet,JL,JL)

## regularized HS
tuning = 10000
Tb = 400
burn = 20000
nmc = 10000
thin = 5
set.seed(2023)
out.RHS <- RHS(y,x,tuning,Tb,burn,nmc,thin)
matrix(out.RHS$BetaHat,JL,JL)

## HS
out.HS <- horseshoe(y,x,method.tau = "halfCauchy",method.sigma = "Jeffreys",burn,nmc,thin)
matrix(out.HS$BetaHat,JL,JL)
