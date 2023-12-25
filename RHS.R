RHS <- function (y,X,tuning,Tb,burn,nmc,thin,alpha_sig = 0.05,a0=0.1,b0=0.1) 
{
  N = burn + nmc + tuning # the length of sampling
  effsamp = nmc/thin # ultimate chain length we get
  n = nrow(X) # sample size
  p = ncol(X) # dimension
  
  ### initial values of the chain
  Beta = rep(0.25, p)
  lambda = rep(1, p) 
  sigmasq = 1
  csq = 1
  tau = 1
  lambda_tilde = sqrt(csq*lambda^2/(csq + tau^2*lambda^2))
  
  ### storing the chain
  betaout = matrix(0, p, effsamp) 
  lambdaout = matrix(0, p, effsamp) 
  tauout = rep(0, effsamp)
  csqout = rep(0, effsamp)
  sigmasqout = rep(0, effsamp)
  
  ### defining variables for tuning
  step.lambda = rep(0.25,p)
  step.tau = 0.25
  
  acc.lambda = rep(0,p)
  acc.tau = 0
  
  ### some one vectors and pre-calculation
  I_n = diag(n)
  l0 = rep(0, p)
  l1 = rep(1, n)
  l2 = rep(1, p)
  
  ### sampling
  for (i in 1:N) {
    
    ## sampling beta
    lambda_star = tau * lambda_tilde
    U = as.numeric(lambda_star^2) * t(X)
    u = stats::rnorm(l2, l0, sqrt(lambda_star))
    v = X %*% u + stats::rnorm(n)
    v_star = solve((X %*% U + I_n), ((y/sqrt(sigmasq)) - v))
    Beta = sqrt(sigmasq) * (u + U %*% v_star)
    
    ## sampling sigmasq
    E_1 = max(t(y - X %*% Beta) %*% (y - X %*% Beta), (1e-10))
    E_2 = max(sum(Beta^2/lambda_star^2), (1e-10))
    sigmasq = 1/stats::rgamma(1, a0+(n + p)/2, scale = 2/(E_1 + E_2 + b0))
    
    ## sampling csq
    csq = 1/stats::rgamma(1, 5+p/2, scale = 2/(sum(Beta^2) + 1))
    
    ## sampling lambdasq
    lambda.new = rep(0,p)
    for (il in 1:p) {
      lambda.star = exp(log(lambda[il]) + rnorm(1,0,step.lambda[il]))
      lambda_tilde_star = sqrt(csq*lambda.star^2/(csq + lambda.star^2*tau^2))
      T0 = log(lambda.star)-log(lambda[il])  
      T1 = log(1 + lambda[il]^2) - log(1 + lambda.star^2)
      T2 = log(abs(lambda_tilde[il])) - log(abs(lambda_tilde_star)) + (0.5*(1/lambda_tilde[il]^2 - 1/lambda_tilde_star^2)*(Beta[il])^2/tau^2) # likelihood
      MH = min(1,exp(T0+T1+T2))
      if(is.nan(MH)){print(c(T0,T1,T2,lambda.star,lambda[il], Beta[il], tau))}
      Ul = runif(1,0,1)
      if(Ul<=MH){lambda.new[il] = lambda.star;acc.lambda[il]=acc.lambda[il]+1;lambda_tilde[il]=lambda_tilde_star}
      else{lambda.new[il] = lambda[il];acc.lambda[il]=acc.lambda[il]}
    }
    lambda = lambda.new
    
    ## sampling tau
    tau.new = 0
    tau.star = exp(log(tau) + rnorm(1,0,step.tau))
    T0 = log(tau.star) - log(tau)
    T1 = log(1 + tau^2) - log(1 + tau.star^2)
    T2 = n*(log(abs(tau))-log(abs(tau.star))) + 0.5*(1/tau^2 - 1/tau.star^2)*sum(Beta^2/(lambda^2^2))
    MH = min(1,exp(T0+T1+T2))
    Ut = runif(1,0,1)
    if(Ut<=MH){tau.new = tau.star;acc.tau = acc.tau + 1}
    else{tau.new = tau;acc.tau = acc.tau}
    tau = tau.new
    
    ## print the progress
    if (i%%1000 == 0) {print(i)}
    
    ## tuning
    if (i<=tuning & i%%Tb == 0){
      ar.lambda = acc.lambda/Tb
      ar.tau = acc.tau/Tb
      for (il in 1:p) {
        if(ar.lambda[il]>0.5){step.lambda[il]=step.lambda[il]*1.1}
        if(ar.lambda[il]<0.3){step.lambda[il]=step.lambda[il]*0.91}
      }
      if(ar.tau>0.5){step.tau = step.tau*1.1}
      if(ar.tau<0.3){step.tau = step.tau*0.91}
      acc.lambda = rep(0,p)
      acc.tau = 0
    }
    
    ## store the thinned chain
    if (i > (burn+tuning) && (i-burn-tuning)%%thin == 0) {
      betaout[, (i - burn - tuning)/thin] = Beta
      lambdaout[, (i - burn - tuning)/thin] = lambda
      tauout[(i - burn - tuning)/thin] = tau
      csqout[(i - burn - tuning)/thin] = csq
      sigmasqout[(i - burn - tuning)/thin] = sigmasq
    }
  }
  
  ## inference
  pMean = apply(betaout, 1, mean)
  pMedian = apply(betaout, 1, stats::median)
  pSigma = mean(sigmasqout)
  pTau = mean(tauout)
  pcsq = mean(csqout)
  
  left <- floor(alpha_sig * effsamp/2)
  right <- ceiling((1 - alpha_sig/2) * effsamp)
  BetaSort <- apply(betaout, 1, sort, decreasing = F)
  left.points <- BetaSort[left, ]
  right.points <- BetaSort[right, ]
  
  ar.lambda = acc.lambda/(nmc + burn)
  ar.tau = acc.tau/(nmc + burn)
  
  ## return results
  result = list(BetaHat = pMean, LeftCI = left.points, RightCI = right.points, 
                BetaMedian = pMedian, Sigma2Hat = pSigma, TauHat = pTau, C2Hat = pcsq,
                BetaSamples = betaout, TauSamples = tauout, Sigma2Samples = sigmasqout, 
                LambdaSamples = lambdaout, C2Samples = csqout,
                ARLambda2 = ar.lambda, ARTau = ar.tau)
  return(result)
}