id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))


library(foreach)
library(doParallel)

#library(rSPDE)
library(cubature)

library(proxy)
library(Matrix)
library(pracma)
Adam = function(f, grad_f, params_init, lr = 0.01, beta1 = 0.9, beta2 = 0.999,
                epsilon = 1e-8, max_iter = 1000, tol = 1e-10,alpha=FALSE,alpha2=FALSE,
                print_grad=FALSE) {
  
  params = log(params_init)
  m = numeric(length(params))
  v = numeric(length(params))
  
  
  iter = 0
  
  
  while (iter < max_iter) {
    iter = iter + 1
    
    grad = grad_f(exp(params))*exp(params)
    m = beta1 * m + (1 - beta1) * grad
    
    v = beta2 * v + (1 - beta2) * (grad^2)
    
    m_hat = m / (1 - beta1^iter)
    v_hat = v / (1 - beta2^iter)
    params_new = params - lr * m_hat / (sqrt(v_hat) + epsilon)
    
    params_new = pmax(params_new, log(tol))
    
    
    if (max(abs(params_new - params)) < tol) {
      break
    }
    
    if(alpha){
      if((exp(params_new)[6]<0)|(exp(params_new)[6]>1)){params_new[6]=params[6]}}
    if(alpha2){
      if((exp(params_new)[10]<0)|(exp(params_new)[10]>1)){params_new[10]=params[10]}}
    
    params = params_new
    if(print_grad){print(grad)}
    #print(c(f(exp(params))))
  }
  #print("iter= ")
  #print(iter)
  return(exp(params))
}

eq_cov_mat = function(x1, x2, l1, sigma1, l2, sigma2, maxEval = 1000,
                      nu1=0,nu2=nu1,fundamental=FALSE) {
  if (fundamental){
    if(nu1==0){
      cov = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                   nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
      
      for(i in 1:(length(as.matrix(x1))/2)){
        #if(i%%100==0){print(i)}
        theta1=ifelse(length(as.matrix(x1))==2,atan2(x1[2],x1[1]),
                      atan2(x1[i,2],x1[i,1]))
        for (j in 1:(length(as.matrix(x2))/2)){
          r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                    sum((x1[i,])^2)^.5)
          r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                    sum((x2[j,])^2)^.5)
          #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
          theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                        atan2(x2[j,2],x2[j,1]))
          
          results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
            diag(c(sigma1^2*exp(-.5*(r1-r2)^2*(l1^2)),sigma2^2*exp(-.5*(r1-r2)^2*(l2^2))))%*%
            cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
          
          cov[i, j] = results[1,1]
          cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
              ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
          cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
          cov[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
        }
      }
      return(cov)
    }else{
      cov = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                   nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
      
      for(i in 1:(length(as.matrix(x1))/2)){
        #if(i%%100==0){print(i)}
        theta1=atan2(x1[i,2],x1[i,1])
        for (j in 1:(length(as.matrix(x2))/2)){
          r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                    sum((x1[i,])^2)^.5)
          r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                    sum((x2[j,])^2)^.5)
          dist=(r1-r2)
          theta2=atan2(x2[j,2],x2[j,1])
          results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
            diag(c(matern.covariance(dist,l1,nu1,sigma1),matern.covariance(dist,l2,nu2,sigma2)))%*%
            cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
          
          cov[i, j] = results[1,1]
          cov[nrow(x1) + i, nrow(x2) + j] = results[2,2]
          cov[nrow(x1) + i, j] = results[2,1]
          cov[i, nrow(x2) + j] = results[1,2]
        }
      }
      return(cov)
      
    }
    
    
  }else{
    
    if(nu1==0){
      cov = matrix(0, ncol = 2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                   nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
      
      for(i in 1:(length(as.matrix(x1))/2)){
        if(i%%100==0){print(i)}
        for (j in 1:(length(as.matrix(x2))/2)){
          if(length(as.matrix(x1)) == 2){
            repr1 = function(theta1,theta2) {
              sigma1^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1)-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2))^2)/(l1^2))
            }
            
            repr2 = function(theta1,theta2) {
              sigma2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1)-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2))^2)/(l2^2))
            }
          }else{
            repr1 = function(theta1,theta2) {
              sigma1^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2)/(l1^2))
            }
            
            repr2 = function(theta1,theta2) {
              sigma2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2)/(l2^2))
            }
          }
          
          integrand1 = function(theta) {
            theta1=theta[1]
            theta2=theta[2]
            repr1(theta1,theta2)*cos(theta1)*cos(theta2)+ repr2(theta1,theta2)*sin(theta1)*sin(theta2)
          }
          
          integrand2 = function(theta) {
            
            theta1=theta[1]
            theta2=theta[2]
            -repr1(theta1,theta2)*cos(theta1)*sin(theta2)+ repr2(theta1,theta2)*sin(theta1)*cos(theta2)
          }
          
          integrand3 = function(theta) {
            theta1=theta[1]
            theta2=theta[2]
            -repr1(theta1,theta2)*sin(theta1)*cos(theta2)+ repr2(theta1,theta2)*sin(theta2)*cos(theta1)
          }
          
          integrand4 = function(theta) {
            theta1=theta[1]
            theta2=theta[2]
            repr2(theta1,theta2)*cos(theta1)*cos(theta2)+ repr1(theta1,theta2)*sin(theta1)*sin(theta2)
          }
          fun=list(integrand1,integrand2,integrand3,integrand4)
          
          results=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                          .export = c("l1", "l2", "sigma1", "sigma2",
                                      "repr1", "repr2", "integrand1",
                                      "integrand2", "integrand3", "integrand4","nu1","nu2"))%dopar%{
                                        adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
                                      }
          
          
          cov[i, j] = results[1]
          cov[nrow(x1) + i, nrow(x2) + j] = results[4]
          cov[nrow(x1) + i, j] = results[3]
          cov[i, nrow(x2) + j] = results[2]
        }
      }
      
      
      return(cov)
    }else{
      cov = matrix(0, ncol = 2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                   nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
      
      for(i in 1:(length(as.matrix(x1))/2)){
        #if(i%%100==0){print(i)}
        for (j in 1:(length(as.matrix(x2))/2)){
          if(length(as.matrix(x1)) == 2){
            repr1 = function(theta1,theta2) {
              matern.covariance(sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1)-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2))^2),
                                l1,nu1,sigma1)
            }
            repr2 = function(theta1,theta2) {
              matern.covariance(sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1)-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2))^2),
                                l2,nu2,sigma2)
            }
          }else{
            repr1 = function(theta1,theta2) {
              matern.covariance(sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2),
                                l1,nu1,sigma1)
            }
            
            repr2 = function(theta1,theta2) {
              matern.covariance(sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2),
                                l2,nu2,sigma2)
            }
          }
          
          
          integrand1 = function(theta) {
            theta1=theta[1]
            theta2=theta[2]
            repr1(theta1,theta2)*cos(theta1)*cos(theta2)+ repr2(theta1,theta2)*sin(theta1)*sin(theta2)
          }
          
          integrand2 = function(theta) {
            
            theta1=theta[1]
            theta2=theta[2]
            -repr1(theta1,theta2)*cos(theta1)*sin(theta2)+ repr2(theta1,theta2)*sin(theta1)*cos(theta2)
          }
          
          integrand3 = function(theta) {
            theta1=theta[1]
            theta2=theta[2]
            -repr1(theta1,theta2)*sin(theta1)*cos(theta2)+ repr2(theta1,theta2)*sin(theta2)*cos(theta1)
          }
          
          integrand4 = function(theta) {
            theta1=theta[1]
            theta2=theta[2]
            repr2(theta1,theta2)*cos(theta1)*cos(theta2)+ repr1(theta1,theta2)*sin(theta1)*sin(theta2)
          }
          fun=list(integrand1,integrand2,integrand3,integrand4)
          
          results=foreach(l=1:4,.packages = c("cubature","rSPDE"), .combine = cbind,
                          .export = c("l1", "l2", "sigma1", "sigma2",
                                      "repr1", "repr2", "integrand1",
                                      "integrand2", "integrand3", "integrand4","nu1","nu2"))%dopar%{
                                        adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
                                      }
          
          
          cov[i, j] = results[1]
          cov[nrow(x1) + i, nrow(x2) + j] = results[4]
          cov[nrow(x1) + i, j] = results[3]
          cov[i, nrow(x2) + j] = results[2]
        }
      }
      
      
      return(cov)
    }
    
  }
  
}

grad_eq=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs_2,maxEval=1000){
  
  
  K=eq_cov_mat(xtr,xtr,l1_2,sigma1_2,l2_2,sigma2_2,fundamental=0,maxEval = maxEval)+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  
  
  K_l1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_l2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  for(i in 1:(length(as.matrix(xtr))/2)){
    #if(i%%100==0){print(i)}
    for (j in 1:(length(as.matrix(xtr))/2)){
      
      repr1 = function(theta1,theta2) {
        sigma1_2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), 
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)),
                                         c(-sin(theta2), cos(theta2)))%*%
                                   as.numeric(xtr[j,]))^2)/(l1_2^2))
      }
      
      repr2 = function(theta1,theta2){
        sigma2_2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), 
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)),
                                         c(-sin(theta2), cos(theta2)))%*%
                                   as.numeric(xtr[j,]))^2)/(l2_2^2))
        
      } 
      
      
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*cos(theta1)*cos(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(l1_2^3)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*cos(theta1)*sin(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(l1_2^3)
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*sin(theta1)*cos(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(l1_2^3)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*sin(theta1)*sin(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(l1_2^3)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results1=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                       .export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0),
                                                    upper=c(2 * pi,2*pi), 
                                                    maxEval = maxEval)$integral
                                   }
      
      repr1 = function(theta1,theta2) {
        2*sigma1_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)),
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)),
                                         c(-sin(theta2), cos(theta2)))%*%
                                   as.numeric(xtr[j,]))^2)/(l1_2^2))
      }
      
      
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*cos(theta1)*cos(theta2)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*cos(theta1)*sin(theta2)
        
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*sin(theta1)*cos(theta2)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*sin(theta1)*sin(theta2)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results2=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                       .export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0),
                                                    upper=c(2 * pi,2*pi), 
                                                    maxEval =maxEval)$integral
                                   }
      
      
      
      repr1 = function(theta1,theta2) {
        sigma1_2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), 
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)),
                                         c(-sin(theta2), cos(theta2)))%*%
                                   as.numeric(xtr[j,]))^2)/(l1_2^2))
      }
      
      repr2 = function(theta1,theta2){
        sigma2_2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)),
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)),
                                         c(-sin(theta2), cos(theta2)))%*%
                                   as.numeric(xtr[j,]))^2)/(l2_2^2))
        
      } 
      
      
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*sin(theta1)*sin(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(l2_2^3)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta2)*sin(theta1)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(l2_2^3)
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*sin(theta2)*cos(theta1)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(l2_2^3)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta1)*cos(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(l2_2^3)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results3=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                       .export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0), 
                                                    upper=c(2 * pi,2*pi), 
                                                    maxEval = maxEval)$integral
                                   }
      
      
      
      repr2 = function(theta1,theta2) {
        2*sigma2_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), 
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)),
                                         c(-sin(theta2),cos(theta2)))%*%
                                   as.numeric(xtr[j,]))^2)/(l2_2^2))
      }
      
      
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*sin(theta1)*sin(theta2)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*sin(theta1)*cos(theta2)
        
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta1)*sin(theta2)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta1)*cos(theta2)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results4=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                       .export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0),
                                                    upper=c(2 * pi,2*pi),
                                                    maxEval = maxEval)$integral
                                   }
      
      K_l1_2[i, j] = results1[1]
      K_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[4]
      K_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[3]
      K_l1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2]
      
      
      
      K_sigma1_2[i, j] = results2[1]
      K_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[4]
      K_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[3]
      K_sigma1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2]
      
      
      
      K_l2_2[i, j] = results3[1]
      K_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[4]
      K_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[3]
      K_l2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2]
      
      
      K_sigma2_2[i, j] = results4[1]
      K_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[4]
      K_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[3]
      K_sigma2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2]
      
    }
    
    
  }
  inv_K=solve(K)
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l1_2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_l1_2),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma1_2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_sigma1_2),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l2_2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_l2_2),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma2_2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_sigma2_2),
             -2*sigma_obs_2*t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  
  return(grad)
}

log_likelihood_grad=function(ytr,xtr,l1,sigma1,l2,sigma2,sigma_obs,
                             equivariant=FALSE,fundamental=FALSE,nu1=0,nu2=nu1){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  if(equivariant){
    if(fundamental){
      Ktr=eq_cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,fundamental = 1,nu1=nu1,nu2=nu2)+
        sigma_obs*diag(nrow=2*nrow(xtr))
    }else{
      Ktr=eq_cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,nu1=nu1,nu2=nu2)+
        sigma_obs*diag(nrow=2*nrow(xtr))
    }
  }else{
    Ktr=cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,nu1=nu1,nu2=nu2)+
      sigma_obs*diag(nrow=2*nrow(xtr))
  }
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}

if(id==2){
  set.seed(123)
  bd=0.5
  nx=seq(-2,2,l=20)
  ny=nx
  
  Xte=expand.grid(nx,ny)
  
  Yte=cbind(Xte[,1]/(bd+(Xte[,1]^2+Xte[,2]^2)^2),Xte[,2]/(bd+(Xte[,1]^2+Xte[,2]^2)^2))
  
  
  #Xtr=cbind(runif(10,-2,2),runif(10,-2,2))
  Xtr=cbind(
    c(-0.35,0.88,-0.06,0.105,0.27,1.35,-0.03,-0.015,-0.3,-0.11)*2/1.4,
    2*c(0.95,0.63,0.2,0.065,0,0,-0.095,-0.2,-0.85,-0.95)
  )
  
  Ytr=cbind(Xtr[,1]/(bd+(Xtr[,1]^2+Xtr[,2]^2)^2),Xtr[,2]/(bd+(Xtr[,1]^2+Xtr[,2]^2)^2))
  
  
  sigma_obs=.15
  noise=rnorm(length(Xtr),0,sigma_obs)
  Ytr=Ytr+cbind(noise[1:(length(Xtr)/2)],noise[(length(Xtr)/2+1):length(Xtr)])
  
  lr=0.2
  max_iter=20
  initial_par=rep(1,5)
}else{
  set.seed(123)
  Xtr=cbind(
    -c(1.45,1.1,1.5,1.5,1.28,1.15,1.53,1.27)/1.6,
    1/1.2*c(1.13,1,0.88,0.32,0.02,-0.45,-0.88,-1.07)
  )
  sigma_obs=0.1
  noise=rnorm(16,0,sigma_obs)
  Ytr=cbind(-Xtr[,2],Xtr[,1])+cbind(noise[1:8],noise[9:16])
  
  
  
  ngrid=17
  
  ny=seq(-1,1,l=ngrid)
  nx=seq(-1,1,l=ngrid)
  
  Xte=expand.grid(nx,ny)
  Yte=cbind(-Xte[,2],Xte[,1])
  max_iter=20
  lr=1/7
  initial_par=c(1,1,2.7,0.368,.135^2)
}


cl=makeCluster(4)
registerDoParallel(cl)

opt_par = Adam(
  function(x) {
    log_likelihood(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],equivariant = 1)
  },
  function(x) {
    grad_eq(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr,print_grad = 0
)^.5



Ktetr_eq=eq_cov_mat(Xte,Xtr,opt_par[1],
                    opt_par[2],
                    opt_par[3],
                    opt_par[4])

Ktetr_trtr_inv_eq=Ktetr_eq%*%solve(eq_cov_mat(Xtr,Xtr,opt_par[1],
                                              opt_par[2],
                                              opt_par[3],
                                              opt_par[4])+
                                     diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_eq=Ktetr_trtr_inv_eq%*%c(Ytr[,1],Ytr[,2])



(rmse_eq=mean(apply((cbind(posterior_mean_eq[1:(0.5*length(posterior_mean_eq))],
                           posterior_mean_eq[(0.5*length(posterior_mean_eq)+1):
                                               (length(posterior_mean_eq))])
                     -Yte)^2,1,sum))^.5)

posterior_cov_eq= eq_cov_mat(Xte,Xte,opt_par[1],
                             opt_par[2],
                             opt_par[3],
                             opt_par[4])-Ktetr_trtr_inv_eq%*%t(Ktetr_eq)

LogS_eq=abs((t(c(Yte[,1],Yte[,2])-posterior_mean_eq)%*%solve(posterior_cov_eq+
                                                               diag(1e-8,2*nrow(Xte)))%*%
                     (c(Yte[,1],Yte[,2])-posterior_mean_eq)+determinant(posterior_cov_eq+
                                                                                diag(1e-8,2*nrow(Xte))
                     )$modulus[1]+
                     log(2*pi)*nrow(Xte))/nrow(Xte))




res=list(opt_par,Ytr,posterior_mean_eq,posterior_cov_eq,LogS_eq)


save(list = "res", file = paste0("rmses_eq_", id, ".rda"))
