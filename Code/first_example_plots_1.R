library(foreach)
library(doParallel)

library(nloptr)
library(cubature)

library(proxy)
library(Matrix)
library(pracma)

lr=0.01;max_iter=1000

test=0

if(test){
  max_iter=2
}

Adam = function(f, grad_f, params_init, lr = .01, beta1 = 0.9, beta2 = 0.999,
                epsilon = 1e-8, max_iter = 100, tol = 1e-10,
                lb=rep(0,length(params_init)),
                ub=rep(Inf,length(params_init))) {
    # Number of iterations (for example)
   
    params = params_init
    
    m = numeric(length(params))
    v = numeric(length(params))
    
    
    iter = 0
    
    
    while (iter < max_iter) {
        iter = iter + 1
        
        grad = grad_f(params)
        m = beta1 * m + (1 - beta1) * grad
        
        v = beta2 * v + (1 - beta2) * (grad^2)
        
        m_hat = m / (1 - beta1^iter)
        v_hat = v / (1 - beta2^iter)
        params_new = params - lr * m_hat / (sqrt(v_hat) + epsilon)
        
        params_new[is.na(params_new)]=params[is.na(params_new)]
        
        params_new[params_new<=lb]=params[params_new<=lb]
        params_new[params_new>ub]=params[params_new>ub]
        
        #params_new = pmax(params_new, tol)
        
        if (max(abs(params_new - params)) < tol) {
            break
        }
        
        
        params = params_new
        #print(params)
        #print(c(f(params)))
        
        
    }
    
    
    #print("iter= ")
    #print(iter)
    return(params)
}

eq_cov_mat = function(x1, x2, l1, sigma1, l2, sigma2, maxEval = 1000) {
  cov = matrix(0, ncol = 2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                   nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
      
      for(i in 1:(length(as.matrix(x1))/2)){
        #if(i%%100==0){print(i)}
        for (j in 1:(length(as.matrix(x2))/2)){
          if(length(as.matrix(x1)) == 2){
            repr1 = function(theta1,theta2) {
              sigma1^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1)-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2))^2)*(l1^2))
            }
            
            repr2 = function(theta1,theta2) {
              sigma2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1)-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2))^2)*(l2^2))
            }
          }else{
            repr1 = function(theta1,theta2) {
              sigma1^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2)*(l1^2))
            }
            
            repr2 = function(theta1,theta2) {
              sigma2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2)*(l2^2))
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
                                      "integrand2", "integrand3", "integrand4"))%dopar%{
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

log_likelihood_eq=function(ytr,xtr,l1,sigma1,l2,sigma2,sigma_obs,
                             maxEval=1000){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  Ktr=eq_cov_mat(xtr,xtr,l1,sigma1,l2,sigma2,maxEval = maxEval)+
        sigma_obs^2*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}
grad_eq=function(xtr,ytr,l1,sigma1,l2,sigma2,sigma_obs,maxEval=1000){
  
  
  K=eq_cov_mat(xtr,xtr,l1,sigma1,l2,sigma2,maxEval = maxEval)+
    diag(sigma_obs^2,nrow=2*nrow(xtr))
  
  
  
  K_l1 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_l2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma1 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  for(i in 1:(length(as.matrix(xtr))/2)){
    #if(i%%100==0){print(i)}
    for (j in 1:(length(as.matrix(xtr))/2)){
      
      repr1 = function(theta1,theta2) {
        sigma1^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), 
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                                 cbind(c(cos(theta2), sin(theta2)),
                                       c(-sin(theta2), cos(theta2)))%*%
                                   as.numeric(xtr[j,]))^2)*(l1^2))
      }
      
      repr2 = function(theta1,theta2){
        sigma2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), 
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                                 cbind(c(cos(theta2), sin(theta2)),
                                       c(-sin(theta2), cos(theta2)))%*%
                                   as.numeric(xtr[j,]))^2)*(l2^2))
        
      } 
      
      
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*cos(theta1)*cos(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)*(-l1)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*cos(theta1)*sin(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)*(-l1)
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*sin(theta1)*cos(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)*(-l1)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*sin(theta1)*sin(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)*(-l1)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results1=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                       .export = c("l1", "l2", "sigma1", "sigma2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0),
                                                    upper=c(2 * pi,2*pi), 
                                                    maxEval = maxEval)$integral
                                   }
      
      repr1 = function(theta1,theta2) {
        2*sigma1*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)),
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                        cbind(c(cos(theta2), sin(theta2)),
                              c(-sin(theta2), cos(theta2)))%*%
                          as.numeric(xtr[j,]))^2)*(l1^2))
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
                       .export = c("l1", "l2", "sigma1", "sigma2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0),
                                                    upper=c(2 * pi,2*pi), 
                                                    maxEval =maxEval)$integral
                                   }
      
      
      
      repr1 = function(theta1,theta2) {
        sigma1^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), 
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                                 cbind(c(cos(theta2), sin(theta2)),
                                       c(-sin(theta2), cos(theta2)))%*%
                                   as.numeric(xtr[j,]))^2)*(l1^2))
      }
      
      repr2 = function(theta1,theta2){
        sigma2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)),
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                                 cbind(c(cos(theta2), sin(theta2)),
                                       c(-sin(theta2), cos(theta2)))%*%
                                   as.numeric(xtr[j,]))^2)*(l2^2))
        
      } 
      
      
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*sin(theta1)*sin(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)*(-l2)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta2)*sin(theta1)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)*(-l2)
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*sin(theta2)*cos(theta1)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)*(-l2)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta1)*cos(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)*(-l2)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results3=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                       .export = c("l1", "l2", "sigma1", "sigma2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0), 
                                                    upper=c(2 * pi,2*pi), 
                                                    maxEval = maxEval)$integral
                                   }
      
      
      
      repr2 = function(theta1,theta2) {
        2*sigma2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), 
                                       c(-sin(theta1), cos(theta1)))%*%
                                   as.numeric(xtr[i,])-
                        cbind(c(cos(theta2), sin(theta2)),
                              c(-sin(theta2),cos(theta2)))%*%
                          as.numeric(xtr[j,]))^2)*(l2^2))
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
                       .export = c("l1", "l2", "sigma1", "sigma2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0),
                                                    upper=c(2 * pi,2*pi),
                                                    maxEval = maxEval)$integral
                                   }
      
      K_l1[i, j] = results1[1]
      K_l1[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[4]
      K_l1[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[3]
      K_l1[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2]
      
      
      
      K_sigma1[i, j] = results2[1]
      K_sigma1[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[4]
      K_sigma1[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[3]
      K_sigma1[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2]
      
      
      
      K_l2[i, j] = results3[1]
      K_l2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[4]
      K_l2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[3]
      K_l2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2]
      
      
      K_sigma2[i, j] = results4[1]
      K_sigma2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[4]
      K_sigma2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[3]
      K_sigma2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2]
      
    }
    
    
  }
  inv_K=solve(K)
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_l1),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_sigma1),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_l2),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_sigma2),
             -2*sigma_obs*t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs*inv_K))
  
  return(grad)
}


Xtr=cbind(
  -c(1.45,1.1,1.5,1.5,1.28,1.15,1.53,1.27)/1.6,
  1/1.2*c(1.13,1,0.88,0.32,0.02,-0.45,-0.88,-1.07)
)
set.seed(12)
sigma_obs=0.15
noise=rnorm(16,0,sigma_obs)
Ytr=cbind(-Xtr[,2],Xtr[,1])+cbind(noise[1:8],noise[9:16])



ngrid=17

ny=seq(-1,1,l=ngrid)
nx=seq(-1,1,l=ngrid)

Xte=expand.grid(nx,ny)
Yte=cbind(-Xte[,2],Xte[,1])

if(test){
 Xtr=Xtr[1:2,];Ytr=Ytr[1:2,];Xte=Xte[1:2,];Yte=Yte[1:2,]
}

cl=makeCluster(4)
registerDoParallel(cl)


initial_sigma1=initial_sigma2=initial_sigma3=initial_sigma4=initial_sigma5=
  initial_sigma6=initial_sigma7=initial_sigma8=initial_sigma9=initial_sigma10=
  initial_sigma11=initial_sigma12=initial_l1=initial_l2=initial_l3=initial_l4=
  initial_l5=initial_l6=initial_range1=initial_l_range1=initial_range2=initial_l_range2 =1

initial_sigma_obs=.1

opt_par=Adam(function(x) {
  log_likelihood_eq(
    Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
},function(x) {
  grad_eq(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
},c(initial_sigma1,
                 initial_sigma2,
                 initial_l1,
                 initial_l2,initial_sigma_obs),lr=lr,max_iter = max_iter)

Ktetr_eq=eq_cov_mat(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])

Ktetr_trtr_inv_eq=Ktetr_eq%*%inv(eq_cov_mat(Xtr,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])
                                               
                                               +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_eq=Ktetr_trtr_inv_eq%*%c(Ytr[,1],Ytr[,2])



(rmse_eq=mean(apply((cbind(
  posterior_mean_eq[1:(0.5*length(posterior_mean_eq))],
                                 posterior_mean_eq
  [(0.5*length(posterior_mean_eq)+1):(length(posterior_mean_eq))])
                           -Yte)^2,1,sum))^.5)



posterior_cov_eq= eq_cov_mat(Xte,Xte,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])-Ktetr_trtr_inv_eq%*%t(Ktetr_eq)


(LogS_eq=(t(c(Yte[,1],Yte[,2])-posterior_mean_eq)%*%
                  inv(posterior_cov_eq+diag(opt_par[5]^2,2*nrow(Xte)))%*%
                  (c(Yte[,1],Yte[,2])-posterior_mean_eq)
                +determinant(posterior_cov_eq+diag(opt_par[5]^2,2*nrow(Xte))
                             )$modulus[1]+
                  log(2*pi)*nrow(Xte))/nrow(Xte))

res=list(cbind(
  posterior_mean_eq[1:(0.5*length(posterior_mean_eq))],
                                 posterior_mean_eq
  [(0.5*length(posterior_mean_eq)+1):(length(posterior_mean_eq))]),
posterior_cov_eq,opt_par)


save(list = "res", file = paste0("Equivariant_results_first_examples_1.rda"))
