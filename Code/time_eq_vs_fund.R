id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

library(foreach)
library(doParallel)

library(nloptr)
library(cubature)

library(proxy)
library(Matrix)
library(pracma)

set.seed(12)
norm=function(x){sum(x^2)^(.5)}

n_train=c(5,10,50,100)
n_test=c(10,50,100,250,500)

grid1=expand.grid(n_train,n_test)
n_train=grid1[id,1];n_test=grid1[id,2]


cov_mat=function(x1,x2,sigma1,sigma2,l1,l2){
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  x1=as.matrix(x1);x2=as.matrix(x2)
  if(n1==1){
    if(n2==1){
      dist_mat=norm(x1-x2)
    }else{dist_mat=distmat(t(x1),x2)
    }
  }else{
    dist_mat=distmat(x1,x2)
  }
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  cov[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))
  cov[(n1+1):(2*n1),1:n2]=cov[1:n1,(n2+1):(2*n2)]=0
  
  
  return(cov)
}

cov_mat_helm=function(x1,x2,l1,sigma1,l2,sigma2){
  x1=as.matrix(x1);x2=as.matrix(x2)
  dist_mat=distmat(x1,x2)
  diff_mat=outer(x1[, 1], x2[, 1], "-") * outer(x1[, 2], x2[, 2], "-")
  dist_mat_x1=distmat(cbind(x1[,1],0),cbind(x2[,1],0))
  dist_mat_x2=distmat(cbind(x1[,2],0),cbind(x2[,2],0))
  
  n1=nrow(x1)
  n2=nrow(x2)
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  
  cov[1:n1,1:n2]=sigma1^2/((l1)^4)*exp(-.5*dist_mat^2/(l1^2))*(l1^2-dist_mat_x1^2)+
    sigma2^2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2))*(l2^2-dist_mat_x2^2)
  
  
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma1^2/((l1)^4)*exp(-.5*dist_mat^2/(l1^2))*
    (l1^2-dist_mat_x2^2)+
    sigma2^2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2))*(l2^2-dist_mat_x1^2)
  
  cov[(n1+1):(2*n1),1:n2]=diff_mat*(-exp(-.5*dist_mat^2/(l1^2))*sigma1^2/((l1)^4)+
                                      sigma2^2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2)) )
  
  
  cov[1:n1,(n2+1):(2*n2)]=cov[(n1+1):(2*n1),1:n2]
  
  
  return(cov)
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

fund_cov_mat=fund_cov_mat_LMC = function(x1, x2, sigma1,sigma2,l1,l2) {
  x1=as.matrix(x1);x2=as.matrix(x2)
  cov = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1,nrow(x2)),
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
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma1,sigma2,l1,l2)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard %*%cbind(c(cos(theta2),-sin(theta2)),
                              c(sin(theta2),cos(theta2)))
      
      cov[i, j] = results[1,1]
      cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
    }
  }
  return(cov)
  
  
}

#set.seed(123)

Xtr=matrix(runif(2*n_train,-1,1),ncol=2)
Xte=matrix(runif(2*n_test,-1,1),ncol=2)
sigma_obs=0.15



noise=matrix(rnorm(2*n_train,0,sigma_obs),ncol=2)
Ytr=cbind(-Xtr[,2],Xtr[,1])+noise
Yte=cbind(-Xte[,2],Xte[,1])

opt_par=c(rep(1,4),0.1)

cl=makeCluster(4)
registerDoParallel(cl)

tic()
Ktetr_eq=eq_cov_mat(Xte,Xtr,opt_par[1],
                    opt_par[2],
                    opt_par[3],
                    opt_par[4])
time_te_tr_eq=toc()

tic()
Ktrtr_eq=eq_cov_mat(Xtr,Xtr,opt_par[1],
                    opt_par[2],
                    opt_par[3],
                    opt_par[4])
time_tr_tr_eq=toc()
Ktetr_trtr_inv_eq=Ktetr_eq%*%inv(Ktrtr_eq
                                 
                                 +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))


posterior_mean_eq=Ktetr_trtr_inv_eq%*%c(Ytr[,1],Ytr[,2])



(rmse_eq=mean(apply((cbind(
  posterior_mean_eq[1:(0.5*length(posterior_mean_eq))],
  posterior_mean_eq
  [(0.5*length(posterior_mean_eq)+1):(length(posterior_mean_eq))])
  -Yte)^2,1,sum))^.5)


tic()
posterior_cov_eq= eq_cov_mat(Xte,Xte,opt_par[1],
                             opt_par[2],
                             opt_par[3],
                             opt_par[4])-Ktetr_trtr_inv_eq%*%t(Ktetr_eq)
time_te_te_eq=toc()

total_time_eq=time_te_te_eq+time_te_tr_eq+time_tr_tr_eq




############
tic()
Ktetr_fund=fund_cov_mat(Xte,Xtr,opt_par[1],
                    opt_par[2],
                    opt_par[3],
                    opt_par[4])
time_te_tr_fund=toc()

tic()
Ktrtr_fund=fund_cov_mat(Xtr,Xtr,opt_par[1],
                    opt_par[2],
                    opt_par[3],
                    opt_par[4])
time_tr_tr_fund=toc()
Ktetr_trtr_inv_fund=Ktetr_fund%*%inv(Ktrtr_fund
                                 
                                 +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))


posterior_mean_fund=Ktetr_trtr_inv_fund%*%c(Ytr[,1],Ytr[,2])



(rmse_fund=mean(apply((cbind(
  posterior_mean_fund[1:(0.5*length(posterior_mean_fund))],
  posterior_mean_fund
  [(0.5*length(posterior_mean_fund)+1):(length(posterior_mean_fund))])
  -Yte)^2,1,sum))^.5)


tic()
posterior_cov_fund= fund_cov_mat(Xte,Xte,opt_par[1],
                             opt_par[2],
                             opt_par[3],
                             opt_par[4])-Ktetr_trtr_inv_fund%*%t(Ktetr_fund)
time_te_te_fund=toc()

total_time_fund=time_te_te_fund+time_te_tr_fund+time_tr_tr_fund

res=list(n_train,n_test,c(time_te_te_fund,time_te_tr_fund,time_tr_tr_fund,total_time_fund),
         c(time_te_te_eq,time_te_tr_eq,time_tr_tr_eq,total_time_eq))

save(list = "res", file = paste0("time_fund_eq", id, ".rda"))



