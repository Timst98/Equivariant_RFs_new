id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#id=26
lr=c(1e-5,1e-4,1e-3,1e-2,1/(20:3))
max_iter=c(seq(10,100,by=5),seq(100,1000,by=10))
grid=expand.grid(lr,max_iter)

lr=grid[id,1]
max_iter=grid[id,2]


size = 30

x1 = seq(-4, 4, length.out = size)
x2 = seq(-4, 4, length.out = size)
L = expand.grid(x1, x2)

# Define functions
F_tilde = function(x) {
  c(x[2], (x[1] - 0.1 * (x[1])^3) * (1 + 0.1 * cos(50 * pi * x[1])))
}

D = function(x, bd, Rd) {
  c((x[1] - (-3)) * bd / Rd, (x[2] - 0) * bd / Rd)
}

C = function(x, bc, Rc) {
  c(-(x[1] - 3) * bc / Rc, -(x[2] - 0) * bc / Rc)
}

# Compute velocities
bd = .5 # Parameter governing size of area of divergence
bc = .5  # Parameter governing size of area of convergence

Yte= t(apply(L, 1, function(x) {
  Rd = sqrt((x[1] - (-3))^2 + (x[2] - 0)^2)
  Rc = sqrt((x[1] - 3)^2 + (x[2] - 0)^2)
  
  F_tilde_val = F_tilde(x)
  D_val = D(x, bd, Rd)
  C_val = C(x, bc, Rc)
  
  F_tilde_val + D_val + C_val
}))

train_points=train_points=list(c(-3.8,3.5),c(-2,3.5),c(-3.4,1),c(-3.5,1.3),
                               c(-3.01,0.1),c(-3,0.05),c(-3.1,-0.5),c(-2.9,-0.7),
                               c(-2.85,-0.75),c(-2,-1.2),c(-2,1),c(-3.25,-0.7),
                               c(-3.3,-0.6),c(-2.2,-2),c(-3.3,-3),c(3,-0.25),
                               c(2.9,-0.1),c(2.9,0.25),c(3.05,0.22),c(3.3,-0.05),
                               c(2,-0.25),c(2.75,1.05),c(3.05,1),c(3.5,0.3),
                               c(3.2,-1.75),c(3.3,-3),c(1.6,-2.8))
Xtr=matrix(ncol=2,nrow=length(train_points))
for(i in 1:length(train_points)){
  Xtr[i,]=train_points[[i]]
}

Ytr=t(apply(Xtr, 1, function(x) {
  Rd = sqrt((x[1] - (-3))^2 + (x[2] - 0)^2)
  Rc = sqrt((x[1] - 3)^2 + (x[2] - 0)^2)
  
  F_tilde_val = F_tilde(x)
  D_val = D(x, bd, Rd)
  C_val = C(x, bc, Rc)
  
  F_tilde_val + D_val + C_val
}))

Xte=L

library(proxy)
library(Matrix)
library(pracma)

Adam = function(f, grad_f, params_init, lr = lr, beta1 = 0.9, beta2 = 0.999,
                epsilon = 1e-8, max_iter = max_iter, tol = 1e-10,alpha=FALSE,alpha2=FALSE) {
  
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
    
    params_new[params_new<=0]=params[params_new<=0]
    #params_new = pmax(params_new, tol)
    
    
    if (max(abs(params_new - params)) < tol) {
      break
    }
    
    if(alpha){
      if((params_new[6]<0)|(params_new[6]>1)){params_new[6]=params[6]}}
    if(alpha2){
      if((params_new[10]<0)|(params_new[10]>1)){params_new[10]=params[10]}}
    
    params = params_new
    print(params)
    print(c(f(params)))
  }
  #print("iter= ")
  #print(iter)
  return(params)
}


cov_mat=function(x1,x2,l1,sigma1,l2,sigma2,nu1=0,nu2=0){
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=dist(t(x1),x2)
  }else{
    dist_mat=dist(x1,x2)
  }
  if(nu1==0){
    cov=matrix(0,nrow=2*n1,ncol=2*n2)
    cov[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2/(l1^2))
    cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*exp(-.5*dist_mat^2/(l2^2))
    
  }else{
    cov=matrix(0,nrow=2*n1,ncol=2*n2)
    cov[1:n1,1:n2]=matern.covariance(dist_mat,l1,nu1,sigma1)
    cov[(n1+1):(2*n1),(n2+1):(2*n2)]=matern.covariance(dist_mat,l2,nu2,sigma2)
    
  }
  return(cov)
}


eq_cov_mat = function(x1, x2, l1, sigma1, l2, sigma2, maxEval = mxevl,
                      nu1=0,nu2=nu1,fundamental=FALSE) {
  
  
  
  if (fundamental){
    if(nu1==0){
      cov = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                   nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
      
      for(i in 1:(length(as.matrix(x1))/2)){
        if(i%%100==0){print(i)}
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
            diag(c(sigma1^2*exp(-.5*(r1-r2)^2/(l1^2)),sigma2^2*exp(-.5*(r1-r2)^2/(l2^2))))%*%
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
        if(i%%100==0){print(i)}
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
        #if(i%%100==0){print(i)}
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
        if(i%%100==0){print(i)}
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



cov_mat_range3=function(x1, x2, l1, sigma1, l2, sigma2,l3, sigma3, l4, sigma4,
                        l5, sigma5, l6, sigma6, range1,range2){
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=dist(t(x1),x2)
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  max_dist_center1=matrix(nrow=n1,ncol=n2)
  for (i in 1:(n1)){
    for(j in 1:(n2)){
      max_dist_center1[i,j]=max(distmat_center1_x1[i],distmat_center1_x2[j])
    }
  }
  sigmoid_max_dist_center1=sigmoid(max_dist_center1-range1)
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  max_dist_center2=matrix(nrow=n1,ncol=n2)
  for (i in 1:(n1)){
    for(j in 1:(n2)){
      max_dist_center2[i,j]=max(distmat_center2_x1[i],distmat_center2_x2[j])
    }
  }
  sigmoid_max_dist_center2=sigmoid(max_dist_center2-range2)
  
  sigmoid_max_dist_center1=cbind(rbind(sigmoid_max_dist_center1,sigmoid_max_dist_center1),
                                 rbind(sigmoid_max_dist_center1,sigmoid_max_dist_center1))
  
  sigmoid_max_dist_center2=cbind(rbind(sigmoid_max_dist_center2,sigmoid_max_dist_center2),
                                 rbind(sigmoid_max_dist_center2,sigmoid_max_dist_center2))
  
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  cov[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2/(l1^2))
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*exp(-.5*dist_mat^2/(l2^2))
  
  cov=cov*(sigmoid_max_dist_center2+sigmoid_max_dist_center1)
  
  
  
  cov_fund1 = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                     nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  cov_fund2=cov_fund1
  for(i in 1:(length(as.matrix(x1))/2)){
    # if(i%%100==0){print(i)}
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
        diag(c(sigma3^2*exp(-.5*(r1-r2)^2/(l3^2)),sigma4^2*exp(-.5*(r1-r2)^2/(l4^2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund1[i, j] = results[1,1]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund1[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
    
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(sigma5^2*exp(-.5*(r1-r2)^2/(l5^2)),sigma6^2*exp(-.5*(r1-r2)^2/(l6^2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund2[i, j] = results[1,1]
      cov_fund2[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund2[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund2[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
    }
  }
  cov_fund1=cov_fund1*(1-sigmoid_max_dist_center1)
  cov_fund2=cov_fund2*(1-sigmoid_max_dist_center2)
  return(cov+cov_fund2+cov_fund1)
  
}

grad_range3=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,l3_2,sigma3_2,l4_2,sigma4_2,
                     l5_2,sigma5_2,l6_2,sigma6_2,
                     range1,range2,sigma_obs_2){
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  
  max_dist_center1=matrix(nrow=n1,ncol=n2)
  del_range1=max_dist_center1
 
  for (i in 1:(n1)){
    for(j in 1:(n2)){
      max_dist_center1[i,j]=max(distmat_center1_x1[i],distmat_center1_x2[j])
      del_range1[i,j]=sigmoid((max_dist_center1[i,j]-range1))^2*(-exp((-max_dist_center1[i,j]+range1)))
     
    }
  }
  sigmoid_max_dist_center1=sigmoid((max_dist_center1-range1))
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  
  max_dist_center2=matrix(nrow=n1,ncol=n2)
  del_range2=max_dist_center2
  del_l_range2=max_dist_center2
  for (i in 1:(n1)){
    for(j in 1:(n2)){
      max_dist_center2[i,j]=max(distmat_center2_x1[i],distmat_center2_x2[j])
      del_range2[i,j]=sigmoid((max_dist_center2[i,j]-range2))^2*(-exp((-max_dist_center2[i,j]+range2)))
      
    }
  }
  
  sigmoid_max_dist_center2=sigmoid((max_dist_center2-range2))
  
  sigmoid_max_dist_center1=cbind(rbind(sigmoid_max_dist_center1,sigmoid_max_dist_center1),
                                 rbind(sigmoid_max_dist_center1,sigmoid_max_dist_center1))
  
  sigmoid_max_dist_center2=cbind(rbind(sigmoid_max_dist_center2,sigmoid_max_dist_center2),
                                 rbind(sigmoid_max_dist_center2,sigmoid_max_dist_center2))
  del_range1=cbind(rbind(del_range1,del_range1),
                   rbind(del_range1,del_range1))
  
  del_range2=cbind(rbind(del_range2,del_range2),
                   rbind(del_range2,del_range2))
 
  
  
  
  
  sigmoid_cov=sigmoid_max_dist_center1+sigmoid_max_dist_center2
  
  K1_l3_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K1_l4_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K1_sigma3_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K1_sigma4_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K2_l5_2=K2_sigma5_2=K2_l6_2=K2_sigma6_2=K1_l3_2
  
  K=cov_mat_range3(xtr,xtr,l1_2,sigma1_2,l2_2,sigma2_2,
                   l3_2,sigma3_2,l4_2,sigma4_2,
                   l5_2,sigma5_2,l6_2,sigma6_2,range1,range2)+
    diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  
  for(i in 1:n1){
    if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[i,])^2)^.5)
      r2=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[j,])^2)^.5)
      #dist=abs(xtr[i,1]-xtr[j,1])+abs(xtr[i,2]-xtr[j,2])
      theta2=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                    atan2(xtr[j,2],xtr[j,1]))
      
      
      
      results1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(sigma3_2^2*exp(-.5*(r1-r2)^2/(l3_2^2))*(r1-r2)^2/(l3_2^3),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,sigma4_2^2*exp(-.5*(r1-r2)^2/(l4_2^2))*(r1-r2)^2/(l4_2^3)))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(2*sigma3_2*exp(-.5*(r1-r2)^2/(l3_2^2)),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,2*sigma4_2*exp(-.5*(r1-r2)^2/(l4_2^2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      
      
      
      K1_l3_2[i, j] = results1[1,1]
      K1_l3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2,2]
      K1_l3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[2,1]
      K1_l3_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[1,2]
      
      
      
      K1_l4_2[i, j] = results2[1,1]
      K1_l4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2,2]
      K1_l4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[2,1]
      K1_l4_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[1,2]
      
      
      
      K1_sigma3_2[i, j] = results3[1,1]
      K1_sigma3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2,2]
      K1_sigma3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[2,1]
      K1_sigma3_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[1,2]
      
      
      
      K1_sigma4_2[i, j] = results4[1,1]
      K1_sigma4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2,2]
      K1_sigma4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[2,1]
      K1_sigma4_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[1,2]
      
      results1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(sigma5_2^2*exp(-.5*(r1-r2)^2/(l3_2^2))*(r1-r2)^2/(l5_2^3),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,sigma6_2^2*exp(-.5*(r1-r2)^2/(l6_2^2))*(r1-r2)^2/(l6_2^3)))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(2*sigma5_2*exp(-.5*(r1-r2)^2/(l5_2^2)),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,2*sigma6_2*exp(-.5*(r1-r2)^2/(l6_2^2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      
      K2_l5_2[i, j] = results1[1,1]
      K2_l5_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2,2]
      K2_l5_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[2,1]
      K2_l5_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[1,2]
      
      
      
      K2_l6_2[i, j] = results2[1,1]
      K2_l6_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2,2]
      K2_l6_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[2,1]
      K2_l6_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[1,2]
      
      
      
      K2_sigma5_2[i, j] = results3[1,1]
      K2_sigma5_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2,2]
      K2_sigma5_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[2,1]
      K2_sigma5_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[1,2]
      
      
      
      K2_sigma6_2[i, j] = results4[1,1]
      K2_sigma6_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2,2]
      K2_sigma6_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[2,1]
      K2_sigma6_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[1,2]
      
      
    }
  }
  
  K1_l3_2=K1_l3_2*(1-sigmoid_max_dist_center1)
  K1_l4_2=K1_l4_2*(1-sigmoid_max_dist_center1)
  K1_sigma3_2=K1_sigma3_2*(1-sigmoid_max_dist_center1)
  K1_sigma4_2=K1_sigma4_2*(1-sigmoid_max_dist_center1)
  
  
  
  K2_l5_2=K2_l5_2*(1-sigmoid_max_dist_center2)
  K2_l6_2=K2_l6_2*(1-sigmoid_max_dist_center2)
  K2_sigma5_2=K2_sigma5_2*(1-sigmoid_max_dist_center2)
  K2_sigma6_2=K2_sigma6_2*(1-sigmoid_max_dist_center2)
  
  if(n1==1){
    dist_mat=dist(t(xtr),xtr)
  }else{
    dist_mat=dist(xtr,xtr)
  }
  del_cov1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov1[1:n1,1:n2]=sigma1_2^2*exp(-.5*dist_mat^2/(l1_2^2))*dist_mat^2/(l1_2^3)
  del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov2[1:n1,1:n2]=0
  del_cov2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2_2^2*exp(-.5*dist_mat^2/(l2_2^2))*dist_mat^2/(l2_2^3)
  
  del_cov3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov3[1:n1,1:n2]=2*sigma1_2*exp(-.5*dist_mat^2/(l1_2^2))
  del_cov3[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov4[1:n1,1:n2]=0
  del_cov4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma2_2*exp(-.5*dist_mat^2/(l2_2^2))
  
  del_cov1=sigmoid_cov*del_cov1
  del_cov2=sigmoid_cov*del_cov2
  del_cov3=sigmoid_cov*del_cov3
  del_cov4=sigmoid_cov*del_cov4
  
  K_fund1=eq_cov_mat(x1,x1,l3_2,sigma3_2,l4_2,sigma4_2,fundamental = 1)
  K_fund2=eq_cov_mat(x1,x1,l5_2,sigma5_2,l6_2,sigma6_2,fundamental = 1)
  
  K_cov=cov_mat(x1,x1,l1_2,sigma1_2,l2_2,sigma2_2)
  
  del_range1=K_fund1*(-del_range1)+K_cov*del_range1
  del_range2=K_fund2*(-del_range2)+K_cov*del_range2
  
  
  inv_K=solve(K)
  b=1
  alpha=1
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_cov1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_cov1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_cov3)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_cov3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_cov2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_cov2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_cov4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_cov4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l3_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_l3_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma3_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma3_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l4_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l4_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma4_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma4_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_l5_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K2_l5_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma5_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K2_sigma5_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_l6_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K2_l6_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma6_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+sum(diag(inv_K%*%(alpha*K2_sigma6_2))),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range2)),
            -t(c(ytr[,1],ytr[,2]))%*%(inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K))
  #print(grad)
  return(grad)
  
}


log_likelihood_grad_range3=function(ytr,xtr,l1,sigma1,l2,sigma2,l3,sigma3,l4,sigma4,
                                    l5,sigma5,l6,sigma6,
                                    range1,range2,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=cov_mat_range3(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,l3^.5,sigma3^.5,l4^.5,sigma4^.5,
                     l5^.5,sigma5^.5,l6^.5,sigma6^.5,range1,range2)+
    sigma_obs*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  print(c(ll))
  return(ll)
  
  
}

log_likelihood_range3=function(ytr,xtr,l1,sigma1,l2,sigma2,l3,sigma3,l4,sigma4,
                                    l5,sigma5,l6,sigma6,
                                    range1,range2,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=cov_mat_range3(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3,l4,sigma4,
                     l5,sigma5,l6,sigma6,range1,range2)+
    sigma_obs*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}



cov_mat_range4=function(x1, x2, l1, sigma1, l2, sigma2,l3, sigma3, l4, sigma4,
                        l5, sigma5, l6, sigma6, range1,range2,l_range1,l_range2){
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=dist(t(x1),x2)
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  max_dist_center1=matrix(nrow=n1,ncol=n2)
  for (i in 1:(n1)){
    for(j in 1:(n2)){
      max_dist_center1[i,j]=max(distmat_center1_x1[i],distmat_center1_x2[j])
    }
  }
  sigmoid_max_dist_center1=sigmoid((max_dist_center1-range1)/l_range1)
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  max_dist_center2=matrix(nrow=n1,ncol=n2)
  for (i in 1:(n1)){
    for(j in 1:(n2)){
      max_dist_center2[i,j]=max(distmat_center2_x1[i],distmat_center2_x2[j])
    }
  }
  sigmoid_max_dist_center2=sigmoid((max_dist_center2-range2)/l_range2)
  
  sigmoid_max_dist_center1=cbind(rbind(sigmoid_max_dist_center1,sigmoid_max_dist_center1),
                                 rbind(sigmoid_max_dist_center1,sigmoid_max_dist_center1))
  
  sigmoid_max_dist_center2=cbind(rbind(sigmoid_max_dist_center2,sigmoid_max_dist_center2),
                                 rbind(sigmoid_max_dist_center2,sigmoid_max_dist_center2))
  
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  cov[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2/(l1^2))
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*exp(-.5*dist_mat^2/(l2^2))
  
  cov=cov*(sigmoid_max_dist_center2+sigmoid_max_dist_center1)
  
  
  
  cov_fund1 = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                     nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  cov_fund2=cov_fund1
  for(i in 1:n1){
    # if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(x1))==2,atan2(x1[2],x1[1]),
                  atan2(x1[i,2],x1[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                sum((x1[i,])^2)^.5)
      r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                sum((x2[j,])^2)^.5)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(sigma3^2*exp(-.5*(r1-r2)^2/(l3^2)),sigma4^2*exp(-.5*(r1-r2)^2/(l4^2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund1[i, j] = results[1,1]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund1[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(sigma5^2*exp(-.5*(r1-r2)^2/(l5^2)),sigma6^2*exp(-.5*(r1-r2)^2/(l6^2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund2[i, j] = results[1,1]
      cov_fund2[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund2[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund2[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
    }
  }
  cov_fund1=cov_fund1*(1-sigmoid_max_dist_center1)
  cov_fund2=cov_fund2*(1-sigmoid_max_dist_center2)
  return(cov+cov_fund2+cov_fund1)
  
}


grad_range4=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,l3_2,sigma3_2,l4_2,sigma4_2,
                     l5_2,sigma5_2,l6_2,sigma6_2,
                     range1,range2,l_range1,l_range2,sigma_obs_2){
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  
  max_dist_center1=matrix(nrow=n1,ncol=n2)
  del_range1=max_dist_center1
  del_l_range1=max_dist_center1
  for (i in 1:(n1)){
    for(j in 1:(n2)){
      max_dist_center1[i,j]=max(distmat_center1_x1[i],distmat_center1_x2[j])
      del_range1[i,j]=sigmoid((max_dist_center1[i,j]-range1)/l_range1)^2*(-exp((-max_dist_center1[i,j]+range1)/l_range1)/l_range1)
      del_l_range1[i,j]=sigmoid((max_dist_center1[i,j]-range1)/l_range1)^2*(exp((-max_dist_center1[i,j]+range1)/l_range1)*
                                                                              (-max_dist_center1[i,j]+range1)/l_range1^2)
      
    }
  }
  sigmoid_max_dist_center1=sigmoid((max_dist_center1-range1)/l_range1)
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  
  max_dist_center2=matrix(nrow=n1,ncol=n2)
  del_range2=max_dist_center2
  del_l_range2=max_dist_center2
  for (i in 1:(n1)){
    for(j in 1:(n2)){
      max_dist_center2[i,j]=max(distmat_center2_x1[i],distmat_center2_x2[j])
      del_range2[i,j]=sigmoid((max_dist_center2[i,j]-range2)/l_range2)^2*(-exp((-max_dist_center2[i,j]+range2)/l_range2)/l_range2)
      del_l_range2[i,j]=sigmoid((max_dist_center2[i,j]-range2)/l_range2)^2*(exp((-max_dist_center2[i,j]+range2)/l_range2)*
                                                                              (-max_dist_center2[i,j]+range2)/l_range2^2)
      
    }
  }
  
  sigmoid_max_dist_center2=sigmoid((max_dist_center2-range2)/l_range2)
  
  sigmoid_max_dist_center1=cbind(rbind(sigmoid_max_dist_center1,sigmoid_max_dist_center1),
                                 rbind(sigmoid_max_dist_center1,sigmoid_max_dist_center1))
  
  sigmoid_max_dist_center2=cbind(rbind(sigmoid_max_dist_center2,sigmoid_max_dist_center2),
                                 rbind(sigmoid_max_dist_center2,sigmoid_max_dist_center2))
  del_range1=cbind(rbind(del_range1,del_range1),
                   rbind(del_range1,del_range1))
  
  del_range2=cbind(rbind(del_range2,del_range2),
                   rbind(del_range2,del_range2))
  
  del_l_range1=cbind(rbind(del_l_range1,del_l_range1),
                     rbind(del_l_range1,del_l_range1))
  
  del_l_range2=cbind(rbind(del_l_range2,del_l_range2),
                     rbind(del_l_range2,del_l_range2))
  
  
  
  
  
  sigmoid_cov=sigmoid_max_dist_center1+sigmoid_max_dist_center2
  
  K1_l3_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K1_l4_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K1_sigma3_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K1_sigma4_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K2_l5_2=K2_sigma5_2=K2_l6_2=K2_sigma6_2=K1_l3_2
  
  K=cov_mat_range4(xtr,xtr,l1_2,sigma1_2,l2_2,sigma2_2,
                   l3_2,sigma3_2,l4_2,sigma4_2,
                   l5_2,sigma5_2,l6_2,sigma6_2,range1,range2,l_range1,l_range2)+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  
  for(i in 1:n1){
    if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[i,])^2)^.5)
      r2=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[j,])^2)^.5)
      #dist=abs(xtr[i,1]-xtr[j,1])+abs(xtr[i,2]-xtr[j,2])
      theta2=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                    atan2(xtr[j,2],xtr[j,1]))
      
      
      
      results1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(sigma3_2^2*exp(-.5*(r1-r2)^2/(l3_2^2))*(r1-r2)^2/(l3_2^3),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,sigma4_2^2*exp(-.5*(r1-r2)^2/(l4_2^2))*(r1-r2)^2/(l4_2^3)))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(2*sigma3_2*exp(-.5*(r1-r2)^2/(l3_2^2)),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,2*sigma4_2*exp(-.5*(r1-r2)^2/(l4_2^2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      
      
      
      K1_l3_2[i, j] = results1[1,1]
      K1_l3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2,2]
      K1_l3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[2,1]
      K1_l3_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[1,2]
      
      
      
      K1_l4_2[i, j] = results2[1,1]
      K1_l4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2,2]
      K1_l4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[2,1]
      K1_l4_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[1,2]
      
      
      
      K1_sigma3_2[i, j] = results3[1,1]
      K1_sigma3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2,2]
      K1_sigma3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[2,1]
      K1_sigma3_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[1,2]
      
      
      
      K1_sigma4_2[i, j] = results4[1,1]
      K1_sigma4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2,2]
      K1_sigma4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[2,1]
      K1_sigma4_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[1,2]
     
      results1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(sigma5_2^2*exp(-.5*(r1-r2)^2/(l3_2^2))*(r1-r2)^2/(l5_2^3),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,sigma6_2^2*exp(-.5*(r1-r2)^2/(l6_2^2))*(r1-r2)^2/(l6_2^3)))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(2*sigma5_2*exp(-.5*(r1-r2)^2/(l5_2^2)),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,2*sigma6_2*exp(-.5*(r1-r2)^2/(l6_2^2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      
      K2_l5_2[i, j] = results1[1,1]
      K2_l5_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2,2]
      K2_l5_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[2,1]
      K2_l5_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[1,2]
      
      
      
      K2_l6_2[i, j] = results2[1,1]
      K2_l6_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2,2]
      K2_l6_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[2,1]
      K2_l6_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[1,2]
      
      
      
      K2_sigma5_2[i, j] = results3[1,1]
      K2_sigma5_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2,2]
      K2_sigma5_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[2,1]
      K2_sigma5_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[1,2]
      
      
      
      K2_sigma6_2[i, j] = results4[1,1]
      K2_sigma6_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2,2]
      K2_sigma6_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[2,1]
      K2_sigma6_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[1,2]
     
      
    }
  }
  
  K1_l3_2=K1_l3_2*(1-sigmoid_max_dist_center1)
  K1_l4_2=K1_l4_2*(1-sigmoid_max_dist_center1)
  K1_sigma3_2=K1_sigma3_2*(1-sigmoid_max_dist_center1)
  K1_sigma4_2=K1_sigma4_2*(1-sigmoid_max_dist_center1)
  
  
  
  K2_l5_2=K2_l5_2*(1-sigmoid_max_dist_center2)
  K2_l6_2=K2_l6_2*(1-sigmoid_max_dist_center2)
  K2_sigma5_2=K2_sigma5_2*(1-sigmoid_max_dist_center2)
  K2_sigma6_2=K2_sigma6_2*(1-sigmoid_max_dist_center2)
  
  if(n1==1){
    dist_mat=dist(t(xtr),xtr)
  }else{
    dist_mat=dist(xtr,xtr)
  }
  del_cov1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov1[1:n1,1:n2]=sigma1_2^2*exp(-.5*dist_mat^2/(l1_2^2))*dist_mat^2/(l1_2^3)
  del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov2[1:n1,1:n2]=0
  del_cov2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2_2^2*exp(-.5*dist_mat^2/(l2_2^2))*dist_mat^2/(l2_2^3)
  
  del_cov3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov3[1:n1,1:n2]=2*sigma1_2*exp(-.5*dist_mat^2/(l1_2^2))
  del_cov3[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov4[1:n1,1:n2]=0
  del_cov4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma2_2*exp(-.5*dist_mat^2/(l2_2^2))
  
  del_cov1=sigmoid_cov*del_cov1
  del_cov2=sigmoid_cov*del_cov2
  del_cov3=sigmoid_cov*del_cov3
  del_cov4=sigmoid_cov*del_cov4
  
  K_fund1=eq_cov_mat(x1,x1,l3_2,sigma3_2,l4_2,sigma4_2,fundamental = 1)
  K_fund2=eq_cov_mat(x1,x1,l5_2,sigma5_2,l6_2,sigma6_2,fundamental = 1)
  
  K_cov=cov_mat(x1,x1,l1_2,sigma1_2,l2_2,sigma2_2)
  
  del_range1=K_fund1*(-del_range1)+K_cov*del_range1
  del_range2=K_fund2*(-del_range2)+K_cov*del_range2
  
  del_l_range1=K_fund1*(-del_l_range1)+K_cov*del_l_range1
  del_l_range2=K_fund2*(-del_l_range2)+K_cov*del_l_range2
  
  
  inv_K=solve(K)
  b=1
  alpha=1
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_cov1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_cov1)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_cov3)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_cov3)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_cov2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_cov2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_cov4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_cov4)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l3_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_l3_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma3_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma3_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l4_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l4_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma4_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma4_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_l5_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K2_l5_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma5_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K2_sigma5_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_l6_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K2_l6_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma6_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+sum(diag(inv_K%*%(alpha*K2_sigma6_2))),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range1)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range1)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range2)),
         
         -t(c(ytr[,1],ytr[,2]))%*%(2*sigma_obs_2*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  #print(grad)
  return(grad)
  
}


log_likelihood_range4=function(ytr,xtr,l1,sigma1,l2,sigma2,l3,sigma3,l4,sigma4,
                                    l5,sigma5,l6,sigma6,
                                    range1,range2,l_range1,l_range2,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=cov_mat_range4(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3,l4,sigma4,
                     l5,sigma5,l6,sigma6,range1,range2,l_range1,l_range2)+
    sigma_obs^2*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*t(ytr)%*%solve(Ktr)%*%ytr-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}

initial_par=c(1,1,2.7,0.368,1,1,2.7,0.368,1,1,2.7,0.368,3,3,0.1)

(opt_par = Adam(
  function(x) {
    log_likelihood_range3(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7],x[8], x[9],x[10],x[11],
      x[12], x[13],x[14],x[15])
  },
  function(x) {
    grad_range3(Xtr, Ytr,x[1], x[2], x[3], x[4], x[5], x[6], x[7],x[8], x[9],x[10],x[11],
                x[12], x[13],x[14],x[15])
  },
  initial_par[1:15],
  max_iter =max_iter,lr=lr
))

posterior_mean_range=cov_mat_range3(Xte,Xtr,opt_par[1],
                                    opt_par[2],
                                    opt_par[3],
                                    opt_par[4],
                                    opt_par[5],
                                    opt_par[6],
                                    opt_par[7],opt_par[8],
                                    opt_par[9],
                                    opt_par[10],
                                    opt_par[11],opt_par[12],
                                    opt_par[13],opt_par[14])%*%
  solve(cov_mat_range3(Xtr,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4],
                       opt_par[5],
                       opt_par[6],
                       opt_par[7],opt_par[8],
                       opt_par[9],
                       opt_par[10],
                       opt_par[11],opt_par[12],
                       opt_par[13],opt_par[14])+
          diag(opt_par[15],nrow=2*nrow(Xtr)),c(Ytr[,1],Ytr[,2]))

posterior_mean_range=cbind(posterior_mean_range[1:(0.5*length(posterior_mean_range))],
                           posterior_mean_range[(0.5*length(posterior_mean_range)+1):(length(posterior_mean_range))])




(rmse_range3=mean(apply((posterior_mean_range-Yte)^2,1,sum))^.5)

initial_par=c(1,1,2.7,0.368,1,1,2.7,0.368,1,1,2.7,0.368,3,3,1,1,0.1)

(opt_par = Adam(
  function(x) {
    log_likelihood_range4(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7],x[8], x[9],x[10],x[11],
      x[12], x[13],x[14],x[15],x[16],x[17])
  },
  function(x) {
    grad_range4(Xtr, Ytr,x[1], x[2], x[3], x[4], x[5], x[6], x[7],x[8], x[9],x[10],x[11],
                x[12], x[13],x[14],x[15],x[16],x[17])
  },
  initial_par,
  max_iter =max_iter,lr=lr
))



posterior_mean_range=cov_mat_range4(Xte,Xtr,opt_par[1],
                                    opt_par[2],
                                    opt_par[3],
                                    opt_par[4],
                                    opt_par[5],
                                    opt_par[6],
                                    opt_par[7],opt_par[8],
                                    opt_par[9],
                                    opt_par[10],
                                    opt_par[11],opt_par[12],
                                    opt_par[13],opt_par[14],
                                    opt_par[15],opt_par[16])%*%
  solve(cov_mat_range4(Xtr,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4],
                       opt_par[5],
                       opt_par[6],
                       opt_par[7],opt_par[8],
                       opt_par[9],
                       opt_par[10],
                       opt_par[11],opt_par[12],
                       opt_par[13],opt_par[14],opt_par[15],opt_par[16])+
          diag(opt_par[17]^2,nrow=2*nrow(Xtr)),c(Ytr[,1],Ytr[,2]))

posterior_mean_range=cbind(posterior_mean_range[1:(0.5*length(posterior_mean_range))],
                           posterior_mean_range[(0.5*length(posterior_mean_range)+1):(length(posterior_mean_range))])


(rmse_range4=mean(apply((posterior_mean_range-Yte)^2,1,sum))^.5)

res=c(lr,max_iter,rmse_range3,rmse_range4)

print(res)

save(list = "res", file = paste0("rmses_LR_MAXITER_range_initial2_", id, ".rda"))


