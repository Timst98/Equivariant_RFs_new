id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#id=2000
lr=c(1e-5,1e-4,1e-3,1e-2,1/(20:3))
max_iter=c(seq(10,100,by=5),seq(100,1000,by=10))
grid=expand.grid(lr,max_iter)

#packageVersion("Matrix")
#packageVersion("pracma")
#packageVersion("proxy")

lr=grid[id,1]
max_iter=grid[id,2]
initial_par=c(rep(1,9),0.55)

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


library(foreach)
library(doParallel)

#library(rSPDE)
library(cubature)

library(proxy)
library(Matrix)
library(pracma)



mxevl=1000


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
    cov[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))
    cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))
    
  }else{
    cov=matrix(0,nrow=2*n1,ncol=2*n2)
    cov[1:n1,1:n2]=matern.covariance(dist_mat,l1,nu1,sigma1)
    cov[(n1+1):(2*n1),(n2+1):(2*n2)]=matern.covariance(dist_mat,l2,nu2,sigma2)
    
  }
  return(cov)
}

cov_mat_helm=function(x1,x2,l1,sigma1,l2,sigma2){
  
  dist_mat=dist(x1,x2)
  diff_mat=outer(x1[, 1], x2[, 1], "-") * outer(x1[, 2], x2[, 2], "-")
  dist_mat_x1=dist(cbind(x1[,1],0),cbind(x2[,1],0))
  dist_mat_x2=dist(cbind(x1[,2],0),cbind(x2[,2],0))
  
  n1=nrow(x1)
  n2=nrow(x2)
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  
  cov[1:n1,1:n2]=sigma1^2/((l1)^4)*exp(-.5*dist_mat^2/(l1^2))*(l1^2-dist_mat_x1^2)+
    sigma2^2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2))*(l2^2-dist_mat_x2^2)
  
  
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma1^2/((l1)^4)*exp(-.5*dist_mat^2/(l1^2))*(l1^2-dist_mat_x2^2)+
    sigma2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2))*(l2^2-dist_mat_x1^2)
  
  cov[(n1+1):(2*n1),1:n2]=diff_mat*(-exp(-.5*dist_mat^2/(l1^2))*sigma1^2/((l1)^4)+
                                      sigma2^2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2)) )
  
  
  cov[1:n1,(n2+1):(2*n2)]=cov[(n1+1):(2*n1),1:n2]
  
  
  return(cov)
}
grad_helm=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs_2){
  
  x1=xtr
  x2=xtr
  
  
  dist_mat=dist(x1,x2)
  diff_mat=outer(x1[, 1], x2[, 1], "-") * outer(x1[, 2], x2[, 2], "-")
  dist_mat_x1=dist(cbind(x1[,1],0),cbind(x2[,1],0))
  dist_mat_x2=dist(cbind(x1[,2],0),cbind(x2[,2],0))
  
  n1=nrow(x1)
  n2=nrow(x2)
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  
  cov[1:n1,1:n2]=sigma1_2/((l1_2)^2)*exp(-.5*dist_mat^2/(l1_2))*(l1_2-dist_mat_x1^2)+
    sigma2_2/((l2_2)^2)*exp(-.5*dist_mat^2/(l2_2))*(l2_2-dist_mat_x2^2)
  
  
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma1_2/((l1_2)^2)*exp(-.5*dist_mat^2/(l1_2))*(l1_2-dist_mat_x2^2)+
    sigma2_2/((l2_2)^2)*exp(-.5*dist_mat^2/(l2_2))*(l2_2-dist_mat_x1^2)
  
  cov[(n1+1):(2*n1),1:n2]=diff_mat*(-exp(-.5*dist_mat^2/(l1_2))*sigma1_2/((l1_2)^2)+
                                      sigma2_2/((l2_2)^2)*exp(-.5*dist_mat^2/(l2_2)) )
  
  
  cov[1:n1,(n2+1):(2*n2)]=cov[(n1+1):(2*n1),1:n2]
  
  
  K=cov+diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  
  dist_mat=dist(x1,x2)
  diff_mat=outer(x1[, 1], x2[, 1], "-") * outer(x1[, 2], x2[, 2], "-")
  dist_mat_x1=dist(cbind(x1[,1],0),cbind(x2[,1],0))
  dist_mat_x2=dist(cbind(x1[,2],0),cbind(x2[,2],0))
  
  n1=nrow(x1)
  n2=nrow(x2)
  
  
  
  del_cov1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov1[1:n1,1:n2]=sigma1_2*exp(-.5*dist_mat^2/(l1_2))/(l1_2^3)*(0.5*dist_mat^2-dist_mat^2*dist_mat_x1^2/(2*l1_2)+2*dist_mat_x1^2-l1_2)
  
  del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma1_2*exp(-.5*dist_mat^2/(l1_2))/(l1_2^3)*(0.5*dist_mat^2-dist_mat^2*dist_mat_x2^2/(2*l1_2)+2*dist_mat_x2^2-l1_2)
  
  del_cov1[(n1+1):(2*n1),1:n2]=sigma1_2*diff_mat*exp(-.5*dist_mat^2/(l1_2))/(l1_2^3)*(2-dist_mat^2/(2*l1_2))
  del_cov1[1:n1,(n2+1):(2*n2)]=del_cov1[(n1+1):(2*n1),1:n2]
  
  
  del_cov2=matrix(0,nrow=2*n1,ncol=2*n2)
  
  del_cov2[1:n1,1:n2]=sigma2_2*exp(-.5*dist_mat^2/(l2_2))/(l2_2^3)*(0.5*dist_mat^2-dist_mat^2*dist_mat_x2^2/(2*l2_2)+2*dist_mat_x2^2-l2_2)
  
  del_cov2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2_2*exp(-.5*dist_mat^2/(l2_2))/(l2_2^3)*(0.5*dist_mat^2-dist_mat^2*dist_mat_x1^2/(2*l2_2)+2*dist_mat_x1^2-l2_2)
  
  del_cov2[(n1+1):(2*n1),1:n2]=sigma2_2*diff_mat*exp(-.5*dist_mat^2/(l2_2))/(l2_2^3)*(dist_mat^2/(2*l2_2)-2)
  
  del_cov2[1:n1,(n2+1):(2*n2)]=del_cov2[(n1+1):(2*n1),1:n2]
  
  
  del_cov3=matrix(0,nrow=2*n1,ncol=2*n2)
  
  del_cov3[1:n1,1:n2]=1/((l1_2)^2)*exp(-.5*dist_mat^2/(l1_2))*(l1_2-dist_mat_x1^2)
  
  del_cov3[(n1+1):(2*n1),(n2+1):(2*n2)]=1/((l1_2)^2)*exp(-.5*dist_mat^2/(l1_2))*(l1_2-dist_mat_x2^2)
  del_cov3[(n1+1):(2*n1),1:n2]=diff_mat/((l1_2)^2)*(-exp(-.5*dist_mat^2/(l1_2)))
  del_cov3[1:n1,(n2+1):(2*n2)]=del_cov3[(n1+1):(2*n1),1:n2]
  
  
  del_cov4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov4[1:n1,1:n2]=1/((l2_2)^2)*exp(-.5*dist_mat^2/(l2_2))*(l2_2-dist_mat_x2^2)
  del_cov4[(n1+1):(2*n1),(n2+1):(2*n2)]=1/((l2_2)^2)*exp(-.5*dist_mat^2/(l2_2))*(l2_2-dist_mat_x1^2)
  del_cov4[(n1+1):(2*n1),1:n2]=diff_mat/((l2_2)^2)*exp(-.5*dist_mat^2/(l2_2)) 
  del_cov4[1:n1,(n2+1):(2*n2)]=del_cov4[(n1+1):(2*n1),1:n2]
  
  
  inv_K=solve(K)
  
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_cov1),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov3%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_cov3),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_cov2),
         
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov4%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_cov4),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K))
  
  return(grad)
  
}


log_likelihood_helm_grad=function(ytr,xtr,l1,sigma1,l2,sigma2,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  Ktr=cov_mat_helm(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5)+sigma_obs*diag(nrow=2*nrow(xtr))
  (ll=-.5*sum((ytr)*solve(Ktr,ytr))-.5*log(det(Ktr))-nrow(xtr)*log(2*pi))
  ll=ifelse(is.na(ll),1e6,-ll)
  #print(ll)
  return(ll)
}


grad_standard=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs_2){
  
  K=cov_mat(xtr,xtr,l1_2^.5,l2_2^.5,sigma1_2^.5,sigma2_2^.5)+
    diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  if(n1==1){
    dist_mat=dist(t(xtr),xtr)
  }else{
    dist_mat=dist(xtr,xtr)
  }
  del_cov1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov1[1:n1,1:n2]=-sigma1_2*exp(-.5*dist_mat^2*(l1_2))*dist_mat^2/2
  del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov2[1:n1,1:n2]=0
  del_cov2[(n1+1):(2*n1),(n2+1):(2*n2)]=-sigma2_2*exp(-.5*dist_mat^2*(l2_2))*dist_mat^2/(2)
  
  del_cov3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov3[1:n1,1:n2]=exp(-.5*dist_mat^2*(l1_2))
  del_cov3[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov4[1:n1,1:n2]=0
  del_cov4[(n1+1):(2*n1),(n2+1):(2*n2)]=exp(-.5*dist_mat^2*(l2_2))
  
  
  inv_K=solve(K)
  
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_cov1),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov3%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_cov3),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_cov2),
         
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov4%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_cov4),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K))
  
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

decomposition_mat2= function(x1, x2, l1, sigma1, l2, sigma2,alpha, maxEval = mxevl) {
  
  eq_cov=eq_cov_mat(x1,x2,l1, sigma1, l2, sigma2)
  cov=cov_mat(x1,x2,l1, sigma1, l2, sigma2)
  return(alpha*cov+(1-2*alpha)*eq_cov)
}

grad_decomp2=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs2,alpha,maxEval=1000){
  
  K_se=cov_mat(xtr,xtr,l1_2^.5,l2_2^.5,sigma1_2^.5,sigma2_2^.5)
  K_eq=eq_cov_mat(xtr,xtr,l1_2^.5,l2_2^.5,sigma1_2^.5,sigma2_2^.5,maxEval = maxEval)
  K= alpha*K_se+(1-2*alpha)*K_eq +diag(sigma_obs2,nrow=2*nrow(xtr))
  inv_K=solve(K)
  del_alpha=K_se-2*K_eq
  
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  if(n1==1){
    dist_mat=dist(t(xtr),xtr)
  }else{
    dist_mat=dist(xtr,xtr)
  }
  del_cov1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov1[1:n1,1:n2]=-sigma1_2*exp(-.5*dist_mat^2*(l1_2))*dist_mat^2/(2)
  del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov2[1:n1,1:n2]=0
  del_cov2[(n1+1):(2*n1),(n2+1):(2*n2)]=-sigma2_2*exp(-.5*dist_mat^2*(l2_2))*dist_mat^2/(2)
  
  del_cov3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov3[1:n1,1:n2]=exp(-.5*dist_mat^2*(l1_2))
  del_cov3[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov4[1:n1,1:n2]=0
  del_cov4[(n1+1):(2*n1),(n2+1):(2*n2)]=exp(-.5*dist_mat^2*(l2_2))
  
  
  
  K_l1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_l2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  for(i in 1:(length(as.matrix(xtr))/2)){
    if(i%%100==0){print(i)}
    for (j in 1:(length(as.matrix(xtr))/2)){
      
      results1=foreach(l=1:4,.packages = c("cubature"), .combine = cbind
                       ,.export = c("maxEval")
                       #c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                       #           "repr1", "repr2", "integrand1","xtr","ytr",
                       #          "integrand2", "integrand3", "integrand4")
                       
      )%dopar%{
        
        repr1 = function(theta1,theta2) {
          sigma1_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)*(l1_2))
        }
        
        repr2 = function(theta1,theta2){
          sigma2_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)*(l2_2))
          
        }
        
        
        integrand1 = function(theta) {
          theta1=theta[1]
          theta2=theta[2]
          -repr1(theta1,theta2)*cos(theta1)*cos(theta2)*
            sum((cbind(c(cos(theta1), sin(theta1)),
                       c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                   as.numeric(xtr[j,]))^2)/(2)
        }
        
        integrand2 = function(theta) {
          
          theta1=theta[1]
          theta2=theta[2]
          repr1(theta1,theta2)*cos(theta1)*sin(theta2)*
            sum((cbind(c(cos(theta1), sin(theta1)),
                       c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                   as.numeric(xtr[j,]))^2)/(2)
        }
        
        integrand3 = function(theta) {
          theta1=theta[1]
          theta2=theta[2]
          repr1(theta1,theta2)*sin(theta1)*cos(theta2)*
            sum((cbind(c(cos(theta1), sin(theta1)),
                       c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                   as.numeric(xtr[j,]))^2)/(2)
        }
        
        integrand4 = function(theta) {
          theta1=theta[1]
          theta2=theta[2]
          -repr1(theta1,theta2)*sin(theta1)*sin(theta2)*
            sum((cbind(c(cos(theta1), sin(theta1)),
                       c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                   as.numeric(xtr[j,]))^2)/(2)
        }
        fun=list(integrand1,integrand2,integrand3,integrand4)
        
        adaptIntegrate(fun[[l]], lower=c(0,0),
                       upper=c(2 * pi,2*pi),
                       maxEval = maxEval)$integral
      }
      
      results2=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,.export = c("maxEval")
                       #  ,.export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                       #             "repr1", "repr2", "integrand1","xtr","ytr",
                       #            "integrand2", "integrand3", "integrand4")
      )%dopar%{
        
        repr1 = function(theta1,theta2) {
          exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                          cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)*(l1_2))
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
        
        
        adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
      }
      
      
      
      results3=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,.export = c("maxEval")
                       #, .export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                       #           "repr1", "repr2", "integrand1","xtr","ytr",
                       #          "integrand2", "integrand3", "integrand4")
                       
      )%dopar%{
        
        repr1 = function(theta1,theta2) {
          sigma1_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)*(l1_2))
        }
        
        repr2 = function(theta1,theta2){
          sigma2_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)*(l2_2))
          
        }
        
        
        integrand1 = function(theta) {
          theta1=theta[1]
          theta2=theta[2]
          -repr2(theta1,theta2)*sin(theta1)*sin(theta2)*
            sum((cbind(c(cos(theta1), sin(theta1)),
                       c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                   as.numeric(xtr[j,]))^2)/(2)
        }
        
        integrand2 = function(theta) {
          
          theta1=theta[1]
          theta2=theta[2]
          -repr2(theta1,theta2)*cos(theta2)*sin(theta1)*
            sum((cbind(c(cos(theta1), sin(theta1)),
                       c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                   as.numeric(xtr[j,]))^2)/(2)
        }
        
        integrand3 = function(theta) {
          theta1=theta[1]
          theta2=theta[2]
          repr2(theta1,theta2)*sin(theta2)*cos(theta1)*
            sum((cbind(c(cos(theta1), sin(theta1)),
                       c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                   as.numeric(xtr[j,]))^2)/(2)
        }
        
        integrand4 = function(theta) {
          theta1=theta[1]
          theta2=theta[2]
          -repr2(theta1,theta2)*cos(theta1)*cos(theta2)*
            sum((cbind(c(cos(theta1), sin(theta1)),
                       c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                   as.numeric(xtr[j,]))^2)/(2)
        }
        fun=list(integrand1,integrand2,integrand3,integrand4)
        
        
        adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
      }
      
      
      
      results4=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,.export = c("maxEval")
                       #,  .export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                       #             "repr1", "repr2", "integrand1","xtr","ytr",
                       #            "integrand2", "integrand3", "integrand4")
      )%dopar%{
        
        
        
        repr2 = function(theta1,theta2) {
          exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                          cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)*(l2_2))
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
        
        
        adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
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
  
  beta=1-2*alpha
  
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(beta*K_l1_2+alpha*del_cov1)%*%inv_K%*%c(ytr[,1],ytr[,2])+
           
           Trace(inv_K%*%(beta*K_l1_2+alpha*del_cov1)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(beta*K_sigma1_2+alpha*del_cov3)%*%inv_K%*%c(ytr[,1],ytr[,2])+
           Trace(inv_K%*%(beta*K_sigma1_2+alpha*del_cov3)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(beta*K_l2_2+alpha*del_cov2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+
           Trace(inv_K%*%(beta*K_l2_2+alpha*del_cov2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(beta*K_sigma2_2+alpha*del_cov4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+
           Trace(inv_K%*%(beta*K_sigma2_2+alpha*del_cov4)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_alpha%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_alpha))
  
  return(grad)
  
  
}

grad_fund=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs_2){
  
  K=eq_cov_mat(xtr,xtr,l1_2^.5,sigma1_2^.5,l2_2^.5,sigma2_2^.5,fundamental=1)+
    diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  K_l1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_l2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  
  for(i in 1:(length(as.matrix(xtr))/2)){
    if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:(length(as.matrix(xtr))/2)){
      r1=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[i,])^2)^.5)
      r2=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[j,])^2)^.5)
      #dist=abs(xtr[i,1]-xtr[j,1])+abs(xtr[i,2]-xtr[j,2])
      theta2=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                    atan2(xtr[j,2],xtr[j,1]))
      
      
      results1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(-sigma1_2*exp(-.5*(r1-r2)^2*(l1_2))*(r1-r2)^2/(2),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,-sigma2_2*exp(-.5*(r1-r2)^2*(l2_2))*(r1-r2)^2/(2)))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(exp(-.5*(r1-r2)^2*(l1_2)),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,exp(-.5*(r1-r2)^2*(l2_2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      
      K_l1_2[i, j] = results1[1,1]
      K_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2,2]
      K_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[2,1]
      K_l1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[1,2]
      
      
      
      K_l2_2[i, j] = results2[1,1]
      K_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2,2]
      K_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[2,1]
      K_l2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[1,2]
      
      
      
      K_sigma1_2[i, j] = results3[1,1]
      K_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2,2]
      K_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[2,1]
      K_sigma1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[1,2]
      
      
      
      K_sigma2_2[i, j] = results4[1,1]
      K_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2,2]
      K_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[2,1]
      K_sigma2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[1,2]
      
      
      
      
    }
  }
  inv_K=solve(K)
  
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l1_2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_l1_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma1_2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_sigma1_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l2_2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_l2_2),
         
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma2_2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_sigma2_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K))
  print(grad)
  return(grad)
  
}


grad_mix=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs_2,alpha){
  
  K1=eq_cov_mat(xtr,xtr,l1_2^.5,sigma1_2^.5,l2_2^.5,sigma2_2^.5,fundamental=1)
  
  K1_l1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K1_l2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K1_sigma1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K1_sigma2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K2=cov_mat(xtr,xtr,l1_2^.5,sigma1_2^.5,l2_2^.5,sigma2_2^.5)
  
  K=(1-2*alpha)*K1+alpha*K2+diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  
  for(i in 1:(length(as.matrix(xtr))/2)){
    if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:(length(as.matrix(xtr))/2)){
      r1=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[i,])^2)^.5)
      r2=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[j,])^2)^.5)
      #dist=abs(xtr[i,1]-xtr[j,1])+abs(xtr[i,2]-xtr[j,2])
      theta2=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                    atan2(xtr[j,2],xtr[j,1]))
      
      
      results1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(-sigma1_2*exp(-.5*(r1-r2)^2*(l1_2))*(r1-r2)^2/(2),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,-sigma2_2*exp(-.5*(r1-r2)^2*(l2_2))*(r1-r2)^2/(2)))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(exp(-.5*(r1-r2)^2*(l1_2)),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,exp(-.5*(r1-r2)^2*(l2_2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      
      K1_l1_2[i, j] = results1[1,1]
      K1_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2,2]
      K1_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[2,1]
      K1_l1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[1,2]
      
      
      
      K1_l2_2[i, j] = results2[1,1]
      K1_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2,2]
      K1_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[2,1]
      K1_l2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[1,2]
      
      
      
      K1_sigma1_2[i, j] = results3[1,1]
      K1_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2,2]
      K1_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[2,1]
      K1_sigma1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[1,2]
      
      
      
      K1_sigma2_2[i, j] = results4[1,1]
      K1_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2,2]
      K1_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[2,1]
      K1_sigma2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[1,2]
      
      
      
      
    }
  }
  
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  if(n1==1){
    dist_mat=dist(t(xtr),xtr)
  }else{
    dist_mat=dist(xtr,xtr)
  }
  del_cov1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov1[1:n1,1:n2]=-sigma1_2*exp(-.5*dist_mat^2*(l1_2))*dist_mat^2/(2)
  del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov2[1:n1,1:n2]=0
  del_cov2[(n1+1):(2*n1),(n2+1):(2*n2)]=-sigma2_2*exp(-.5*dist_mat^2*(l2_2))*dist_mat^2/(2)
  
  del_cov3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov3[1:n1,1:n2]=exp(-.5*dist_mat^2*(l1_2))
  del_cov3[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov4[1:n1,1:n2]=0
  del_cov4[(n1+1):(2*n1),(n2+1):(2*n2)]=exp(-.5*dist_mat^2*(l2_2))
  
  del_alpha=K2-2*K1
  
  inv_K=solve(K)
  b=1-2*alpha
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*K1_l1_2+alpha*del_cov1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(b*K1_l1_2+alpha*del_cov1)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*K1_sigma1_2+alpha*del_cov3)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(b*K1_sigma1_2+alpha*del_cov3)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*K1_l2_2+alpha*del_cov2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(b*K1_l2_2+alpha*del_cov2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*K1_sigma2_2+alpha*del_cov4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(b*K1_sigma2_2+alpha*del_cov4)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_alpha%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_alpha))
  print(grad)
  return(grad)
  
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


log_likelihood_grad_decomp2=function(ytr,xtr,l1,sigma1,l2,sigma2,sigma_obs,alpha){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=decomposition_mat2(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,alpha)+
    sigma_obs*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  print(c(ll))
  return(ll)
  
  
}


log_likelihood_decomp2=function(ytr,xtr,l1,sigma1,l2,sigma2,sigma_obs,alpha){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=decomposition_mat2(xtr,xtr,l1,sigma1,l2,sigma2,alpha)+
    sigma_obs^2*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  print(c(ll))
  return(ll)
  
  
}


log_likelihood_mixture_grad=function(ytr,xtr,l1,sigma1,l2,sigma2,sigma_obs,alpha){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=(1-2*alpha)*eq_cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,fundamental = 1)+
    alpha*cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5)+
    sigma_obs*diag(nrow=2*nrow(xtr))
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}

log_likelihood_mixture_grad2=function(ytr,xtr,l1,sigma1,l2,sigma2,l1_2,sigma1_2,
                                      l2_2,sigma2_2,sigma_obs,alpha){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=(1-alpha)*eq_cov_mat(xtr,xtr,l1_2^.5,sigma1_2^.5,l2_2^.5,sigma2_2^.5,fundamental = 1)+
    alpha*(cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5)-
             eq_cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,fundamental = 1))+
    sigma_obs*diag(nrow=2*nrow(xtr))
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}


mix_cov_mat=function(x1, x2, l1, sigma1, l2, sigma2,alpha){
  K1=eq_cov_mat(x1, x2,l1,sigma1,l2,sigma2,fundamental=1)
  
  K2=cov_mat(x1, x2,l1,sigma1,l2,sigma2)
  
  (1-2*alpha)*K1+alpha*K2
}


grad_mix2=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,l3_2,sigma3_2,l4_2,sigma4_2,sigma_obs_2,alpha){
  
  K1=eq_cov_mat(xtr,xtr,l1_2^.5,sigma1_2^.5,l2_2^.5,sigma2_2^.5,fundamental=1)
  
  K1_l1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K1_l2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K1_sigma1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K1_sigma2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  K2_l3_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K2_l4_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K2_sigma3_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K2_sigma4_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  K2=eq_cov_mat(xtr,xtr,l3_2^.5,sigma3_2^.5,l4_2^.5,sigma4_2^.5,fundamental=1)
  
  K3=cov_mat(xtr,xtr,l1_2^.5,sigma1_2^.5,l2_2^.5,sigma2_2^.5)
  
  K=(1-alpha)*K2+alpha*(K3-K1)+diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  
  for(i in 1:(length(as.matrix(xtr))/2)){
    if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:(length(as.matrix(xtr))/2)){
      r1=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[i,])^2)^.5)
      r2=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[j,])^2)^.5)
      #dist=abs(xtr[i,1]-xtr[j,1])+abs(xtr[i,2]-xtr[j,2])
      theta2=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                    atan2(xtr[j,2],xtr[j,1]))
      
      
      results1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(-sigma1_2*exp(-.5*(r1-r2)^2*(l1_2))*(r1-r2)^2/2,0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,-sigma2_2*exp(-.5*(r1-r2)^2*(l2_2))*(r1-r2)^2/2))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(exp(-.5*(r1-r2)^2*(l1_2)),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,exp(-.5*(r1-r2)^2*(l2_2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(-sigma3_2*exp(-.5*(r1-r2)^2*(l3_2))*(r1-r2)^2/2,0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,-sigma4_2*exp(-.5*(r1-r2)^2*(l4_2))*(r1-r2)^2/2))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results7=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(exp(-.5*(r1-r2)^2*(l3_2)),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results8=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,exp(-.5*(r1-r2)^2*(l2_2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      K1_l1_2[i, j] = results1[1,1]
      K1_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2,2]
      K1_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[2,1]
      K1_l1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[1,2]
      
      
      
      K1_l2_2[i, j] = results2[1,1]
      K1_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2,2]
      K1_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[2,1]
      K1_l2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[1,2]
      
      
      
      K1_sigma1_2[i, j] = results3[1,1]
      K1_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2,2]
      K1_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[2,1]
      K1_sigma1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[1,2]
      
      
      
      K1_sigma2_2[i, j] = results4[1,1]
      K1_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2,2]
      K1_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[2,1]
      K1_sigma2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[1,2]
      
      
      K2_l3_2[i, j] = results5[1,1]
      K2_l3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results5[2,2]
      K2_l3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results5[2,1]
      K2_l3_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results5[1,2]
      
      
      
      K2_l4_2[i, j] = results6[1,1]
      K2_l4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
              ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results6[2,2]
      K2_l4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results6[2,1]
      K2_l4_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results6[1,2]
      
      
      
      K2_sigma3_2[i, j] = results7[1,1]
      K2_sigma3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results7[2,2]
      K2_sigma3_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results7[2,1]
      K2_sigma3_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results7[1,2]
      
      
      
      K2_sigma4_2[i, j] = results8[1,1]
      K2_sigma4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                  ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results8[2,2]
      K2_sigma4_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results8[2,1]
      K2_sigma4_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results8[1,2]
      
      
      
    }
  }
  
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  if(n1==1){
    dist_mat=dist(t(xtr),xtr)
  }else{
    dist_mat=dist(xtr,xtr)
  }
  del_cov1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov1[1:n1,1:n2]=-sigma1_2*exp(-.5*dist_mat^2*(l1_2))*dist_mat^2/2
  del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov2[1:n1,1:n2]=0
  del_cov2[(n1+1):(2*n1),(n2+1):(2*n2)]=-sigma2_2*exp(-.5*dist_mat^2*(l2_2))*dist_mat^2/2
  
  del_cov3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov3[1:n1,1:n2]=exp(-.5*dist_mat^2*(l1_2))
  del_cov3[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov4[1:n1,1:n2]=0
  del_cov4[(n1+1):(2*n1),(n2+1):(2*n2)]=exp(-.5*dist_mat^2*(l2_2))
  
  del_alpha=-K2+K3-K1
  
  #K=(1-alpha)*K2+alpha*(K3-K1)+diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  inv_K=solve(K)
  b=1-alpha
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(alpha*(del_cov1-K1_l1_2))%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*(del_cov1-K1_l1_2))),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(alpha*(del_cov3-K1_sigma1_2))%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*(del_cov3-K1_sigma1_2))),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(alpha*(del_cov2-K1_l2_2))%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*(del_cov2-K1_l2_2))),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(alpha*(del_cov4-K1_sigma2_2))%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*(del_cov4-K1_sigma2_2))),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*K2_l3_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(b*K2_l3_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*K2_sigma3_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(b*K2_sigma3_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*K2_l4_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(b*K2_l4_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*K2_sigma4_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(b*K2_sigma4_2)),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_alpha%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_alpha))
  
  return(grad)
  
}

mix_cov_mat2=function(x1, x2, l1, sigma1, l2, sigma2, l3,sigma3,l4,sigma4,alpha){
  K1=eq_cov_mat(x1,x2,l1,sigma1,l2,sigma2,fundamental=1)
  
  K2=eq_cov_mat(x1,x2,l3,sigma3,l4,sigma4,fundamental=1)
  
  K3=cov_mat(x1,x2,l1,sigma1,l2,sigma2)
  
  K=(1-alpha)*K2+alpha*(K3-K1)
  
  K
}


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
    params_new[is.na(params_new)]=params[is.na(params_new)]
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



#initial_par=c(1,1,2.7,0.368,0.135^.5,1,1,2.7,0.368,0.135^.5,0.46)
#initial_par=c(rep(1.5,9),0.55)



opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  function(x) {
    grad_mix(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  initial_par[c(1:5,10)],
  max_iter = max_iter,lr=lr,alpha=1,tol=1e-2
)^.5

Ktetr_mix=mix_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],
                      opt_par[6]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat(Xtr,Xtr,opt_par[1],
                                                 opt_par[2],
                                                 opt_par[3],
                                                 opt_par[4],
                                                 opt_par[6]^2)
                                     
                                     +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                            posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                 (length(posterior_mean_mix))])
                      -Yte)^2,1,sum))^.5)
opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad2(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  function(x) {
    grad_mix2(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  initial_par,
  max_iter = max_iter,lr=lr,alpha2=1,tol=1e-8
)^.5


Ktetr_mix=mix_cov_mat2(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4],opt_par[5], opt_par[6], opt_par[7], opt_par[8],
                       opt_par[10]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat2(Xtr,Xtr,opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],
                                                  opt_par[5],
                                                  opt_par[6],
                                                  opt_par[7], 
                                                  opt_par[8],
                                                  opt_par[10]^2)
                                     
                                     +diag(opt_par[9]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix2=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                             posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                  (length(posterior_mean_mix))])
                       -Yte)^2,1,sum))^.5)





opt_par = Adam(
  function(x) {
    log_likelihood_helm_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_helm(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

# lΦ = 1.1131, σΦ = 0.0342, lΨ = 1.5142, σΨ = 0.8884, σ2 = 0.1597
#opt_par=c(1.1131,0.0342,1.5142,0.8884,0.1597)
Ktetr_helm=cov_mat_helm(Xte,Xtr,opt_par[1],
                        opt_par[2],
                        opt_par[3],
                        opt_par[4])

Ktetr_trtr_inv_helm=Ktetr_helm%*%solve(cov_mat_helm(Xtr,Xtr,opt_par[1],
                                                    opt_par[2],
                                                    opt_par[3],
                                                    opt_par[4])
                                       
                                       +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_helm=Ktetr_trtr_inv_helm%*%c(Ytr[,1],Ytr[,2])



(rmse_helm=mean(apply((cbind(posterior_mean_helm[1:(0.5*length(posterior_mean_helm))],
                             posterior_mean_helm[(0.5*length(posterior_mean_helm)+1):
                                                   (length(posterior_mean_helm))])
                       -Yte)^2,1,sum))^.5)

opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5
# l1 = 1.6191, σ1 = 0.9710, l2 = 2.7183, σ2 = 0.5811, σ2 obs 0.1759.

Ktetr_standard=cov_mat(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])

Ktetr_trtr_inv_standard=Ktetr_standard%*%solve(cov_mat(Xtr,Xtr,opt_par[1],
                                                       opt_par[2],
                                                       opt_par[3],
                                                       opt_par[4])
                                               
                                               +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_standard=Ktetr_trtr_inv_standard%*%c(Ytr[,1],Ytr[,2])



(rmse_standard=mean(apply((cbind(posterior_mean_standard[1:(0.5*length(posterior_mean_standard))],
                                 posterior_mean_standard[(0.5*length(posterior_mean_standard)+1):
                                                           (length(posterior_mean_standard))])
                           -Yte)^2,1,sum))^.5)



opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],0)
  },
  function(x) {
    grad_fund(Xtr, Ytr,x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

Ktetr_fund=eq_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],fundamental=1)

Ktetr_trtr_inv_fund=Ktetr_fund%*%solve(eq_cov_mat(Xtr,Xtr,
                                                  opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],fundamental=1)+
                                         diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_fund=Ktetr_trtr_inv_fund%*%c(Ytr[,1],Ytr[,2])

(rmse_fund=mean(apply((cbind(posterior_mean_fund[1:(0.5*length(posterior_mean_fund))],
                             posterior_mean_fund[(0.5*length(posterior_mean_fund)+1):
                                                   (length(posterior_mean_fund))])
                       -Yte)^2,1,sum))^.5)


res1=c(initial_par[1],lr,max_iter,rmse_standard,rmse_helm,rmse_mix,rmse_mix2,rmse_fund)
print(res1)


initial_par=c(rep(1.5,9),0.55)



opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  function(x) {
    grad_mix(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  initial_par[c(1:5,10)],
  max_iter = max_iter,lr=lr,alpha=1,tol=1e-2
)^.5

Ktetr_mix=mix_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],
                      opt_par[6]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat(Xtr,Xtr,opt_par[1],
                                                 opt_par[2],
                                                 opt_par[3],
                                                 opt_par[4],
                                                 opt_par[6]^2)
                                     
                                     +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                            posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                 (length(posterior_mean_mix))])
                      -Yte)^2,1,sum))^.5)
opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad2(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  function(x) {
    grad_mix2(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  initial_par,
  max_iter = max_iter,lr=lr,alpha2=1,tol=1e-8
)^.5


Ktetr_mix=mix_cov_mat2(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4],opt_par[5], opt_par[6], opt_par[7], opt_par[8],
                       opt_par[10]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat2(Xtr,Xtr,opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],
                                                  opt_par[5],
                                                  opt_par[6],
                                                  opt_par[7], 
                                                  opt_par[8],
                                                  opt_par[10]^2)
                                     
                                     +diag(opt_par[9]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix2=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                             posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                  (length(posterior_mean_mix))])
                       -Yte)^2,1,sum))^.5)





opt_par = Adam(
  function(x) {
    log_likelihood_helm_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_helm(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

# lΦ = 1.1131, σΦ = 0.0342, lΨ = 1.5142, σΨ = 0.8884, σ2 = 0.1597
#opt_par=c(1.1131,0.0342,1.5142,0.8884,0.1597)
Ktetr_helm=cov_mat_helm(Xte,Xtr,opt_par[1],
                        opt_par[2],
                        opt_par[3],
                        opt_par[4])

Ktetr_trtr_inv_helm=Ktetr_helm%*%solve(cov_mat_helm(Xtr,Xtr,opt_par[1],
                                                    opt_par[2],
                                                    opt_par[3],
                                                    opt_par[4])
                                       
                                       +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_helm=Ktetr_trtr_inv_helm%*%c(Ytr[,1],Ytr[,2])



(rmse_helm=mean(apply((cbind(posterior_mean_helm[1:(0.5*length(posterior_mean_helm))],
                             posterior_mean_helm[(0.5*length(posterior_mean_helm)+1):
                                                   (length(posterior_mean_helm))])
                       -Yte)^2,1,sum))^.5)

opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5
# l1 = 1.6191, σ1 = 0.9710, l2 = 2.7183, σ2 = 0.5811, σ2 obs 0.1759.

Ktetr_standard=cov_mat(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])

Ktetr_trtr_inv_standard=Ktetr_standard%*%solve(cov_mat(Xtr,Xtr,opt_par[1],
                                                       opt_par[2],
                                                       opt_par[3],
                                                       opt_par[4])
                                               
                                               +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_standard=Ktetr_trtr_inv_standard%*%c(Ytr[,1],Ytr[,2])



(rmse_standard=mean(apply((cbind(posterior_mean_standard[1:(0.5*length(posterior_mean_standard))],
                                 posterior_mean_standard[(0.5*length(posterior_mean_standard)+1):
                                                           (length(posterior_mean_standard))])
                           -Yte)^2,1,sum))^.5)



opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],0)
  },
  function(x) {
    grad_fund(Xtr, Ytr,x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

Ktetr_fund=eq_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],fundamental=1)

Ktetr_trtr_inv_fund=Ktetr_fund%*%solve(eq_cov_mat(Xtr,Xtr,
                                                  opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],fundamental=1)+
                                         diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_fund=Ktetr_trtr_inv_fund%*%c(Ytr[,1],Ytr[,2])

(rmse_fund=mean(apply((cbind(posterior_mean_fund[1:(0.5*length(posterior_mean_fund))],
                             posterior_mean_fund[(0.5*length(posterior_mean_fund)+1):
                                                   (length(posterior_mean_fund))])
                       -Yte)^2,1,sum))^.5)


res15=c(initial_par[1],lr,max_iter,rmse_standard,rmse_helm,rmse_mix,rmse_mix2,rmse_fund)
print(res15)




initial_par=c(rep(2,9),0.55)



opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  function(x) {
    grad_mix(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  initial_par[c(1:5,10)],
  max_iter = max_iter,lr=lr,alpha=1,tol=1e-2
)^.5

Ktetr_mix=mix_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],
                      opt_par[6]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat(Xtr,Xtr,opt_par[1],
                                                 opt_par[2],
                                                 opt_par[3],
                                                 opt_par[4],
                                                 opt_par[6]^2)
                                     
                                     +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                            posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                 (length(posterior_mean_mix))])
                      -Yte)^2,1,sum))^.5)
opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad2(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  function(x) {
    grad_mix2(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  initial_par,
  max_iter = max_iter,lr=lr,alpha2=1,tol=1e-8
)^.5


Ktetr_mix=mix_cov_mat2(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4],opt_par[5], opt_par[6], opt_par[7], opt_par[8],
                       opt_par[10]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat2(Xtr,Xtr,opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],
                                                  opt_par[5],
                                                  opt_par[6],
                                                  opt_par[7], 
                                                  opt_par[8],
                                                  opt_par[10]^2)
                                     
                                     +diag(opt_par[9]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix2=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                             posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                  (length(posterior_mean_mix))])
                       -Yte)^2,1,sum))^.5)





opt_par = Adam(
  function(x) {
    log_likelihood_helm_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_helm(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

# lΦ = 1.1131, σΦ = 0.0342, lΨ = 1.5142, σΨ = 0.8884, σ2 = 0.1597
#opt_par=c(1.1131,0.0342,1.5142,0.8884,0.1597)
Ktetr_helm=cov_mat_helm(Xte,Xtr,opt_par[1],
                        opt_par[2],
                        opt_par[3],
                        opt_par[4])

Ktetr_trtr_inv_helm=Ktetr_helm%*%solve(cov_mat_helm(Xtr,Xtr,opt_par[1],
                                                    opt_par[2],
                                                    opt_par[3],
                                                    opt_par[4])
                                       
                                       +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_helm=Ktetr_trtr_inv_helm%*%c(Ytr[,1],Ytr[,2])



(rmse_helm=mean(apply((cbind(posterior_mean_helm[1:(0.5*length(posterior_mean_helm))],
                             posterior_mean_helm[(0.5*length(posterior_mean_helm)+1):
                                                   (length(posterior_mean_helm))])
                       -Yte)^2,1,sum))^.5)

opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5
# l1 = 1.6191, σ1 = 0.9710, l2 = 2.7183, σ2 = 0.5811, σ2 obs 0.1759.

Ktetr_standard=cov_mat(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])

Ktetr_trtr_inv_standard=Ktetr_standard%*%solve(cov_mat(Xtr,Xtr,opt_par[1],
                                                       opt_par[2],
                                                       opt_par[3],
                                                       opt_par[4])
                                               
                                               +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_standard=Ktetr_trtr_inv_standard%*%c(Ytr[,1],Ytr[,2])



(rmse_standard=mean(apply((cbind(posterior_mean_standard[1:(0.5*length(posterior_mean_standard))],
                                 posterior_mean_standard[(0.5*length(posterior_mean_standard)+1):
                                                           (length(posterior_mean_standard))])
                           -Yte)^2,1,sum))^.5)



opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],0)
  },
  function(x) {
    grad_fund(Xtr, Ytr,x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

Ktetr_fund=eq_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],fundamental=1)

Ktetr_trtr_inv_fund=Ktetr_fund%*%solve(eq_cov_mat(Xtr,Xtr,
                                                  opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],fundamental=1)+
                                         diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_fund=Ktetr_trtr_inv_fund%*%c(Ytr[,1],Ytr[,2])

(rmse_fund=mean(apply((cbind(posterior_mean_fund[1:(0.5*length(posterior_mean_fund))],
                             posterior_mean_fund[(0.5*length(posterior_mean_fund)+1):
                                                   (length(posterior_mean_fund))])
                       -Yte)^2,1,sum))^.5)


res2=c(initial_par[1],lr,max_iter,rmse_standard,rmse_helm,rmse_mix,rmse_mix2,rmse_fund)
print(res2)


initial_par=c(rep(0.5,9),0.55)



opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  function(x) {
    grad_mix(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  initial_par[c(1:5,10)],
  max_iter = max_iter,lr=lr,alpha=1,tol=1e-2
)^.5

Ktetr_mix=mix_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],
                      opt_par[6]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat(Xtr,Xtr,opt_par[1],
                                                 opt_par[2],
                                                 opt_par[3],
                                                 opt_par[4],
                                                 opt_par[6]^2)
                                     
                                     +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                            posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                 (length(posterior_mean_mix))])
                      -Yte)^2,1,sum))^.5)
opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad2(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  function(x) {
    grad_mix2(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  initial_par,
  max_iter = max_iter,lr=lr,alpha2=1,tol=1e-8
)^.5


Ktetr_mix=mix_cov_mat2(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4],opt_par[5], opt_par[6], opt_par[7], opt_par[8],
                       opt_par[10]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat2(Xtr,Xtr,opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],
                                                  opt_par[5],
                                                  opt_par[6],
                                                  opt_par[7], 
                                                  opt_par[8],
                                                  opt_par[10]^2)
                                     
                                     +diag(opt_par[9]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix2=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                             posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                  (length(posterior_mean_mix))])
                       -Yte)^2,1,sum))^.5)





opt_par = Adam(
  function(x) {
    log_likelihood_helm_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_helm(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

# lΦ = 1.1131, σΦ = 0.0342, lΨ = 1.5142, σΨ = 0.8884, σ2 = 0.1597
#opt_par=c(1.1131,0.0342,1.5142,0.8884,0.1597)
Ktetr_helm=cov_mat_helm(Xte,Xtr,opt_par[1],
                        opt_par[2],
                        opt_par[3],
                        opt_par[4])

Ktetr_trtr_inv_helm=Ktetr_helm%*%solve(cov_mat_helm(Xtr,Xtr,opt_par[1],
                                                    opt_par[2],
                                                    opt_par[3],
                                                    opt_par[4])
                                       
                                       +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_helm=Ktetr_trtr_inv_helm%*%c(Ytr[,1],Ytr[,2])



(rmse_helm=mean(apply((cbind(posterior_mean_helm[1:(0.5*length(posterior_mean_helm))],
                             posterior_mean_helm[(0.5*length(posterior_mean_helm)+1):
                                                   (length(posterior_mean_helm))])
                       -Yte)^2,1,sum))^.5)

opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5
# l1 = 1.6191, σ1 = 0.9710, l2 = 2.7183, σ2 = 0.5811, σ2 obs 0.1759.

Ktetr_standard=cov_mat(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])

Ktetr_trtr_inv_standard=Ktetr_standard%*%solve(cov_mat(Xtr,Xtr,opt_par[1],
                                                       opt_par[2],
                                                       opt_par[3],
                                                       opt_par[4])
                                               
                                               +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_standard=Ktetr_trtr_inv_standard%*%c(Ytr[,1],Ytr[,2])



(rmse_standard=mean(apply((cbind(posterior_mean_standard[1:(0.5*length(posterior_mean_standard))],
                                 posterior_mean_standard[(0.5*length(posterior_mean_standard)+1):
                                                           (length(posterior_mean_standard))])
                           -Yte)^2,1,sum))^.5)



opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],0)
  },
  function(x) {
    grad_fund(Xtr, Ytr,x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

Ktetr_fund=eq_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],fundamental=1)

Ktetr_trtr_inv_fund=Ktetr_fund%*%solve(eq_cov_mat(Xtr,Xtr,
                                                  opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],fundamental=1)+
                                         diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_fund=Ktetr_trtr_inv_fund%*%c(Ytr[,1],Ytr[,2])

(rmse_fund=mean(apply((cbind(posterior_mean_fund[1:(0.5*length(posterior_mean_fund))],
                             posterior_mean_fund[(0.5*length(posterior_mean_fund)+1):
                                                   (length(posterior_mean_fund))])
                       -Yte)^2,1,sum))^.5)


res05=c(initial_par[1],lr,max_iter,rmse_standard,rmse_helm,rmse_mix,rmse_mix2,rmse_fund)
print(res05)

initial_par=c(rep(2.5,9),0.55)



opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  function(x) {
    grad_mix(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  initial_par[c(1:5,10)],
  max_iter = max_iter,lr=lr,alpha=1,tol=1e-2
)^.5

Ktetr_mix=mix_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],
                      opt_par[6]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat(Xtr,Xtr,opt_par[1],
                                                 opt_par[2],
                                                 opt_par[3],
                                                 opt_par[4],
                                                 opt_par[6]^2)
                                     
                                     +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                            posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                 (length(posterior_mean_mix))])
                      -Yte)^2,1,sum))^.5)
opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad2(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  function(x) {
    grad_mix2(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  initial_par,
  max_iter = max_iter,lr=lr,alpha2=1,tol=1e-8
)^.5


Ktetr_mix=mix_cov_mat2(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4],opt_par[5], opt_par[6], opt_par[7], opt_par[8],
                       opt_par[10]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat2(Xtr,Xtr,opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],
                                                  opt_par[5],
                                                  opt_par[6],
                                                  opt_par[7], 
                                                  opt_par[8],
                                                  opt_par[10]^2)
                                     
                                     +diag(opt_par[9]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix2=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                             posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                  (length(posterior_mean_mix))])
                       -Yte)^2,1,sum))^.5)





opt_par = Adam(
  function(x) {
    log_likelihood_helm_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_helm(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

# lΦ = 1.1131, σΦ = 0.0342, lΨ = 1.5142, σΨ = 0.8884, σ2 = 0.1597
#opt_par=c(1.1131,0.0342,1.5142,0.8884,0.1597)
Ktetr_helm=cov_mat_helm(Xte,Xtr,opt_par[1],
                        opt_par[2],
                        opt_par[3],
                        opt_par[4])

Ktetr_trtr_inv_helm=Ktetr_helm%*%solve(cov_mat_helm(Xtr,Xtr,opt_par[1],
                                                    opt_par[2],
                                                    opt_par[3],
                                                    opt_par[4])
                                       
                                       +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_helm=Ktetr_trtr_inv_helm%*%c(Ytr[,1],Ytr[,2])



(rmse_helm=mean(apply((cbind(posterior_mean_helm[1:(0.5*length(posterior_mean_helm))],
                             posterior_mean_helm[(0.5*length(posterior_mean_helm)+1):
                                                   (length(posterior_mean_helm))])
                       -Yte)^2,1,sum))^.5)

opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5
# l1 = 1.6191, σ1 = 0.9710, l2 = 2.7183, σ2 = 0.5811, σ2 obs 0.1759.

Ktetr_standard=cov_mat(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])

Ktetr_trtr_inv_standard=Ktetr_standard%*%solve(cov_mat(Xtr,Xtr,opt_par[1],
                                                       opt_par[2],
                                                       opt_par[3],
                                                       opt_par[4])
                                               
                                               +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_standard=Ktetr_trtr_inv_standard%*%c(Ytr[,1],Ytr[,2])



(rmse_standard=mean(apply((cbind(posterior_mean_standard[1:(0.5*length(posterior_mean_standard))],
                                 posterior_mean_standard[(0.5*length(posterior_mean_standard)+1):
                                                           (length(posterior_mean_standard))])
                           -Yte)^2,1,sum))^.5)



opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],0)
  },
  function(x) {
    grad_fund(Xtr, Ytr,x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

Ktetr_fund=eq_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],fundamental=1)

Ktetr_trtr_inv_fund=Ktetr_fund%*%solve(eq_cov_mat(Xtr,Xtr,
                                                  opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],fundamental=1)+
                                         diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_fund=Ktetr_trtr_inv_fund%*%c(Ytr[,1],Ytr[,2])

(rmse_fund=mean(apply((cbind(posterior_mean_fund[1:(0.5*length(posterior_mean_fund))],
                             posterior_mean_fund[(0.5*length(posterior_mean_fund)+1):
                                                   (length(posterior_mean_fund))])
                       -Yte)^2,1,sum))^.5)


res25=c(initial_par[1],lr,max_iter,rmse_standard,rmse_helm,rmse_mix,rmse_mix2,rmse_fund)
print(res25)

initial_par=c(rep(4,9),0.55)



opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  function(x) {
    grad_mix(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  initial_par[c(1:5,10)],
  max_iter = max_iter,lr=lr,alpha=1,tol=1e-2
)^.5

Ktetr_mix=mix_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],
                      opt_par[6]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat(Xtr,Xtr,opt_par[1],
                                                 opt_par[2],
                                                 opt_par[3],
                                                 opt_par[4],
                                                 opt_par[6]^2)
                                     
                                     +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                            posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                 (length(posterior_mean_mix))])
                      -Yte)^2,1,sum))^.5)
opt_par = Adam(
  function(x) {
    log_likelihood_mixture_grad2(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  function(x) {
    grad_mix2(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
  },
  initial_par,
  max_iter = max_iter,lr=lr,alpha2=1,tol=1e-8
)^.5


Ktetr_mix=mix_cov_mat2(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4],opt_par[5], opt_par[6], opt_par[7], opt_par[8],
                       opt_par[10]^2)

Ktetr_trtr_inv_mix=Ktetr_mix%*%solve(mix_cov_mat2(Xtr,Xtr,opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],
                                                  opt_par[5],
                                                  opt_par[6],
                                                  opt_par[7], 
                                                  opt_par[8],
                                                  opt_par[10]^2)
                                     
                                     +diag(opt_par[9]^2,nrow=2*nrow(Xtr)))

posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])



(rmse_mix2=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                             posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                  (length(posterior_mean_mix))])
                       -Yte)^2,1,sum))^.5)





opt_par = Adam(
  function(x) {
    log_likelihood_helm_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_helm(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

# lΦ = 1.1131, σΦ = 0.0342, lΨ = 1.5142, σΨ = 0.8884, σ2 = 0.1597
#opt_par=c(1.1131,0.0342,1.5142,0.8884,0.1597)
Ktetr_helm=cov_mat_helm(Xte,Xtr,opt_par[1],
                        opt_par[2],
                        opt_par[3],
                        opt_par[4])

Ktetr_trtr_inv_helm=Ktetr_helm%*%solve(cov_mat_helm(Xtr,Xtr,opt_par[1],
                                                    opt_par[2],
                                                    opt_par[3],
                                                    opt_par[4])
                                       
                                       +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_helm=Ktetr_trtr_inv_helm%*%c(Ytr[,1],Ytr[,2])



(rmse_helm=mean(apply((cbind(posterior_mean_helm[1:(0.5*length(posterior_mean_helm))],
                             posterior_mean_helm[(0.5*length(posterior_mean_helm)+1):
                                                   (length(posterior_mean_helm))])
                       -Yte)^2,1,sum))^.5)

opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
  },
  function(x) {
    grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5
# l1 = 1.6191, σ1 = 0.9710, l2 = 2.7183, σ2 = 0.5811, σ2 obs 0.1759.

Ktetr_standard=cov_mat(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])

Ktetr_trtr_inv_standard=Ktetr_standard%*%solve(cov_mat(Xtr,Xtr,opt_par[1],
                                                       opt_par[2],
                                                       opt_par[3],
                                                       opt_par[4])
                                               
                                               +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_standard=Ktetr_trtr_inv_standard%*%c(Ytr[,1],Ytr[,2])



(rmse_standard=mean(apply((cbind(posterior_mean_standard[1:(0.5*length(posterior_mean_standard))],
                                 posterior_mean_standard[(0.5*length(posterior_mean_standard)+1):
                                                           (length(posterior_mean_standard))])
                           -Yte)^2,1,sum))^.5)



opt_par = Adam(
  function(x) {
    log_likelihood_grad(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],0)
  },
  function(x) {
    grad_fund(Xtr, Ytr,x[1], x[2], x[3], x[4], x[5])
  },
  initial_par[1:5],
  max_iter = max_iter,lr=lr
)^.5

Ktetr_fund=eq_cov_mat(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],fundamental=1)

Ktetr_trtr_inv_fund=Ktetr_fund%*%solve(eq_cov_mat(Xtr,Xtr,
                                                  opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],fundamental=1)+
                                         diag(opt_par[5]^2,nrow=2*nrow(Xtr)))

posterior_mean_fund=Ktetr_trtr_inv_fund%*%c(Ytr[,1],Ytr[,2])

(rmse_fund=mean(apply((cbind(posterior_mean_fund[1:(0.5*length(posterior_mean_fund))],
                             posterior_mean_fund[(0.5*length(posterior_mean_fund)+1):
                                                   (length(posterior_mean_fund))])
                       -Yte)^2,1,sum))^.5)


res4=c(initial_par[1],lr,max_iter,rmse_standard,rmse_helm,rmse_mix,rmse_mix2,rmse_fund)
print(res4)

res=rbind(res05,res1,res15,res2,res25,res4)
print(res)
save(list = "res", file = paste0("rmses_inital_REPARAM_2", id, ".rda"))
