id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

library(progress)
library(doParallel)
library(pracma)
library(proxy)
library(nloptr)
library(RColorBrewer)
library(rgl)
notcluster=is.na(id)
lr=.01;max_iter=1000
QUICK_RUN=0
jit=1e-6

right_bound=5
discontinuities=1000
discontinuity_rate=1/discontinuities*right_bound

A=function(r){
  for(i in 1:discontinuities){
    if(((i-1)*discontinuity_rate<=r)&(r<(i)*discontinuity_rate)){
      if(i%%2==0){
        return(1)
      }else{
        return(-1)
      }
    }
  }
  if(r>=right_bound){
    return(ifelse(discontinuities%%2==0,1,-1))
  }
}


cov_mat=function(x1,x2,sigma1,sigma2,l1,l2){
  norm=function(x){sum(x^2)^(.5)}
  
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
#Helmholtz kernel

cov_mat_helm=function(x1,x2,l1,sigma1,l2,sigma2){
  norm=function(x){sum(x^2)^(.5)}
  
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

# Fundamental domain equivariant kernel

fund_cov_mat= function(x1, x2, sigma1,sigma2,l1,l2) {
  norm=function(x){sum(x^2)^(.5)}
  
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
      #print(c(A(r1),A(r2)))
      theta1=theta1+pi*ifelse(A(r1)==-1,1,0)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))+pi*ifelse(A(r2)==-1,1,0)
      cov_standard=cov_mat(c(A(r1),0),c(A(r2),0),
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
# Double integration equivariant kernel 


log_likelihood_fund=function(ytr,xtr,sigma1,sigma2,tau,l1,l2){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  Ktr=fund_cov_mat(xtr,xtr,sigma1,sigma2,l1,l2)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}

norm=function(x){sum(x^2)^.5}
del_l_1=function(x1,x2,sigma1,sigma2,l1,l2){
  
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  if(n1==1){
    if(n2==1){
      dist_mat=norm(x1-x2)
    }else{dist_mat=distmat(t(x1),x2)
    }
  }else{
    dist_mat=distmat(x1,x2) #needs pracma library
  }
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=0
  
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  return(del_l1)
}

del_l_2=function(x1,x2,sigma1,sigma2,l1,l2){
  
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  if(n1==1){
    if(n2==1){
      dist_mat=norm(x1-x2)
    }else{dist_mat=distmat(t(x1),x2)
    }
  }else{
    dist_mat=distmat(x1,x2) #needs pracma library
  }
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=0
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=0
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  return(del_l2)
}

del_sig1=function(x1,x2,sigma1,sigma2,l1,l2){
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  if(n1==1){
    if(n2==1){
      dist_mat=norm(x1-x2)
    }else{dist_mat=distmat(t(x1),x2)
    }
  }else{
    dist_mat=distmat(x1,x2) #needs pracma library
  }
  
  del_sigma1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),1:n2]=0
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  return(del_sigma1)
}

del_sig2=function(x1,x2,sigma1,sigma2,l1,l2){
  
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  if(n1==1){
    if(n2==1){
      dist_mat=norm(x1-x2)
    }else{dist_mat=distmat(t(x1),x2)
    }
  }else{
    dist_mat=distmat(x1,x2) #needs pracma library
  }
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=0
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=0
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  return(del_sigma2)
}


grad_fund=function(xtr,ytr,sigma1,sigma2,tau,l1,l2){
  x1=x2=xtr
  K=fund_cov_mat(xtr,xtr,sigma1,sigma2,l1,l2)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  cov = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1,nrow(x2)),
               nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  
  del_sigma1=del_sigma2=del_l1=del_l2=cov
  
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
      theta1=theta1+pi*ifelse(A(r1)==-1,1,0)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))+pi*ifelse(A(r2)==-1,1,0)
      
      cov_standard=cov_mat(c(A(r1),0),c(A(r2),0),
                            sigma1,sigma2,l1,l2)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard %*%cbind(c(cos(theta2),-sin(theta2)),
                              c(sin(theta2),cos(theta2)))
      
      cov[i, j] = results[1,1]
      cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      
      results_sigma1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(A(r1),0),c(A(r2),0),sigma1,sigma2,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(A(r1),0),c(A(r2),0),sigma1,sigma2,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      results_l1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(A(r1),0),c(A(r2),0),sigma1,sigma2,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(A(r1),0),c(A(r2),0),sigma1,sigma2,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      del_sigma1[i, j] = results_sigma1[1,1]
      del_sigma1[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma1[2,2]
      del_sigma1[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma1[2,1]
      del_sigma1[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma1[1,2]
      
      
      
      del_sigma2[i, j] = results_sigma2[1,1]
      del_sigma2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma2[2,2]
      del_sigma2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma2[2,1]
      del_sigma2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma2[1,2]
      
      
      
      del_l1[i, j] = results_l1[1,1]
      del_l1[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l1[2,2]
      del_l1[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l1[2,1]
      del_l1[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l1[1,2]
      
      del_l2[i, j] = results_l2[1,1]
      del_l2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l2[2,2]
      del_l2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l2[2,1]
      del_l2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l2[1,2]
      
      
      
    }
  }
  inv_K=inv(K)
  alpha=1
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma2)),
            -t(c(ytr[,1],ytr[,2]))%*%(2*tau*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*tau*inv_K),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l2))
  )
  #print(grad)
  return(grad)
  
  
} 

grad_standard=function(xtr,ytr,sigma1,sigma2,tau,l1,l2){
  x1=x2=xtr
  K=cov_mat(xtr,xtr,sigma1,sigma2,l1,l2)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  if(n1==1){
    if(n2==1){
      dist_mat=norm(x1-x2)
    }else{dist_mat=distmat(t(x1),x2)
    }
  }else{
    dist_mat=distmat(x1,x2) #needs pracma library
  }
  
  
  
  del_sigma1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=0
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=0
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=0
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  
  
  
  del_tau=2*tau*diag(nrow=2*nrow(xtr))
  
  
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=0
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=0
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=0
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  
  inv_K=inv(K)
  
  #sigma1,sigma2,sigma3,sigma4,nu1=1,nu2=1,tau,a1,a2
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma1),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma2),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_tau%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_tau),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_l1),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_l2)
  )
  
  return(grad)
  
}

grad_helm=function(x1,ytr,l1,sigma1,l2,sigma2,sigma_obs_2){
  x1=x2=as.matrix(x1)
  dist_mat=distmat(x1,x2)
  diff_mat=outer(x1[, 1], x2[, 1], "-") * outer(x1[, 2], x2[, 2], "-")
  dist_mat_x1=distmat(cbind(x1[,1],0),cbind(x2[,1],0))
  dist_mat_x2=distmat(cbind(x1[,2],0),cbind(x2[,2],0))
  
  n1=nrow(x1)
  n2=nrow(x2)
  
  del_cov1=del_cov2=del_cov3=del_cov4=cov=matrix(0,nrow=2*n1,ncol=2*n2)
  
  
  del_cov1[1:n1,1:n2]=-4*sigma1^2/((l1)^5)*exp(-.5*dist_mat^2/(l1^2))*(l1^2-dist_mat_x1^2)+
    sigma1^2/((l1)^4)*exp(-.5*dist_mat^2/(l1^2))*((l1^2-dist_mat_x1^2)*dist_mat^2/(l1^3)+2*l1)
  
  del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=-4*sigma1^2/((l1)^5)*exp(-.5*dist_mat^2/(l1^2))*(l1^2-dist_mat_x2^2)+
    sigma1^2/((l1)^4)*exp(-.5*dist_mat^2/(l1^2))*((l1^2-dist_mat_x2^2)*dist_mat^2/(l1^3)+2*l1)
  
  del_cov1[1:n1,(n2+1):(2*n2)]=del_cov1[(n1+1):(2*n1),1:n2]=
    -diff_mat*exp(-.5*dist_mat^2/(l1^2))*(dist_mat^2/(l1^3)*sigma1^2/((l1)^4)-4*sigma1^2/((l1)^5))
  
  del_cov2[1:n1,1:n2]=-4*sigma2^2/((l2)^5)*exp(-.5*dist_mat^2/(l2^2))*(l2^2-dist_mat_x2^2)+
    sigma2^2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2))*((l2^2-dist_mat_x2^2)*dist_mat^2/(l2^3)+2*l2)
  
  del_cov2[(n1+1):(2*n1),(n2+1):(2*n2)]=-4*sigma2^2/((l2)^5)*exp(-.5*dist_mat^2/(l2^2))*(l2^2-dist_mat_x1^2)+
    sigma2^2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2))*((l2^2-dist_mat_x1^2)*dist_mat^2/(l2^3)+2*l2)
  
  del_cov2[1:n1,(n2+1):(2*n2)]=del_cov2[(n1+1):(2*n1),1:n2]=
    diff_mat*exp(-.5*dist_mat^2/(l2^2))*(dist_mat^2/(l2^3)*sigma2^2/((l2)^4)-4*sigma2^2/((l2)^5))
  
  
  del_cov3[1:n1,1:n2]=2*sigma1/((l1)^4)*exp(-.5*dist_mat^2/(l1^2))*(l1^2-dist_mat_x1^2)
  del_cov3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma1/((l1)^4)*exp(-.5*dist_mat^2/(l1^2))*
    (l1^2-dist_mat_x2^2)
  del_cov3[1:n1,(n2+1):(2*n2)]=del_cov3[(n1+1):(2*n1),1:n2]=
    diff_mat*(-exp(-.5*dist_mat^2/(l1^2))*2*sigma1/((l1)^4))
  
  del_cov4[1:n1,1:n2]=2*sigma2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2))*(l2^2-dist_mat_x2^2)
  del_cov4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma2/((l2)^4)*exp(-.5*dist_mat^2/
                                                                (l2^2))*(l2^2-dist_mat_x1^2)
  del_cov4[1:n1,(n2+1):(2*n2)]=del_cov4[(n1+1):(2*n1),1:n2]=
    diff_mat*(2*sigma2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2)) )
  
  cov[1:n1,1:n2]=sigma1^2/((l1)^4)*exp(-.5*dist_mat^2/(l1^2))*(l1^2-dist_mat_x1^2)+
    sigma2^2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2))*(l2^2-dist_mat_x2^2)
  
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma1^2/((l1)^4)*exp(-.5*dist_mat^2/(l1^2))*
    (l1^2-dist_mat_x2^2)+
    sigma2^2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2))*(l2^2-dist_mat_x1^2)
  
  cov[(n1+1):(2*n1),1:n2]=diff_mat*(-exp(-.5*dist_mat^2/(l1^2))*sigma1^2/((l1)^4)+
                                      sigma2^2/((l2)^4)*exp(-.5*dist_mat^2/(l2^2)) )
  
  
  cov[1:n1,(n2+1):(2*n2)]=cov[(n1+1):(2*n1),1:n2]
  
  K=cov+diag(sigma_obs_2^2,nrow=2*nrow(x1))
  inv_K=inv(K)
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_cov1),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov3%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_cov3),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_cov2),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov4%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_cov4),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%diag(2*sigma_obs_2,nrow=2*nrow(x1))%*%inv_K%*%c(ytr[,1],ytr[,2])+
              Trace(inv_K%*%diag(2*sigma_obs_2,nrow=2*nrow(x1))))
  
  
  return(grad)
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
Adam = function(f, grad_f, params_init, lr = 0.01, beta1 = 0.9, beta2 = 0.999,
                epsilon = 1e-8, max_iter = 100, tol = 1e-10,
                lb = rep(0, length(params_init)),
                ub = rep(Inf, length(params_init))) {
  
  # Initialize progress bar if not in cluster
  if (notcluster) {
    pb = progress_bar$new(
      format = "  Progress [:bar] :percent in :elapsed",
      total = max_iter,
      clear = FALSE,
      width = 60
    )
  }
  
  params = params_init
  m = numeric(length(params))
  v = numeric(length(params))
  iter = 0
  
  while (iter < max_iter) {
    iter = iter + 1
    
    # Compute gradient
    grad = grad_f(params)
    
    # Update moments
    m = beta1 * m + (1 - beta1) * grad
    v = beta2 * v + (1 - beta2) * (grad^2)
    
    # Compute bias-corrected moments
    m_hat = m / (1 - beta1^iter)
    v_hat = v / (1 - beta2^iter)
    
    # Update parameters
    params_new = params - lr * m_hat / (sqrt(v_hat) + epsilon)
    
    # Handle boundary constraints
    params_new[is.na(params_new)] = params[is.na(params_new)]
    params_new[params_new <= lb] = params[params_new <= lb]
    params_new[params_new > ub] = params[params_new > ub]
    
    # Check convergence
    if (max(abs(params_new - params)) < tol) {
      break
    }
    
    params = params_new
    
    
    if (notcluster) {
      pb$tick()
    }
  }
  
  return(params)
}



initial_par=rep(1,4)
initial_sigma_obs=0.1


Xtr=matrix(runif(16,-1,1),ncol=2)
sigma_obs=0.15
set.seed(12)
noise=rnorm(16,0,sigma_obs)
Ytr=cbind(-Xtr[,2],Xtr[,1])+cbind(noise[1:8],noise[9:16])



ngrid=17


Xte=matrix(runif(2*17^2,-1,1),ncol=2)
Yte=cbind(-Xte[,2],Xte[,1])

Xte1=Xte;Yte1=Yte;Xtr1=Xtr;Ytr1=Ytr

initial_sigma_obs=.1
initial_sigma1=initial_sigma2=initial_l1=initial_l2=1



opt_par=Adam(function(x) {
  log_likelihood_fund(
    Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
},function(x) {
  grad_fund(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
},c(initial_sigma1,
    initial_sigma2,  initial_sigma_obs,
    initial_l1,
    initial_l2),lr=lr,max_iter = max_iter)

Ktetr_fund=fund_cov_mat(Xte,Xtr,opt_par[1],
                        opt_par[2],opt_par[4],
                        opt_par[5])

Ktetr_trtr_inv_fund=Ktetr_fund%*%inv(fund_cov_mat(Xtr,Xtr,opt_par[1],
                                                  opt_par[2],opt_par[4],
                                                  opt_par[5])+diag(opt_par[3]^2,nrow=2*nrow(Xtr)))

posterior_mean_fund1=Ktetr_trtr_inv_fund%*%c(Ytr[,1],Ytr[,2])



(rmse_fund1=mean(apply((cbind(
  posterior_mean_fund1[1:(0.5*length(posterior_mean_fund1))],
  posterior_mean_fund1
  [(0.5*length(posterior_mean_fund1)+1):(length(posterior_mean_fund1))])
  -Yte)^2,1,sum))^.5)



posterior_cov_fund1= fund_cov_mat(Xte,Xte,opt_par[1],
                                  opt_par[2], opt_par[4],
                                  opt_par[5])-Ktetr_trtr_inv_fund%*%t(Ktetr_fund)


(LogS_fund1=(t(c(Yte[,1],Yte[,2])-posterior_mean_fund1)%*%
               inv(posterior_cov_fund1+diag(opt_par[3]^2,2*nrow(Xte)))%*%
               (c(Yte[,1],Yte[,2])-posterior_mean_fund1)
             +determinant(posterior_cov_fund1+diag(opt_par[3]^2,2*nrow(Xte))
             )$modulus[1]+
               log(2*pi)*nrow(Xte))/nrow(Xte))




bd=0.5

Xte=cbind(runif(20^2,-2,2),runif(20^2,-2,2))

Yte=cbind(Xte[,1]/(bd+(Xte[,1]^2+Xte[,2]^2)^2),Xte[,2]/(bd+(Xte[,1]^2+Xte[,2]^2)^2))


#Xtr=cbind(runif(10,-2,2),runif(10,-2,2))
Xtr=cbind(runif(10,-2,2),runif(10,-2,2))

Ytr=cbind(Xtr[,1]/(bd+(Xtr[,1]^2+Xtr[,2]^2)^2),Xtr[,2]/(bd+(Xtr[,1]^2+Xtr[,2]^2)^2))

sigma_obs=.1
noise=rnorm(length(Xtr),0,sigma_obs)
Ytr=Ytr+cbind(noise[1:(length(Xtr)/2)],noise[(length(Xtr)/2+1):length(Xtr)])




opt_par=Adam(function(x) {
  log_likelihood_fund(
    Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
},function(x) {
  grad_fund(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
},c(initial_sigma1,
    initial_sigma2,
    
    initial_sigma_obs,
    initial_l1,
    initial_l2),lr=lr,max_iter = max_iter)

Ktetr_fund=fund_cov_mat(Xte,Xtr,opt_par[1],
                        opt_par[2],
                        
                        opt_par[4],
                        opt_par[5])

Ktetr_trtr_inv_fund=Ktetr_fund%*%inv(fund_cov_mat(Xtr,Xtr,opt_par[1],
                                                  opt_par[2],opt_par[4],
                                                  opt_par[5]) +diag(opt_par[3]^2,nrow=2*nrow(Xtr)))

posterior_mean_fund=Ktetr_trtr_inv_fund%*%c(Ytr[,1],Ytr[,2])



(rmse_fund=mean(apply((cbind(
  posterior_mean_fund[1:(0.5*length(posterior_mean_fund))],
  posterior_mean_fund
  [(0.5*length(posterior_mean_fund)+1):(length(posterior_mean_fund))])
  -Yte)^2,1,sum))^.5)



posterior_cov_fund= fund_cov_mat(Xte,Xte,opt_par[1],
                                 opt_par[2], opt_par[4],
                                 opt_par[5])-Ktetr_trtr_inv_fund%*%t(Ktetr_fund)


(LogS_fund=(t(c(Yte[,1],Yte[,2])-posterior_mean_fund)%*%
              inv(posterior_cov_fund+diag(opt_par[3]^2,2*nrow(Xte)))%*%
              (c(Yte[,1],Yte[,2])-posterior_mean_fund)
            +determinant(posterior_cov_fund+diag(opt_par[3]^2,2*nrow(Xte))
            )$modulus[1]+
              log(2*pi)*nrow(Xte))/nrow(Xte))

RMSES=c(rmse_fund1,rmse_fund)
LogS=c(LogS_fund1,LogS_fund)

res=list(RMSES,LogS)
save(list = "res", file = paste0("UQ_first_example_missspec1000_", id, ".rda"))
