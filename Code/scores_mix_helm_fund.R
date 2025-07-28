id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
test=0
library(pracma)
library(progress)
reduced_sample=FALSE
notcluster=ifelse(is.na(id),1,0)
scale_data=1
#data=1 #6hours
data="Gulf" 

lr=.01;max_iter=1000
bd=.5

ALPHA=seq(0,1,by=.2)

train_ind=sample(564,10)
#test_ind=1:4
train_size=10
N_train=10

norm=function(x){sum(x^2)^.5}

gulf_data_train = read.csv(
  "https://raw.githubusercontent.com/JaxGaussianProcesses/static/main/data/gulfdata_train.csv"
)

xtr=cbind(gulf_data_train$lon,gulf_data_train$lat)
ytr=cbind(gulf_data_train$ubar,gulf_data_train$vbar)

gulf_data_test = read.csv(
  "https://raw.githubusercontent.com/JaxGaussianProcesses/static/main/data/gulfdata_test.csv"
)

xte=cbind(gulf_data_test$lon,gulf_data_test$lat)
yte=cbind(gulf_data_test$ubar,gulf_data_test$vbar)

X=rbind(xtr,xte)
Y=rbind(ytr,yte)

Xtr=xtr
Ytr=ytr


Xte=xte
Yte=yte


if(scale_data){
  X=scale(rbind(xtr,xte))
  Y=scale(rbind(ytr,yte))
}

center=apply(X,2,function(x){mean(range(x))})
F=function(x){(x-center)/(bd+(norm(x-center)^2))}
Eq_noise=scale(t(apply(X,1,F)))


initial_alpha=.5
initial_sigma_obs=.1
initial_par=c(rep(1,4),initial_sigma_obs)
initial_sigma1=initial_sigma2=initial_sigma3=initial_sigma4=initial_sigma5=
  initial_sigma6=initial_sigma7=initial_sigma8=initial_sigma9=initial_sigma10=
  initial_sigma11=initial_sigma12=initial_l1=initial_l2=initial_l3=initial_l4=
  initial_l5=initial_l6=initial_range1=initial_l_range1=initial_range2=initial_l_range2 =1
#Xtr=Xtr[1:3,];Ytr=Ytr[1:3,];Xte=Xte[1:3,];Yte=Yte[1:3,];lb_1=-4;lb_2=-4;ub_1=4;ub_2=4;m1=0;m2=0
noise=.01



#Xtr=matrix(runif(200,-4,4),ncol=2);N_train=1


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


fund_cov_mat= function(x1, x2, sigma1,sigma2,l1,l2) {
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



fund_cov_mat_helm= function(x1, x2, sigma1,sigma2,l1,l2) {
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
      
      cov1=matrix(0,2,2)
      
      cov1[1,1]=sigma1^2/((l1)^4)*exp(-.5*(r1-r2)^2/(l1^2))*(l1^2-(r1-r2)^2)+
        sigma2^2/((l2)^2)*exp(-.5*(r1-r2)^2/(l2^2))
      
      
      cov1[2,2]=sigma1^2/((l1)^2)*exp(-.5*(r1-r2)^2/(l1^2))+ 
        sigma2^2/((l2)^4)*exp(-.5*(r1-r2)^2/(l2^2))*(l2^2-(r1-r2)^2)
      
      cov1[1,2]=cov1[2,1]=0
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov1%*%cbind(c(cos(theta2),-sin(theta2)),
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




mix_cov_mat3=function(x1, x2, sigma1,sigma2,l1,l2,alpha){
  x1=as.matrix(x1);x2=as.matrix(x2)
  K1=fund_cov_mat(x1, x2, sigma1,sigma2,l1,l2)
  
  K2=cov_mat_helm(x1, x2,l1,sigma1,l2,sigma2)
  
  (1-alpha)^2*K1+alpha^2*K2
}



log_likelihood_mix3=function(ytr,xtr,sigma1,sigma2,tau,l1,l2,alpha){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=mix_cov_mat3(xtr,xtr,sigma1,sigma2,l1,l2,alpha)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
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
      
      
      results_sigma1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma1,sigma2,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma1,sigma2,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      results_l1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma1,sigma2,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma1,sigma2,l1,l2)%*%
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
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%diag(2*sigma_obs_2,nrow=2*nrow(x1))%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%diag(2*sigma_obs_2,nrow=2*nrow(x1))))
  
  
  return(grad)
}

grad_fund_helm=function(xtr,ytr,sigma1,sigma2,tau,l1,l2){
  x1=x2=xtr
  K=fund_cov_mat_helm(xtr,xtr,sigma1,sigma2,l1,l2)+
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
      cov1=results_sigma1=results_sigma2=results_l1=results_l2=matrix(0,2,2)
      
      cov1[1,1]=sigma1^2/((l1)^4)*exp(-.5*(r1-r2)^2/(l1^2))*(l1^2-(r1-r2)^2)+
        sigma2^2/((l2)^2)*exp(-.5*(r1-r2)^2/(l2^2))
      
      
      cov1[2,2]=sigma1^2/((l1)^2)*exp(-.5*(r1-r2)^2/(l1^2))
      +sigma2^2/((l2)^4)*exp(-.5*(r1-r2)^2/(l2^2))*(l2^2-(r1-r2)^2)
      
      cov1[1,2]=cov1[2,1]=0
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov1%*%cbind(c(cos(theta2),-sin(theta2)),
                     c(sin(theta2),cos(theta2)))
      
      
      cov[i, j] = results[1,1]
      cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      results_sigma1[1,1]=2*sigma1/((l1)^4)*exp(-.5*(r1-r2)^2/(l1^2))*(l1^2-(r1-r2)^2)
      results_sigma1[2,2]=2*sigma1/((l1)^2)*exp(-.5*(r1-r2)^2/(l1^2))
      
      results_sigma1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        results_sigma1%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma2[1,1]=2*sigma2/((l2)^2)*exp(-.5*(r1-r2)^2/(l2^2))
      results_sigma2[2,2]=2*sigma2/((l2)^4)*exp(-.5*(r1-r2)^2/(l2^2))*(l2^2-(r1-r2)^2)
      
      
      results_sigma2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        results_sigma2%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l1[1,1]=-4*sigma1^2/((l1)^5)*exp(-.5*(r1-r2)^2/(l1^2))*(l1^2-(r1-r2)^2)+
        sigma1^2/((l1)^4)*exp(-.5*(r1-r2)^2/(l1^2))*((r1-r2)^2/(l1^3)*(l1^2-(r1-r2)^2)+2*l1)
      results_l1[2,2]=sigma1^2/((l1)^2)*exp(-.5*(r1-r2)^2/(l1^2))*(-2/l1+(r1-r2)^2/(l1^3))
      
      results_l1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        results_l1%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l2[1,1]=sigma2^2/((l2)^2)*exp(-.5*(r1-r2)^2/(l2^2))*(-2/l2+(r1-r2)^2/(l2^3))
      
      results_l2[2,2]=sigma2^2/((l2)^4)*exp(-.5*(r1-r2)^2/(l2^2))*
        (-4/l2*(l2^2-(r1-r2)^2)+ (r1-r2)^2/(l2^3)*(l2^2-(r1-r2)^2)+2*l2)
      
      
      results_l2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        results_l2%*%
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

grad_mix3=function(xtr,ytr,sigma1,sigma2,tau,l1,l2,alpha){
  x1=x2=xtr
  b=(1-alpha)^2
  
  K1=fund_cov_mat(x1, x2, sigma1,sigma2,l1,l2)
  
  K2=cov_mat_helm(x1, x2,l1,sigma1,l2,sigma2)
  
  K=(1-alpha)^2*K1+alpha^2*K2+
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
      cov1=results_sigma1=results_sigma2=results_l1=results_l2=matrix(0,2,2)
      
      cov1[1,1]=sigma1^2/((l1)^4)*exp(-.5*(r1-r2)^2/(l1^2))*(l1^2-(r1-r2)^2)+
        sigma2^2/((l2)^2)*exp(-.5*(r1-r2)^2/(l2^2))
      
      
      cov1[2,2]=sigma1^2/((l1)^2)*exp(-.5*(r1-r2)^2/(l1^2))
      +sigma2^2/((l2)^4)*exp(-.5*(r1-r2)^2/(l2^2))*(l2^2-(r1-r2)^2)
      
      cov1[1,2]=cov1[2,1]=0
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov1%*%cbind(c(cos(theta2),-sin(theta2)),
                     c(sin(theta2),cos(theta2)))
      
      
      cov[i, j] = results[1,1]
      cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      
      
      results_sigma1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma1,sigma2,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma1,sigma2,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      results_l1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma1,sigma2,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma1,sigma2,l1,l2)%*%
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
  x1=x2=as.matrix(xtr)
  dist_mat=distmat(x1,x2)
  diff_mat=outer(x1[, 1], x2[, 1], "-") * outer(x1[, 2], x2[, 2], "-")
  dist_mat_x1=distmat(cbind(x1[,1],0),cbind(x2[,1],0))
  dist_mat_x2=distmat(cbind(x1[,2],0),cbind(x2[,2],0))
  
  n1=nrow(x1)
  n2=nrow(x2)
  
  del_cov1=del_cov2=del_cov3=del_cov4=cov=matrix(0,nrow=2*n1,ncol=2*n2)
  
  
  
  del_cov1[1:n1,1:n2]=exp(-.5*dist_mat^2/(l1^2))*(
    dist_mat^2/(l1^3)*sigma1^2/((l1)^4)*(l1^2-dist_mat_x1^2)-
      4*sigma1^2/((l1)^5)*(l1^2-dist_mat_x1^2)+
      2*sigma1^2/((l1)^3))
  del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=exp(-.5*dist_mat^2/(l1^2))*(
    dist_mat^2/(l1^3)*sigma1^2/((l1)^4)*(l1^2-dist_mat_x2^2)-
      4*sigma1^2/((l1)^5)*(l1^2-dist_mat_x2^2)+
      2*sigma1^2/((l1)^3))
  
  del_cov1[1:n1,(n2+1):(2*n2)]=del_cov1[(n1+1):(2*n1),1:n2]=
    diff_mat*(-exp(-.5*dist_mat^2/(l1^2)))*(dist_mat^2/(l1^3)*sigma1^2/((l1)^4)-
                                             4*sigma1^2/((l1)^5))
  
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
  
  del_sigma1=b*del_sigma1+alpha^2*del_cov3
  del_sigma2=b*del_sigma2+alpha^2*del_cov4
 # del_11=b*del_l1+alpha^2*del_cov1
  #del_12=b*del_l2+alpha^2*del_cov2
  
  
  inv_K=inv(K)
  del_alpha= -2*(1-alpha)*K1+2*alpha*K2
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(del_sigma1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(del_sigma2)),
            -t(c(ytr[,1],ytr[,2]))%*%(2*tau*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*tau*inv_K),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*del_l1+alpha^2*del_cov1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(b*del_l1+alpha^2*del_cov1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*del_l2+alpha^2*del_cov2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(b*del_l2+alpha^2*del_cov2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_alpha)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(del_alpha))
  )
  #print(grad)
  return(grad)
  
  
}


Adam = function(f, grad_f, params_init, lr = .01, beta1 = 0.9, beta2 = 0.999,
                epsilon = 1e-8, max_iter = 100, tol = 1e-10,
                lb=rep(0,length(params_init)),
                ub=rep(Inf,length(params_init))) {
  # Number of iterations (for example)
  if(notcluster){pb <- progress_bar$new(
    format = "  Progress [:bar] :percent in :elapsed",
    total = max_iter,
    clear = FALSE,
    width = 60)}
  
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
    if(notcluster){pb$tick()}
    
  }
  
  
  #print("iter= ")
  #print(iter)
  return(params)
}



LogS=RMSES=opt_alpha=opt_alpha_nlopt2=opt_alpha_helm_nlopt=opt_alpha_helm_nlopt2=ALPHA

train_ind=sample(564,100)
test_ind=-train_ind
#print(max_iter)
#train_ind=1:100;test_ind=30:40;max_iter=10

for(i in 1:length(ALPHA)){
  Xtr=X[train_ind,]
  Ytr=Y[train_ind,];F_tr=Eq_noise[train_ind,]
  
  
  Xte=X[test_ind,]
  Yte=Y[test_ind,];F_te=Eq_noise[test_ind,]
  
  
  
  Ytr=(1-ALPHA[i])*Ytr+ALPHA[i]*F_tr
  Yte=(1-ALPHA[i])*Yte+ALPHA[i]*F_te
  
  #test_ind=-train_ind
  
  opt_par=Adam(function(x) {log_likelihood_mix3(
    Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6])
  },function(x) {
    grad_mix3(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6])
  },
  c(initial_sigma1,
    initial_sigma2,
    initial_sigma_obs,
    initial_l1,
    initial_l2,
    initial_alpha),
  lr=lr,max_iter=max_iter,
  lb=c(rep(0.1,5),0.01),ub=c(rep(Inf,5),0.99))
  
  opt_alpha[i]=opt_par[6]
  
  
  Ktetr_mix=mix_cov_mat3(Xte,Xtr,opt_par[1],
                          opt_par[2],
                          
                          opt_par[4],
                          opt_par[5],opt_par[6])
  
  Ktetr_trtr_inv_mix=Ktetr_mix%*%inv(mix_cov_mat3(Xtr,Xtr,opt_par[1],
                                                    opt_par[2],opt_par[4],
                                                    opt_par[5],opt_par[6]) +diag(opt_par[3]^2,nrow=2*nrow(Xtr)))
  
  posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])
  
  
  
  (rmse_mix=mean(apply((cbind(
    posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
    posterior_mean_mix
    [(0.5*length(posterior_mean_mix)+1):(length(posterior_mean_mix))])
    -Yte)^2,1,sum))^.5)
  
  
  
  posterior_cov_mix= mix_cov_mat3(Xte,Xte,opt_par[1],
                                   opt_par[2], opt_par[4],
                                   opt_par[5],opt_par[6])-Ktetr_trtr_inv_mix%*%t(Ktetr_mix)
  
  
  (LogS_mix=(t(c(Yte[,1],Yte[,2])-posterior_mean_mix)%*%
                inv(posterior_cov_mix+diag(opt_par[3]^2,2*nrow(Xte)))%*%
                (c(Yte[,1],Yte[,2])-posterior_mean_mix)
              +determinant(posterior_cov_mix+diag(opt_par[3]^2,2*nrow(Xte))
              )$modulus[1]+
                log(2*pi)*nrow(Xte))/nrow(Xte))
  LogS[i]=LogS_mix
  RMSES[i]=rmse_mix
  
}

res=list(opt_alpha,RMSES,LogS)


save(list = "res", file = paste0("Scores_mix_helm_fund", id, ".rda"))





