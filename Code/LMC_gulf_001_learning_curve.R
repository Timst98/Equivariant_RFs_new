id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

library(foreach)
library(doParallel)

library(nloptr)
library(cubature)

library(proxy)
library(Matrix)
library(pracma)

reduced_sample=FALSE
#data=1 #6hours
data="Gulf" 

lr=.01;max_iter=100
initial_par=rep(1,10)
initial_alpha=.5
initial_sigma_obs=.1
initial_sigma1=initial_sigma2=initial_sigma3=initial_sigma4=initial_sigma5=
  initial_sigma6=initial_sigma7=initial_sigma8=initial_sigma9=initial_sigma10=
  initial_sigma11=initial_sigma12=initial_l1=initial_l2=initial_l3=initial_l4=
  initial_l5=initial_l6=initial_range1=initial_l_range1=initial_range2=initial_l_range2 =1
#Xtr=Xtr[1:3,];Ytr=Ytr[1:3,];Xte=Xte[1:3,];Yte=Yte[1:3,];lb_1=-4;lb_2=-4;ub_1=4;ub_2=4;m1=0;m2=0
noise=.01
N_train=20
tic()

#Xtr=matrix(runif(200,-4,4),ncol=2);N_train=1

norm=function(x){sum(x^2)^(.5)}
cov_mat=cov_mat_LMC=function(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2){
  
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
  cov[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))+
    sigma2^2*exp(-.5*dist_mat^2*(l2^2))
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=
    sigma3^2*exp(-.5*dist_mat^2*(l1^2))+sigma4^2*exp(-.5*dist_mat^2*(l2^2))
  cov[(n1+1):(2*n1),1:n2]=
    cov[1:n1,(n2+1):(2*n2)]=
    sigma1*sigma3*exp(-.5*dist_mat^2*(l1^2))+
    sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  return(cov)
}




cov_mat_helm=function(x1,x2,l1,sigma1,l2,sigma2){
  x1=as.matrix(x1);x2=as.matrix(x2)
  dist_mat=distmat(x1,x2)
  diff_mat=outer(x1[, 1], x2[, 1], "-") * outer(x1[, 2], x2[, 2], "-")
  dist_mat_x1=dist(cbind(x1[,1],0),cbind(x2[,1],0))
  dist_mat_x2=dist(cbind(x1[,2],0),cbind(x2[,2],0))
  
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


fund_cov_mat_LMC = function(x1, x2, sigma1,sigma2,sigma3,sigma4,l1,l2) {
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
                           sigma1,sigma2,sigma3,sigma4,l1,l2)
      
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

mix_cov_mat=function(x1, x2, sigma1,sigma2,sigma3,sigma4,l1,l2,alpha){
  x1=as.matrix(x1);x2=as.matrix(x2)
  K1=fund_cov_mat_LMC(x1, x2, sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  K2=cov_mat(x1, x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  (1-2*alpha)*K1+alpha*K2
}


sigmoid=function(x){1/(1+exp(-x))}

cov_mat_range=function(x1, x2, sigma1,sigma2,sigma3,sigma4,l1,l2,
                       sigma5,sigma6,sigma7,sigma8,l3,l4,
                       sigma9,sigma10,sigma11,sigma12,l5,l6,
                       range1,range2,l_range1,l_range2){
  x1=as.matrix(x1);x2=as.matrix(x2)
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=distmat(t(x1),x2)
  }else{
    dist_mat=distmat(x1,x2)
  }
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  
  sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  sigmoid_dist_center2=(1-sigmoid_21)%*%t(1-sigmoid_22)
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  sigmoid_dist_center2=cbind(rbind(sigmoid_dist_center2,sigmoid_dist_center2),
                             rbind(sigmoid_dist_center2,sigmoid_dist_center2))
  
  sigmoid_se=(sigmoid_11*sigmoid_21%*%t(sigmoid_12*sigmoid_22))
  
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  weight_sum=sigmoid_se+sigmoid_dist_center1+sigmoid_dist_center2
  
  
  cov= cov_mat(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  cov=cov*sigmoid_se
  
  
  
  cov_fund1 = fund_cov_mat_LMC(x1,x2,sigma5,sigma6,sigma7,sigma8,l3,l4)
  cov_fund2 = fund_cov_mat_LMC(x1,x2,sigma9,sigma10,sigma11,sigma12,l5,l6)
  
  cov_fund1=cov_fund1*sigmoid_dist_center1
  cov_fund2=cov_fund2*sigmoid_dist_center2
  return(1/(weight_sum)*(cov+cov_fund2+cov_fund1))
  
}

cov_mat_range_reduced=function(x1, x2, sigma1,sigma2,sigma3,sigma4,l1,l2,
                               sigma5,sigma6,sigma7,sigma8,l3,l4,
                               range1,range2,l_range1,l_range2){
  x1=as.matrix(x1);x2=as.matrix(x2)
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=distmat(t(x1),x2)
  }else{
    dist_mat=distmat(x1,x2)
  }
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  
  sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  sigmoid_dist_center2=(1-sigmoid_21)%*%t(1-sigmoid_22)
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  sigmoid_dist_center2=cbind(rbind(sigmoid_dist_center2,sigmoid_dist_center2),
                             rbind(sigmoid_dist_center2,sigmoid_dist_center2))
  
  sigmoid_se=(sigmoid_11*sigmoid_21%*%t(sigmoid_12*sigmoid_22))
  
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  weight_sum=sigmoid_se+sigmoid_dist_center1+sigmoid_dist_center2
  
  
  cov= cov_mat(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  cov=cov*sigmoid_se
  
  
  
  cov_fund1 = cov_fund2=fund_cov_mat_LMC(x1,x2,sigma5,
                                         sigma6,sigma7,
                                         sigma8,l3,l4)
  
  
  cov_fund1=cov_fund1*sigmoid_dist_center1
  cov_fund2=cov_fund2*sigmoid_dist_center2
  return(1/(weight_sum)*(cov+cov_fund2+cov_fund1))
  
}
cov_mat_range_reduced2=function(x1, x2, sigma1,sigma2,sigma3,sigma4,l1,l2,
                                range1,range2,l_range1,l_range2){
  
  x1=as.matrix(x1);x2=as.matrix(x2)
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=distmat(t(x1),x2)
  }else{
    dist_mat=distmat(x1,x2)
  }
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  
  sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  sigmoid_dist_center2=(1-sigmoid_21)%*%t(1-sigmoid_22)
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  sigmoid_dist_center2=cbind(rbind(sigmoid_dist_center2,sigmoid_dist_center2),
                             rbind(sigmoid_dist_center2,sigmoid_dist_center2))
  
  sigmoid_se=(sigmoid_11*sigmoid_21%*%t(sigmoid_12*sigmoid_22))
  
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  weight_sum=sigmoid_se+sigmoid_dist_center1+sigmoid_dist_center2
  
  
  cov= cov_mat(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  cov=cov*sigmoid_se
  
  
  
  cov_fund1 = cov_fund2=fund_cov_mat_LMC(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  
  cov_fund1=cov_fund1*sigmoid_dist_center1
  cov_fund2=cov_fund2*sigmoid_dist_center2
  return(1/(weight_sum)*(cov+cov_fund2+cov_fund1))
  
}




#No range argument (cov_mat_range 5)
cov_mat_no_range=function(x1, x2, sigma1,sigma2,sigma3,sigma4,l1,l2,
                          sigma5,sigma6,sigma7,sigma8,l3,l4,
                          sigma9,sigma10,sigma11,sigma12,l5,l6,
                          l_range1,l_range2){
  x1=as.matrix(x1);x2=as.matrix(x2)
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    if(n2==1){
      dist_mat=sum((x1-x2)^2)^.5
    }else{dist_mat=distmat(t(x1),x2)}
  }else{
    dist_mat=distmat(x1,x2)
  }
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  
  sigmoid_21=sigmoid((distmat_center2_x1)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2)/l_range2)
  
  sigmoid_dist_center2=(1-sigmoid_21)%*%t(1-sigmoid_22)
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  sigmoid_dist_center2=cbind(rbind(sigmoid_dist_center2,sigmoid_dist_center2),
                             rbind(sigmoid_dist_center2,sigmoid_dist_center2))
  
  sigmoid_se=(sigmoid_11*sigmoid_21%*%t(sigmoid_12*sigmoid_22))
  
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  weight_sum=sigmoid_se+sigmoid_dist_center1+sigmoid_dist_center2
  
  
  
  
  cov= cov_mat(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  cov=cov*sigmoid_se
  
  
  
  cov_fund1 = fund_cov_mat_LMC(x1,x2,sigma5,sigma6,sigma7,sigma8,l3,l4)
  cov_fund2 = fund_cov_mat_LMC(x1,x2,sigma9,sigma10,sigma11,sigma12,l5,l6)
  
  cov_fund1=cov_fund1*sigmoid_dist_center1
  cov_fund2=cov_fund2*sigmoid_dist_center2
  return(1/(weight_sum)*(cov+cov_fund2+cov_fund1))
}

# Range covariance with one center (for GULF)


cov_mat_range_one_center=function(x1, x2, 
                                  sigma1,sigma2,sigma3,sigma4,l1,l2,
                                  sigma5,sigma6,sigma7,sigma8,l3,l4,
                                  range1,l_range1){
  x1=as.matrix(x1);x2=as.matrix(x2)
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=dist(t(x1),x2)
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(-87.1,n1),25.3))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(-87.1,n2),25.3))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  
  sigmoid_se=(sigmoid_11%*%t(sigmoid_12))
  
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  weight_sum=sigmoid_se+sigmoid_dist_center1
  
  
  cov= cov_mat(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  cov=cov*sigmoid_se
  
  
  
  cov_fund1 = fund_cov_mat_LMC(x1,x2,sigma5,sigma6,sigma7,sigma8,l3,l4)
  
  
  cov_fund1=cov_fund1*sigmoid_dist_center1
  
  
  return(1/(weight_sum)*(cov+cov_fund1))
  
}
cov_mat_range_one_center_reduced=function(x1, x2, 
                                          sigma1,sigma2,sigma3,sigma4,l1,l2,
                                          
                                          range1,l_range1){
  x1=as.matrix(x1);x2=as.matrix(x2)
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=dist(t(x1),x2)
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(-87.1,n1),25.3))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(-87.1,n2),25.3))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  
  sigmoid_se=(sigmoid_11%*%t(sigmoid_12))
  
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  weight_sum=sigmoid_se+sigmoid_dist_center1
  
  
  cov= cov_mat(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  cov=cov*sigmoid_se
  
  
  
  cov_fund1 = fund_cov_mat_LMC(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  
  cov_fund1=cov_fund1*sigmoid_dist_center1
  
  
  return(1/(weight_sum)*(cov+cov_fund1))
  
}


# Range covariance with one center, no range parameter

cov_mat_range_one_center_norange=function(x1, x2, 
                                          sigma1,sigma2,sigma3,sigma4,l1,l2,
                                          sigma5,sigma6,sigma7,sigma8,l3,l4,
                                          l_range1){
  x1=as.matrix(x1);x2=as.matrix(x2)
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=dist(t(x1),x2)
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(-87.1,n1),25.3))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(-87.1,n2),25.3))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  
  sigmoid_se=(sigmoid_11%*%t(sigmoid_12))
  
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  weight_sum=sigmoid_se+sigmoid_dist_center1
  
  
  cov= cov_mat(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  cov=cov*sigmoid_se
  
  
  
  cov_fund1 = fund_cov_mat_LMC(x1,x2,sigma5,sigma6,sigma7,sigma8,l3,l4)
  
  
  cov_fund1=cov_fund1*sigmoid_dist_center1
  
  return(1/(weight_sum)*(cov+cov_fund1))
  
}

cov_mat_range_one_center_norange_reduced=function(x1, x2, 
                                                  sigma1,sigma2,sigma3,sigma4,l1,l2,
                                                  
                                                  l_range1){
  x1=as.matrix(x1);x2=as.matrix(x2)
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=dist(t(x1),x2)
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(-87.1,n1),25.3))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(-87.1,n2),25.3))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  
  sigmoid_se=(sigmoid_11%*%t(sigmoid_12))
  
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  weight_sum=sigmoid_se+sigmoid_dist_center1
  
  
  cov= cov_mat(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  cov=cov*sigmoid_se
  
  
  
  cov_fund1 = fund_cov_mat_LMC(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  
  cov_fund1=cov_fund1*sigmoid_dist_center1
  
  return(1/(weight_sum)*(cov+cov_fund1))
  
}






log_likelihood=function(ytr,xtr,sigma1,sigma2,sigma3,sigma4,tau,l1,l2){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  Ktr=cov_mat(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}

log_likelihood_fund=function(ytr,xtr,sigma1,sigma2,sigma3,sigma4,tau,l1,l2){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  Ktr=fund_cov_mat_LMC(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}

log_likelihood_mixture=function(ytr,xtr,sigma1,sigma2,sigma3,sigma4,tau,l1,l2,alpha){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=mix_cov_mat(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,alpha)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}
#sigma1=sigma2=sigma3=sigma4=l1=l2=sigma5=sigma6=sigma7=sigma8=l3=l4=sigma9=
#  sigma10=sigma11=sigma12=l5=l6=range1=range2=l_range1=l_range2=sigma_obs=1

log_likelihood_range=function(ytr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                              sigma5,sigma6,sigma7,sigma8,l3,l4,
                              sigma9,sigma10,sigma11,sigma12,l5,l6,
                              range1,range2,l_range1,l_range2,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=cov_mat_range(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                    sigma5,sigma6,sigma7,sigma8,l3,l4,
                    sigma9,sigma10,sigma11,sigma12,l5,l6,
                    range1,range2,l_range1,l_range2)+
    sigma_obs^2*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}
log_likelihood_range_reduced=function(ytr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                                      sigma5,sigma6,sigma7,sigma8,l3,l4,
                                      range1,range2,l_range1,l_range2,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=cov_mat_range_reduced(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                            sigma5,sigma6,sigma7,sigma8,l3,l4,
                            
                            range1,range2,l_range1,l_range2)+
    sigma_obs^2*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}

log_likelihood_range_reduced2=function(ytr,xtr,sigma1,sigma2,sigma3,sigma4,
                                       l1,l2,range1,range2,l_range1,
                                       l_range2,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=cov_mat_range_reduced2(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                             range1,range2,l_range1,l_range2)+
    sigma_obs^2*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}


log_likelihood_helm=function(ytr,xtr,l1,sigma1,l2,sigma2,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  Ktr=cov_mat_helm(xtr,xtr,l1,sigma1,l2,sigma2)+sigma_obs^2*
    diag(nrow=2*nrow(xtr))
  (ll=-.5*((ytr)%*%inv(Ktr)%*%ytr)-.5*determinant(Ktr,logarithm=1)$modulus-
      nrow(xtr)*log(2*pi))
  ll=ifelse(is.na(ll),1e6,-ll)
  #print(ll)
  return(ll)
}



log_likelihood_no_range=function(ytr,xtr,
                                 sigma1,sigma2,sigma3,sigma4,l1,l2,
                                 sigma5,sigma6,sigma7,sigma8,l3,l4,
                                 sigma9,sigma10,sigma11,sigma12,l5,l6,
                                 l_range1,l_range2,sigma_obs){
  if(length(ytr)>2){
    if(ncol(ytr)==2){
      ytr=c(ytr[,1],ytr[,2])
    }
  }
  ntr=ifelse(length(xtr)==2,1,nrow(xtr))
  
  Ktr=cov_mat_no_range(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                       sigma5,sigma6,sigma7,sigma8,l3,l4,
                       sigma9,sigma10,sigma11,sigma12,l5,l6,l_range1,l_range2)+
    sigma_obs^2*diag(nrow=2*ntr)
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    ntr*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}
log_likelihood_range_one_center=function(ytr,xtr,
                                         sigma1,sigma2,sigma3,sigma4,l1,l2,
                                         sigma5,sigma6,sigma7,sigma8,l3,l4,
                                         range1,l_range1,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=cov_mat_range_one_center(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                               sigma5,sigma6,sigma7,sigma8,l3,l4,
                               range1,l_range1)+
    sigma_obs^2*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}

log_likelihood_range_one_center_norange=function(ytr,xtr,
                                                 sigma1,sigma2,sigma3,sigma4,l1,l2,
                                                 sigma5,sigma6,sigma7,sigma8,l3,l4,
                                                 l_range1,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=cov_mat_range_one_center_norange(xtr,xtr,
                                       sigma1,sigma2,sigma3,sigma4,l1,l2,
                                       sigma5,sigma6,sigma7,sigma8,l3,l4,
                                       
                                       l_range1)+
    sigma_obs^2*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}

log_likelihood_range_one_center_reduced=function(ytr,xtr,
                                                 sigma1,sigma2,sigma3,sigma4,l1,l2,
                                                 
                                                 range1,l_range1,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=cov_mat_range_one_center_reduced(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                                       
                                       range1,l_range1)+
    sigma_obs^2*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}

log_likelihood_range_one_center_norange_reduced=function(ytr,xtr,
                                                         sigma1,sigma2,sigma3,sigma4,l1,l2,
                                                         
                                                         l_range1,sigma_obs){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  Ktr=cov_mat_range_one_center_norange_reduced(xtr,xtr,
                                               sigma1,sigma2,sigma3,sigma4,l1,l2,
                                               
                                               
                                               l_range1)+
    sigma_obs^2*diag(nrow=2*nrow(xtr))
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm=1)$modulus-
    nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}



del_sig1=function(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2){
  
  
  
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
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  return(del_sigma1)
}

del_sig2=function(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2){
  
  
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
  del_sigma2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  return(del_sigma2)
}
del_sig3=function(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2){
  
  
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
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  return(del_sigma3)
}

del_sig4=function(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2){
  
  
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
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  return(del_sigma4)
}

del_l_1=function(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2){
  
  
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
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  return(del_l1)
}

del_l_2=function(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2){
  
  
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
  del_l2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  return(del_l2)
}

grad_fund=function(xtr,ytr,sigma1,sigma2,sigma3,sigma4,tau,l1,l2){
  x1=x2=xtr
  K=fund_cov_mat_LMC(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  cov = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1,nrow(x2)),
               nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  
  del_sigma1=del_sigma2=del_sigma3=del_sigma4=del_l1=del_l2=cov
  
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
                           sigma1,sigma2,sigma3,sigma4,l1,l2)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard %*%cbind(c(cos(theta2),-sin(theta2)),
                              c(sin(theta2),cos(theta2)))
      
      cov[i, j] = results[1,1]
      cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      
      results_sigma1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
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
      
      
      
      del_sigma3[i, j] = results_sigma3[1,1]
      del_sigma3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma3[2,2]
      del_sigma3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma3[2,1]
      del_sigma3[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma3[1,2]
      
      
      
      del_sigma4[i, j] = results_sigma4[1,1]
      del_sigma4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma4[2,2]
      del_sigma4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma4[2,1]
      del_sigma4[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma4[1,2]
      
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
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma4)),
            -t(c(ytr[,1],ytr[,2]))%*%(2*tau*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*tau*inv_K),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l2))
  )
  #print(grad)
  return(grad)
  
  
}

grad_standard=function(xtr,ytr,sigma1,sigma2,sigma3,sigma4,tau,l1,l2){
  x1=x2=xtr
  K=cov_mat(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2)+
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
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  del_tau=2*tau*diag(nrow=2*nrow(xtr))
  
  
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  
  inv_K=inv(K)
  
  #sigma1,sigma2,sigma3,sigma4,nu1=1,nu2=1,tau,a1,a2
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma1),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma2),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma3%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma3),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma4%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma4),
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
  dist_mat_x1=dist(cbind(x1[,1],0),cbind(x2[,1],0))
  dist_mat_x2=dist(cbind(x1[,2],0),cbind(x2[,2],0))
  
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

grad_no_range=function(xtr,ytr,
                       sigma1,sigma2,sigma3,sigma4,l1,l2,
                       sigma5,sigma6,sigma7,sigma8,l3,l4,
                       sigma9,sigma10,sigma11,sigma12,l5,l6,
                       l_range1,l_range2,sigma_obs_2){
  
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  del_l_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1)/l_range1)*
                    (distmat_center1_x1)/l_range1^2)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2)/l_range1)*
        (distmat_center1_x2)/l_range1^2)
  
  
  del_l_range1_1=cbind(rbind(del_l_range1_1,del_l_range1_1),
                       rbind(del_l_range1_1,del_l_range1_1))
  
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  
  sigmoid_21=sigmoid((distmat_center2_x1)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2)/l_range2)
  
  sigmoid_dist_center2=(1-sigmoid_21)%*%t(1-sigmoid_22)
  
  
  del_l_range2_1=((sigmoid_21)^2*exp(-(distmat_center2_x1)/l_range2)*
                    (distmat_center2_x1)/l_range2^2)%*%t(1-sigmoid_22)+
    (1-sigmoid_21)%*%t(
      (sigmoid_22)^2*exp(-(distmat_center2_x2)/l_range2)*
        (distmat_center2_x2)/l_range2^2)
  
  del_l_range2_1=cbind(rbind(del_l_range2_1,del_l_range2_1),
                       rbind(del_l_range2_1,del_l_range2_1))
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  sigmoid_dist_center2=cbind(rbind(sigmoid_dist_center2,sigmoid_dist_center2),
                             rbind(sigmoid_dist_center2,sigmoid_dist_center2))
  
  sigmoid_se=(sigmoid_11*sigmoid_21%*%t(sigmoid_12*sigmoid_22))
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  
  
  
  del_l_range1_2=((sigmoid_11^2*exp(-(distmat_center1_x1)/l_range1)*
                     (-distmat_center1_x1)/(l_range1^2))*sigmoid_21)%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t((sigmoid_12^2*exp(-(distmat_center1_x2)/l_range1)*
                                   (-distmat_center1_x2)/(l_range1^2))*sigmoid_22)
  
  
  del_l_range1_2=cbind(rbind(del_l_range1_2,del_l_range1_2),
                       rbind(del_l_range1_2,del_l_range1_2))
  
  
  
  del_l_range2_2=(sigmoid_11*(sigmoid_21^2*exp(-(distmat_center2_x1)/l_range2)*
                                (-distmat_center2_x1)/l_range2^2))%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t(sigmoid_12*(sigmoid_22^2*exp(-(distmat_center2_x2)/l_range2)*
                                              (-distmat_center2_x2)/l_range2^2))
  
  del_l_range2_2=cbind(rbind(del_l_range2_2,del_l_range2_2),
                       rbind(del_l_range2_2,del_l_range2_2))
  
  
  sigmoid_cov=sigmoid_se
  weight_sum=sigmoid_se+sigmoid_dist_center1+sigmoid_dist_center2
  
  K1_sigma5=K1_sigma6=K1_sigma7=K1_sigma8=K1_l3=K1_l4=
    matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
           nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K2_sigma9=K2_sigma10=K2_sigma11=K2_sigma12=K2_l5=K2_l6=K1_sigma5
  
  K=cov_mat_no_range(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                     sigma5,sigma6,sigma7,sigma8,l3,l4,
                     sigma9,sigma10,sigma11,sigma12,l5,l6,
                     l_range1,l_range2)+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  cov_fund1 = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                     nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  cov_fund2=cov_fund1
  
  
  for(i in 1:n1){
    #if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                sum((x1[i,])^2)^.5)
      r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                sum((x2[j,])^2)^.5)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))
      
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma5,sigma6,sigma7,sigma8,l3,l4)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund1[i, j] = results[1,1]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund1[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma9,sigma10,sigma11,sigma12,l5,l6)
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund2[i, j] = results[1,1]
      cov_fund2[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund2[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund2[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      results_sigma5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma7=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma8=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      #####
      results_sigma9=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma10=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma11=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma12=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      #####
      
      
      K1_sigma5[i, j] = results_sigma5[1,1]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[2,2]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma5[2,1]
      K1_sigma5[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[1,2]
      
      
      
      K1_sigma6[i, j] = results_sigma6[1,1]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[2,2]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma6[2,1]
      K1_sigma6[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[1,2]
      
      
      
      K1_sigma7[i, j] = results_sigma7[1,1]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[2,2]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma7[2,1]
      K1_sigma7[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[1,2]
      
      
      
      K1_sigma8[i, j] = results_sigma8[1,1]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[2,2]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma8[2,1]
      K1_sigma8[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[1,2]
      
      K1_l3[i, j] = results_l3[1,1]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[2,2]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l3[2,1]
      K1_l3[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[1,2]
      
      K1_l4[i, j] = results_l4[1,1]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[2,2]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l4[2,1]
      K1_l4[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[1,2]
      
      
      K2_sigma9[i, j] = results_sigma9[1,1]
      K2_sigma9[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma9[2,2]
      K2_sigma9[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma9[2,1]
      K2_sigma9[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma9[1,2]
      
      
      
      K2_sigma10[i, j] = results_sigma10[1,1]
      K2_sigma10[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma10[2,2]
      K2_sigma10[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma10[2,1]
      K2_sigma10[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma10[1,2]
      
      
      
      K2_sigma11[i, j] = results_sigma11[1,1]
      K2_sigma11[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma11[2,2]
      K2_sigma11[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma11[2,1]
      K2_sigma11[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma11[1,2]
      
      
      
      K2_sigma12[i, j] = results_sigma12[1,1]
      K2_sigma12[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma12[2,2]
      K2_sigma12[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma12[2,1]
      K2_sigma12[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma12[1,2]
      
      K2_l5[i, j] = results_l5[1,1]
      K2_l5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l5[2,2]
      K2_l5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l5[2,1]
      K2_l5[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l5[1,2]
      
      K2_l6[i, j] = results_l6[1,1]
      K2_l6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l6[2,2]
      K2_l6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l6[2,1]
      K2_l6[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l6[1,2]
      
    }
  }
  
  K1_sigma5=K1_sigma5*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma6=K1_sigma6*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma7=K1_sigma7*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma8=K1_sigma8*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l3=K1_l3*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l4=K1_l4*(sigmoid_dist_center1)*1/(weight_sum)
  
  K2_sigma9=K2_sigma9*(sigmoid_dist_center2)*1/(weight_sum)
  K2_sigma10=K2_sigma10*(sigmoid_dist_center2)*1/(weight_sum)
  K2_sigma11=K2_sigma11*(sigmoid_dist_center2)*1/(weight_sum)
  K2_sigma12=K2_sigma12*(sigmoid_dist_center2)*1/(weight_sum)
  K2_l5=K2_l5*(sigmoid_dist_center2)*1/(weight_sum)
  K2_l6=K2_l6*(sigmoid_dist_center2)*1/(weight_sum)
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    if(n2==1){
      dist_mat=sum((x1-x2)^2)^.5
    }else{dist_mat=dist(t(x1),x2)}
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  
  del_sigma1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  
  
  del_sigma1=sigmoid_cov*del_sigma1*1/(weight_sum)
  del_sigma2=sigmoid_cov*del_sigma2*1/(weight_sum)
  del_sigma3=sigmoid_cov*del_sigma3*1/(weight_sum)
  del_sigma4=sigmoid_cov*del_sigma4*1/(weight_sum)
  
  del_l1=sigmoid_cov*del_l1*1/(weight_sum)
  del_l2=sigmoid_cov*del_l2*1/(weight_sum)
  
  K_fund1=cov_fund1
  K_fund2=cov_fund2
  
  K_cov=cov_mat(x1,x1,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  del_l_range1=(K_fund1*del_l_range1_1+K_cov*del_l_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_l_range1_2+del_l_range1_1)
  
  del_l_range2=(K_fund2*del_l_range2_1+K_cov*del_l_range2_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_l_range2_2+del_l_range2_1)
  
  
  
  inv_K=inv(K)
  b=1
  alpha=1
  
  
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma5)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma5)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma6)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma6)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma7)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma7)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma8)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma8)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma9)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K2_sigma9)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma10)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K2_sigma10)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma11)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K2_sigma11)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma12)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K2_sigma12)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_l5)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K2_l5)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_l6)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K2_l6)),
            
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range2)),
            
            -t(c(ytr[,1],ytr[,2]))%*%(2*sigma_obs_2*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))  #print(grad)
  return(grad)
  
}

grad_range=function(xtr,ytr,
                    sigma1,sigma2,sigma3,sigma4,l1,l2,
                    sigma5,sigma6,sigma7,sigma8,l3,l4,
                    sigma9,sigma10,sigma11,sigma12,l5,l6,         
                    range1,range2,l_range1,l_range2,sigma_obs_2){
  
  
  
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  del_l_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                    (distmat_center1_x1-range1)/l_range1^2)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)*
        (distmat_center1_x2-range1)/l_range1^2)
  
  
  del_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)/l_range1)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)/l_range1)
  
  
  del_l_range1_1=cbind(rbind(del_l_range1_1,del_l_range1_1),
                       rbind(del_l_range1_1,del_l_range1_1))
  
  
  del_range1_1=cbind(rbind(del_range1_1,del_range1_1),
                     rbind(del_range1_1,del_range1_1))
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  
  sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  sigmoid_dist_center2=(1-sigmoid_21)%*%t(1-sigmoid_22)
  
  
  del_l_range2_1=((sigmoid_21)^2*exp(-(distmat_center2_x1-range2)/l_range2)*
                    (distmat_center2_x1-range2)/l_range2^2)%*%t(1-sigmoid_22)+
    (1-sigmoid_21)%*%t(
      (sigmoid_22)^2*exp(-(distmat_center2_x2-range2)/l_range2)*
        (distmat_center2_x2-range2)/l_range2^2)
  
  del_l_range2_1=cbind(rbind(del_l_range2_1,del_l_range2_1),
                       rbind(del_l_range2_1,del_l_range2_1))
  
  
  del_range2_1=((sigmoid_21)^2*exp(-(distmat_center2_x1-range2)/l_range2)/l_range2)%*%t(1-sigmoid_22)+
    (1-sigmoid_21)%*%t(
      (sigmoid_22)^2*exp(-(distmat_center2_x2-range2)/l_range2)/l_range2)
  
  
  del_range2_1=cbind(rbind(del_range2_1,del_range2_1),
                     rbind(del_range2_1,del_range2_1))
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  
  sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  sigmoid_dist_center2=cbind(rbind(sigmoid_dist_center2,sigmoid_dist_center2),
                             rbind(sigmoid_dist_center2,sigmoid_dist_center2))
  
  sigmoid_se=(sigmoid_11*sigmoid_21%*%t(sigmoid_12*sigmoid_22))
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  
  
  
  del_l_range1_2=((sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                     (-distmat_center1_x1+range1)/(l_range1^2))*sigmoid_21)%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t((sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)*
                                   (-distmat_center1_x2+range1)/(l_range1^2))*sigmoid_22)
  
  
  del_l_range1_2=cbind(rbind(del_l_range1_2,del_l_range1_2),
                       rbind(del_l_range1_2,del_l_range1_2))
  
  del_range1_2=((-sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)/l_range1)*sigmoid_21)%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t((-sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)/l_range1)*sigmoid_22)
  
  
  del_range1_2=cbind(rbind(del_range1_2,del_range1_2),
                     rbind(del_range1_2,del_range1_2))
  
  
  
  
  del_l_range2_2=(sigmoid_11*(sigmoid_21^2*exp(-(distmat_center2_x1-range2)/l_range2)*
                                (-distmat_center2_x1+range2)/l_range2^2))%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t(sigmoid_12*(sigmoid_22^2*exp(-(distmat_center2_x2-range2)/l_range2)*
                                              (-distmat_center2_x2+range2)/l_range2^2))
  
  del_l_range2_2=cbind(rbind(del_l_range2_2,del_l_range2_2),
                       rbind(del_l_range2_2,del_l_range2_2))
  
  del_range2_2=(sigmoid_11*(-sigmoid_21^2*exp(-(distmat_center2_x1-range2)/l_range2)/l_range2))%*%
    t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t(sigmoid_12*(-sigmoid_22^2*exp(-(distmat_center2_x2-range2)/l_range2)/l_range2))
  
  del_range2_2=cbind(rbind(del_range2_2,del_range2_2),
                     rbind(del_range2_2,del_range2_2))
  
  
  sigmoid_cov=sigmoid_se
  weight_sum=sigmoid_se+sigmoid_dist_center1+sigmoid_dist_center2
  
  K1_sigma5=K1_sigma6=K1_sigma7=K1_sigma8=K1_l3=K1_l4=
    matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
           nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K2_sigma9=K2_sigma10=K2_sigma11=K2_sigma12=K2_l5=K2_l6=K1_sigma5
  
  K=cov_mat_range(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                  sigma5,sigma6,sigma7,sigma8,l3,l4,
                  sigma9,sigma10,sigma11,sigma12,l5,l6,
                  range1,range2,l_range1,l_range2)+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  cov_fund1 = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                     nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  cov_fund2=cov_fund1
  
  
  for(i in 1:n1){
    #if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                sum((x1[i,])^2)^.5)
      r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                sum((x2[j,])^2)^.5)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))
      
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma5,sigma6,sigma7,sigma8,l3,l4)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund1[i, j] = results[1,1]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund1[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma9,sigma10,sigma11,sigma12,l5,l6)
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund2[i, j] = results[1,1]
      cov_fund2[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund2[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund2[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      results_sigma5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma7=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma8=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      #####
      results_sigma9=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma10=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma11=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma12=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma9,sigma10,sigma11,sigma12,l5,l6)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      #####
      
      
      K1_sigma5[i, j] = results_sigma5[1,1]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[2,2]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma5[2,1]
      K1_sigma5[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[1,2]
      
      
      
      K1_sigma6[i, j] = results_sigma6[1,1]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[2,2]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma6[2,1]
      K1_sigma6[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[1,2]
      
      
      
      K1_sigma7[i, j] = results_sigma7[1,1]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[2,2]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma7[2,1]
      K1_sigma7[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[1,2]
      
      
      
      K1_sigma8[i, j] = results_sigma8[1,1]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[2,2]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma8[2,1]
      K1_sigma8[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[1,2]
      
      K1_l3[i, j] = results_l3[1,1]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[2,2]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l3[2,1]
      K1_l3[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[1,2]
      
      K1_l4[i, j] = results_l4[1,1]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[2,2]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l4[2,1]
      K1_l4[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[1,2]
      
      
      K2_sigma9[i, j] = results_sigma9[1,1]
      K2_sigma9[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma9[2,2]
      K2_sigma9[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma9[2,1]
      K2_sigma9[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma9[1,2]
      
      
      
      K2_sigma10[i, j] = results_sigma10[1,1]
      K2_sigma10[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma10[2,2]
      K2_sigma10[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma10[2,1]
      K2_sigma10[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma10[1,2]
      
      
      
      K2_sigma11[i, j] = results_sigma11[1,1]
      K2_sigma11[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma11[2,2]
      K2_sigma11[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma11[2,1]
      K2_sigma11[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma11[1,2]
      
      
      
      K2_sigma12[i, j] = results_sigma12[1,1]
      K2_sigma12[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma12[2,2]
      K2_sigma12[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma12[2,1]
      K2_sigma12[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma12[1,2]
      
      K2_l5[i, j] = results_l5[1,1]
      K2_l5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l5[2,2]
      K2_l5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l5[2,1]
      K2_l5[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l5[1,2]
      
      K2_l6[i, j] = results_l6[1,1]
      K2_l6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l6[2,2]
      K2_l6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l6[2,1]
      K2_l6[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l6[1,2]
      
    }
  }
  
  K1_sigma5=K1_sigma5*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma6=K1_sigma6*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma7=K1_sigma7*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma8=K1_sigma8*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l3=K1_l3*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l4=K1_l4*(sigmoid_dist_center1)*1/(weight_sum)
  
  K2_sigma9=K2_sigma9*(sigmoid_dist_center2)*1/(weight_sum)
  K2_sigma10=K2_sigma10*(sigmoid_dist_center2)*1/(weight_sum)
  K2_sigma11=K2_sigma11*(sigmoid_dist_center2)*1/(weight_sum)
  K2_sigma12=K2_sigma12*(sigmoid_dist_center2)*1/(weight_sum)
  K2_l5=K2_l5*(sigmoid_dist_center2)*1/(weight_sum)
  K2_l6=K2_l6*(sigmoid_dist_center2)*1/(weight_sum)
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    if(n2==1){
      dist_mat=sum((x1-x2)^2)^.5
    }else{dist_mat=dist(t(x1),x2)}
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  
  del_sigma1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  
  
  del_sigma1=sigmoid_cov*del_sigma1*1/(weight_sum)
  del_sigma2=sigmoid_cov*del_sigma2*1/(weight_sum)
  del_sigma3=sigmoid_cov*del_sigma3*1/(weight_sum)
  del_sigma4=sigmoid_cov*del_sigma4*1/(weight_sum)
  
  del_l1=sigmoid_cov*del_l1*1/(weight_sum)
  del_l2=sigmoid_cov*del_l2*1/(weight_sum)
  
  K_fund1=cov_fund1
  K_fund2=cov_fund2
  
  K_cov=cov_mat(x1,x1,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  del_l_range1=(K_fund1*del_l_range1_1+K_cov*del_l_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_l_range1_2+del_l_range1_1)
  
  del_l_range2=(K_fund2*del_l_range2_1+K_cov*del_l_range2_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_l_range2_2+del_l_range2_1)
  
  del_range1=(K_fund1*(del_range1_1)+K_cov*del_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_range1_2+del_range1_1)
  
  
  del_range2=(K_fund2*(del_range2_1)+K_cov*del_range2_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_range2_2+del_range2_1)
  
  
  inv_K=inv(K)
  b=1
  alpha=1
  
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma5)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma5)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma6)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma6)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma7)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma7)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma8)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma8)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma9)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K2_sigma9)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma10)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K2_sigma10)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma11)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K2_sigma11)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_sigma12)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K2_sigma12)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_l5)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K2_l5)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K2_l6)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K2_l6)),
            
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range2)),
            
            -t(c(ytr[,1],ytr[,2]))%*%(2*sigma_obs_2*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  #print(grad)
  return(grad)
  
}





grad_mix=function(xtr,ytr,sigma1,sigma2,sigma3,sigma4,tau,l1,l2,alpha){
  
  x1=x2=xtr
  K1=fund_cov_mat_LMC(x1, x2, sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  K2=cov_mat(x1, x2,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  sigma_obs_2=tau
  
  K=(1-2*alpha)*K1+alpha*K2+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  cov = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1,nrow(x2)),
               nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  
  del_sigma1=del_sigma2=del_sigma3=del_sigma4=del_l1=del_l2=cov
  
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
      
      
      results_sigma1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
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
      
      
      
      del_sigma3[i, j] = results_sigma3[1,1]
      del_sigma3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma3[2,2]
      del_sigma3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma3[2,1]
      del_sigma3[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma3[1,2]
      
      
      
      del_sigma4[i, j] = results_sigma4[1,1]
      del_sigma4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma4[2,2]
      del_sigma4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma4[2,1]
      del_sigma4[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma4[1,2]
      
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
  
  
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  if(n1==1){
    dist_mat=distmat(t(xtr),xtr)
  }else{
    dist_mat=distmat(xtr,xtr)
  }
  b=1-2*alpha
  
  del_sigma1_2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1_2[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1_2[(n1+1):(2*n1),1:n2]=del_sigma1_2[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1_2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2_2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2_2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2_2[(n1+1):(2*n1),1:n2]=del_sigma2_2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2_2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  
  del_sigma3_2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3_2[1:n1,1:n2]=0
  del_sigma3_2[(n1+1):(2*n1),1:n2]=del_sigma3_2[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3_2[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  
  del_sigma4_2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4_2[1:n1,1:n2]=0
  del_sigma4_2[(n1+1):(2*n1),1:n2]=del_sigma4_2[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4_2[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  del_tau=2*tau*diag(nrow=2*nrow(xtr))
  
  
  
  del_l1_2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1_2[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1_2[(n1+1):(2*n1),1:n2]=del_l1_2[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1_2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  
  del_l2_2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2_2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2_2[(n1+1):(2*n1),1:n2]=del_l2_2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2_2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  del_alpha=K2-2*K1
  
  inv_K=inv(K)
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*del_sigma1+alpha*del_sigma1_2)%*%inv_K%*%c(ytr[,1],ytr[,2])+
              Trace(inv_K%*%((b*del_sigma1+alpha*del_sigma1_2))),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%((b*del_sigma2+alpha*del_sigma2_2))%*%inv_K%*%c(ytr[,1],ytr[,2])+
              Trace(inv_K%*%(b*del_sigma2+alpha*del_sigma2_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*del_sigma3+alpha*del_sigma3_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+
              Trace(inv_K%*%(b*del_sigma3+alpha*del_sigma3_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*del_sigma4+alpha*del_sigma4_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+
              Trace(inv_K%*%(b*del_sigma4+alpha*del_sigma4_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_tau%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_tau),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*del_l1+alpha*del_l1_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+
              Trace(inv_K%*%(b*del_l1+alpha*del_l1_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(b*del_l2+alpha*del_l2_2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+
              Trace(inv_K%*%(b*del_l2+alpha*del_l2_2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_alpha%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_alpha))
  #print(grad)
  return(grad)
  
}


grad_range_one_center=function(xtr,ytr,
                               sigma1,sigma2,sigma3,sigma4,l1,l2,
                               sigma5,sigma6,sigma7,sigma8,l3,l4,
                               range1,l_range1,sigma_obs_2){
  
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(-87.1,n1),25.3))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(-87.1,n2),25.3))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  del_l_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                    (distmat_center1_x1-range1)/l_range1^2)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)*
        (distmat_center1_x2-range1)/l_range1^2)
  
  
  del_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)/l_range1)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)/l_range1)
  
  
  del_l_range1_1=cbind(rbind(del_l_range1_1,del_l_range1_1),
                       rbind(del_l_range1_1,del_l_range1_1))
  
  
  del_range1_1=cbind(rbind(del_range1_1,del_range1_1),
                     rbind(del_range1_1,del_range1_1))
  
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  
  sigmoid_se=(sigmoid_11%*%t(sigmoid_12))
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  
  
  
  del_l_range1_2=((sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                     (-distmat_center1_x1+range1)/(l_range1^2)))%*%t(sigmoid_12)+
    (sigmoid_11)%*%t((sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)*
                        (-distmat_center1_x2+range1)/(l_range1^2)))
  
  
  del_l_range1_2=cbind(rbind(del_l_range1_2,del_l_range1_2),
                       rbind(del_l_range1_2,del_l_range1_2))
  
  del_range1_2=((-sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)/l_range1))%*%t(sigmoid_12)+
    (sigmoid_11)%*%t((-sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)/l_range1))
  
  
  del_range1_2=cbind(rbind(del_range1_2,del_range1_2),
                     rbind(del_range1_2,del_range1_2))
  
  
  
  
  
  sigmoid_cov=sigmoid_se
  weight_sum=sigmoid_se+sigmoid_dist_center1
  
  
  
  
  K=cov_mat_range_one_center(xtr,xtr,
                             sigma1,sigma2,sigma3,sigma4,l1,l2,
                             sigma5,sigma6,sigma7,sigma8,l3,l4
                             ,range1,l_range1)+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  
  K1_sigma5=K1_sigma6=K1_sigma7=K1_sigma8=K1_l3=K1_l4=
    matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
           nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  cov_fund1 = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                     nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  
  
  
  for(i in 1:n1){
    #if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                sum((x1[i,])^2)^.5)
      r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                sum((x2[j,])^2)^.5)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))
      
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma5,sigma6,sigma7,sigma8,l3,l4)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund1[i, j] = results[1,1]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund1[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      
      results_sigma5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma7=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma8=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      K1_sigma5[i, j] = results_sigma5[1,1]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[2,2]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma5[2,1]
      K1_sigma5[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[1,2]
      
      
      
      K1_sigma6[i, j] = results_sigma6[1,1]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[2,2]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma6[2,1]
      K1_sigma6[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[1,2]
      
      
      
      K1_sigma7[i, j] = results_sigma7[1,1]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[2,2]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma7[2,1]
      K1_sigma7[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[1,2]
      
      
      
      K1_sigma8[i, j] = results_sigma8[1,1]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[2,2]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma8[2,1]
      K1_sigma8[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[1,2]
      
      K1_l3[i, j] = results_l3[1,1]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[2,2]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l3[2,1]
      K1_l3[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[1,2]
      
      K1_l4[i, j] = results_l4[1,1]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[2,2]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l4[2,1]
      K1_l4[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[1,2]
      
      
    }
  }
  
  K1_sigma5=K1_sigma5*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma6=K1_sigma6*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma7=K1_sigma7*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma8=K1_sigma8*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l3=K1_l3*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l4=K1_l4*(sigmoid_dist_center1)*1/(weight_sum)
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    if(n2==1){
      dist_mat=sum((x1-x2)^2)^.5
    }else{dist_mat=distmat(t(x1),x2)}
  }else{
    dist_mat=distmat(x1,x2)
  }
  
  
  
  del_sigma1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  
  
  del_sigma1=sigmoid_cov*del_sigma1*1/(weight_sum)
  del_sigma2=sigmoid_cov*del_sigma2*1/(weight_sum)
  del_sigma3=sigmoid_cov*del_sigma3*1/(weight_sum)
  del_sigma4=sigmoid_cov*del_sigma4*1/(weight_sum)
  
  del_l1=sigmoid_cov*del_l1*1/(weight_sum)
  del_l2=sigmoid_cov*del_l2*1/(weight_sum)
  
  K_fund1=cov_fund1
  
  
  K_cov=cov_mat(x1,x1,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  K_fund1=cov_fund1
  
  
  
  
  del_l_range1=(K_fund1*del_l_range1_1+K_cov*del_l_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se)/(weight_sum^2)*
    (del_l_range1_2+del_l_range1_1)
  
  
  del_range1=(K_fund1*(del_range1_1)+K_cov*del_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se)/(weight_sum^2)*
    (del_range1_2+del_range1_1)
  
  
  
  inv_K=inv(K)
  b=1
  alpha=1
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l2)),
            
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma5)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma5)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma6)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma6)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma7)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma7)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma8)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma8)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l4)),
            
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range1)),
            
            -t(c(ytr[,1],ytr[,2]))%*%(2*sigma_obs_2*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  #print(grad)
  return(grad)
  
}
grad_range_one_center_norange=function(xtr,ytr,
                                       sigma1,sigma2,sigma3,sigma4,l1,l2,
                                       sigma5,sigma6,sigma7,sigma8,l3,l4,
                                       l_range1,sigma_obs_2){
  
  
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(-87.1,n1),25.3))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(-87.1,n2),25.3))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  del_l_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1)/l_range1)*
                    (distmat_center1_x1)/l_range1^2)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2)/l_range1)*
        (distmat_center1_x2)/l_range1^2)
  
  
  
  
  
  del_l_range1_1=cbind(rbind(del_l_range1_1,del_l_range1_1),
                       rbind(del_l_range1_1,del_l_range1_1))
  
  
  
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  
  sigmoid_se=(sigmoid_11%*%t(sigmoid_12))
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  
  
  
  del_l_range1_2=((sigmoid_11^2*exp(-(distmat_center1_x1)/l_range1)*
                     (-distmat_center1_x1)/(l_range1^2)))%*%t(sigmoid_12)+
    (sigmoid_11)%*%t((sigmoid_12^2*exp(-(distmat_center1_x2)/l_range1)*
                        (-distmat_center1_x2)/(l_range1^2)))
  
  
  del_l_range1_2=cbind(rbind(del_l_range1_2,del_l_range1_2),
                       rbind(del_l_range1_2,del_l_range1_2))
  
  
  
  
  
  
  sigmoid_cov=sigmoid_se
  weight_sum=sigmoid_se+sigmoid_dist_center1
  
  
  
  K=cov_mat_range_one_center_norange(xtr,xtr,
                                     sigma1,sigma2,sigma3,sigma4,l1,l2,
                                     sigma5,sigma6,sigma7,sigma8,l3,l4
                                     ,l_range1)+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  
  K1_sigma5=K1_sigma6=K1_sigma7=K1_sigma8=K1_l3=K1_l4=
    matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
           nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  cov_fund1 = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                     nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  
  
  
  for(i in 1:n1){
    #if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                sum((x1[i,])^2)^.5)
      r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                sum((x2[j,])^2)^.5)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))
      
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma5,sigma6,sigma7,sigma8,l3,l4)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund1[i, j] = results[1,1]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund1[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      
      results_sigma5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma7=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma8=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      K1_sigma5[i, j] = results_sigma5[1,1]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[2,2]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma5[2,1]
      K1_sigma5[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[1,2]
      
      
      
      K1_sigma6[i, j] = results_sigma6[1,1]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[2,2]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma6[2,1]
      K1_sigma6[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[1,2]
      
      
      
      K1_sigma7[i, j] = results_sigma7[1,1]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[2,2]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma7[2,1]
      K1_sigma7[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[1,2]
      
      
      
      K1_sigma8[i, j] = results_sigma8[1,1]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[2,2]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma8[2,1]
      K1_sigma8[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[1,2]
      
      K1_l3[i, j] = results_l3[1,1]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[2,2]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l3[2,1]
      K1_l3[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[1,2]
      
      K1_l4[i, j] = results_l4[1,1]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[2,2]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l4[2,1]
      K1_l4[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[1,2]
      
      
    }
  }
  
  K1_sigma5=K1_sigma5*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma6=K1_sigma6*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma7=K1_sigma7*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma8=K1_sigma8*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l3=K1_l3*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l4=K1_l4*(sigmoid_dist_center1)*1/(weight_sum)
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    if(n2==1){
      dist_mat=sum((x1-x2)^2)^.5
    }else{dist_mat=distmat(t(x1),x2)}
  }else{
    dist_mat=distmat(x1,x2)
  }
  
  
  
  del_sigma1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  
  
  del_sigma1=sigmoid_cov*del_sigma1*1/(weight_sum)
  del_sigma2=sigmoid_cov*del_sigma2*1/(weight_sum)
  del_sigma3=sigmoid_cov*del_sigma3*1/(weight_sum)
  del_sigma4=sigmoid_cov*del_sigma4*1/(weight_sum)
  
  del_l1=sigmoid_cov*del_l1*1/(weight_sum)
  del_l2=sigmoid_cov*del_l2*1/(weight_sum)
  
  K_fund1=cov_fund1
  
  
  K_cov=cov_mat(x1,x1,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  
  
  
  del_l_range1=(K_fund1*del_l_range1_1+K_cov*del_l_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se)/(weight_sum^2)*
    (del_l_range1_2+del_l_range1_1)
  
  
  
  
  
  inv_K=inv(K)
  b=1
  alpha=1
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l2)),
            
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma5)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma5)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma6)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma6)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma7)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma7)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma8)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma8)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l4)),
            
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range1)),
            
            -t(c(ytr[,1],ytr[,2]))%*%(2*sigma_obs_2*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  #print(grad)
  return(grad)
  
}

grad_range_reduced=function(xtr,ytr,
                            sigma1,sigma2,sigma3,sigma4,l1,l2,
                            sigma5,sigma6,sigma7,sigma8,l3,l4,
                            range1,range2,l_range1,l_range2,sigma_obs_2){
  
  
  
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  del_l_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                    (distmat_center1_x1-range1)/l_range1^2)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)*
        (distmat_center1_x2-range1)/l_range1^2)
  
  
  del_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)/l_range1)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)/l_range1)
  
  
  del_l_range1_1=cbind(rbind(del_l_range1_1,del_l_range1_1),
                       rbind(del_l_range1_1,del_l_range1_1))
  
  
  del_range1_1=cbind(rbind(del_range1_1,del_range1_1),
                     rbind(del_range1_1,del_range1_1))
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  
  sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  sigmoid_dist_center2=(1-sigmoid_21)%*%t(1-sigmoid_22)
  
  
  del_l_range2_1=((sigmoid_21)^2*exp(-(distmat_center2_x1-range2)/l_range2)*
                    (distmat_center2_x1-range2)/l_range2^2)%*%t(1-sigmoid_22)+
    (1-sigmoid_21)%*%t(
      (sigmoid_22)^2*exp(-(distmat_center2_x2-range2)/l_range2)*
        (distmat_center2_x2-range2)/l_range2^2)
  
  del_l_range2_1=cbind(rbind(del_l_range2_1,del_l_range2_1),
                       rbind(del_l_range2_1,del_l_range2_1))
  
  
  del_range2_1=((sigmoid_21)^2*exp(-(distmat_center2_x1-range2)/l_range2)/l_range2)%*%t(1-sigmoid_22)+
    (1-sigmoid_21)%*%t(
      (sigmoid_22)^2*exp(-(distmat_center2_x2-range2)/l_range2)/l_range2)
  
  
  del_range2_1=cbind(rbind(del_range2_1,del_range2_1),
                     rbind(del_range2_1,del_range2_1))
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  
  sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  sigmoid_dist_center2=cbind(rbind(sigmoid_dist_center2,sigmoid_dist_center2),
                             rbind(sigmoid_dist_center2,sigmoid_dist_center2))
  
  sigmoid_se=(sigmoid_11*sigmoid_21%*%t(sigmoid_12*sigmoid_22))
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  
  
  
  del_l_range1_2=((sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                     (-distmat_center1_x1+range1)/(l_range1^2))*sigmoid_21)%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t((sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)*
                                   (-distmat_center1_x2+range1)/(l_range1^2))*sigmoid_22)
  
  
  del_l_range1_2=cbind(rbind(del_l_range1_2,del_l_range1_2),
                       rbind(del_l_range1_2,del_l_range1_2))
  
  del_range1_2=((-sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)/l_range1)*sigmoid_21)%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t((-sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)/l_range1)*sigmoid_22)
  
  
  del_range1_2=cbind(rbind(del_range1_2,del_range1_2),
                     rbind(del_range1_2,del_range1_2))
  
  
  
  
  del_l_range2_2=(sigmoid_11*(sigmoid_21^2*exp(-(distmat_center2_x1-range2)/l_range2)*
                                (-distmat_center2_x1+range2)/l_range2^2))%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t(sigmoid_12*(sigmoid_22^2*exp(-(distmat_center2_x2-range2)/l_range2)*
                                              (-distmat_center2_x2+range2)/l_range2^2))
  
  del_l_range2_2=cbind(rbind(del_l_range2_2,del_l_range2_2),
                       rbind(del_l_range2_2,del_l_range2_2))
  
  del_range2_2=(sigmoid_11*(-sigmoid_21^2*exp(-(distmat_center2_x1-range2)/l_range2)/l_range2))%*%
    t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t(sigmoid_12*(-sigmoid_22^2*exp(-(distmat_center2_x2-range2)/l_range2)/l_range2))
  
  del_range2_2=cbind(rbind(del_range2_2,del_range2_2),
                     rbind(del_range2_2,del_range2_2))
  
  
  sigmoid_cov=sigmoid_se
  weight_sum=sigmoid_se+sigmoid_dist_center1+sigmoid_dist_center2
  
  K1_sigma5=K1_sigma6=K1_sigma7=K1_sigma8=K1_l3=K1_l4=
    matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
           nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  
  K=cov_mat_range_reduced(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                          sigma5,sigma6,sigma7,sigma8,l3,l4,
                          range1,range2,l_range1,l_range2)+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  cov_fund1 = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                     nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  cov_fund2=cov_fund1
  
  
  for(i in 1:n1){
    #if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                sum((x1[i,])^2)^.5)
      r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                sum((x2[j,])^2)^.5)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))
      
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma5,sigma6,sigma7,sigma8,l3,l4)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund1[i, j] = results[1,1]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund1[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      
      
      results_sigma5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma7=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma8=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma5,sigma6,sigma7,sigma8,l3,l4)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      #####
      
      
      K1_sigma5[i, j] = results_sigma5[1,1]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[2,2]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma5[2,1]
      K1_sigma5[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[1,2]
      
      
      
      K1_sigma6[i, j] = results_sigma6[1,1]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[2,2]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma6[2,1]
      K1_sigma6[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[1,2]
      
      
      
      K1_sigma7[i, j] = results_sigma7[1,1]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[2,2]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma7[2,1]
      K1_sigma7[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[1,2]
      
      
      
      K1_sigma8[i, j] = results_sigma8[1,1]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[2,2]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma8[2,1]
      K1_sigma8[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[1,2]
      
      K1_l3[i, j] = results_l3[1,1]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[2,2]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l3[2,1]
      K1_l3[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[1,2]
      
      K1_l4[i, j] = results_l4[1,1]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[2,2]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l4[2,1]
      K1_l4[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[1,2]
      
      
    }
  }
  
  K1_sigma5=K1_sigma5*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  K1_sigma6=K1_sigma6*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  K1_sigma7=K1_sigma7*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  K1_sigma8=K1_sigma8*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  K1_l3=K1_l3*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  K1_l4=K1_l4*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    if(n2==1){
      dist_mat=sum((x1-x2)^2)^.5
    }else{dist_mat=dist(t(x1),x2)}
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  
  del_sigma1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  
  
  del_sigma1=sigmoid_cov*del_sigma1*1/(weight_sum)
  del_sigma2=sigmoid_cov*del_sigma2*1/(weight_sum)
  del_sigma3=sigmoid_cov*del_sigma3*1/(weight_sum)
  del_sigma4=sigmoid_cov*del_sigma4*1/(weight_sum)
  
  del_l1=sigmoid_cov*del_l1*1/(weight_sum)
  del_l2=sigmoid_cov*del_l2*1/(weight_sum)
  
  K_fund1=K_fund2=cov_fund1
  #K_fund2=cov_fund2
  
  K_cov=cov_mat(x1,x1,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  
  del_l_range1=(K_fund1*del_l_range1_1+K_cov*del_l_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_l_range1_2+del_l_range1_1)
  
  del_l_range2=(K_fund2*del_l_range2_1+K_cov*del_l_range2_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_l_range2_2+del_l_range2_1)
  
  del_range1=(K_fund1*(del_range1_1)+K_cov*del_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_range1_2+del_range1_1)
  
  
  del_range2=(K_fund2*(del_range2_1)+K_cov*del_range2_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_range2_2+del_range2_1)
  
  
  inv_K=inv(K)
  b=1
  alpha=1
  
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma5)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma5)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma6)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*K1_sigma6)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma7)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma7)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_sigma8)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_sigma8)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(K1_l4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*K1_l4)),
            
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range2)),
            
            -t(c(ytr[,1],ytr[,2]))%*%(2*sigma_obs_2*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  #print(grad)
  return(grad)
  
}


grad_range_reduced2=function(xtr,ytr,
                             sigma1,sigma2,sigma3,sigma4,l1,l2,
                             range1,range2,l_range1,l_range2,sigma_obs_2){
  
  
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(3,n1),0))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(3,n2),0))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  del_l_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                    (distmat_center1_x1-range1)/l_range1^2)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)*
        (distmat_center1_x2-range1)/l_range1^2)
  
  
  del_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)/l_range1)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)/l_range1)
  
  
  del_l_range1_1=cbind(rbind(del_l_range1_1,del_l_range1_1),
                       rbind(del_l_range1_1,del_l_range1_1))
  
  
  del_range1_1=cbind(rbind(del_range1_1,del_range1_1),
                     rbind(del_range1_1,del_range1_1))
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(-3,n1),0))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(-3,n2),0))^2,1,sum)^.5)
  
  sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  sigmoid_dist_center2=(1-sigmoid_21)%*%t(1-sigmoid_22)
  
  
  del_l_range2_1=((sigmoid_21)^2*exp(-(distmat_center2_x1-range2)/l_range2)*
                    (distmat_center2_x1-range2)/l_range2^2)%*%t(1-sigmoid_22)+
    (1-sigmoid_21)%*%t(
      (sigmoid_22)^2*exp(-(distmat_center2_x2-range2)/l_range2)*
        (distmat_center2_x2-range2)/l_range2^2)
  
  del_l_range2_1=cbind(rbind(del_l_range2_1,del_l_range2_1),
                       rbind(del_l_range2_1,del_l_range2_1))
  
  
  del_range2_1=((sigmoid_21)^2*exp(-(distmat_center2_x1-range2)/l_range2)/l_range2)%*%t(1-sigmoid_22)+
    (1-sigmoid_21)%*%t(
      (sigmoid_22)^2*exp(-(distmat_center2_x2-range2)/l_range2)/l_range2)
  
  
  del_range2_1=cbind(rbind(del_range2_1,del_range2_1),
                     rbind(del_range2_1,del_range2_1))
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  
  sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  sigmoid_dist_center2=cbind(rbind(sigmoid_dist_center2,sigmoid_dist_center2),
                             rbind(sigmoid_dist_center2,sigmoid_dist_center2))
  
  sigmoid_se=(sigmoid_11*sigmoid_21%*%t(sigmoid_12*sigmoid_22))
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  
  
  
  del_l_range1_2=((sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                     (-distmat_center1_x1+range1)/(l_range1^2))*sigmoid_21)%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t((sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)*
                                   (-distmat_center1_x2+range1)/(l_range1^2))*sigmoid_22)
  
  
  del_l_range1_2=cbind(rbind(del_l_range1_2,del_l_range1_2),
                       rbind(del_l_range1_2,del_l_range1_2))
  
  del_range1_2=((-sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)/l_range1)*sigmoid_21)%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t((-sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)/l_range1)*sigmoid_22)
  
  
  del_range1_2=cbind(rbind(del_range1_2,del_range1_2),
                     rbind(del_range1_2,del_range1_2))
  
  
  
  
  del_l_range2_2=(sigmoid_11*(sigmoid_21^2*exp(-(distmat_center2_x1-range2)/l_range2)*
                                (-distmat_center2_x1+range2)/l_range2^2))%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t(sigmoid_12*(sigmoid_22^2*exp(-(distmat_center2_x2-range2)/l_range2)*
                                              (-distmat_center2_x2+range2)/l_range2^2))
  
  del_l_range2_2=cbind(rbind(del_l_range2_2,del_l_range2_2),
                       rbind(del_l_range2_2,del_l_range2_2))
  
  del_range2_2=(sigmoid_11*(-sigmoid_21^2*exp(-(distmat_center2_x1-range2)/l_range2)/l_range2))%*%
    t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t(sigmoid_12*(-sigmoid_22^2*exp(-(distmat_center2_x2-range2)/l_range2)/l_range2))
  
  del_range2_2=cbind(rbind(del_range2_2,del_range2_2),
                     rbind(del_range2_2,del_range2_2))
  
  
  sigmoid_cov=sigmoid_se
  weight_sum=sigmoid_se+sigmoid_dist_center1+sigmoid_dist_center2
  
  K1_sigma5=K1_sigma6=K1_sigma7=K1_sigma8=K1_l3=K1_l4=
    matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
           nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  
  K=cov_mat_range_reduced2(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2,
                           range1,range2,l_range1,l_range2)+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  cov_fund1 = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                     nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  cov_fund2=cov_fund1
  
  
  for(i in 1:n1){
    #if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                sum((x1[i,])^2)^.5)
      r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                sum((x2[j,])^2)^.5)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))
      
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma1,sigma2,sigma3,sigma4,l1,l2)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund1[i, j] = results[1,1]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund1[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      
      
      results_sigma5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma7=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma8=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      #####
      
      
      K1_sigma5[i, j] = results_sigma5[1,1]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[2,2]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma5[2,1]
      K1_sigma5[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[1,2]
      
      
      
      K1_sigma6[i, j] = results_sigma6[1,1]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[2,2]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma6[2,1]
      K1_sigma6[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[1,2]
      
      
      
      K1_sigma7[i, j] = results_sigma7[1,1]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[2,2]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma7[2,1]
      K1_sigma7[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[1,2]
      
      
      
      K1_sigma8[i, j] = results_sigma8[1,1]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[2,2]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma8[2,1]
      K1_sigma8[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[1,2]
      
      K1_l3[i, j] = results_l3[1,1]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[2,2]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l3[2,1]
      K1_l3[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[1,2]
      
      K1_l4[i, j] = results_l4[1,1]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[2,2]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l4[2,1]
      K1_l4[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[1,2]
      
      
    }
  }
  
  K1_sigma5=K1_sigma5*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  K1_sigma6=K1_sigma6*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  K1_sigma7=K1_sigma7*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  K1_sigma8=K1_sigma8*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  K1_l3=K1_l3*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  K1_l4=K1_l4*(sigmoid_dist_center1+sigmoid_dist_center2)*1/(weight_sum)
  
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    if(n2==1){
      dist_mat=sum((x1-x2)^2)^.5
    }else{dist_mat=dist(t(x1),x2)}
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  
  del_sigma1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  
  
  del_sigma1=sigmoid_cov*del_sigma1*1/(weight_sum)+K1_sigma5
  del_sigma2=sigmoid_cov*del_sigma2*1/(weight_sum)+K1_sigma6
  del_sigma3=sigmoid_cov*del_sigma3*1/(weight_sum)+K1_sigma7
  del_sigma4=sigmoid_cov*del_sigma4*1/(weight_sum)+K1_sigma8
  
  del_l1=sigmoid_cov*del_l1*1/(weight_sum)+K1_l3
  del_l2=sigmoid_cov*del_l2*1/(weight_sum)+K1_l4
  
  K_fund1=K_fund2=cov_fund1
  #K_fund2=cov_fund2
  
  K_cov=cov_mat(x1,x1,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  
  del_l_range1=(K_fund1*del_l_range1_1+K_cov*del_l_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_l_range1_2+del_l_range1_1)
  
  del_l_range2=(K_fund2*del_l_range2_1+K_cov*del_l_range2_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_l_range2_2+del_l_range2_1)
  
  del_range1=(K_fund1*(del_range1_1)+K_cov*del_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_range1_2+del_range1_1)
  
  
  del_range2=(K_fund2*(del_range2_1)+K_cov*del_range2_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_range2_2+del_range2_1)
  
  
  inv_K=inv(K)
  b=1
  alpha=1
  
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l2)),
            
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range2)),
            
            -t(c(ytr[,1],ytr[,2]))%*%(2*sigma_obs_2*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  #print(grad)
  return(grad)
  
}

grad_range_one_center_reduced=function(xtr,ytr,
                                       sigma1,sigma2,sigma3,sigma4,l1,l2,
                                       range1,l_range1,sigma_obs_2){
  
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(-87.1,n1),25.3))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(-87.1,n2),25.3))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  del_l_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                    (distmat_center1_x1-range1)/l_range1^2)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)*
        (distmat_center1_x2-range1)/l_range1^2)
  
  
  del_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)/l_range1)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)/l_range1)
  
  
  del_l_range1_1=cbind(rbind(del_l_range1_1,del_l_range1_1),
                       rbind(del_l_range1_1,del_l_range1_1))
  
  
  del_range1_1=cbind(rbind(del_range1_1,del_range1_1),
                     rbind(del_range1_1,del_range1_1))
  
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  
  sigmoid_se=(sigmoid_11%*%t(sigmoid_12))
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  
  
  
  del_l_range1_2=((sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                     (-distmat_center1_x1+range1)/(l_range1^2)))%*%t(sigmoid_12)+
    (sigmoid_11)%*%t((sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)*
                        (-distmat_center1_x2+range1)/(l_range1^2)))
  
  
  del_l_range1_2=cbind(rbind(del_l_range1_2,del_l_range1_2),
                       rbind(del_l_range1_2,del_l_range1_2))
  
  del_range1_2=((-sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)/l_range1))%*%t(sigmoid_12)+
    (sigmoid_11)%*%t((-sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)/l_range1))
  
  
  del_range1_2=cbind(rbind(del_range1_2,del_range1_2),
                     rbind(del_range1_2,del_range1_2))
  
  
  
  
  
  sigmoid_cov=sigmoid_se
  weight_sum=sigmoid_se+sigmoid_dist_center1
  
  
  
  
  K=cov_mat_range_one_center_reduced(xtr,xtr,
                                     sigma1,sigma2,sigma3,sigma4,l1,l2
                                     
                                     ,range1,l_range1)+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  
  K1_sigma5=K1_sigma6=K1_sigma7=K1_sigma8=K1_l3=K1_l4=
    matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
           nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  cov_fund1 = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                     nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  
  
  
  for(i in 1:n1){
    #if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                sum((x1[i,])^2)^.5)
      r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                sum((x2[j,])^2)^.5)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))
      
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma1,sigma2,sigma3,sigma4,l1,l2)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund1[i, j] = results[1,1]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund1[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      
      results_sigma5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma7=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma8=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      K1_sigma5[i, j] = results_sigma5[1,1]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[2,2]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma5[2,1]
      K1_sigma5[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[1,2]
      
      
      
      K1_sigma6[i, j] = results_sigma6[1,1]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[2,2]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma6[2,1]
      K1_sigma6[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[1,2]
      
      
      
      K1_sigma7[i, j] = results_sigma7[1,1]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[2,2]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma7[2,1]
      K1_sigma7[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[1,2]
      
      
      
      K1_sigma8[i, j] = results_sigma8[1,1]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[2,2]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma8[2,1]
      K1_sigma8[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[1,2]
      
      K1_l3[i, j] = results_l3[1,1]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[2,2]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l3[2,1]
      K1_l3[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[1,2]
      
      K1_l4[i, j] = results_l4[1,1]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[2,2]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l4[2,1]
      K1_l4[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[1,2]
      
      
    }
  }
  
  K1_sigma5=K1_sigma5*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma6=K1_sigma6*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma7=K1_sigma7*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma8=K1_sigma8*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l3=K1_l3*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l4=K1_l4*(sigmoid_dist_center1)*1/(weight_sum)
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    if(n2==1){
      dist_mat=sum((x1-x2)^2)^.5
    }else{dist_mat=distmat(t(x1),x2)}
  }else{
    dist_mat=distmat(x1,x2)
  }
  
  
  
  del_sigma1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  
  
  del_sigma1=sigmoid_cov*del_sigma1*1/(weight_sum)+K1_sigma5
  del_sigma2=sigmoid_cov*del_sigma2*1/(weight_sum)+K1_sigma6
  del_sigma3=sigmoid_cov*del_sigma3*1/(weight_sum)+K1_sigma7
  del_sigma4=sigmoid_cov*del_sigma4*1/(weight_sum)+K1_sigma8
  
  del_l1=sigmoid_cov*del_l1*1/(weight_sum)+K1_l3
  del_l2=sigmoid_cov*del_l2*1/(weight_sum)+K1_l4
  
  K_fund1=cov_fund1
  
  
  K_cov=cov_mat(x1,x1,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  K_fund1=cov_fund1
  
  
  
  
  del_l_range1=(K_fund1*del_l_range1_1+K_cov*del_l_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se)/(weight_sum^2)*
    (del_l_range1_2+del_l_range1_1)
  
  
  del_range1=(K_fund1*(del_range1_1)+K_cov*del_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se)/(weight_sum^2)*
    (del_range1_2+del_range1_1)
  
  
  
  inv_K=inv(K)
  b=1
  alpha=1
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_range1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range1)),
            
            -t(c(ytr[,1],ytr[,2]))%*%(2*sigma_obs_2*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  #print(grad)
  return(grad)
  
}
grad_range_one_center_norange_reduced=function(xtr,ytr,
                                               sigma1,sigma2,sigma3,sigma4,l1,l2,
                                               
                                               l_range1,sigma_obs_2){
  
  
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(-87.1,n1),25.3))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(-87.1,n2),25.3))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  del_l_range1_1=((sigmoid_11)^2*exp(-(distmat_center1_x1)/l_range1)*
                    (distmat_center1_x1)/l_range1^2)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2)/l_range1)*
        (distmat_center1_x2)/l_range1^2)
  
  
  
  
  
  del_l_range1_1=cbind(rbind(del_l_range1_1,del_l_range1_1),
                       rbind(del_l_range1_1,del_l_range1_1))
  
  
  
  
  
  sigmoid_dist_center1=cbind(rbind(sigmoid_dist_center1,sigmoid_dist_center1),
                             rbind(sigmoid_dist_center1,sigmoid_dist_center1))
  
  
  sigmoid_se=(sigmoid_11%*%t(sigmoid_12))
  sigmoid_se=cbind(rbind(sigmoid_se,sigmoid_se),rbind(sigmoid_se,sigmoid_se))
  
  
  
  
  del_l_range1_2=((sigmoid_11^2*exp(-(distmat_center1_x1)/l_range1)*
                     (-distmat_center1_x1)/(l_range1^2)))%*%t(sigmoid_12)+
    (sigmoid_11)%*%t((sigmoid_12^2*exp(-(distmat_center1_x2)/l_range1)*
                        (-distmat_center1_x2)/(l_range1^2)))
  
  
  del_l_range1_2=cbind(rbind(del_l_range1_2,del_l_range1_2),
                       rbind(del_l_range1_2,del_l_range1_2))
  
  
  
  
  
  
  sigmoid_cov=sigmoid_se
  weight_sum=sigmoid_se+sigmoid_dist_center1
  
  
  
  K=cov_mat_range_one_center_norange_reduced(xtr,xtr,
                                             sigma1,sigma2,sigma3,sigma4,l1,l2,
                                             
                                             l_range1)+
    diag(sigma_obs_2^2,nrow=2*nrow(xtr))
  
  
  K1_sigma5=K1_sigma6=K1_sigma7=K1_sigma8=K1_l3=K1_l4=
    matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
           nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  cov_fund1 = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                     nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  
  
  
  for(i in 1:n1){
    #if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:n2){
      r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                sum((x1[i,])^2)^.5)
      r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                sum((x2[j,])^2)^.5)
      #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
      theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                    atan2(x2[j,2],x2[j,1]))
      
      cov_standard=cov_mat(c((r1),0),c((r2),0),
                           sigma1,sigma2,sigma3,sigma4,l1,l2)
      
      results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        cov_standard%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      cov_fund1[i, j] = results[1,1]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
                ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
      cov_fund1[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
      cov_fund1[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
      
      
      results_sigma5=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig1(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_sigma6=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig2(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma7=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig3(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_sigma8=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_sig4(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results_l3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_1(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      results_l4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        del_l_2(c(r1,0),c(r2,0),sigma1,sigma2,sigma3,sigma4,l1,l2)%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      K1_sigma5[i, j] = results_sigma5[1,1]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[2,2]
      K1_sigma5[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma5[2,1]
      K1_sigma5[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma5[1,2]
      
      
      
      K1_sigma6[i, j] = results_sigma6[1,1]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[2,2]
      K1_sigma6[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma6[2,1]
      K1_sigma6[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma6[1,2]
      
      
      
      K1_sigma7[i, j] = results_sigma7[1,1]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[2,2]
      K1_sigma7[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma7[2,1]
      K1_sigma7[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma7[1,2]
      
      
      
      K1_sigma8[i, j] = results_sigma8[1,1]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[2,2]
      K1_sigma8[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_sigma8[2,1]
      K1_sigma8[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_sigma8[1,2]
      
      K1_l3[i, j] = results_l3[1,1]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[2,2]
      K1_l3[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l3[2,1]
      K1_l3[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l3[1,2]
      
      K1_l4[i, j] = results_l4[1,1]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
            ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[2,2]
      K1_l4[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results_l4[2,1]
      K1_l4[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results_l4[1,2]
      
      
    }
  }
  
  K1_sigma5=K1_sigma5*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma6=K1_sigma6*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma7=K1_sigma7*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma8=K1_sigma8*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l3=K1_l3*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l4=K1_l4*(sigmoid_dist_center1)*1/(weight_sum)
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    if(n2==1){
      dist_mat=sum((x1-x2)^2)^.5
    }else{dist_mat=distmat(t(x1),x2)}
  }else{
    dist_mat=distmat(x1,x2)
  }
  
  
  
  del_sigma1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma1[1:n1,1:n2]=2*sigma1*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*exp(-.5*dist_mat^2*(l1^2))
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=2*sigma2*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*exp(-.5*dist_mat^2*(l2^2))
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    exp(-.5*dist_mat^2*(l1^2))
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*exp(-.5*dist_mat^2*(l1^2))
  
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    exp(-.5*dist_mat^2*(l2^2))
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  
  
  del_l1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l1[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),1:n2]=del_l1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*
    exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  del_l1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))*(-l1*dist_mat^2)
  
  del_l2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_l2[1:n1,1:n2]=sigma2^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),1:n2]=del_l2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  del_l2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*exp(-.5*dist_mat^2*(l2^2))*(-l2*dist_mat^2)
  
  
  
  del_sigma1=sigmoid_cov*del_sigma1*1/(weight_sum)+K1_sigma5
  del_sigma2=sigmoid_cov*del_sigma2*1/(weight_sum)+K1_sigma6
  del_sigma3=sigmoid_cov*del_sigma3*1/(weight_sum)+K1_sigma7
  del_sigma4=sigmoid_cov*del_sigma4*1/(weight_sum)+K1_sigma8
  
  del_l1=sigmoid_cov*del_l1*1/(weight_sum)+K1_l3
  del_l2=sigmoid_cov*del_l2*1/(weight_sum)++K1_l4
  
  K_fund1=cov_fund1
  
  
  K_cov=cov_mat(x1,x1,sigma1,sigma2,sigma3,sigma4,l1,l2)
  
  
  
  
  del_l_range1=(K_fund1*del_l_range1_1+K_cov*del_l_range1_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se)/(weight_sum^2)*
    (del_l_range1_2+del_l_range1_1)
  
  
  
  
  
  inv_K=inv(K)
  b=1
  alpha=1
  
  grad=.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma1)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma2)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%(alpha*del_sigma2)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma3)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma3)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_sigma4)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_sigma4)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l1)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l1)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%(del_l2)%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l2)),
            
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_l_range1%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_l_range1)),
            
            -t(c(ytr[,1],ytr[,2]))%*%(2*sigma_obs_2*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  #print(grad)
  return(grad)
  
}







# Set up a for loop with a progress bar

Adam = function(f, grad_f, params_init, lr = .01, beta1 = 0.9, beta2 = 0.999,
                epsilon = 1e-8, max_iter = 100, tol = 1e-10,
                lb=rep(0,length(params_init)),
                ub=rep(Inf,length(params_init))) {
  
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



if(data!="Gulf"){
  RMSES=LOGS=matrix(nrow=N_train,ncol=7)
  Xtr=c()
  for( i in 1:N_train){
    
    size = 30
    
    x1 = seq(-4, 4, length.out = size)
    x2 = seq(-4, 4, length.out = size)
    L = expand.grid(x1, x2)
    
    F_tilde = function(x) {
      c(x[2], (x[1] - 0.1 * (x[1])^3) * (1 + 0.1 * cos(50 * pi * x[1])))
    }
    
    D = function(x, bd, Rd) {
      c((x[1] - (-3)) * bd / Rd, (x[2] - 0) * bd / Rd)
    }
    
    C = function(x, bc, Rc) {
      c(-(x[1] - 3) * bc / Rc, -(x[2] - 0) * bc / Rc)
    }
    
    
    bd = .5 
    bc = .5  
    
    Yte= t(apply(L, 1, function(x) {
      Rd = sqrt((x[1] - (-3))^2 + (x[2] - 0)^2)
      Rc = sqrt((x[1] - 3)^2 + (x[2] - 0)^2)
      
      F_tilde_val = F_tilde(x)
      D_val = D(x, bd, Rd)
      C_val = C(x, bc, Rc)
      
      F_tilde_val + D_val + C_val
    }))
    
    
    Xtr=rbind(Xtr,matrix(runif(10,min=-4,max=4),ncol=2))
    
    Ytr=noise+t(apply(Xtr, 1, function(x) {
      Rd = sqrt((x[1] - (-3))^2 + (x[2] - 0)^2)
      Rc = sqrt((x[1] - 3)^2 + (x[2] - 0)^2)
      
      F_tilde_val = F_tilde(x)
      D_val = D(x, bd, Rd)
      C_val = C(x, bc, Rc)
      
      F_tilde_val + D_val + C_val
    }))
    
    Xte=L
    
    
    opt_par=Adam(function(x) {
      log_likelihood(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7])
    },function(x) {
      grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],x[6], x[7])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_sigma_obs,
        initial_l1,
        initial_l2),lb=rep(0.001,7),lr=lr,max_iter=max_iter)
    
    Ktetr_standard=cov_mat(Xte,Xtr,opt_par[1],
                           opt_par[2],
                           opt_par[3],
                           opt_par[4],
                           opt_par[6],
                           opt_par[7])
    
    Ktetr_trtr_inv_standard=Ktetr_standard%*%inv(cov_mat(Xtr,Xtr,opt_par[1],
                                                         opt_par[2],
                                                         opt_par[3],
                                                         opt_par[4],
                                                         opt_par[6],
                                                         opt_par[7])
                                                 
                                                 +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_standard=Ktetr_trtr_inv_standard%*%c(Ytr[,1],Ytr[,2])
    
    
    
    (RMSES[i,1]=rmse_standard=mean(apply((cbind(
      posterior_mean_standard[1:(0.5*length(posterior_mean_standard))],
      posterior_mean_standard
      [(0.5*length(posterior_mean_standard)+1):(length(posterior_mean_standard))])
      -Yte)^2,1,sum))^.5)
    
    
    
    posterior_cov_standard= cov_mat(Xte,Xte,opt_par[1],
                                    opt_par[2],
                                    opt_par[3],
                                    opt_par[4],
                                    opt_par[6],
                                    opt_par[7])-Ktetr_trtr_inv_standard%*%t(Ktetr_standard)
    
    
    (LOGS[i,1]=LogS_standard=(t(c(Yte[,1],Yte[,2])-posterior_mean_standard)%*%
                                inv(posterior_cov_standard+diag(opt_par[5]^2,2*nrow(Xte)))%*%
                                (c(Yte[,1],Yte[,2])-posterior_mean_standard)
                              +determinant(posterior_cov_standard+diag(opt_par[5]^2,2*nrow(Xte))
                              )$modulus[1]+
                                log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    
    
    opt_par=Adam(function(x) {log_likelihood_helm(
      Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
    },function(x) {
      grad_helm(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
    },c(initial_par[1:4],initial_sigma_obs),lb=rep(0.001,5),lr=lr,max_iter=max_iter)
    
    
    Ktetr_helm=cov_mat_helm(Xte,Xtr,opt_par[1],
                            opt_par[2],
                            opt_par[3],
                            opt_par[4])
    
    Ktetr_trtr_inv_helm=Ktetr_helm%*%inv(cov_mat_helm(Xtr,Xtr,opt_par[1],
                                                      opt_par[2],
                                                      opt_par[3],
                                                      opt_par[4])
                                         
                                         +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_helm=Ktetr_trtr_inv_helm%*%c(Ytr[,1],Ytr[,2])
    
    
    
    (RMSES[i,2]=rmse_helm=mean(apply((cbind(posterior_mean_helm[1:(0.5*length(posterior_mean_helm))],
                                            posterior_mean_helm[(0.5*length(posterior_mean_helm)+1):
                                                                  (length(posterior_mean_helm))])
                                      -Yte)^2,1,sum))^.5)
    
    posterior_cov_helm= cov_mat_helm(Xte,Xte,opt_par[1],
                                     opt_par[2],
                                     opt_par[3],
                                     opt_par[4])-Ktetr_trtr_inv_helm%*%t(Ktetr_helm)
    
    
    (LOGS[i,2]=LogS_helm=(t(c(Yte[,1],Yte[,2])-posterior_mean_helm)%*%
                            inv(posterior_cov_helm+diag(opt_par[5]^2,2*nrow(Xte)))%*%
                            (c(Yte[,1],Yte[,2])-posterior_mean_helm)
                          +determinant(posterior_cov_helm+diag(opt_par[5]^2,2*nrow(Xte)))$modulus[1]+
                            log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    
    
    opt_par=Adam(function(x) {
      log_likelihood_mixture(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6],x[7],x[8])
    },function(x) {
      grad_mix(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6],x[7],x[8])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_sigma_obs,
        initial_l1,
        initial_l2,
        initial_alpha),lb=c(rep(0.001,7),.2),ub=c(rep(Inf,7),.8),
    lr=lr,max_iter=max_iter)
    
    
    
    Ktetr_mix=mix_cov_mat(Xte,Xtr,opt_par[1],
                          opt_par[2],
                          opt_par[3],
                          opt_par[4],
                          opt_par[6],
                          opt_par[7],
                          opt_par[8])
    
    Ktetr_trtr_inv_mix=Ktetr_mix%*%inv(mix_cov_mat(Xtr,Xtr,opt_par[1],
                                                   opt_par[2],
                                                   opt_par[3],
                                                   opt_par[4],
                                                   opt_par[6],
                                                   opt_par[7],
                                                   opt_par[8])
                                       
                                       +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])
    
    
    
    (RMSES[i,3]=rmse_mix=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                                           posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                                (length(posterior_mean_mix))])
                                     -Yte)^2,1,sum))^.5)
    
    posterior_cov_mix= mix_cov_mat(Xte,Xte,opt_par[1],
                                   opt_par[2],
                                   opt_par[3],
                                   opt_par[4],opt_par[6],
                                   opt_par[7],opt_par[8])-Ktetr_trtr_inv_mix%*%
      t(Ktetr_mix) 
    
    (LOGS[i,3]=LogS_mix=(t(c(Yte[,1],Yte[,2])-posterior_mean_mix)%*%
                           inv(posterior_cov_mix+diag(opt_par[5]^2,2*nrow(Xte)))%*%
                           (c(Yte[,1],Yte[,2])-posterior_mean_mix)
                         +determinant(posterior_cov_mix+diag(opt_par[5]^2,2*nrow(Xte)))$modulus[1]+
                           log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    
    opt_par=Adam(function(x) {
      log_likelihood_range(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7], x[8], x[9], x[10],
        x[11], x[12], x[13],x[14], x[15], x[16], x[17],
        x[18], x[19], x[20], x[21],x[22],x[23])
    },function(x) {
      grad_range(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                 x[6], x[7], x[8], x[9], x[10],
                 x[11], x[12], x[13],x[14], x[15], x[16], x[17],
                 x[18], x[19], x[20], x[21],x[22],x[23])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_l1,
        initial_l2,
        initial_sigma5,
        initial_sigma6,
        initial_sigma7,
        initial_sigma8,
        initial_l3,
        initial_l4,
        initial_sigma9,
        initial_sigma10,
        initial_sigma11,
        initial_sigma12,
        initial_l5,
        initial_l6,
        initial_range1,
        initial_l_range1,
        initial_range2,
        initial_l_range2,
        initial_sigma_obs),lr=lr,max_iter=max_iter)
    
    Ktetr_range=cov_mat_range(Xte,Xtr,opt_par[1],
                              opt_par[2],
                              opt_par[3],
                              opt_par[4],
                              opt_par[5],
                              opt_par[6],
                              opt_par[7],
                              opt_par[8],
                              opt_par[9],
                              opt_par[10],
                              opt_par[11],
                              opt_par[12],
                              opt_par[13],
                              opt_par[14],
                              opt_par[15],
                              opt_par[16],
                              opt_par[17],
                              opt_par[18],
                              opt_par[19],
                              opt_par[20],
                              opt_par[21],
                              opt_par[22])
    
    Ktetr_trtr_inv_range=Ktetr_range%*%inv(cov_mat_range(Xtr,Xtr,opt_par[1],
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
                                                         opt_par[15],opt_par[16],
                                                         opt_par[17],
                                                         opt_par[18],
                                                         opt_par[19],
                                                         opt_par[20],
                                                         opt_par[21],
                                                         opt_par[22])
                                           
                                           +diag(opt_par[23]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_range=Ktetr_trtr_inv_range%*%c(Ytr[,1],Ytr[,2])
    
    
    
    
    (RMSES[i,4]=rmse_range4=mean(apply(
      (cbind(posterior_mean_range[1:(0.5*length(posterior_mean_range))],
             posterior_mean_range[
               (0.5*length(posterior_mean_range)+1):
                 (length(posterior_mean_range))])
       -Yte)^2,1,sum))^.5)
    
    
    
    posterior_cov_range= cov_mat_range(Xte,Xte,opt_par[1],
                                       opt_par[2],
                                       opt_par[3],
                                       opt_par[4],opt_par[5],
                                       opt_par[6],
                                       opt_par[7],opt_par[8],
                                       opt_par[9],
                                       opt_par[10],
                                       opt_par[11],opt_par[12],
                                       opt_par[13],opt_par[14],
                                       opt_par[15],opt_par[16],
                                       opt_par[17],
                                       opt_par[18],
                                       opt_par[19],
                                       opt_par[20],
                                       opt_par[21],
                                       opt_par[22])-
      Ktetr_trtr_inv_range%*%t(Ktetr_range)
    
    
    (LOGS[i,4]=LogS_range4=(t(c(Yte[,1],Yte[,2])-posterior_mean_range)%*%
                              inv(posterior_cov_range+diag(opt_par[23]^2,2*nrow(Xte)))%*%
                              (c(Yte[,1],Yte[,2])-posterior_mean_range)
                            +determinant(posterior_cov_range+diag(opt_par[23]^2,2*nrow(Xte)))$modulus[1]+
                              log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    opt_par=Adam(function(x) {
      log_likelihood_range_reduced(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7], x[8], x[9], x[10],
        x[11], x[12], x[13],x[14], x[15], x[16], x[17])
    },function(x) {
      grad_range_reduced(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                         x[6], x[7], x[8], x[9], x[10],
                         x[11], x[12], x[13],x[14], x[15], x[16], x[17])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_l1,
        initial_l2,
        initial_sigma5,
        initial_sigma6,
        initial_sigma7,
        initial_sigma8,
        initial_l3,
        initial_l4,
        initial_range1,
        initial_l_range1,
        initial_range2,
        initial_l_range2,
        initial_sigma_obs),lr=lr,max_iter=max_iter)
    
    Ktetr_range_reduced=cov_mat_range_reduced(Xte,Xtr,opt_par[1],
                                              opt_par[2],
                                              opt_par[3],
                                              opt_par[4],
                                              opt_par[5],
                                              opt_par[6],
                                              opt_par[7],
                                              opt_par[8],
                                              opt_par[9],
                                              opt_par[10],
                                              opt_par[11],
                                              opt_par[12],
                                              opt_par[13],
                                              opt_par[14],
                                              opt_par[15],
                                              opt_par[16])
    
    Ktetr_trtr_inv_range_reduced=Ktetr_range_reduced%*%inv(cov_mat_range_reduced(Xtr,Xtr,opt_par[1],
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
                                                                                 opt_par[15],opt_par[16])
                                                           
                                                           +diag(opt_par[17]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_range_reduced=Ktetr_trtr_inv_range_reduced%*%c(Ytr[,1],Ytr[,2])
    
    
    
    
    (RMSES[i,5]=rmse_range_reduced4=mean(apply(
      (cbind(posterior_mean_range_reduced[1:(0.5*length(posterior_mean_range_reduced))],
             posterior_mean_range_reduced[
               (0.5*length(posterior_mean_range_reduced)+1):
                 (length(posterior_mean_range_reduced))])
       -Yte)^2,1,sum))^.5)
    
    
    
    posterior_cov_range_reduced= cov_mat_range_reduced(Xte,Xte,opt_par[1],
                                                       opt_par[2],
                                                       opt_par[3],
                                                       opt_par[4],opt_par[5],
                                                       opt_par[6],
                                                       opt_par[7],opt_par[8],
                                                       opt_par[9],
                                                       opt_par[10],
                                                       opt_par[11],opt_par[12],
                                                       opt_par[13],opt_par[14],
                                                       opt_par[15],opt_par[16])-
      Ktetr_trtr_inv_range_reduced%*%t(Ktetr_range_reduced)
    
    
    (LOGS[i,5]=LogS_range_reduced4_reduced=(t(c(Yte[,1],Yte[,2])-posterior_mean_range_reduced)%*%
                                              inv(posterior_cov_range_reduced+diag(opt_par[17]^2,2*nrow(Xte)))%*%
                                              (c(Yte[,1],Yte[,2])-posterior_mean_range_reduced)
                                            +determinant(posterior_cov_range_reduced+diag(opt_par[17]^2,2*nrow(Xte)))$modulus[1]+
                                              log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    opt_par=Adam(function(x) {
      log_likelihood_range_reduced2(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7], x[8], x[9], x[10],
        x[11])
    },function(x) {
      grad_range_reduced2(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                          x[6], x[7], x[8], x[9], x[10],
                          x[11])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_l1,
        initial_l2,
        initial_range1,
        initial_l_range1,
        initial_range2,
        initial_l_range2,
        initial_sigma_obs),lr=.1,max_iter = 20)
    
    Ktetr_range_reduced2=cov_mat_range_reduced2(Xte,Xtr,opt_par[1],
                                                opt_par[2],
                                                opt_par[3],
                                                opt_par[4],
                                                opt_par[5],
                                                opt_par[6],
                                                opt_par[7],
                                                opt_par[8],
                                                opt_par[9],
                                                opt_par[10])
    
    Ktetr_trtr_inv_range_reduced2=Ktetr_range_reduced2%*%inv(cov_mat_range_reduced2(Xtr,Xtr,opt_par[1],
                                                                                    opt_par[2],
                                                                                    opt_par[3],
                                                                                    opt_par[4],
                                                                                    opt_par[5],
                                                                                    opt_par[6],
                                                                                    opt_par[7],opt_par[8],
                                                                                    opt_par[9],
                                                                                    opt_par[10])
                                                             
                                                             +diag(opt_par[11]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_range_reduced2=Ktetr_trtr_inv_range_reduced2%*%c(Ytr[,1],Ytr[,2])
    
    
    
    
    (RMSES[i,6]=rmse_range_reduced24=mean(apply(
      (cbind(posterior_mean_range_reduced2[1:(0.5*length(posterior_mean_range_reduced2))],
             posterior_mean_range_reduced2[
               (0.5*length(posterior_mean_range_reduced2)+1):
                 (length(posterior_mean_range_reduced2))])
       -Yte)^2,1,sum))^.5)
    
    
    
    posterior_cov_range_reduced2= cov_mat_range_reduced2(Xte,Xte,opt_par[1],
                                                         opt_par[2],
                                                         opt_par[3],
                                                         opt_par[4],opt_par[5],
                                                         opt_par[6],
                                                         opt_par[7],opt_par[8],
                                                         opt_par[9],
                                                         opt_par[10])-
      Ktetr_trtr_inv_range_reduced2%*%t(Ktetr_range_reduced2)
    
    
    (LOGS[i,6]=LogS_range_reduced24_reduced2=(t(c(Yte[,1],Yte[,2])-posterior_mean_range_reduced2)%*%
                                                inv(posterior_cov_range_reduced2+diag(opt_par[11]^2,2*nrow(Xte)))%*%
                                                (c(Yte[,1],Yte[,2])-posterior_mean_range_reduced2)
                                              +determinant(posterior_cov_range_reduced2+diag(opt_par[11]^2,2*nrow(Xte)))$modulus[1]+
                                                log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    
    
    opt_par=Adam(function(x) {
      log_likelihood_no_range (
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7], x[8], x[9], x[10],
        x[11], x[12], x[13],x[14], x[15], x[16], x[17],x[18], x[19], x[20], x[21])
    },function(x) {
      grad_no_range(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                    x[6], x[7], x[8], x[9], x[10],
                    x[11], x[12], x[13],x[14], x[15],
                    x[16], x[17],x[18], x[19], x[20], x[21])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_l1,
        initial_l2,
        initial_sigma5,
        initial_sigma6,
        initial_sigma7,
        initial_sigma8,
        initial_l3,
        initial_l4,
        initial_sigma9,
        initial_sigma10,
        initial_sigma11,
        initial_sigma12,
        initial_l5,
        initial_l6,
        initial_l_range1,
        initial_l_range2,
        initial_sigma_obs),
    lr=lr,max_iter=max_iter)
    
    Ktetr_range=cov_mat_no_range(Xte,Xtr,opt_par[1],
                                 opt_par[2],
                                 opt_par[3],
                                 opt_par[4],
                                 opt_par[5],
                                 opt_par[6],
                                 opt_par[7],
                                 opt_par[8],
                                 opt_par[9],
                                 opt_par[10],
                                 opt_par[11],
                                 opt_par[12],
                                 opt_par[13],
                                 opt_par[14],
                                 opt_par[15],
                                 opt_par[16],
                                 opt_par[17],
                                 opt_par[18],
                                 opt_par[19],
                                 opt_par[20]
    )
    
    Ktetr_trtr_inv_range=Ktetr_range%*%inv(cov_mat_no_range(Xtr,Xtr,
                                                            opt_par[1],
                                                            opt_par[2],
                                                            opt_par[3],
                                                            opt_par[4],
                                                            opt_par[5],
                                                            opt_par[6],
                                                            opt_par[7],
                                                            opt_par[8],
                                                            opt_par[9],
                                                            opt_par[10],
                                                            opt_par[11],
                                                            opt_par[12],
                                                            opt_par[13],
                                                            opt_par[14],
                                                            opt_par[15],
                                                            opt_par[16],
                                                            opt_par[17],
                                                            opt_par[18],
                                                            opt_par[19],
                                                            opt_par[20])+
                                             diag(opt_par[21]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_range=Ktetr_trtr_inv_range%*%c(Ytr[,1],Ytr[,2])
    
    
    
    
    (RMSES[i,7]=rmse_range5=mean(apply((cbind(
      posterior_mean_range[1:(0.5*length(posterior_mean_range))],
      posterior_mean_range[(0.5*length(posterior_mean_range)+1):
                             (length(posterior_mean_range))])
      -Yte)^2,1,sum))^.5)
    
    
    
    posterior_cov_range= cov_mat_no_range(Xte,Xte,opt_par[1],
                                          opt_par[2],
                                          opt_par[3],
                                          opt_par[4],opt_par[5],
                                          opt_par[6],
                                          opt_par[7],opt_par[8],
                                          opt_par[9],
                                          opt_par[10],
                                          opt_par[11],opt_par[12],
                                          opt_par[13],opt_par[14],
                                          opt_par[15],
                                          opt_par[16],
                                          opt_par[17],
                                          opt_par[18],
                                          opt_par[19],
                                          opt_par[20])-
      Ktetr_trtr_inv_range%*%t(Ktetr_range)
    
    
    (LOGS[i,7]=LogS_range5=(t(c(Yte[,1],Yte[,2])-posterior_mean_range)%*%
                              inv(posterior_cov_range+diag(opt_par[21]^2,2*nrow(Xte)))%*%
                              (c(Yte[,1],Yte[,2])-posterior_mean_range)
                            +determinant(posterior_cov_range+
                                           diag(opt_par[21]^2,2*nrow(Xte)))$modulus[1]+
                              log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    
  }
  
}else{
  RMSES=LOGS=matrix(nrow=N_train,ncol=7)  
  for( i in 1:N_train){
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
    train_ind=sample(564,5*i)
    Xtr=X[train_ind,];Xte=X[-train_ind,]
    Ytr=Y[train_ind,];Yte=Y[-train_ind,]
    
    opt_par=Adam(function(x) {
      log_likelihood(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7])
    },function(x) {
      grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],x[6], x[7])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_sigma_obs,
        initial_l1,
        initial_l2),lr=lr,max_iter = max_iter)
    
    Ktetr_standard=cov_mat(Xte,Xtr,opt_par[1],
                           opt_par[2],
                           opt_par[3],
                           opt_par[4],
                           opt_par[6],
                           opt_par[7])
    
    Ktetr_trtr_inv_standard=Ktetr_standard%*%inv(cov_mat(Xtr,Xtr,opt_par[1],
                                                         opt_par[2],
                                                         opt_par[3],
                                                         opt_par[4],
                                                         opt_par[6],
                                                         opt_par[7])
                                                 
                                                 +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_standard=Ktetr_trtr_inv_standard%*%c(Ytr[,1],Ytr[,2])
    
    
    
    (RMSES[i,1]=rmse_standard=mean(apply((cbind(
      posterior_mean_standard[1:(0.5*length(posterior_mean_standard))],
      posterior_mean_standard
      [(0.5*length(posterior_mean_standard)+1):(length(posterior_mean_standard))])
      -Yte)^2,1,sum))^.5)
    
    
    
    posterior_cov_standard= cov_mat(Xte,Xte,opt_par[1],
                                    opt_par[2],
                                    opt_par[3],
                                    opt_par[4],
                                    opt_par[6],
                                    opt_par[7])-Ktetr_trtr_inv_standard%*%t(Ktetr_standard)
    
    
    (LOGS[i,1]=LogS_standard=(t(c(Yte[,1],Yte[,2])-posterior_mean_standard)%*%
                      inv(posterior_cov_standard+diag(opt_par[5]^2,2*nrow(Xte)))%*%
                      (c(Yte[,1],Yte[,2])-posterior_mean_standard)
                    +determinant(posterior_cov_standard+diag(opt_par[5]^2,2*nrow(Xte))
                    )$modulus[1]+
                      log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    
    opt_par=Adam(function(x) {
      log_likelihood_helm(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
    },function(x) {
      grad_helm(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
    },c(initial_par[1:4],initial_sigma_obs),lr=lr,max_iter=max_iter)
    
    Ktetr_helm=cov_mat_helm(Xte,Xtr,opt_par[1],
                            opt_par[2],
                            opt_par[3],
                            opt_par[4])
    
    Ktetr_trtr_inv_helm=Ktetr_helm%*%inv(cov_mat_helm(Xtr,Xtr,opt_par[1],
                                                      opt_par[2],
                                                      opt_par[3],
                                                      opt_par[4])
                                         
                                         +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_helm=Ktetr_trtr_inv_helm%*%c(Ytr[,1],Ytr[,2])
    
    
    
    (RMSES[i,2]=rmse_helm=mean(apply((cbind(posterior_mean_helm[1:(0.5*length(posterior_mean_helm))],
                                 posterior_mean_helm[(0.5*length(posterior_mean_helm)+1):
                                                       (length(posterior_mean_helm))])
                           -Yte)^2,1,sum))^.5)
    
    posterior_cov_helm= cov_mat_helm(Xte,Xte,opt_par[1],
                                     opt_par[2],
                                     opt_par[3],
                                     opt_par[4])-Ktetr_trtr_inv_helm%*%t(Ktetr_helm)
    
    
    (LOGS[i,2]=LogS_helm=(t(c(Yte[,1],Yte[,2])-posterior_mean_helm)%*%
                  inv(posterior_cov_helm+diag(opt_par[5]^2,2*nrow(Xte)))%*%
                  (c(Yte[,1],Yte[,2])-posterior_mean_helm)
                +determinant(posterior_cov_helm+diag(opt_par[5]^2,2*nrow(Xte)))$modulus[1]+
                  log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    
    
    opt_par=Adam(function(x) {
      log_likelihood_mixture(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5], x[6],x[7],x[8])
    },function(x) {
      grad_mix(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5], x[6],x[7],x[8])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_sigma_obs,
        initial_l1,
        initial_l2,
        initial_alpha),lb=c(rep(0.1,7),0.2),ub=c(rep(Inf,7),0.8),lr=lr,
    max_iter=max_iter
    )
    
    
    Ktetr_mix=mix_cov_mat(Xte,Xtr,opt_par[1],
                          opt_par[2],
                          opt_par[3],
                          opt_par[4],
                          opt_par[6],
                          opt_par[7],
                          opt_par[8])
    
    Ktetr_trtr_inv_mix=Ktetr_mix%*%inv(mix_cov_mat(Xtr,Xtr,opt_par[1],
                                                   opt_par[2],
                                                   opt_par[3],
                                                   opt_par[4],
                                                   opt_par[6],
                                                   opt_par[7],
                                                   opt_par[8])
                                       
                                       +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_mix=Ktetr_trtr_inv_mix%*%c(Ytr[,1],Ytr[,2])
    
    
    
    (RMSES[i,3]=rmse_mix=mean(apply((cbind(posterior_mean_mix[1:(0.5*length(posterior_mean_mix))],
                                posterior_mean_mix[(0.5*length(posterior_mean_mix)+1):
                                                     (length(posterior_mean_mix))])
                          -Yte)^2,1,sum))^.5)
    
    posterior_cov_mix= mix_cov_mat(Xte,Xte,opt_par[1],
                                   opt_par[2],
                                   opt_par[3],
                                   opt_par[4],opt_par[6],
                                   opt_par[7],opt_par[8])-Ktetr_trtr_inv_mix%*%
      t(Ktetr_mix) 
    
    (LOGS[i,3]=LogS_mix=(t(c(Yte[,1],Yte[,2])-posterior_mean_mix)%*%
                 inv(posterior_cov_mix+diag(opt_par[5]^2,2*nrow(Xte)))%*%
                 (c(Yte[,1],Yte[,2])-posterior_mean_mix)
               +determinant(posterior_cov_mix+diag(opt_par[5]^2,2*nrow(Xte)))$modulus[1]+
                 log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    opt_par=Adam(function(x) {
      log_likelihood_range_one_center(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7], x[8], x[9], x[10],
        x[11], x[12], x[13],x[14], x[15])
    },function(x) {
      grad_range_one_center(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                            x[6], x[7], x[8], x[9], x[10],
                            x[11], x[12], x[13],x[14], x[15])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_l1,
        initial_l2,
        initial_sigma5,
        initial_sigma6,
        initial_sigma7,
        initial_sigma8,
        initial_l3,
        initial_l4,
        initial_range1,
        initial_l_range1,
        initial_sigma_obs),max_iter = max_iter,lr=lr)
    
    Ktetr_range=cov_mat_range_one_center(Xte,Xtr,opt_par[1],
                                         opt_par[2],
                                         opt_par[3],
                                         opt_par[4],
                                         opt_par[5],
                                         opt_par[6],
                                         opt_par[7],opt_par[8],
                                         opt_par[9],
                                         opt_par[10],
                                         opt_par[11],opt_par[12],
                                         opt_par[13],opt_par[14])
    
    Ktetr_trtr_inv_range=Ktetr_range%*%inv(cov_mat_range_one_center(Xtr,Xtr,
                                                                    opt_par[1],
                                                                    opt_par[2],
                                                                    opt_par[3],
                                                                    opt_par[4],
                                                                    opt_par[5],
                                                                    opt_par[6],
                                                                    opt_par[7],opt_par[8],
                                                                    opt_par[9],
                                                                    opt_par[10],
                                                                    opt_par[11],opt_par[12],
                                                                    opt_par[13],opt_par[14])
                                           
                                           +diag(opt_par[15]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_range=Ktetr_trtr_inv_range%*%c(Ytr[,1],Ytr[,2])
    
    
    
    
    (RMSES[i,4]=rmse_range4=mean(apply((cbind(posterior_mean_range[
      1:(0.5*length(posterior_mean_range))],posterior_mean_range[(0.5*length(
        posterior_mean_range)+1):(length(posterior_mean_range))])
      -Yte)^2,1,sum))^.5)
    
    
    
    posterior_cov_range= cov_mat_range_one_center(Xte,Xte,opt_par[1],
                                                  opt_par[2],
                                                  opt_par[3],
                                                  opt_par[4],opt_par[5],
                                                  opt_par[6],
                                                  opt_par[7],opt_par[8],
                                                  opt_par[9],
                                                  opt_par[10],
                                                  opt_par[11],opt_par[12],
                                                  opt_par[13],opt_par[14]
    )-Ktetr_trtr_inv_range%*%
      t(Ktetr_range)
    
    
    (LOGS[i,4]=LogS_range4=(t(c(Yte[,1],Yte[,2])-posterior_mean_range)%*%
                    inv(posterior_cov_range+diag(opt_par[15]^2,2*nrow(Xte)))%*%
                    (c(Yte[,1],Yte[,2])-posterior_mean_range)
                  +determinant(posterior_cov_range+diag(opt_par[15]^2,2*nrow(Xte)))$modulus[1]+
                    log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    
    opt_par=Adam(function(x) {
      log_likelihood_range_one_center_reduced(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7], x[8], x[9])
    },function(x) {
      grad_range_one_center_reduced(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                                    x[6], x[7], x[8], x[9])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_l1,
        initial_l2,
        
        initial_range1,
        initial_l_range1,
        initial_sigma_obs),max_iter = max_iter,lr=lr)
    
    Ktetr_range=cov_mat_range_one_center_reduced(Xte,Xtr,opt_par[1],
                                                 opt_par[2],
                                                 opt_par[3],
                                                 opt_par[4],
                                                 opt_par[5],
                                                 opt_par[6],
                                                 opt_par[7],opt_par[8])
    
    Ktetr_trtr_inv_range=Ktetr_range%*%inv(cov_mat_range_one_center_reduced(Xtr,Xtr,
                                                                            opt_par[1],
                                                                            opt_par[2],
                                                                            opt_par[3],
                                                                            opt_par[4],
                                                                            opt_par[5],
                                                                            opt_par[6],
                                                                            opt_par[7],opt_par[8])
                                           
                                           +diag(opt_par[9]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_range=Ktetr_trtr_inv_range%*%c(Ytr[,1],Ytr[,2])
    
    
    
    
    (RMSES[i,5]=rmse_range4=mean(apply((cbind(posterior_mean_range[
      1:(0.5*length(posterior_mean_range))],posterior_mean_range[(0.5*length(
        posterior_mean_range)+1):(length(posterior_mean_range))])
      -Yte)^2,1,sum))^.5)
    
    
    
    posterior_cov_range= cov_mat_range_one_center_reduced(Xte,Xte,opt_par[1],
                                                          opt_par[2],
                                                          opt_par[3],
                                                          opt_par[4],opt_par[5],
                                                          opt_par[6],
                                                          opt_par[7],opt_par[8])-Ktetr_trtr_inv_range%*%
      t(Ktetr_range)
    
    
    (LOGS[i,5]=LogS_range4=(t(c(Yte[,1],Yte[,2])-posterior_mean_range)%*%
                    inv(posterior_cov_range+diag(opt_par[9]^2,2*nrow(Xte)))%*%
                    (c(Yte[,1],Yte[,2])-posterior_mean_range)
                  +determinant(posterior_cov_range+diag(opt_par[9]^2,2*nrow(Xte)))$modulus[1]+
                    log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    opt_par=Adam(function(x) {
      log_likelihood_range_one_center_norange(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7], x[8], x[9], x[10],
        x[11], x[12], x[13],x[14])
    },function(x) {
      grad_range_one_center_norange(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                                    x[6], x[7], x[8], x[9], x[10],
                                    x[11], x[12], x[13],x[14])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_l1,
        initial_l2,
        initial_sigma5,
        initial_sigma6,
        initial_sigma7,
        initial_sigma8,
        initial_l3,
        initial_l4,
        initial_l_range1,
        initial_sigma_obs),lr=lr,max_iter = max_iter)
    
    Ktetr_range=cov_mat_range_one_center_norange(Xte,Xtr,opt_par[1],
                                                 opt_par[2],
                                                 opt_par[3],
                                                 opt_par[4],
                                                 opt_par[5],
                                                 opt_par[6],
                                                 opt_par[7],
                                                 opt_par[8],
                                                 opt_par[9],
                                                 opt_par[10],
                                                 opt_par[11],
                                                 opt_par[12],
                                                 opt_par[13]
    )
    
    Ktetr_trtr_inv_range=Ktetr_range%*%inv(cov_mat_range_one_center_norange(Xtr,Xtr,
                                                                            opt_par[1],
                                                                            opt_par[2],
                                                                            opt_par[3],
                                                                            opt_par[4],
                                                                            opt_par[5],
                                                                            opt_par[6],
                                                                            opt_par[7],
                                                                            opt_par[8],
                                                                            opt_par[9],
                                                                            opt_par[10],
                                                                            opt_par[11],
                                                                            opt_par[12],
                                                                            opt_par[13])+
                                             diag(opt_par[14]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_range=Ktetr_trtr_inv_range%*%c(Ytr[,1],Ytr[,2])
    
    
    
    
    (RMSES[i,6]=rmse_range5=mean(
      apply((cbind(posterior_mean_range[1:(0.5*length(posterior_mean_range))],
                   posterior_mean_range[(0.5*length(posterior_mean_range)+1):
                                          (length(posterior_mean_range))])
             -Yte)^2,1,sum))^.5)
    
    
    
    posterior_cov_range= cov_mat_range_one_center_norange(Xte,Xte,opt_par[1],
                                                          opt_par[2],
                                                          opt_par[3],
                                                          opt_par[4],opt_par[5],
                                                          opt_par[6],
                                                          opt_par[7],opt_par[8],
                                                          opt_par[9],
                                                          opt_par[10],
                                                          opt_par[11],opt_par[12],
                                                          opt_par[13])-
      Ktetr_trtr_inv_range%*%t(Ktetr_range)
    
    
    (LOGS[i,6]=LogS_range5=(t(c(Yte[,1],Yte[,2])-posterior_mean_range)%*%
                    inv(posterior_cov_range+diag(opt_par[14]^2,2*nrow(Xte)))%*%
                    (c(Yte[,1],Yte[,2])-posterior_mean_range)
                  +determinant(posterior_cov_range+diag(opt_par[14]^2,2*nrow(Xte)))$modulus[1]+
                    log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    opt_par=Adam(function(x) {
      log_likelihood_range_one_center_norange_reduced(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7], x[8])
    },function(x) {
      grad_range_one_center_norange_reduced(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                                            x[6], x[7], x[8])
    },c(initial_sigma1,
        initial_sigma2,
        initial_sigma3,
        initial_sigma4,
        initial_l1,
        initial_l2,
        
        initial_l_range1,
        initial_sigma_obs),lr=lr,max_iter = max_iter)
    
    Ktetr_range=cov_mat_range_one_center_norange_reduced(Xte,Xtr,opt_par[1],
                                                         opt_par[2],
                                                         opt_par[3],
                                                         opt_par[4],
                                                         opt_par[5],
                                                         opt_par[6],
                                                         opt_par[7])
    
    Ktetr_trtr_inv_range=Ktetr_range%*%inv(cov_mat_range_one_center_norange_reduced(Xtr,Xtr,
                                                                                    opt_par[1],
                                                                                    opt_par[2],
                                                                                    opt_par[3],
                                                                                    opt_par[4],
                                                                                    opt_par[5],
                                                                                    opt_par[6],
                                                                                    opt_par[7])+
                                             diag(opt_par[8]^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_range=Ktetr_trtr_inv_range%*%c(Ytr[,1],Ytr[,2])
    
    
    
    
    (RMSES[i,7]=rmse_range5=mean(
      apply((cbind(posterior_mean_range[1:(0.5*length(posterior_mean_range))],
                   posterior_mean_range[(0.5*length(posterior_mean_range)+1):
                                          (length(posterior_mean_range))])
             -Yte)^2,1,sum))^.5)
    
    
    
    posterior_cov_range= cov_mat_range_one_center_norange_reduced(Xte,Xte,opt_par[1],
                                                                  opt_par[2],
                                                                  opt_par[3],
                                                                  opt_par[4],opt_par[5],
                                                                  opt_par[6],
                                                                  opt_par[7])-
      Ktetr_trtr_inv_range%*%t(Ktetr_range)
    
    
    (LOGS[i,7]=LogS_range5=(t(c(Yte[,1],Yte[,2])-posterior_mean_range)%*%
                    inv(posterior_cov_range+diag(opt_par[8]^2,2*nrow(Xte)))%*%
                    (c(Yte[,1],Yte[,2])-posterior_mean_range)
                  +determinant(posterior_cov_range+diag(opt_par[8]^2,2*nrow(Xte)))$modulus[1]+
                    log(2*pi)*nrow(Xte))/nrow(Xte))
    
  }
  
}

res=list(RMSES,LOGS)

save(list = "res", file = paste0("Learning_curves_gulf_001_", id, ".rda"))

