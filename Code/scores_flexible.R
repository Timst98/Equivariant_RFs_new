id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

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

train_points=list(c(-3.8,3.5),c(-2,3.5),c(-3.5,1.3)
                               ,c(-3.1,-0.5),
                               c(-2.85,-0.75),c(-2,-1.2),c(-2,1),c(-3.25,-0.7),
                               c(-3.3,-0.6),c(-2.2,-2),c(-3.3,-3),
                               c(2.9,-0.1),c(3.3,-0.05),
                               c(2,-0.25),c(2.75,1.05),c(3.5,0.3),
                               c(3.2,-1.75),c(3.3,-3),c(1.6,-2.8),c(0.4,3.5),
                               c(-1,-3),c(1,0.5),c(2.5,2.5),c(-1.2,1.5),
                               c(0.2,2.5),c(0.03,-2)
)### NEW points

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

nsim=500

cov_mat_range_flexible=function(x1, x2, l1, sigma1, l2, sigma2,l3, sigma3, l4, sigma4,
                        l5, sigma5, l6, sigma6, range1,range2,l_range1,l_range2,
                        c11,c12,c21,c22){
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    if(n2==1){
      dist_mat=sum((x1-x2)^2)^.5
    }else{dist_mat=dist(t(x1),x2)}
  }else{
    dist_mat=dist(x1,x2)
  }
  
  
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(c11,n1),c12))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(c11,n2),c12))^2,1,sum)^.5)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(c21,n1),c22))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(c21,n2),c22))^2,1,sum)^.5)
  
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
  
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  cov[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2/(l1^2))
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*exp(-.5*dist_mat^2/(l2^2))
  
  cov=cov*sigmoid_se
  
  
  
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
  cov_fund1=cov_fund1*sigmoid_dist_center1
  cov_fund2=cov_fund2*sigmoid_dist_center2
  return(1/(weight_sum)*(cov+cov_fund2+cov_fund1))
  
}


grad_range_flexible=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,l3_2,sigma3_2,l4_2,sigma4_2,
                     l5_2,sigma5_2,l6_2,sigma6_2,
                     range1,range2,l_range1,l_range2,
                     c11,c12,c21,c22,sigma_obs_2){
  
  l1=l1_2;sigma1=sigma1_2;l2=l2_2;sigma2=sigma2_2;l3=l3_2;
  sigma3=sigma3_2;l4=l4_2;sigma4=sigma4_2;l5=l5_2;sigma5=sigma5_2;l6=l6_2;sigma6=sigma6_2
  
  x1=xtr
  x2=x1
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  distmat_center1_x1=as.numeric(apply((x1-cbind(rep(c11,n1),c12))^2,1,sum)^.5)
  distmat_center1_x2=as.numeric(apply((x2-cbind(rep(c11,n2),c12))^2,1,sum)^.5)
  
  del_distmat_center1_x1_c11=as.numeric(apply((x1-cbind(rep(c11,n1),c12))^2,1,sum)^(-.5))*(-x1[,1]+c11)
  del_distmat_center1_x1_c12=as.numeric(apply((x1-cbind(rep(c11,n1),c12))^2,1,sum)^(-.5))*(-x1[,2]+c12)
  
  del_distmat_center1_x2_c11=as.numeric(apply((x2-cbind(rep(c11,n2),c12))^2,1,sum)^(-.5))*(-x2[,1]+c11)
  del_distmat_center1_x2_c12=as.numeric(apply((x2-cbind(rep(c11,n2),c12))^2,1,sum)^(-.5))*(-x2[,2]+c12)
  
  
  sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  sigmoid_dist_center1=(1-sigmoid_11)%*%t(1-sigmoid_12)
  
  del_c11_1=(-(sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)*
             (del_distmat_center1_x1_c11)/l_range1)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)*
        (-del_distmat_center1_x2_c11)/l_range1)
  
  del_c12_1=(-(sigmoid_11)^2*exp(-(distmat_center1_x1-range1)/l_range1)*
               (del_distmat_center1_x1_c12)/l_range1)%*%t(1-sigmoid_12)+
    (1-sigmoid_11)%*%t(
      (sigmoid_12)^2*exp(-(distmat_center1_x2-range1)/l_range1)*
        (-del_distmat_center1_x2_c12)/l_range1)
  
  
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
  
  del_c11_1=cbind(rbind(del_c11_1,del_c11_1),
                     rbind(del_c11_1,del_c11_1))
  del_c12_1=cbind(rbind(del_c12_1,del_c12_1),
                  rbind(del_c12_1,del_c12_1))
  
  
  distmat_center2_x1=as.numeric(apply((x1-cbind(rep(c21,n1),c22))^2,1,sum)^.5)
  distmat_center2_x2=as.numeric(apply((x2-cbind(rep(c21,n2),c22))^2,1,sum)^.5)
  
  
  del_distmat_center2_x1_c21=as.numeric(apply((x1-cbind(rep(c21,n1),c22))^2,1,sum)^(-.5))*(-x1[,1]+c21)
  del_distmat_center2_x1_c22=as.numeric(apply((x1-cbind(rep(c21,n1),c22))^2,1,sum)^(-.5))*(-x1[,2]+c22)
  
  del_distmat_center2_x2_c21=as.numeric(apply((x2-cbind(rep(c21,n2),c22))^2,1,sum)^(-.5))*(-x2[,1]+c21)
  del_distmat_center2_x2_c22=as.numeric(apply((x2-cbind(rep(c21,n2),c22))^2,1,sum)^(-.5))*(-x2[,2]+c22)
  
  
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
  
  #sigmoid_dist_center2=(1-sigmoid_21)%*%t(1-sigmoid_22)
  
  
  del_c21_1=(-(sigmoid_21)^2*exp(-(distmat_center2_x1-range2)/l_range2)*
               (del_distmat_center2_x1_c21)/l_range2)%*%t(1-sigmoid_22)+
    (1-sigmoid_21)%*%t(
      (sigmoid_22)^2*exp(-(distmat_center2_x2-range2)/l_range2)*
        (-del_distmat_center2_x2_c21)/l_range2)
  
  del_c22_1=(-(sigmoid_21)^2*exp(-(distmat_center2_x1-range2)/l_range2)*
               (del_distmat_center2_x1_c22)/l_range2)%*%t(1-sigmoid_22)+
    (1-sigmoid_21)%*%t(
      (sigmoid_22)^2*exp(-(distmat_center2_x2-range2)/l_range2)*
        (-del_distmat_center2_x2_c22)/l_range2)
  
  del_c21_1=cbind(rbind(del_c21_1,del_c21_1),
                  rbind(del_c21_1,del_c21_1))
  del_c22_1=cbind(rbind(del_c22_1,del_c22_1),
                  rbind(del_c22_1,del_c22_1))
  
 
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
  
  #sigmoid_se=(sigmoid_11*sigmoid_21%*%t(sigmoid_12*sigmoid_22))
  
  #sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  #sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  #sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  #sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  
  #distmat_center1_x1=as.numeric(apply((x1-cbind(rep(c11,n1),c12))^2,1,sum)^.5)
  #distmat_center1_x2=as.numeric(apply((x2-cbind(rep(c11,n2),c12))^2,1,sum)^.5)
  
  #distmat_center2_x1=as.numeric(apply((x1-cbind(rep(c21,n1),c22))^2,1,sum)^.5)
  #distmat_center2_x2=as.numeric(apply((x2-cbind(rep(c21,n2),c22))^2,1,sum)^.5)
  
  #del_distmat_center1_x1_c11=as.numeric(apply((x1-cbind(rep(c11,n1),c12))^2,1,sum)^(-.5))*(-x1[,1]+c11)
  #del_distmat_center1_x1_c12=as.numeric(apply((x1-cbind(rep(c11,n1),c12))^2,1,sum)^(-.5))*(-x1[,2]+c12)
  
  #del_distmat_center1_x2_c11=as.numeric(apply((x2-cbind(rep(c11,n2),c12))^2,1,sum)^(-.5))*(-x2[,1]+c11)
  #del_distmat_center1_x2_c12=as.numeric(apply((x2-cbind(rep(c11,n2),c12))^2,1,sum)^(-.5))*(-x2[,2]+c12)
  
  #del_distmat_center2_x1_c21=as.numeric(apply((x1-cbind(rep(c21,n1),c22))^2,1,sum)^(-.5))*(-x1[,1]+c21)
  #del_distmat_center2_x1_c22=as.numeric(apply((x1-cbind(rep(c21,n1),c22))^2,1,sum)^(-.5))*(-x1[,2]+c22)
  
  #del_distmat_center2_x2_c21=as.numeric(apply((x2-cbind(rep(c21,n2),c22))^2,1,sum)^(-.5))*(-x2[,1]+c21)
  #del_distmat_center2_x2_c22=as.numeric(apply((x2-cbind(rep(c21,n2),c22))^2,1,sum)^(-.5))*(-x2[,2]+c22)
  
  
  #sigmoid_se=(sigmoid_11*sigmoid_21%*%t(sigmoid_12*sigmoid_22))
  
  #sigmoid_11=sigmoid((distmat_center1_x1-range1)/l_range1)
  #sigmoid_12=sigmoid((distmat_center1_x2-range1)/l_range1)
  
  #sigmoid_21=sigmoid((distmat_center2_x1-range2)/l_range2)
  #sigmoid_22=sigmoid((distmat_center2_x2-range2)/l_range2)
  
  del_c11_2=-(((sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                (-del_distmat_center1_x1_c11)/(l_range1))*sigmoid_21)%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t((sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)*
                                   (-del_distmat_center1_x2_c11)/(l_range1))*sigmoid_22))
  
  del_c12_2=-(((sigmoid_11^2*exp(-(distmat_center1_x1-range1)/l_range1)*
                (-del_distmat_center1_x1_c12)/(l_range1))*sigmoid_21)%*%t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t((sigmoid_12^2*exp(-(distmat_center1_x2-range1)/l_range1)*
                                   (-del_distmat_center1_x2_c12)/(l_range1))*sigmoid_22))
  
  del_c21_2=-(sigmoid_11*(sigmoid_21^2*exp(-(distmat_center2_x1-range2)/l_range2)*(-(del_distmat_center2_x1_c21))/l_range2))%*%
    t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t(sigmoid_12*(-sigmoid_22^2*exp(-(distmat_center2_x2-range2)/l_range2)*(-(del_distmat_center2_x2_c21))/l_range2))
  
  del_c22_2=-(sigmoid_11*(sigmoid_21^2*exp(-(distmat_center2_x1-range2)/l_range2)*(-(del_distmat_center2_x1_c22))/l_range2))%*%
    t(sigmoid_12*sigmoid_22)+
    (sigmoid_11*sigmoid_21)%*%t(sigmoid_12*(-sigmoid_22^2*exp(-(distmat_center2_x2-range2)/l_range2)*(-(del_distmat_center2_x2_c22))/l_range2))
  
  del_c21_2=cbind(rbind(del_c21_2,del_c21_2),
                  rbind(del_c21_2,del_c21_2))
  
  del_c22_2=cbind(rbind(del_c22_2,del_c22_2),
                  rbind(del_c22_2,del_c22_2))
  
  del_c12_2=cbind(rbind(del_c12_2,del_c12_2),
                  rbind(del_c12_2,del_c12_2))
  
  del_c11_2=cbind(rbind(del_c11_2,del_c11_2),
                  rbind(del_c11_2,del_c11_2))
  

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
  
  K1_l3_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K1_l4_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                   nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K1_sigma3_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K1_sigma4_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                       nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  K2_l5_2=K2_sigma5_2=K2_l6_2=K2_sigma6_2=K1_l3_2
  
  K=cov_mat_range_flexible(xtr,xtr,l1_2,sigma1_2,l2_2,sigma2_2,
                   l3_2,sigma3_2,l4_2,sigma4_2,
                   l5_2,sigma5_2,l6_2,sigma6_2,range1,range2,l_range1,l_range2,
                   c11,c12,c21,c22)+
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
        diag(c(sigma5_2^2*exp(-.5*(r1-r2)^2/(l5_2^2))*(r1-r2)^2/(l5_2^3),0))%*%
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
  
  K1_l3_2=K1_l3_2*(sigmoid_dist_center1)*1/(weight_sum)
  K1_l4_2=K1_l4_2*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma3_2=K1_sigma3_2*(sigmoid_dist_center1)*1/(weight_sum)
  K1_sigma4_2=K1_sigma4_2*(sigmoid_dist_center1)*1/(weight_sum)
  
  
  
  K2_l5_2=K2_l5_2*(sigmoid_dist_center2)*1/(weight_sum)
  K2_l6_2=K2_l6_2*(sigmoid_dist_center2)*1/(weight_sum)
  K2_sigma5_2=K2_sigma5_2*(sigmoid_dist_center2)*1/(weight_sum)
  K2_sigma6_2=K2_sigma6_2*(sigmoid_dist_center2)*1/(weight_sum)
  
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
  
  del_cov1=sigmoid_cov*del_cov1*1/(weight_sum)
  del_cov2=sigmoid_cov*del_cov2*1/(weight_sum)
  del_cov3=sigmoid_cov*del_cov3*1/(weight_sum)
  del_cov4=sigmoid_cov*del_cov4*1/(weight_sum)
  
  
  K_fund1=cov_fund1
  K_fund2=cov_fund2
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  cov[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2/(l1^2))
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*exp(-.5*dist_mat^2/(l2^2))
  
  K_cov=cov
  
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
  
  
  del_c11=(K_fund1*del_c11_1+K_cov*del_c11_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_c11_1+del_c11_2)
  
  del_c12=(K_fund1*del_c12_1+K_cov*del_c12_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_c12_2+del_c12_1)
  
  del_c21=(K_fund2*(del_c21_1)+K_cov*del_c21_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_c21_2+del_c21_1)
  
  
  del_c22=(K_fund2*(del_c22_1)+K_cov*del_c22_2)*1/(weight_sum)-
    (K_fund1*sigmoid_dist_center1+K_cov*sigmoid_se+K_fund2*sigmoid_dist_center2)/(weight_sum^2)*
    (del_c22_2+del_c22_1)
  
  
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
            
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_c11%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_c11)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_c12%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_c12)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_c21%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_c21)),
            -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_c22%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%(alpha*del_c22)),
            
            -t(c(ytr[,1],ytr[,2]))%*%(2*sigma_obs_2*inv_K)%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  #print(grad)
  return(grad)
  
}


log_likelihood_range_flexible=function(ytr,xtr,l1,sigma1,l2,sigma2,l3,sigma3,l4,sigma4,
                               l5,sigma5,l6,sigma6,
                               range1,range2,l_range1,l_range2,c11,c12,c21,c22,
                               sigma_obs){
  if(length(ytr)>2){
    if(ncol(ytr)==2){
      ytr=c(ytr[,1],ytr[,2])
    }
  }
  ntr=ifelse(length(xtr)==2,1,nrow(xtr))
  
  Ktr=cov_mat_range_flexible(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3,l4,sigma4,
                     l5,sigma5,l6,sigma6,range1,range2,l_range1,l_range2,c11,c12,
                     c21,c22)+
    sigma_obs^2*diag(nrow=2*ntr)
  
  ll=-0.5*t(ytr)%*%solve(Ktr)%*%ytr-0.5*log(det(Ktr))-ntr*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}



library(nloptr)

nloptr::nloptr(runif(21,0.5,20),function(x) {
  log_likelihood_range_flexible(
    Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7], x[8], x[9], x[10],
    x[11], x[12], x[13],x[14], x[15], x[16], x[17],x[18],x[19], x[20], x[21])
},function(x) {
  grad_range_flexible(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
              x[6], x[7], x[8], x[9], x[10],
              x[11], x[12], x[13],x[14], x[15], x[16], x[17],
              x[18],x[19], x[20], x[21])
},opts=list(algorithm="NLOPT_LD_LBFGS",check_derivatives=TRUE))
lb_1=range(Xte[,1])[1]
ub_1=range(Xte[,1])[2]

lb_2=range(Xte[,2])[1]
ub_2=range(Xte[,2])[2]

ntr=nrow(Xtr)
for(j in 1:nsim){
  
  initial_par=c(runif(16,0.8,3),
                runif(1,lb_1,ub_1),
                runif(1,lb_2,ub_2),
                runif(1,lb_1,ub_1),
                runif(1,lb_2,ub_2),
                runif(1,0.0001,0.9999))
  
  new_score=0
  for (i in 1:ntr) {
    ind_val=i
    opt_par=nloptr(initial_par,eval_f=function(x) {
      
      log_likelihood_range_flexible(
        Ytr[-ind_val,], Xtr[-ind_val,],x[1], x[2], x[3], x[4], x[5],x[6], x[7],x[8], x[9],x[10],x[11],
        x[12], x[13],x[14],x[15],x[16],x[17],
        x[18],x[19], x[20], x[21])
    },
    eval_grad_f=function(x) {
      grad_range_flexible(Ytr[-ind_val,], Xtr[-ind_val,],x[1], x[2], x[3], x[4], x[5],x[6], x[7],x[8], x[9],x[10],x[11],
                  x[12], x[13],x[14],x[15],x[16],x[17],
                  x[18],x[19], x[20], x[21])
    },lb= c(rep(0.8,16),lb_1,lb_2,lb_1,lb_2,0.0001),
    ub=c(rep(4,16),ub_1,ub_2,ub_1,ub_2,1),opts=list(algorithm="NLOPT_LD_LBFGS",
                                                    "max_eval"=100))$solution
        
    
    
    
    
    (new_score=new_score+1/ntr*log_likelihood_range_flexible(Ytr[ind_val,], Xtr[ind_val,],opt_par[1],
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
                                                     opt_par[15],opt_par[16],opt_par[17],
                                                     opt_par[18] ,opt_par[19], 
                                                     opt_par[20], opt_par[21]))
    
  }
  
  
  if(j==1){
    best_initial_par=initial_par
    best_score=new_score
  }else{
    if(best_score>new_score){
      best_initial_par=initial_par
      best_score=new_score
      print(best_score)
    }
  }
  
}
time=toc()

initial_par=best_initial_par



opt_par=nloptr(initial_par,eval_f=function(x) {
  
  log_likelihood_range_flexible(
    Ytr, Xtr,x[1], x[2], x[3], x[4], x[5],x[6], x[7],x[8], x[9],x[10],x[11],
    x[12], x[13],x[14],x[15],x[16],x[17],
    x[18],x[19], x[20], x[21])
},
eval_grad_f=function(x) {
  grad_range_flexible(Ytr, Xtr,x[1], x[2], x[3], x[4], x[5],x[6], x[7],x[8], x[9],x[10],x[11],
                      x[12], x[13],x[14],x[15],x[16],x[17],
                      x[18],x[19], x[20], x[21])
},
lb= c(rep(0.8,16),lb_1,lb_2,lb_1,lb_2,0.0001),
ub=c(rep(4,16),ub_1,ub_2,ub_1,ub_2,1),opts=list(algorithm="NLOPT_LD_LBFGS"))$solution


Ktetr_range=cov_mat_range_flexible(Xte,Xtr,opt_par[1],
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
                           opt_par[15],opt_par[16],opt_par[17],
                           opt_par[18] ,opt_par[19], 
                           opt_par[20])

Ktetr_trtr_inv_range=Ktetr_range%*%solve(cov_mat_range_flexible(Xtr,Xtr,opt_par[1],
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
                                                        opt_par[15],opt_par[16],opt_par[17],
                                                        opt_par[18] ,opt_par[19], 
                                                        opt_par[20])+
                                         
                                         diag(opt_par[21]^2,nrow=2*nrow(Xtr)))

posterior_mean_range=Ktetr_trtr_inv_range%*%c(Ytr[,1],Ytr[,2])




(rmse_range=mean(apply((cbind(posterior_mean_range[1:(0.5*length(posterior_mean_range))],
                               posterior_mean_range[(0.5*length(posterior_mean_range)+1):
                                                      (length(posterior_mean_range))])
                         -Yte)^2,1,sum))^.5)



posterior_cov_range= cov_mat_range_flexible(Xte,Xte,opt_par[1],
                                    opt_par[2],
                                    opt_par[3],
                                    opt_par[4],opt_par[5],
                                    opt_par[6],
                                    opt_par[7],opt_par[8],
                                    opt_par[9],
                                    opt_par[10],
                                    opt_par[11],opt_par[12],
                                    opt_par[13],opt_par[14],opt_par[15],
                                    opt_par[16],opt_par[17],
                                    opt_par[18] ,opt_par[19], 
                                    opt_par[20])-Ktetr_trtr_inv_range%*%t(Ktetr_range)


(LogS_range=(t(c(Yte[,1],Yte[,2])-posterior_mean_range)%*%solve(posterior_cov_range+
                                                                  diag(1e-6,2*nrow(Xte)))%*%
               (c(Yte[,1],Yte[,2])-posterior_mean_range)
             +determinant(posterior_cov_range+diag(1e-6,2*nrow(Xte)))$modulus[1]+
               log(2*pi)*nrow(Xte))/nrow(Xte))


rmses=c(rmse_range)
LogS_scores=c(LogS_range)
centers=c(opt_par[17:20])


res=list(rmses,LogS_scores,centers)

save(list = "res", file = paste0("Scores_flexible_ex3_", id, ".rda"))


