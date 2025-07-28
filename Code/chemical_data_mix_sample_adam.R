id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

ind_train=sample(3641,330)
ind_test=-ind_train


Adam=function(f, grad_f, params_init, lr = lr, beta1 = 0.9, beta2 = 0.999,
              epsilon = 1e-8, max_iter = max_iter, tol = 1e-10) {
  
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


jit=1e-8

fund_cov_mat3d=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  norm=function(x){sum(x^2)^.5}
  theta1=atan2(x1[1],x1[2])
  x2_1=cos(theta1)*x1[2]+sin(theta1)*x1[1]
  phi1=atan2(x2_1,x1[3])
  
  theta2=atan2(x2[1],x2[2])
  x2_2=cos(theta2)*x2[2]+sin(theta2)*x2[1]
  phi2=atan2(x2_2,x2[3])
  
  t(matrix(c(1,0,0,0,cos(phi1),-sin(phi1),0,sin(phi1),cos(phi1)
  ),nrow=3,byrow = 1)%*%matrix(c(cos(theta1),-sin(theta1),0,
                                 sin(theta1),cos(theta1),0,0,0,1),nrow=3,
                               byrow = 1 ))%*%
    (exp(diag(c(-(norm(x1)-norm(x2))^2/(2*l1^2),
                -(norm(x1)-norm(x2))^2/(2*l2^2),
                -(norm(x1)-norm(x2))^2/(2*l3^2))))*
       diag(c(sigma1^2,sigma2^2,sigma3^2)))%*%
    matrix(c(1,0,0,0,cos(phi2),-sin(phi2),0,sin(phi2),cos(phi2)
    ),nrow=3,byrow = 1)%*%matrix(c(cos(theta2),
                                   -sin(theta2),0,
                                   sin(theta2),cos(theta2),0,0,0,1),nrow=3,
                                 byrow = 1 )
}


cross_fund_cov_mat3d=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  norm=function(x){sum(x^2)^.5}
  cov_O=cov_H1=cov_H2 = matrix(0, ncol =3 * ifelse(length(as.matrix(x2)) == 9, 1,nrow(x2)),
                               nrow = 3 * ifelse(length(as.matrix(x1)) == 9, 1, nrow(x1)))
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      if(length(as.matrix(x2)) == 9){
        if(length(as.matrix(x1)) == 9){
          cov1= fund_cov_mat3d(x1[1:3],x2[1:3],l1,sigma1,l2,sigma2,l3,sigma3)
          cov2= fund_cov_mat3d(x1[4:6],x2[4:6],l1,sigma1,l2,sigma2,l3,sigma3)
          cov3= fund_cov_mat3d(x1[7:9],x2[7:9],l1,sigma1,l2,sigma2,l3,sigma3)
          
        }else{
          cov1= fund_cov_mat3d(x1[i,1:3],x2[1:3],l1,sigma1,l2,sigma2,l3,sigma3)
          cov2= fund_cov_mat3d(x1[i,4:6],x2[4:6],l1,sigma1,l2,sigma2,l3,sigma3)
          cov3= fund_cov_mat3d(x1[i,7:9],x2[7:9],l1,sigma1,l2,sigma2,l3,sigma3)
          
        }
      }else if(length(as.matrix(x1)) == 9){
        cov1= fund_cov_mat3d(x1[1:3],x2[j,1:3],l1,sigma1,l2,sigma2,l3,sigma3)
        cov2= fund_cov_mat3d(x1[4:6],x2[j,4:6],l1,sigma1,l2,sigma2,l3,sigma3)
        cov3= fund_cov_mat3d(x1[7:9],x2[j,7:9],l1,sigma1,l2,sigma2,l3,sigma3)
        
      }else{
        cov1= fund_cov_mat3d(x1[i,1:3],x2[j,1:3],l1,sigma1,l2,sigma2,l3,sigma3)
        cov2= fund_cov_mat3d(x1[i,4:6],x2[j,4:6],l1,sigma1,l2,sigma2,l3,sigma3)
        cov3= fund_cov_mat3d(x1[i,7:9],x2[j,7:9],l1,sigma1,l2,sigma2,l3,sigma3)
        
      }
      
      cov_O[i, j] = cov1[1,1]
      cov_O[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
            ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,2]
      
      cov_O[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
            2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,3]
      
      cov_O[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[2,1]
      cov_O[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[3,1]
      
      cov_O[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,2]
      cov_O[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,3]
      
      cov_O[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
            ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,2]
      
      cov_O[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
            2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,3]
      
      
      
      cov_H1[i, j] = cov2[1,1]
      cov_H1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,2]
      
      cov_H1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,3]
      
      cov_H1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[2,1]
      cov_H1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[3,1]
      
      cov_H1[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,2]
      cov_H1[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,3]
      
      cov_H1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,2]
      
      cov_H1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,3]
      
      
      cov_H2[i, j] = cov3[1,1]
      cov_H2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[2,2]
      
      cov_H2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[3,3]
      
      cov_H2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov3[2,1]
      cov_H2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov3[3,1]
      
      cov_H2[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[1,2]
      cov_H2[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[1,3]
      
      cov_H2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[3,2]
      
      cov_H2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[2,3]
      
      
    }
  }
  return(cov_O+cov_H1+cov_H2)
}


cov_mat3d=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  norm=function(x){sum(x^2)^.5}
  (exp(diag(c(-(norm(x1)-norm(x2))^2/(2*l1^2),
              -(norm(x1)-norm(x2))^2/(2*l2^2),
              -(norm(x1)-norm(x2))^2/(2*l3^2))))*
      diag(c(sigma1^2,sigma2^2,sigma3^2)))
}



cross_cov_mat3d=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  norm=function(x){sum(x^2)^.5}
  cov_O=cov_H1=cov_H2 = matrix(0, ncol =3 * ifelse(length(as.matrix(x2)) == 9, 1,nrow(x2)),
                               nrow = 3 * ifelse(length(as.matrix(x1)) == 9, 1, nrow(x1)))
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      if(length(as.matrix(x2)) == 9){
        if(length(as.matrix(x1)) == 9){
          cov1= cov_mat3d(x1[1:3],x2[1:3],l1,sigma1,l2,sigma2,l3,sigma3)
          cov2= cov_mat3d(x1[4:6],x2[4:6],l1,sigma1,l2,sigma2,l3,sigma3)
          cov3= cov_mat3d(x1[7:9],x2[7:9],l1,sigma1,l2,sigma2,l3,sigma3)
          
        }else{
          cov1= cov_mat3d(x1[i,1:3],x2[1:3],l1,sigma1,l2,sigma2,l3,sigma3)
          cov2= cov_mat3d(x1[i,4:6],x2[4:6],l1,sigma1,l2,sigma2,l3,sigma3)
          cov3= cov_mat3d(x1[i,7:9],x2[7:9],l1,sigma1,l2,sigma2,l3,sigma3)
          
        }
      }else if(length(as.matrix(x1)) == 9){
        cov1= cov_mat3d(x1[1:3],x2[j,1:3],l1,sigma1,l2,sigma2,l3,sigma3)
        cov2= cov_mat3d(x1[4:6],x2[j,4:6],l1,sigma1,l2,sigma2,l3,sigma3)
        cov3= cov_mat3d(x1[7:9],x2[j,7:9],l1,sigma1,l2,sigma2,l3,sigma3)
        
      }else{
        cov1= cov_mat3d(x1[i,1:3],x2[j,1:3],l1,sigma1,l2,sigma2,l3,sigma3)
        cov2= cov_mat3d(x1[i,4:6],x2[j,4:6],l1,sigma1,l2,sigma2,l3,sigma3)
        cov3= cov_mat3d(x1[i,7:9],x2[j,7:9],l1,sigma1,l2,sigma2,l3,sigma3)
        
      }
      
      cov_O[i, j] = cov1[1,1]
      cov_O[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
            ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,2]
      
      cov_O[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
            2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,3]
      
      cov_O[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[2,1]
      cov_O[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[3,1]
      
      cov_O[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,2]
      cov_O[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,3]
      
      cov_O[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
            ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,2]
      
      cov_O[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
            2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,3]
      
      
      
      cov_H1[i, j] = cov2[1,1]
      cov_H1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,2]
      
      cov_H1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,3]
      
      cov_H1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[2,1]
      cov_H1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[3,1]
      
      cov_H1[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,2]
      cov_H1[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,3]
      
      cov_H1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,2]
      
      cov_H1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,3]
      
      
      cov_H2[i, j] = cov3[1,1]
      cov_H2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[2,2]
      
      cov_H2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[3,3]
      
      cov_H2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov3[2,1]
      cov_H2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov3[3,1]
      
      cov_H2[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[1,2]
      cov_H2[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[1,3]
      
      cov_H2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[3,2]
      
      cov_H2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
             2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[2,3]
      
      
    }
  }
  return(cov_O+cov_H1+cov_H2)
}


cov_mat3d_2=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  norm=function(x){sum(x^2)^.5}
  (exp(diag(c(-(norm(x1)-norm(x2))^2/(2*l1^2),
              -(norm(x1)-norm(x2))^2/(2*l2^2),
              -(norm(x1)-norm(x2))^2/(2*l3^2))))*
      diag(c(sigma1^2,sigma2^2,sigma3^2)))
}



cross_cov_mat3d_2=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  cov= matrix(0, ncol =3 * ifelse(length(as.matrix(x2)) == 9, 1,nrow(x2)),
              nrow = 3 * ifelse(length(as.matrix(x1)) == 9, 1, nrow(x1)))
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      if(length(as.matrix(x2)) == 9){
        if(length(as.matrix(x1)) == 9){
          cov1= cov_mat3d_2(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3)
          
          
        }else{
          cov1= cov_mat3d_2(x1[i,],x2,l1,sigma1,l2,sigma2,l3,sigma3)
          
        }
      }else if(length(as.matrix(x1)) == 9){
        cov1= cov_mat3d_2(x1,x2[j,],l1,sigma1,l2,sigma2,l3,sigma3)
        
        
      }else{
        cov1= cov_mat3d_2(x1[i,],x2[j,],l1,sigma1,l2,sigma2,l3,sigma3)
        
      }
      
      cov[i, j] = cov1[1,1]
      cov[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,2]
      
      cov[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,3]
      
      cov[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[2,1]
      cov[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[3,1]
      
      cov[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,2]
      cov[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,3]
      
      cov[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,2]
      
      cov[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,3]
      
      
      
      
    }
  }
  return(cov)
}


log_likelihood_fund=function(ytr,xtr,l1,sigma1,l2,sigma2,l3,sigma3){
  if(length(ytr)>3){
    if(ncol(ytr)==3){
      ytr=c(ytr[,1],ytr[,2],ytr[,3])
    }
  }
  
  
  ntr=ifelse(length(xtr)==9,1,nrow(xtr))
  Ktr=cross_fund_cov_mat3d(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*nrow(xtr))
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-ntr*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}



grad_fund=function(xtr,ytr,l1,sigma1,l2,sigma2,l3,sigma3){
  Trace=function(x){sum(diag(x))}
  x1=x2=xtr
  norm=function(x){sum(x^2)^.5}
  
  K=cross_fund_cov_mat3d(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*nrow(xtr))
  
  K_l1 =  K_l2 = K_l3 = K_sigma1 = K_sigma2 = K_sigma3 =
    matrix(0, ncol =3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)),
           nrow = 3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)))
  
  inv_K=solve(K)
  
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      if(length(as.matrix(x1)) == 9){
        cov1= fund_cov_mat3d(x1[1:3],x2[1:3],l1,sigma1,l2,0,l3,0)*
          # exp(-(norm(x1[1:3])-norm(x2[1:3]))^2/(2*l1^2))*
          (norm(x1[1:3])-norm(x2[1:3]))^2/(l1^3)+
          
          fund_cov_mat3d(x1[4:6],x2[4:6],l1,sigma1,l2,0,l3,0)*
          # exp(-(norm(x1[4:6])-norm(x2[4:6]))^2/(2*l1^2))*
          (norm(x1[4:6])-norm(x2[4:6]))^2/(l1^3)+
          
          fund_cov_mat3d(x1[7:9],x2[7:9],l1,sigma1,l2,0,l3,0)*
          # exp(-(norm(x1[7:9])-norm(x2[7:9]))^2/(2*l1^2))*
          (norm(x1[7:9])-norm(x2[7:9]))^2/(l1^3)
        
        cov2= fund_cov_mat3d(x1[4:6],x2[4:6],l1,0,l2,sigma2,l3,0)*
          # exp(-(norm(x1[4:6])-norm(x2[4:6]))^2/(2*l2^2))*
          (norm(x1[4:6])-norm(x2[4:6]))^2/(l2^3)+
          
          fund_cov_mat3d(x1[1:3],x2[1:3],l1,0,l2,sigma2,l3,0)*
          #  exp(-(norm(x1[1:3])-norm(x2[1:3]))^2/(2*l2^2))*
          (norm(x1[1:3])-norm(x2[1:3]))^2/(l2^3)+
          
          fund_cov_mat3d(x1[7:9],x2[7:9],l1,0,l2,sigma2,l3,0)*
          #  exp(-(norm(x1[7:9])-norm(x2[7:9]))^2/(2*l2^2))*
          (norm(x1[7:9])-norm(x2[7:9]))^2/(l2^3)
        
        cov3= fund_cov_mat3d(x1[7:9],x2[7:9],l1,0,l2,0,l3,sigma3)*
          #  exp(-(norm(x1[4:6])-norm(x2[4:6]))^2/(2*l3^2))*
          (norm(x1[7:9])-norm(x2[7:9]))^2/(l3^3)+
          
          fund_cov_mat3d(x1[1:3],x2[1:3],l1,0,l2,0,l3,sigma3)*
          #   exp(-(norm(x1[1:3])-norm(x2[1:3]))^2/(2*l3^2))*
          (norm(x1[1:3])-norm(x2[1:3]))^2/(l3^3)+
          
          
          fund_cov_mat3d(x1[4:6],x2[4:6],l1,0,l2,0,l3,sigma3)*
          #   exp(-(norm(x1[4:6])-norm(x2[4:6]))^2/(2*l3^2))*
          (norm(x1[4:6])-norm(x2[4:6]))^2/(l3^3)
        
        
        cov4= fund_cov_mat3d(x1[1:3],x2[1:3],l1,1,l2,0,l3,0)*2*sigma1+
          fund_cov_mat3d(x1[4:6],x2[4:6],l1,1,l2,0,l3,0)*2*sigma1+
          fund_cov_mat3d(x1[7:9],x2[7:9],l1,1,l2,0,l3,0)*2*sigma1
        
        
        cov5= fund_cov_mat3d(x1[4:6],x2[4:6],l1,0,l2,1,l3,0)*2*sigma2+
          fund_cov_mat3d(x1[1:3],x2[1:3],l1,0,l2,1,l3,0)*2*sigma2+
          fund_cov_mat3d(x1[7:9],x2[7:9],l1,0,l2,1,l3,0)*2*sigma2
        
        cov6= fund_cov_mat3d(x1[7:9],x2[7:9],l1,0,l2,0,l3,1)*2*sigma3+
          fund_cov_mat3d(x1[1:3],x2[1:3],l1,0,l2,0,l3,1)*2*sigma3+
          fund_cov_mat3d(x1[4:6],x2[4:6],l1,0,l2,0,l3,1)*2*sigma3
      }else{
        cov1= fund_cov_mat3d(x1[i,1:3],x2[j,1:3],l1,sigma1,l2,0,l3,0)*
          #   exp(-(norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(2*l1^2))*
          (norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(l1^3)+
          
          fund_cov_mat3d(x1[i,4:6],x2[j,4:6],l1,sigma1,l2,0,l3,0)*
          #   exp(-(norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(2*l1^2))*
          (norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(l1^3)+
          
          fund_cov_mat3d(x1[i,7:9],x2[j,7:9],l1,sigma1,l2,0,l3,0)*
          #    exp(-(norm(x1[i,7:9])-norm(x2[j,7:9]))^2/(2*l1^2))*
          (norm(x1[i,7:9])-norm(x2[j,7:9]))^2/(l1^3)
        
        cov2= fund_cov_mat3d(x1[i,4:6],x2[j,4:6],l1,0,l2,sigma2,l3,0)*
          #    exp(-(norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(2*l2^2))*
          (norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(l2^3)+
          
          fund_cov_mat3d(x1[i,1:3],x2[j,1:3],l1,0,l2,sigma2,l3,0)*
          #     exp(-(norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(2*l2^2))*
          (norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(l2^3)+
          
          fund_cov_mat3d(x1[i,7:9],x2[j,7:9],l1,0,l2,sigma2,l3,0)*
          #   exp(-(norm(x1[i,7:9])-norm(x2[j,7:9]))^2/(2*l2^2))*
          (norm(x1[i,7:9])-norm(x2[j,7:9]))^2/(l2^3)
        
        cov3= fund_cov_mat3d(x1[i,7:9],x2[j,7:9],l1,0,l2,0,l3,sigma3)*
          #    exp(-(norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(2*l3^2))*
          (norm(x1[i,7:9])-norm(x2[j,7:9]))^2/(l3^3)+
          
          fund_cov_mat3d(x1[i,1:3],x2[j,1:3],l1,0,l2,0,l3,sigma3)*
          #    exp(-(norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(2*l3^2))*
          (norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(l3^3)+
          
          
          fund_cov_mat3d(x1[i,4:6],x2[j,4:6],l1,0,l2,0,l3,sigma3)*
          #    exp(-(norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(2*l3^2))*
          (norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(l3^3)
        
        
        cov4= fund_cov_mat3d(x1[i,1:3],x2[j,1:3],l1,1,l2,0,l3,0)*2*sigma1+
          fund_cov_mat3d(x1[i,4:6],x2[j,4:6],l1,1,l2,0,l3,0)*2*sigma1+
          fund_cov_mat3d(x1[i,7:9],x2[j,7:9],l1,1,l2,0,l3,0)*2*sigma1
        
        
        cov5= fund_cov_mat3d(x1[i,4:6],x2[j,4:6],l1,0,l2,1,l3,0)*2*sigma2+
          fund_cov_mat3d(x1[i,1:3],x2[j,1:3],l1,0,l2,1,l3,0)*2*sigma2+
          fund_cov_mat3d(x1[i,7:9],x2[j,7:9],l1,0,l2,1,l3,0)*2*sigma2
        
        cov6= fund_cov_mat3d(x1[i,7:9],x2[j,7:9],l1,0,l2,0,l3,1)*2*sigma3+
          fund_cov_mat3d(x1[i,1:3],x2[j,1:3],l1,0,l2,0,l3,1)*2*sigma3+
          fund_cov_mat3d(x1[i,4:6],x2[j,4:6],l1,0,l2,0,l3,1)*2*sigma3
        
        
      }
      
      K_l1[i, j] = cov1[1,1]
      K_l1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,2]
      
      K_l1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,3]
      
      K_l1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[2,1]
      K_l1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[3,1]
      
      K_l1[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,2]
      K_l1[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,3]
      
      K_l1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,2]
      
      K_l1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,3]
      
      
      
      K_l2[i, j] = cov2[1,1]
      K_l2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,2]
      
      K_l2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,3]
      
      K_l2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[2,1]
      K_l2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[3,1]
      
      K_l2[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,2]
      K_l2[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,3]
      
      K_l2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,2]
      
      K_l2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,3]
      
      
      K_l3[i, j] = cov3[1,1]
      K_l3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[2,2]
      
      K_l3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[3,3]
      
      K_l3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov3[2,1]
      K_l3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov3[3,1]
      
      K_l3[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[1,2]
      K_l3[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[1,3]
      
      K_l3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[3,2]
      
      K_l3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[2,3]
      
      
      
      K_sigma1[i, j] = cov4[1,1]
      K_sigma1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[2,2]
      
      K_sigma1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[3,3]
      
      K_sigma1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov4[2,1]
      K_sigma1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov4[3,1]
      
      K_sigma1[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[1,2]
      K_sigma1[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[1,3]
      
      K_sigma1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[3,2]
      
      K_sigma1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[2,3]
      
      
      
      K_sigma2[i, j] = cov5[1,1]
      K_sigma2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[2,2]
      
      K_sigma2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[3,3]
      
      K_sigma2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov5[2,1]
      K_sigma2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov5[3,1]
      
      K_sigma2[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[1,2]
      K_sigma2[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[1,3]
      
      K_sigma2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[3,2]
      
      K_sigma2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[2,3]
      
      
      
      K_sigma3[i, j] = cov6[1,1]
      K_sigma3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[2,2]
      
      K_sigma3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[3,3]
      
      K_sigma3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov6[2,1]
      K_sigma3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov6[3,1]
      
      K_sigma3[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[1,2]
      K_sigma3[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[1,3]
      
      K_sigma3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[3,2]
      
      K_sigma3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[2,3]
      
    }
  }
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l1%*%inv_K%*%c(ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l1),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma1%*%inv_K%*%c(ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_sigma1),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l2%*%inv_K%*%(c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_l2),
             
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma2%*%inv_K%*%(c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_sigma2),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l3%*%inv_K%*%c(ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l3),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma3%*%inv_K%*%(c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_sigma3)
  )
  
  return(grad)
  
  
}



log_likelihood_standard=function(ytr,xtr,l1,sigma1,l2,sigma2,l3,sigma3){
  if(length(ytr)>3){
    if(ncol(ytr)==3){
      ytr=c(ytr[,1],ytr[,2],ytr[,3])
    }
  }
  
  
  ntr=ifelse(length(xtr)==9,1,nrow(xtr))
  Ktr=cross_cov_mat3d(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*nrow(xtr))
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-ntr*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}



grad_standard=function(xtr,ytr,l1,sigma1,l2,sigma2,l3,sigma3){
  Trace=function(x){sum(diag(x))}
  x1=x2=xtr
  norm=function(x){sum(x^2)^.5}
  
  K=cross_cov_mat3d(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*nrow(xtr))
  
  K_l1 =  K_l2 = K_l3 = K_sigma1 = K_sigma2 = K_sigma3 =
    matrix(0, ncol =3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)),
           nrow = 3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)))
  
  inv_K=solve(K)
  
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      if(length(as.matrix(x1)) == 9){
        cov1= cov_mat3d(x1[1:3],x2[1:3],l1,sigma1,l2,0,l3,0)*
          # exp(-(norm(x1[1:3])-norm(x2[1:3]))^2/(2*l1^2))*
          (norm(x1[1:3])-norm(x2[1:3]))^2/(l1^3)+
          
          cov_mat3d(x1[4:6],x2[4:6],l1,sigma1,l2,0,l3,0)*
          # exp(-(norm(x1[4:6])-norm(x2[4:6]))^2/(2*l1^2))*
          (norm(x1[4:6])-norm(x2[4:6]))^2/(l1^3)+
          
          cov_mat3d(x1[7:9],x2[7:9],l1,sigma1,l2,0,l3,0)*
          # exp(-(norm(x1[7:9])-norm(x2[7:9]))^2/(2*l1^2))*
          (norm(x1[7:9])-norm(x2[7:9]))^2/(l1^3)
        
        cov2= cov_mat3d(x1[4:6],x2[4:6],l1,0,l2,sigma2,l3,0)*
          # exp(-(norm(x1[4:6])-norm(x2[4:6]))^2/(2*l2^2))*
          (norm(x1[4:6])-norm(x2[4:6]))^2/(l2^3)+
          
          cov_mat3d(x1[1:3],x2[1:3],l1,0,l2,sigma2,l3,0)*
          #  exp(-(norm(x1[1:3])-norm(x2[1:3]))^2/(2*l2^2))*
          (norm(x1[1:3])-norm(x2[1:3]))^2/(l2^3)+
          
          cov_mat3d(x1[7:9],x2[7:9],l1,0,l2,sigma2,l3,0)*
          #  exp(-(norm(x1[7:9])-norm(x2[7:9]))^2/(2*l2^2))*
          (norm(x1[7:9])-norm(x2[7:9]))^2/(l2^3)
        
        cov3= cov_mat3d(x1[7:9],x2[7:9],l1,0,l2,0,l3,sigma3)*
          #  exp(-(norm(x1[4:6])-norm(x2[4:6]))^2/(2*l3^2))*
          (norm(x1[7:9])-norm(x2[7:9]))^2/(l3^3)+
          
          cov_mat3d(x1[1:3],x2[1:3],l1,0,l2,0,l3,sigma3)*
          #   exp(-(norm(x1[1:3])-norm(x2[1:3]))^2/(2*l3^2))*
          (norm(x1[1:3])-norm(x2[1:3]))^2/(l3^3)+
          
          
          cov_mat3d(x1[4:6],x2[4:6],l1,0,l2,0,l3,sigma3)*
          #   exp(-(norm(x1[4:6])-norm(x2[4:6]))^2/(2*l3^2))*
          (norm(x1[4:6])-norm(x2[4:6]))^2/(l3^3)
        
        
        cov4= cov_mat3d(x1[1:3],x2[1:3],l1,1,l2,0,l3,0)*2*sigma1+
          cov_mat3d(x1[4:6],x2[4:6],l1,1,l2,0,l3,0)*2*sigma1+
          cov_mat3d(x1[7:9],x2[7:9],l1,1,l2,0,l3,0)*2*sigma1
        
        
        cov5= cov_mat3d(x1[4:6],x2[4:6],l1,0,l2,1,l3,0)*2*sigma2+
          cov_mat3d(x1[1:3],x2[1:3],l1,0,l2,1,l3,0)*2*sigma2+
          cov_mat3d(x1[7:9],x2[7:9],l1,0,l2,1,l3,0)*2*sigma2
        
        cov6= cov_mat3d(x1[7:9],x2[7:9],l1,0,l2,0,l3,1)*2*sigma3+
          cov_mat3d(x1[1:3],x2[1:3],l1,0,l2,0,l3,1)*2*sigma3+
          cov_mat3d(x1[4:6],x2[4:6],l1,0,l2,0,l3,1)*2*sigma3
      }else{
        cov1= cov_mat3d(x1[i,1:3],x2[j,1:3],l1,sigma1,l2,0,l3,0)*
          #   exp(-(norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(2*l1^2))*
          (norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(l1^3)+
          
          cov_mat3d(x1[i,4:6],x2[j,4:6],l1,sigma1,l2,0,l3,0)*
          #   exp(-(norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(2*l1^2))*
          (norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(l1^3)+
          
          cov_mat3d(x1[i,7:9],x2[j,7:9],l1,sigma1,l2,0,l3,0)*
          #    exp(-(norm(x1[i,7:9])-norm(x2[j,7:9]))^2/(2*l1^2))*
          (norm(x1[i,7:9])-norm(x2[j,7:9]))^2/(l1^3)
        
        cov2= cov_mat3d(x1[i,4:6],x2[j,4:6],l1,0,l2,sigma2,l3,0)*
          #    exp(-(norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(2*l2^2))*
          (norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(l2^3)+
          
          cov_mat3d(x1[i,1:3],x2[j,1:3],l1,0,l2,sigma2,l3,0)*
          #     exp(-(norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(2*l2^2))*
          (norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(l2^3)+
          
          cov_mat3d(x1[i,7:9],x2[j,7:9],l1,0,l2,sigma2,l3,0)*
          #   exp(-(norm(x1[i,7:9])-norm(x2[j,7:9]))^2/(2*l2^2))*
          (norm(x1[i,7:9])-norm(x2[j,7:9]))^2/(l2^3)
        
        cov3= cov_mat3d(x1[i,7:9],x2[j,7:9],l1,0,l2,0,l3,sigma3)*
          #    exp(-(norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(2*l3^2))*
          (norm(x1[i,7:9])-norm(x2[j,7:9]))^2/(l3^3)+
          
          cov_mat3d(x1[i,1:3],x2[j,1:3],l1,0,l2,0,l3,sigma3)*
          #    exp(-(norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(2*l3^2))*
          (norm(x1[i,1:3])-norm(x2[j,1:3]))^2/(l3^3)+
          
          
          cov_mat3d(x1[i,4:6],x2[j,4:6],l1,0,l2,0,l3,sigma3)*
          #    exp(-(norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(2*l3^2))*
          (norm(x1[i,4:6])-norm(x2[j,4:6]))^2/(l3^3)
        
        
        cov4= cov_mat3d(x1[i,1:3],x2[j,1:3],l1,1,l2,0,l3,0)*2*sigma1+
          cov_mat3d(x1[i,4:6],x2[j,4:6],l1,1,l2,0,l3,0)*2*sigma1+
          cov_mat3d(x1[i,7:9],x2[j,7:9],l1,1,l2,0,l3,0)*2*sigma1
        
        
        cov5= cov_mat3d(x1[i,4:6],x2[j,4:6],l1,0,l2,1,l3,0)*2*sigma2+
          cov_mat3d(x1[i,1:3],x2[j,1:3],l1,0,l2,1,l3,0)*2*sigma2+
          cov_mat3d(x1[i,7:9],x2[j,7:9],l1,0,l2,1,l3,0)*2*sigma2
        
        cov6= cov_mat3d(x1[i,7:9],x2[j,7:9],l1,0,l2,0,l3,1)*2*sigma3+
          cov_mat3d(x1[i,1:3],x2[j,1:3],l1,0,l2,0,l3,1)*2*sigma3+
          cov_mat3d(x1[i,4:6],x2[j,4:6],l1,0,l2,0,l3,1)*2*sigma3
        
        
      }
      
      K_l1[i, j] = cov1[1,1]
      K_l1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,2]
      
      K_l1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,3]
      
      K_l1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[2,1]
      K_l1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[3,1]
      
      K_l1[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,2]
      K_l1[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,3]
      
      K_l1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,2]
      
      K_l1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,3]
      
      
      
      K_l2[i, j] = cov2[1,1]
      K_l2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,2]
      
      K_l2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,3]
      
      K_l2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[2,1]
      K_l2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[3,1]
      
      K_l2[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,2]
      K_l2[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,3]
      
      K_l2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,2]
      
      K_l2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,3]
      
      
      K_l3[i, j] = cov3[1,1]
      K_l3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[2,2]
      
      K_l3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[3,3]
      
      K_l3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov3[2,1]
      K_l3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov3[3,1]
      
      K_l3[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[1,2]
      K_l3[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[1,3]
      
      K_l3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[3,2]
      
      K_l3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[2,3]
      
      
      
      K_sigma1[i, j] = cov4[1,1]
      K_sigma1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[2,2]
      
      K_sigma1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[3,3]
      
      K_sigma1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov4[2,1]
      K_sigma1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov4[3,1]
      
      K_sigma1[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[1,2]
      K_sigma1[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[1,3]
      
      K_sigma1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[3,2]
      
      K_sigma1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[2,3]
      
      
      
      K_sigma2[i, j] = cov5[1,1]
      K_sigma2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[2,2]
      
      K_sigma2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[3,3]
      
      K_sigma2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov5[2,1]
      K_sigma2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov5[3,1]
      
      K_sigma2[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[1,2]
      K_sigma2[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[1,3]
      
      K_sigma2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[3,2]
      
      K_sigma2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[2,3]
      
      
      
      K_sigma3[i, j] = cov6[1,1]
      K_sigma3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[2,2]
      
      K_sigma3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[3,3]
      
      K_sigma3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov6[2,1]
      K_sigma3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov6[3,1]
      
      K_sigma3[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[1,2]
      K_sigma3[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[1,3]
      
      K_sigma3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[3,2]
      
      K_sigma3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[2,3]
      
    }
  }
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l1%*%inv_K%*%c(ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l1),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma1%*%inv_K%*%c(ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_sigma1),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l2%*%inv_K%*%(c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_l2),
             
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma2%*%inv_K%*%(c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_sigma2),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l3%*%inv_K%*%c(ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l3),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma3%*%inv_K%*%(c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_sigma3)
  )
  
  return(grad)
  
  
}




log_likelihood_standard2=function(ytr,xtr,l1,sigma1,l2,sigma2,l3,sigma3){
  if(length(ytr)>3){
    if(ncol(ytr)==3){
      ytr=c(ytr[,1],ytr[,2],ytr[,3])
    }
  }
  
  
  ntr=ifelse(length(xtr)==9,1,nrow(xtr))
  Ktr=cross_cov_mat3d_2(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*nrow(xtr))
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-ntr*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}



grad_standard2=function(xtr,ytr,l1,sigma1,l2,sigma2,l3,sigma3){
  Trace=function(x){sum(diag(x))}
  x1=x2=xtr
  norm=function(x){sum(x^2)^.5}
  
  K=cross_cov_mat3d_2(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*nrow(xtr))
  
  K_l1 =  K_l2 = K_l3 = K_sigma1 = K_sigma2 = K_sigma3 =
    matrix(0, ncol =3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)),
           nrow = 3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)))
  
  inv_K=solve(K)
  
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      if(length(as.matrix(x1)) == 9){
        cov1= cov_mat3d_2(x1,x2,l1,sigma1,l2,0,l3,0)*
          
          (norm(x1)-norm(x2))^2/(l1^3)
        
        
        cov2= cov_mat3d_2(x1,x2,l1,0,l2,sigma2,l3,0)*
          # exp(-(norm(x1)-norm(x2))^2/(2*l2^2))*
          (norm(x1)-norm(x2))^2/(l2^3)
        
        cov3= cov_mat3d_2(x1,x2,l1,0,l2,0,l3,sigma3)*
          #  exp(-(norm(x1)-norm(x2))^2/(2*l3^2))*
          (norm(x1)-norm(x2))^2/(l3^3)
        
        cov4= cov_mat3d_2(x1,x2,l1,1,l2,0,l3,0)*2*sigma1
        
        
        cov5= cov_mat3d_2(x1,x2,l1,0,l2,1,l3,0)*2*sigma2
        
        cov6= cov_mat3d_2(x1,x2,l1,0,l2,0,l3,1)*2*sigma3
      }else{
        cov1= cov_mat3d_2(x1[i,],x2[j,],l1,sigma1,l2,0,l3,0)*
          #   exp(-(norm(x1[i,])-norm(x2[j,]))^2/(2*l1^2))*
          (norm(x1[i,])-norm(x2[j,]))^2/(l1^3)
        
        cov2= cov_mat3d_2(x1[i,],x2[j,],l1,0,l2,sigma2,l3,0)*
          #    exp(-(norm(x1[i,])-norm(x2[j,]))^2/(2*l2^2))*
          (norm(x1[i,])-norm(x2[j,]))^2/(l2^3)
        
        cov3= cov_mat3d_2(x1[i,],x2[j,],l1,0,l2,0,l3,sigma3)*
          #    exp(-(norm(x1[i,])-norm(x2[j,]))^2/(2*l3^2))*
          (norm(x1[i,])-norm(x2[j,]))^2/(l3^3)
        
        
        cov4= cov_mat3d_2(x1[i,],x2[j,],l1,1,l2,0,l3,0)*2*sigma1
        
        
        cov5= cov_mat3d_2(x1[i,],x2[j,],l1,0,l2,1,l3,0)*2*sigma2
        
        cov6= cov_mat3d_2(x1[i,],x2[j,],l1,0,l2,0,l3,1)*2*sigma3
        
      }
      
      K_l1[i, j] = cov1[1,1]
      K_l1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,2]
      
      K_l1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,3]
      
      K_l1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[2,1]
      K_l1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[3,1]
      
      K_l1[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,2]
      K_l1[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,3]
      
      K_l1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,2]
      
      K_l1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,3]
      
      
      
      K_l2[i, j] = cov2[1,1]
      K_l2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,2]
      
      K_l2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,3]
      
      K_l2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[2,1]
      K_l2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[3,1]
      
      K_l2[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,2]
      K_l2[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,3]
      
      K_l2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,2]
      
      K_l2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,3]
      
      
      K_l3[i, j] = cov3[1,1]
      K_l3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[2,2]
      
      K_l3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[3,3]
      
      K_l3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov3[2,1]
      K_l3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov3[3,1]
      
      K_l3[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[1,2]
      K_l3[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[1,3]
      
      K_l3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[3,2]
      
      K_l3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
           2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov3[2,3]
      
      
      
      K_sigma1[i, j] = cov4[1,1]
      K_sigma1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[2,2]
      
      K_sigma1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[3,3]
      
      K_sigma1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov4[2,1]
      K_sigma1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov4[3,1]
      
      K_sigma1[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[1,2]
      K_sigma1[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[1,3]
      
      K_sigma1[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[3,2]
      
      K_sigma1[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov4[2,3]
      
      
      
      K_sigma2[i, j] = cov5[1,1]
      K_sigma2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[2,2]
      
      K_sigma2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[3,3]
      
      K_sigma2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov5[2,1]
      K_sigma2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov5[3,1]
      
      K_sigma2[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[1,2]
      K_sigma2[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[1,3]
      
      K_sigma2[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[3,2]
      
      K_sigma2[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov5[2,3]
      
      
      
      K_sigma3[i, j] = cov6[1,1]
      K_sigma3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[2,2]
      
      K_sigma3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[3,3]
      
      K_sigma3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov6[2,1]
      K_sigma3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov6[3,1]
      
      K_sigma3[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[1,2]
      K_sigma3[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[1,3]
      
      K_sigma3[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[3,2]
      
      K_sigma3[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
               2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov6[2,3]
      
    }
  }
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l1%*%inv_K%*%c(ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l1),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma1%*%inv_K%*%c(ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_sigma1),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l2%*%inv_K%*%(c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_l2),
             
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma2%*%inv_K%*%(c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_sigma2),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l3%*%inv_K%*%c(ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l3),
             -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma3%*%inv_K%*%(c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_sigma3)
  )
  
  return(grad)
  
  
}



data <- readLines("results/original_struct.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

original_struct=data.frame(Y=Y,X=X)
# Print or use the matrices as needed

data <- readLines("results/original_struct.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

original_struct=data.frame(Y=Y,X=X)

data <- readLines("results/rot_0.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

rot_0=data.frame(Y=Y,X=X)



data <- readLines("results/rot_1.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

rot_1=data.frame(Y=Y,X=X)


data <- readLines("results/rot_2.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

X[1,]
Y[1,]
Y[331,]
X[331,]

rot_2=data.frame(Y=Y,X=X)



data <- readLines("results/rot_3.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

X[1,]
Y[1,]
Y[331,]
X[331,]

rot_4=data.frame(Y=Y,X=X)

data <- readLines("results/rot_4.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

X[1,]
Y[1,]
Y[331,]
X[331,]

rot_4=data.frame(Y=Y,X=X)


data <- readLines("results/rot_5.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

X[1,]
Y[1,]
Y[331,]
X[331,]

rot_5=data.frame(Y=Y,X=X)



data <- readLines("results/rot_6.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

X[1,]
Y[1,]
Y[331,]
X[331,]

rot_6=data.frame(Y=Y,X=X)



data <- readLines("results/rot_7.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

X[1,]
Y[1,]
Y[331,]
X[331,]

rot_7=data.frame(Y=Y,X=X)


data <- readLines("results/rot_8.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

X[1,]
Y[1,]
Y[331,]
X[331,]

rot_8=data.frame(Y=Y,X=X)


data <- readLines("results/rot_9.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)


rot_9=data.frame(Y=Y,X=X)



data <- readLines("results/rot_3.xyz")

X_DIPOLE_MP2 <- c()
Y_DIPOLE_MP2 <- c()
Z_DIPOLE_MP2 <- c()
X_O <- c()
Y_O <- c()
Z_O <- c()
X_H1 <- c()
Y_H1 <- c()
Z_H1 <- c()
X_H2 <- c()
Y_H2 <- c()
Z_H2 <- c()

for (i in 1:length(data)) {
  
  if (grepl("X_DIPOLE_MP2", data[i])) {
    X_DIPOLE_MP2 <- append(X_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i], "=| "))[4]))
    Y_DIPOLE_MP2 <- append(Y_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+1], "=| "))[4]))
    Z_DIPOLE_MP2 <- append(Z_DIPOLE_MP2, as.numeric(unlist(strsplit(data[i+2], "=| "))[4]))
    
    loc1=as.numeric(unlist(strsplit(data[i+3], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+4], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+5], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)



rot_3=data.frame(Y=Y,X=X)

X_names=c("X.X_O","X.Y_O","X.Z_O","X.X_H1","X.Y_H1","X.Z_H1",
          "X.X_H2","X.Y_H2","X.Z_H2")
Y_names=c("Y.X_DIPOLE_MP2","Y.Y_DIPOLE_MP2","Y.Z_DIPOLE_MP2")


X=as.matrix(original_struct[X_names])
Y=as.matrix(original_struct[Y_names])
for(dat in list(rot_0,rot_1,rot_2,rot_3,rot_4,rot_5,rot_6,rot_7,rot_8,rot_9)){
  X=rbind(X,as.matrix(dat[X_names]))
  Y=rbind(Y,as.matrix(dat[Y_names]))
}
Xtr=X[ind_train,]
Ytr=Y[ind_train,]
Xte=X[ind_test,]
Yte=Y[ind_test,]


opt_par=Adam(function(x) {
  log_likelihood_fund(
    Ytr,Xtr, x[1], x[2], x[3], x[4], x[5],x[6])
},function(x) {
  
  g=grad_fund(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
              x[6])
  
  
},rep(1,6),lr=.1,max_iter=20)
#opt_par=rep(1,6)

Ktetr_eq=cross_fund_cov_mat3d(Xte,Xtr,opt_par[1],
                              opt_par[2],
                              opt_par[3],
                              opt_par[4],
                              opt_par[5],
                              opt_par[6])

Ktetr_trtr_inv_eq=Ktetr_eq%*%solve(cross_fund_cov_mat3d(Xtr,Xtr,opt_par[1],
                                                        opt_par[2],
                                                        opt_par[3],
                                                        opt_par[4],
                                                        opt_par[5],
                                                        opt_par[6])+diag(jit,nrow=3*nrow(Xtr)))



posterior_mean_eq=Ktetr_trtr_inv_eq%*%c(Ytr[,1],Ytr[,2],Ytr[,3])


(rmse_eq=mean(apply((cbind(posterior_mean_eq[1:(length(posterior_mean_eq)/3)],
                           posterior_mean_eq[(length(posterior_mean_eq)/3+1):
                                               (2*length(posterior_mean_eq)/3)],
                           posterior_mean_eq[(2*length(posterior_mean_eq)/3+1):
                                               (length(posterior_mean_eq))])-Yte)^2,1,sum)
)^.5)

posterior_cov_eq= cross_fund_cov_mat3d(Xte,Xte,opt_par[1],
                                       opt_par[2],
                                       opt_par[3],
                                       opt_par[4],opt_par[5],
                                       opt_par[6])-
  Ktetr_trtr_inv_eq%*%t(Ktetr_eq)


(LogS_eq=(t(c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean_eq)%*%
            solve(posterior_cov_eq
                  +diag(1e-6,3*nrow(Xte))
            )%*%
            (c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean_eq)
          +determinant(posterior_cov_eq
                       +diag(1e-6,3*nrow(Xte))
          )$modulus[1]+
            log(2*pi)*nrow(Xte))/nrow(Xte))
rm(posterior_cov_eq)

opt_par=Adam(function(x) {
  log_likelihood_standard(
    Ytr,Xtr, x[1], x[2], x[3], x[4], x[5],x[6])
},function(x) {
  grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                x[6])
  
},rep(1,6),lr=.1,max_iter=30)

Ktetr=cross_cov_mat3d(Xte,Xtr,opt_par[1],
                      opt_par[2],
                      opt_par[3],
                      opt_par[4],
                      opt_par[5],
                      opt_par[6])

Ktetr_trtr_inv=Ktetr%*%solve(cross_cov_mat3d(Xtr,Xtr,opt_par[1],
                                             opt_par[2],
                                             opt_par[3],
                                             opt_par[4],
                                             opt_par[5],
                                             opt_par[6])+diag(jit,nrow=3*nrow(Xtr)))



posterior_mean=Ktetr_trtr_inv%*%c(Ytr[,1],Ytr[,2],Ytr[,3])


(rmse=mean(apply((cbind(posterior_mean[1:(length(posterior_mean)/3)],
                        posterior_mean[(length(posterior_mean)/3+1):
                                         (2*length(posterior_mean)/3)],
                        posterior_mean[(2*length(posterior_mean)/3+1):
                                         (length(posterior_mean))])-Yte)^2,1,sum)
)^.5)

posterior_cov= cross_cov_mat3d(Xte,Xte,opt_par[1],
                               opt_par[2],
                               opt_par[3],
                               opt_par[4],opt_par[5],
                               opt_par[6])-
  Ktetr_trtr_inv%*%t(Ktetr)


(LogS=(t(c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean)%*%
         solve(posterior_cov
               +diag(1e-6,3*nrow(Xte))
         )%*%
         (c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean)
       +determinant(posterior_cov
                    +diag(1e-6,3*nrow(Xte))
       )$modulus[1]+
         log(2*pi)*nrow(Xte))/nrow(Xte))
rm(posterior_cov)


opt_par=Adam(function(x) {
  log_likelihood_standard2(
    Ytr,Xtr, x[1], x[2], x[3], x[4], x[5],x[6])
},function(x) {
  grad_standard2(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                 x[6])
},rep(1,6),lr=.1,max_iter=30)

Ktetr2=cross_cov_mat3d_2(Xte,Xtr,opt_par[1],
                         opt_par[2],
                         opt_par[3],
                         opt_par[4],
                         opt_par[5],
                         opt_par[6])

Ktetr_trtr_inv2=Ktetr2%*%solve(cross_cov_mat3d_2(Xtr,Xtr,opt_par[1],
                                                 opt_par[2],
                                                 opt_par[3],
                                                 opt_par[4],
                                                 opt_par[5],
                                                 opt_par[6])+diag(jit,nrow=3*nrow(Xtr)))



posterior_mean2=Ktetr_trtr_inv2%*%c(Ytr[,1],Ytr[,2],Ytr[,3])


(rmse2=mean(apply((cbind(posterior_mean2[1:(length(posterior_mean)/3)],
                         posterior_mean2[(length(posterior_mean)/3+1):
                                           (2*length(posterior_mean)/3)],
                         posterior_mean2[(2*length(posterior_mean)/3+1):
                                           (length(posterior_mean))])-Yte)^2,1,sum)
)^.5)


posterior_cov2= cross_cov_mat3d_2(Xte,Xte,opt_par[1],
                                  opt_par[2],
                                  opt_par[3],
                                  opt_par[4],opt_par[5],
                                  opt_par[6])-
  Ktetr_trtr_inv2%*%t(Ktetr2)


(LogS2=(t(c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean2)%*%
          solve(posterior_cov2
                +diag(1e-6,3*nrow(Xte))
          )%*%
          (c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean2)
        +determinant(posterior_cov2
                     +diag(1e-6,3*nrow(Xte))
        )$modulus[1]+
          log(2*pi)*nrow(Xte))/nrow(Xte))

rm(posterior_cov2)


res=rbind(c(rmse,rmse2,rmse_eq),c(LogS,LogS2,LogS_eq))

save(list = "res", file = paste0("Chemical_data_scores_mix_sample_Adam", id, ".rda"))





