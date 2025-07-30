id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
library(pracma)
library(cubature)
library(doParallel)
library(foreach)
library(progress)
library(future.apply)
plan(multicore, workers = parallel::detectCores() - 1)
#notcluster=is.na(id)

notcluster=1
quant=quant_fund=.9

Packages=c('doParallel','foreach')
Export=c("Packages","cov_mat","cross_cov_mat","grad_standard","fund_cov_mat",
         "cross_fund_cov_mat","grad_fund","Xtr","Ytr","cross_fund_cov_mat","jit","inv")
Xtr=Ytr=c()
data_set_size=10000
test_set_size=1000


lr=0.001;max_iter=100
jit=1e-8

n_ind=10
n_ind_points=seq(20,100,l=n_ind)

jit=1e-8


load('Data/data_NMF.rda')
X=as.matrix(unname(data[[1]]))
Y=as.matrix(unname(data[[2]]))
# Display the resulting matrix
meanY=apply(Y,2,mean)

Y=t(apply(Y,1,function(x){x-meanY}))


ind=sample(nrow(X),data_set_size)

ind_test=sample((1:nrow(X))[-ind],test_set_size)
X=X[ind,];Y=Y[ind,]
Xte=as.matrix(unname(data[[1]]))[ind_test,]
Yte=as.matrix(unname(data[[2]]))[ind_test,]


######## REMOVE DUPLICATES

projection_fund=function(x1){
  norm=function(x){sum(x^2)^.5}
  X1=x1-rep(x1[1:3],9)
  
  x_H_bar_1=X1[4:6];x_H_bar_2=X1[7:9]
  r1=norm(x_H_bar_1)
  
  theta1=atan2(x_H_bar_1[1],x_H_bar_1[2])
  
  x_H_tilde_1=c(0,sin(theta1)*x_H_bar_1[1]+cos(theta1)*x_H_bar_1[2],x_H_bar_1[3])
  phi1=-atan2(x_H_tilde_1[3],x_H_tilde_1[2])
  
  psi_1=matrix(c(cos(theta1),-sin(theta1),0,sin(theta1),cos(theta1),0,0,0,1),
               nrow=3,byrow = 1)
  psi_2=matrix(c(1,0,0,0,cos(phi1),-sin(phi1),0,sin(phi1),cos(phi1)),
               nrow=3,byrow = 1)
  
  x_H_tilde_2=psi_2%*%psi_1%*%x_H_bar_2
  
  phi2=-atan2(x_H_tilde_2[3],x_H_tilde_2[1])
  psi_3=matrix(c(cos(phi2),0,-sin(phi2),0,1,0,sin(phi2),0,cos(phi2)),
               nrow=3,byrow = 1)
  X1=X1[-(1:9)]
  rho=psi_3%*%psi_2%*%psi_1
  
  A1=c(r1,(psi_3%*%x_H_tilde_2)[1],(psi_3%*%x_H_tilde_2)[2],
       foreach(i=1:6,.combine='c')%dopar%{rho%*%(X1[(3*(i-1)+1):(3*i)])})
  
  return(A1)
}

#A=t(apply(X,1,projection_fund))
#A=matrix(rnorm(3*3641,mean(A),sd(A)),ncol=3)

rots=matrix(rnorm(nrow(Y)*3,0,2*pi),ncol=3)
X=t(sapply(1:nrow(rots),function(x){
  
  theta1=rots[x,1];phi1=rots[x,2];phi2=rots[x,3]
  psi_1=matrix(c(cos(theta1),-sin(theta1),0,sin(theta1),cos(theta1),0,0,0,1),
               nrow=3,byrow = 1)
  psi_2=matrix(c(1,0,0,0,cos(phi1),-sin(phi1),0,sin(phi1),cos(phi1)),
               nrow=3,byrow = 1)
  psi_3=matrix(c(cos(phi2),0,-sin(phi2),0,1,0,sin(phi2),0,cos(phi2)),
               nrow=3,byrow = 1)
  as.numeric(psi_3%*%psi_2%*%psi_1%*%matrix(X[x,],ncol=9))
}))


Y=t(sapply(1:nrow(rots),function(x){
  theta1=rots[x,1];phi1=rots[x,2];phi2=rots[x,3]
  psi_1=matrix(c(cos(theta1),-sin(theta1),0,sin(theta1),cos(theta1),0,0,0,1),
               nrow=3,byrow = 1)
  psi_2=matrix(c(1,0,0,0,cos(phi1),-sin(phi1),0,sin(phi1),cos(phi1)),
               nrow=3,byrow = 1)
  psi_3=matrix(c(cos(phi2),0,-sin(phi2),0,1,0,sin(phi2),0,cos(phi2)),
               nrow=3,byrow = 1)
  as.numeric(psi_3%*%psi_2%*%psi_1%*%Y[x,])
}))


####################################################################


fund_cov_mat=function(x1, x2, l, sigma) {
  norm=function(x) sqrt(sum(x^2))
  
  transform_point=function(x) {
    X=x - rep(x[1:3], 9)
    x_H1=X[4:6]
    x_H2=X[7:9]
    r=norm(x_H1)
    
    theta=atan2(x_H1[1], x_H1[2])
    x_H_tilde1=c(0, sin(theta)*x_H1[1] + cos(theta)*x_H1[2], x_H1[3])
    phi1=-atan2(x_H_tilde1[3], x_H_tilde1[2])
    
    psi1=matrix(c(cos(theta), -sin(theta), 0,
                     sin(theta),  cos(theta), 0,
                     0, 0, 1), nrow = 3, byrow = TRUE)
    
    psi2=matrix(c(1, 0, 0,
                     0, cos(phi1), -sin(phi1),
                     0, sin(phi1),  cos(phi1)), nrow = 3, byrow = TRUE)
    
    x_H_tilde2=psi2 %*% psi1 %*% x_H2
    phi2=-atan2(x_H_tilde2[3], x_H_tilde2[1])
    
    psi3=matrix(c(cos(phi2), 0, -sin(phi2),
                     0, 1, 0,
                     sin(phi2), 0, cos(phi2)), nrow = 3, byrow = TRUE)
    
    rho=psi3 %*% psi2 %*% psi1
    X_remain=X[-(1:9)]
    
    A=numeric(3 + 3 * 6)
    A[1]=r
    x_H_rot=psi3 %*% x_H_tilde2
    A[2:3]=x_H_rot[1:2]
    
    for (i in 1:6) {
      idx=(3 * (i - 1) + 1):(3 * i)
      A[3 + (3 * (i - 1) + 1):(3 * i)]=rho %*% X_remain[idx]
    }
    
    list(A = A, rho = rho)
  }
  
  out1=transform_point(x1)
  out2=transform_point(x2)
  
  diff=out1$A - out2$A
  K_scalar=sigma^2 * exp(-sum(diff^2) / (2 * l^2))
  return(K_scalar * t(out1$rho) %*% out2$rho)
}



cross_fund_cov_mat=function(x1, x2, l, sigma) {
  norm=function(x) sqrt(sum(x^2))
  
  n1=nrow(x1)
  n2=nrow(x2)
  
  cov=matrix(0, nrow = 3 * n1, ncol = 3 * n2)
  
  blocks=future_lapply(1:n1, function(i) {
    lapply(1:n2, function(j) {
      fund_cov_mat(x1[i, ], x2[j, ], l, sigma)
    })
  })
  
 
  for (i in 1:n1) {
    for (j in 1:n2) {
      cov1=blocks[[i]][[j]]  # 3x3 matrix
      
      row_idx1=i
      row_idx2=n1 + i
      row_idx3=2 * n1 + i
      
      col_idx1=j
      col_idx2=n2 + j
      col_idx3=2 * n2 + j
      
      # Assign diagonal elements
      cov[row_idx1, col_idx1]=cov1[1, 1]
      cov[row_idx2, col_idx2]=cov1[2, 2]
      cov[row_idx3, col_idx3]=cov1[3, 3]
      
      # Assign off-diagonal elements (following original indexing)
      cov[row_idx2, col_idx1]=cov1[2, 1]
      cov[row_idx3, col_idx1]=cov1[3, 1]
      
      cov[row_idx1, col_idx2]=cov1[1, 2]
      cov[row_idx1, col_idx3]=cov1[1, 3]
      
      cov[row_idx3, col_idx2]=cov1[3, 2]
      cov[row_idx2, col_idx3]=cov1[2, 3]
    }
  }
  
  return(cov)
}


log_likelihood_fund=function(ytr,xtr,l,sigma){
  if(length(ytr)>3){
    if(ncol(ytr)==3){
      ytr=c(ytr[,1],ytr[,2],ytr[,3])
    }
  }
  
  
  ntr=ifelse(length(xtr)==27,1,nrow(xtr))
  Ktr=cross_fund_cov_mat(xtr,xtr,l,sigma)+
    diag(jit,nrow=3*ntr)
  
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*log(det(Ktr))-ntr*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}






# Adam optimizer


Adam = function(f, grad_f, params_init, lr = 0.01, beta1 = 0.9, beta2 = 0.999,
                epsilon = 1e-8, max_iter = 100, tol = 1e-10,
                lb = rep(0, length(params_init)),
                ub = rep(Inf, length(params_init))) {
  
  pb = progress_bar$new(
    format = "  Progress [:bar] :percent in :elapsed",
    total = max_iter,
    clear = FALSE,
    width = 60
  )
  
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
    
    
    pb$tick()
  }
  
  return(params)
}

cov_mat=function(x1,x2,l,sigma){
  norm=function(x){sum(x^2)^.5}
  
  diag(c(sigma^2*exp(-norm(x1-x2)^2/(2*l^2)),
         sigma^2*exp(-norm(x1-x2)^2/(2*l^2)),
         sigma^2*exp(-norm(x1-x2)^2/(2*l^2))))
}



cross_cov_mat=function(x1,x2,l,sigma){
  cov= matrix(0, ncol =3 * ifelse(length(as.matrix(x2)) == 27, 1,nrow(x2)),
              nrow = 3 * ifelse(length(as.matrix(x1)) == 27, 
                                1, nrow(x1)))
  for(i in 1:(length(as.matrix(x1))/27)){
    for (j in 1:(length(as.matrix(x2))/27)){
      
      
      
      if (length(as.matrix(x2)) == 27) {
        if (length(as.matrix(x1)) == 27) {
          cov1 <- cov_mat(x1, x2, l, sigma)
        } else {
          cov1 <- cov_mat(x1[i, ], x2, l, sigma)
        }
      } else if (length(as.matrix(x1)) == 27) {
        cov1 <- cov_mat(x1, x2[j, ], l, sigma)
      } else {
        cov1 <-cov_mat(x1[i, ], x2[j, ], l, sigma)
      }
      
      cov[i, j] = cov1[1,1]
      cov[ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[2,2]
      
      cov[2*ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
          2*ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[3,3]
      
      cov[ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i, j] = cov1[2,1]
      cov[2*ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i, j] = cov1[3,1]
      
      cov[i,ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[1,2]
      cov[i,2*ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[1,3]
      
      cov[2*ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[3,2]
      
      cov[ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
          2*ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[2,3]
      
      
      
      
    }
  }
  return(cov)
}





grad_fund=function(xtr,ytr,l,sigma){
  Trace=function(x){sum(diag(x))}
    x1=x2=xtr
  norm=function(x){sum(x^2)^.5}
  ntr=ifelse(is.null(dim(xtr)),1,nrow(xtr))
  
  K=cross_fund_cov_mat(xtr,xtr,l,sigma)+
    diag(jit,nrow=3*ntr)
  
  K_l=K_sigma=
    matrix(0, ncol =3 * ifelse(length(as.matrix(xtr)) == 27, 1, ntr),
           nrow = 3 * ifelse(length(as.matrix(xtr)) == 27, 1, ntr))
  
  inv_K=inv(K)
  
  N1=N2=ifelse(length(as.matrix(xtr)) == 27, 1, ntr)
  
  results <- foreach(idx = 1:(N1 * N2)) %dopar% {
    i <- ((idx - 1) %% N1) + 1
    j <- ((idx - 1) %/% N1) + 1
    if(N1==1){
      x1 <- xtr
      x2 <- xtr
    }else{
      x1 <- xtr[i,]
      x2 <- xtr[j,]
    }
    
    X1 <- x1 - rep(x1[1:3], 9)
    X2 <- x2 - rep(x2[1:3], 9)
    
    # === First rotation chain ===
    x_H_bar_1 <- X1[4:6]; x_H_bar_2 <- X1[7:9]
    r1 <- norm(x_H_bar_1)
    theta1 <- atan2(x_H_bar_1[1], x_H_bar_1[2])
    x_H_tilde_1 <- c(0, sin(theta1)*x_H_bar_1[1] + cos(theta1)*x_H_bar_1[2], x_H_bar_1[3])
    phi1 <- -atan2(x_H_tilde_1[3], x_H_tilde_1[2])
    
    psi_1 <- matrix(c(cos(theta1), -sin(theta1), 0,
                      sin(theta1), cos(theta1), 0,
                      0, 0, 1), nrow=3, byrow=TRUE)
    psi_2 <- matrix(c(1, 0, 0,
                      0, cos(phi1), -sin(phi1),
                      0, sin(phi1), cos(phi1)), nrow=3, byrow=TRUE)
    
    x_H_tilde_2 <- psi_2 %*% psi_1 %*% x_H_bar_2
    phi2 <- -atan2(x_H_tilde_2[3], x_H_tilde_2[1])
    psi_3 <- matrix(c(cos(phi2), 0, -sin(phi2),
                      0, 1, 0,
                      sin(phi2), 0, cos(phi2)), nrow=3, byrow=TRUE)
    
    X1 <- X1[-(1:9)]
    rho <- psi_3 %*% psi_2 %*% psi_1
    A1 <- c(r1, (psi_3 %*% x_H_tilde_2)[1], (psi_3 %*% x_H_tilde_2)[2],
            unlist(lapply(1:6, function(k) rho %*% X1[(3*(k-1)+1):(3*k)])))
    
    # === Second rotation chain ===
    x_H2_bar_1 <- X2[4:6]; x_H2_bar_2 <- X2[7:9]
    r1_2 <- norm(x_H2_bar_1)
    theta1 <- atan2(x_H2_bar_1[1], x_H2_bar_1[2])
    x_H2_tilde_1 <- c(0, sin(theta1)*x_H2_bar_1[1] + cos(theta1)*x_H2_bar_1[2], x_H2_bar_1[3])
    phi1 <- -atan2(x_H2_tilde_1[3], x_H2_tilde_1[2])
    
    psi_1_2 <- matrix(c(cos(theta1), -sin(theta1), 0,
                        sin(theta1), cos(theta1), 0,
                        0, 0, 1), nrow=3, byrow=TRUE)
    psi_2_2 <- matrix(c(1, 0, 0,
                        0, cos(phi1), -sin(phi1),
                        0, sin(phi1), cos(phi1)), nrow=3, byrow=TRUE)
    
    x_H2_tilde_2 <- psi_2_2 %*% psi_1_2 %*% x_H2_bar_2
    phi2 <- -atan2(x_H2_tilde_2[3], x_H2_tilde_2[1])
    psi_3_2 <- matrix(c(cos(phi2), 0, -sin(phi2),
                        0, 1, 0,
                        sin(phi2), 0, cos(phi2)), nrow=3, byrow=TRUE)
    
    X2 <- X2[-(1:9)]
    rho2 <- psi_3_2 %*% psi_2_2 %*% psi_1_2
    A2 <- c(r1_2, (psi_3_2 %*% x_H2_tilde_2)[1], (psi_3_2 %*% x_H2_tilde_2)[2],
            unlist(lapply(1:6, function(k) rho2 %*% X2[(3*(k-1)+1):(3*k)])))
    
    # === Covariances ===
    cov1 <- sigma^2 * exp(-(norm(A1 - A2))^2 / (2 * l^2)) / (l^3) * (norm(A1 - A2))^2 * t(rho) %*% rho2
    cov2 <- 2 * sigma * exp(-(norm(A1 - A2))^2 / (2 * l^2)) * t(rho) %*% rho2
    
    list(i = i, j = j, cov1 = cov1, cov2 = cov2)
  }
  
  K_l <- matrix(0, 3*N1, 3*N2)
  K_sigma <- matrix(0, 3*N1, 3*N2)
  
  for (res in results) {
    i <- res$i
    j <- res$j
    cov1 <- res$cov1
    cov2 <- res$cov2
    
    offset1 <- ifelse(N1 == 1, 0, N1)
    offset2 <- ifelse(N2 == 1, 0, N2)
    
    K_l[i, j] <- cov1[1,1]
    K_l[offset1 + i, offset2 + j] <- cov1[2,2]
    K_l[2*offset1 + i, 2*offset2 + j] <- cov1[3,3]
    K_l[offset1 + i, j] <- cov1[2,1]
    K_l[2*offset1 + i, j] <- cov1[3,1]
    K_l[i, offset2 + j] <- cov1[1,2]
    K_l[i, 2*offset2 + j] <- cov1[1,3]
    K_l[2*offset1 + i, offset2 + j] <- cov1[3,2]
    K_l[offset1 + i, 2*offset2 + j] <- cov1[2,3]
    
    K_sigma[i, j] <- cov2[1,1]
    K_sigma[offset1 + i, offset2 + j] <- cov2[2,2]
    K_sigma[2*offset1 + i, 2*offset2 + j] <- cov2[3,3]
    K_sigma[offset1 + i, j] <- cov2[2,1]
    K_sigma[2*offset1 + i, j] <- cov2[3,1]
    K_sigma[i, offset2 + j] <- cov2[1,2]
    K_sigma[i, 2*offset2 + j] <- cov2[1,3]
    K_sigma[2*offset1 + i, offset2 + j] <- cov2[3,2]
    K_sigma[offset1 + i, 2*offset2 + j] <- cov2[2,3]
  }
  
  if(N1==1){
    grad=0.5*c(-t(c(ytr))%*%inv_K%*%K_l%*%inv_K%*%c(
      ytr)+Trace(inv_K%*%K_l),
      -t(c(ytr))%*%inv_K%*%K_sigma%*%inv_K%*%c(
        ytr)+Trace(inv_K%*%K_sigma)
    )
    
  }else{
    grad=0.5*c(-t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l%*%inv_K%*%c(
      ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l),
      -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma%*%inv_K%*%c(
        ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_sigma)
    )
    
  }
  
  return(grad)
  
  
}


grad_standard=function(xtr,ytr,l,sigma){
  Trace=function(x){sum(diag(x))}
  x1=x2=xtr
  ntr=ifelse(is.null(dim(xtr)),1,nrow(xtr))
  
  norm=function(x){sum(x^2)^.5}
  
  K=cross_cov_mat(xtr,xtr,l,sigma)+
    diag(jit,nrow=3*ntr)
  
  
  K_l = K_sigma =
    matrix(0, ncol =3 * ifelse(length(as.matrix(xtr)) == 27, 1, ntr),
           nrow = 3 * ifelse(length(as.matrix(xtr)) == 27, 1, ntr))
  
  inv_K=solve(K)
  
  
  for(i in 1:(length(as.matrix(x1))/27)){
    for (j in 1:(length(as.matrix(x2))/27)){
   
      if(!is.null(dim(x1))){
        x1=x1[i,]
      }
      if(!is.null(dim(x2))){
        x2=x2[j,]
      }
      cov1= cov_mat(x1,x2,l,sigma)*(norm(x1-x2))^2/(l^3)
      
      cov2= cov_mat(x1,x2,l,1)*sigma*2
      
      x1=x2=xtr
      
      K_l[i, j] = cov1[1,1]
      K_l[ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[2,2]
      
      K_l[2*ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
          2*ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[3,3]
      
      K_l[ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i, j] = cov1[2,1]
      K_l[2*ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i, j] = cov1[3,1]
      
      K_l[i,ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[1,2]
      K_l[i,2*ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[1,3]
      
      K_l[2*ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[3,2]
      
      K_l[ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
          2*ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov1[2,3]
      
      
      
      K_sigma[i, j] = cov2[1,1]
      K_sigma[ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
              ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov2[2,2]
      
      K_sigma[2*ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
              2*ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov2[3,3]
      
      K_sigma[ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i, j] = cov2[2,1]
      K_sigma[2*ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i, j] = cov2[3,1]
      
      K_sigma[i,ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov2[1,2]
      K_sigma[i,2*ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov2[1,3]
      
      K_sigma[2*ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
              ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov2[3,2]
      
      K_sigma[ifelse(length(as.matrix(x1))==27,1,nrow(x1)) + i,
              2*ifelse(length(as.matrix(x2))==27,1,nrow(x2)) + j] = cov2[2,3]
      
    }
  }
  if(is.null(dim(ytr))){
    grad=0.5*c(-t(ytr)%*%inv_K%*%K_l%*%inv_K%*%c(
      ytr)+Trace(inv_K%*%K_l),
      -t(ytr)%*%inv_K%*%K_sigma%*%inv_K%*%c(
        ytr)+Trace(inv_K%*%K_sigma)
    )
  }else{
    grad=0.5*c(-t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l%*%inv_K%*%c(
      ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l),
      -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma%*%inv_K%*%c(
        ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_sigma)
    )
  }
  
  
  return(grad)
  
  
  
}


log_likelihood_standard=function(ytr,xtr,l,sigma){
  if(length(ytr)>3){
    if(ncol(ytr)==3){
      ytr=c(ytr[,1],ytr[,2],ytr[,3])
    }
  }
  
  
  ntr=ifelse(length(xtr)==27,1,nrow(xtr))
  Ktr=cross_cov_mat(xtr,xtr,l,sigma)+
    diag(jit,nrow=3*ntr)
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*determinant(Ktr)$modulus[1]-ntr*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}


RMSES=LogS=matrix(nrow=n_ind,ncol=2)
colnames(RMSES)<-c("standard","fundamental")
colnames(LogS)<-c("standard","fundamental")

print('Number of parallel workers: \n')
future::nbrOfWorkers()

ind_points_fund=ind_points=sample(data_set_size,10)

opt_par=Adam(function(x) {
  log_likelihood_standard(
    X[ind_points,], Y[ind_points,], x[1], x[2])
},function(x) {
  grad_standard(X[ind_points,], Y[ind_points,], x[1], x[2])
},rep(1,2),lr=lr,max_iter = max_iter)

cat('Opt par: ', opt_par,'\n')


                          
distances = outer(ind_points, ind_points,
                  FUN = Vectorize(function(x, y) {
                    Kxy=cov_mat(X[x, ], X[y, ], opt_par[1], opt_par[2])
                    Kxx=cov_mat(X[x, ], X[x, ], opt_par[1], opt_par[2])
                    Kyy=cov_mat(X[y, ], X[y, ], opt_par[1], opt_par[2])
                    
                    base::norm(
                      -2 * Kxy/sqrt(base::norm(Kxx,type='F')*base::norm(Kyy,type='F')) +
                        Kxx/base::norm(Kxx,type='F') +
                        Kyy/base::norm(Kyy,type='F'),
                      type = 'F'
                    )
                  }))
tau=quantile(distances,quant)
                          
cat('Tau : ', tau,'\n')

opt_par_fund=Adam(function(x) {
  log_likelihood_fund(
    X[ind_points_fund,], Y[ind_points_fund,], x[1], x[2])
},function(x) {
  grad_fund(X[ind_points_fund,], Y[ind_points_fund,], x[1], x[2])
},rep(1,2),lr=lr,max_iter = max_iter)

distances = outer(ind_points_fund, ind_points_fund,
                  FUN = Vectorize(function(x, y) {
                    Kxy=fund_cov_mat(X[x, ], X[y, ], opt_par_fund[1], opt_par_fund[2])
                    Kxx=fund_cov_mat(X[x, ], X[x, ], opt_par_fund[1], opt_par_fund[2])
                    Kyy=fund_cov_mat(X[y, ], X[y, ], opt_par_fund[1], opt_par_fund[2])
                    
                    base::norm(
                      -2 * Kxy/sqrt(base::norm(Kxx,type='F')*base::norm(Kyy,type='F')) +
                        Kxx/base::norm(Kxx,type='F') +
                        Kyy/base::norm(Kyy,type='F'),
                      type = 'F'
                    )
                  }))
tau_fund=quantile(distances,quant_fund)

cat('Opt par fund: ', opt_par_fund,'\n')

cat('Tau fund: ', tau_fund,'\n')
                         

for( i in 1:n_ind){
  Trace=function(Mat){sum(diag(Mat))}
  
  while(length(ind_points)<=n_ind_points[i]){
    
     for(kk in 1:nrow(X)){
      Similarities=numeric(nrow(X))
      if(!(kk %in% ind_points)){
         Similarities[kk]=min(sapply(ind_points,function(k){
        Kxy=cov_mat(X[k, ], X[kk, ], opt_par[1], opt_par[2])
                    Kxx=cov_mat(X[k, ], X[k, ], opt_par[1], opt_par[2])
                    Kyy=cov_mat(X[kk, ], X[kk, ], opt_par[1], opt_par[2])
                    
                    base::norm(
                      -2 * Kxy/sqrt(base::norm(Kxx,type='F')*base::norm(Kyy,type='F')) +
                        Kxx/base::norm(Kxx,type='F') +
                        Kyy/base::norm(Kyy,type='F'),
                      type = 'F'
                    )
      }))
        if(Similarities[kk]>tau){
        ind_points=c(ind_points,kk)
        print(kk)
        break
        }
      }
      if(kk== nrow(X)){
         print('No new points found')
         quant=.9*quant
         ind_points=c(ind_points,which.max(Similarities))
      }      
    }
    
    opt_par=Adam(function(x) {
      log_likelihood_standard(
        X[ind_points,], Y[ind_points,], x[1], x[2])
    },function(x) {
      grad_standard(X[ind_points,], Y[ind_points,], x[1], x[2])
    },opt_par,lr=lr,max_iter = max_iter)
    
    cat('Opt par: ', opt_par,'\n')
    distances = outer(ind_points, ind_points,
                  FUN = Vectorize(function(x, y) {
                    Kxy=cov_mat(X[x, ], X[y, ], opt_par[1], opt_par[2])
                    Kxx=cov_mat(X[x, ], X[x, ], opt_par[1], opt_par[2])
                    Kyy=cov_mat(X[y, ], X[y, ], opt_par[1], opt_par[2])
                    
                    base::norm(
                      -2 * Kxy/sqrt(base::norm(Kxx,type='F')*base::norm(Kyy,type='F')) +
                        Kxx/base::norm(Kxx,type='F') +
                        Kyy/base::norm(Kyy,type='F'),
                      type = 'F'
                    )
                  }))
    tau=quantile(distances,quant)

    cat('Tau : ', tau,'\n')

                          
  }
         
  while(length(ind_points_fund)<=n_ind_points[i]){
    
    for(kk in 1:nrow(X)){
      Similarities=numeric(nrow(X))
      if(!(kk %in% ind_points_fund)){
         Similarities[kk]=min(sapply(ind_points_fund,function(k){
        Kxy=fund_cov_mat(X[k, ], X[kk, ], opt_par_fund[1], opt_par_fund[2])
                    Kxx=fund_cov_mat(X[k, ], X[k, ], opt_par_fund[1], opt_par_fund[2])
                    Kyy=fund_cov_mat(X[kk, ], X[kk, ], opt_par_fund[1], opt_par_fund[2])
                    
                    base::norm(
                      -2 * Kxy/sqrt(base::norm(Kxx,type='F')*base::norm(Kyy,type='F')) +
                        Kxx/base::norm(Kxx,type='F') +
                        Kyy/base::norm(Kyy,type='F'),
                      type = 'F'
                    )
      }))
        if(Similarities[kk]>tau_fund){
        ind_points_fund=c(ind_points_fund,kk)
        print(kk)
        break
        }
      }
      if(kk== nrow(X)){
         print('No new points found')
         quant_fund=.9*quant_fund
         ind_points_fund=c(ind_points_fund,which.max(Similarities))
      }      
    }
    
    opt_par_fund=Adam(function(x) {
      log_likelihood_fund(
        X[ind_points_fund,], Y[ind_point_funds,], x[1], x[2])
    },function(x) {
      grad_fund(X[ind_points_fund,], Y[ind_points_fund,], x[1], x[2])
    },opt_par_fund,lr=lr,max_iter = max_iter)
    
    cat('Opt par fund: ', opt_par_fund,'\n')

    distances = outer(ind_points_fund, ind_points_fund,
                  FUN = Vectorize(function(x, y) {
                    Kxy=fund_cov_mat(X[x, ], X[y, ], opt_par_fund[1], opt_par_fund[2])
                    Kxx=fund_cov_mat(X[x, ], X[x, ], opt_par_fund[1], opt_par_fund[2])
                    Kyy=fund_cov_mat(X[y, ], X[y, ], opt_par_fund[1], opt_par_fund[2])
                    
                    base::norm(
                      -2 * Kxy/sqrt(base::norm(Kxx,type='F')*base::norm(Kyy,type='F')) +
                        Kxx/base::norm(Kxx,type='F') +
                        Kyy/base::norm(Kyy,type='F'),
                      type = 'F'
                    )
                  }))
   tau_fund=quantile(distances,quant_fund)
   cat('Tau fund: ', tau_fund,'\n')

  }

  cat('Ind points: ', ind_points,'\n')
  cat('Ind points fund: ', ind_points_fund,'\n')
  Xtr=X[ind_points,]
  #Xte=X[test_ind,]
  Ytr=Y[ind_points,]
  #Yte=Y[test_ind,]
  
  
  Ktetr2=cross_cov_mat(Xte,Xtr,opt_par[1],
                       opt_par[2])
  
  Ktetr_trtr_inv2=Ktetr2%*%solve(cross_cov_mat(Xtr,Xtr,opt_par[1],
                                               opt_par[2])+diag(jit,nrow=3*nrow(Xtr)))
  
  
  
  posterior_mean=posterior_mean2=Ktetr_trtr_inv2%*%c(Ytr[,1],Ytr[,2],Ytr[,3])
  
  
  (RMSES[i,1]=rmse2=mean(apply((cbind(posterior_mean2[1:(length(posterior_mean)/3)],
                                      posterior_mean2[(length(posterior_mean)/3+1):
                                                        (2*length(posterior_mean)/3)],
                                      posterior_mean2[(2*length(posterior_mean)/3+1):
                                                        (length(posterior_mean))])-Yte)^2,1,sum)
  )^.5)
  
  print(rmse2)
  posterior_cov2= cross_cov_mat(Xte,Xte,opt_par[1],
                                opt_par[2])-
    Ktetr_trtr_inv2%*%t(Ktetr2)
  
  
  (LogS[i,1]=LogS2=(t(c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean2)%*%
                      solve(posterior_cov2
                            +diag(1e-6,3*nrow(Xte))
                      )%*%
                      (c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean2)
                    +determinant(posterior_cov2
                                 +diag(1e-6,3*nrow(Xte))
                    )$modulus[1]+
                      log(2*pi)*nrow(Xte))/nrow(Xte))
  
  print(LogS2)  
  
  opt_par=opt_par_fund
  Xtr=X[ind_points_fund,]
  Ytr=Y[ind_points_fund,]
  
  
  Ktetr_eq=cross_fund_cov_mat(Xte,Xtr,opt_par[1],
                              opt_par[2])
  
  Ktetr_trtr_inv_eq=Ktetr_eq%*%inv(cross_fund_cov_mat(Xtr,Xtr,opt_par[1],
                                                      opt_par[2])+
                                     diag(jit,nrow=3*nrow(Xtr)))
  
  
  
  posterior_mean_eq=Ktetr_trtr_inv_eq%*%c(Ytr[,1],Ytr[,2],Ytr[,3])
  
  
  (RMSES[i,2]=rmse_eq=mean(apply((cbind(posterior_mean_eq[1:(length(posterior_mean_eq)/3)],
                                        posterior_mean_eq[(length(posterior_mean_eq)/3+1):
                                                            (2*length(posterior_mean_eq)/3)],
                                        posterior_mean_eq[(2*length(posterior_mean_eq)/3+1):
                                                            (length(posterior_mean_eq))])-Yte)^2,1,sum)
  )^.5)
  
  print(rmse_eq)
  
  posterior_cov_eq= cross_fund_cov_mat(Xte,Xte,opt_par[1],
                                       opt_par[2])-
    Ktetr_trtr_inv_eq%*%t(Ktetr_eq)
  
  
  (LogS[i,2]=LogS_eq=(t(c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean_eq)%*%
                        solve(posterior_cov_eq
                              +diag(1e-6,3*nrow(Xte))
                        )%*%
                        (c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean_eq)
                      +determinant(posterior_cov_eq
                                   +diag(1e-6,3*nrow(Xte))
                      )$modulus[1]+
                        log(2*pi)*nrow(Xte))/nrow(Xte))
  
  print(LogS_eq)
  
}
stopCluster(cl)
res=list(RMSES,LogS)
save(list = "res", file = paste0("Results/learning_curves_NMF_inducing_points_", id, ".rda"))
