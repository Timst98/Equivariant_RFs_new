id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
library(sp)
library(pracma)

cov_mat=function(x1,x2,sigma1,sigma2,sigma3,sigma4,l1,l2){
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=distmat(t(x1),x2)
  }else{
    dist_mat=distmat(x1,x2)
  }
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  cov[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2*(l1^2))+sigma2^2*exp(-.5*dist_mat^2*(l2^2))
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*exp(-.5*dist_mat^2*(l1^2))+sigma4^2*exp(-.5*dist_mat^2*(l2^2))
  cov[(n1+1):(2*n1),1:n2]=
    cov[1:n1,(n2+1):(2*n2)]=
    sigma1*sigma3*exp(-.5*dist_mat^2*(l1^2))+
    sigma2*sigma4*exp(-.5*dist_mat^2*(l2^2))
  

  return(cov)
}

log_likelihood=function(ytr,xtr,sigma1,sigma2,sigma3,sigma4,tau,l1,l2){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  Ktr=cov_mat(xtr,xtr,sigma1,sigma2,sigma3,sigma4,l1,l2)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  
  ll=-0.5*t(ytr)%*%solve(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm = 1)$modulus-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
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



norm=function(x){sum(x^2)^.5}

M=function(H,n,a){
  k=0:n
  l=length(H)
  ncH=ncol(H)
  matrix(sapply(H,function(h){exp(-a*h)*sum(factorial(n+k)/factorial(2*n)*
                                              (factorial(n)/(factorial(k)*factorial(n-k)))*
                                              (2*a*h)^(n-k))}),ncol=ifelse(l>1,ncH,1),
         byrow=1)
  
}


matern_cov_mat=function(x1,x2,sigma1,sigma2,sigma3,sigma4,nu1,nu2,a1,a2){
 
  x1=as.matrix(x1);x2=as.matrix(x2)
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
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  cov[1:n1,1:n2]=sigma1^2*M(dist_mat,nu1,a1)+sigma2^2*M(dist_mat,nu2,a2)
  cov[(n1+1):(2*n1),1:n2]=cov[(1):n1,(n2+1):(2*n2)]=
   sigma1*sigma3*M(dist_mat,nu1,a1)+sigma2*sigma4*M(dist_mat,nu2,a2)
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*M(dist_mat,nu1,a1)+sigma4^2*M(dist_mat,nu2,a2)
  return(cov)
  
}




log_likelihood_matern=function(ytr,xtr,sigma1,sigma2,sigma3,sigma4,nu1,nu2,tau,a1,a2){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  Ktr=matern_cov_mat(xtr,xtr,sigma1,sigma2,sigma3,sigma4,nu1,nu2,a1,a2)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  
  ll=-0.5*t(ytr)%*%solve(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm = 1)$modulus-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}


grad_matern=function(xtr,ytr,sigma1,sigma2,sigma3,sigma4,nu1,nu2,tau,a1,a2){
  x1=x2=xtr
  K=matern_cov_mat(xtr,xtr,sigma1,sigma2,sigma3,sigma4,nu1,nu2,a1,a2)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  M_a=function(H,n,a){
    k=0:n
    l=length(H)
    ncH=ncol(H)
    matrix(sapply(H,function(h){exp(-a*h)*(-h*sum(factorial(n+k)/factorial(2*n)*
                                                    (factorial(n)/(factorial(k)*factorial(n-k)))*
                                                    (2*a*h)^(n-k))+
                                             sum(factorial(n+k)/factorial(2*n)*
                                                   (factorial(n)/(factorial(k)*factorial(n-k)))*
                                                   (2*h)^(n-k)*a^(n-k-1)*(n-k)))}),ncol=ifelse(l>1,ncH,1),
           byrow=1)
    
  }
  
  
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
  del_sigma1[1:n1,1:n2]=2*sigma1*M(dist_mat,nu1,a1)
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=sigma3*
    M(dist_mat,nu1,a1)
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=2*sigma2*M(dist_mat,nu2,a2)
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=sigma4*
    M(dist_mat,nu2,a2)
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
 
  del_sigma3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma3[1:n1,1:n2]=0
  del_sigma3[(n1+1):(2*n1),1:n2]=del_sigma3[(1):n1,(n2+1):(2*n2)]=sigma1*
    M(dist_mat,nu1,a1)
  del_sigma3[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma3*M(dist_mat,nu1,a1)
  
  del_sigma4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma4[1:n1,1:n2]=0
  del_sigma4[(n1+1):(2*n1),1:n2]=del_sigma4[(1):n1,(n2+1):(2*n2)]=sigma2*
    M(dist_mat,nu2,a2)
  del_sigma4[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma4*M(dist_mat,nu2,a2)
  
  
  del_tau=2*tau*diag(nrow=2*nrow(xtr))
  
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  cov[1:n1,1:n2]=sigma1^2*M(dist_mat,nu1,a1)+sigma2^2*M(dist_mat,nu2,a2)
  cov[(n1+1):(2*n1),1:n2]=cov[(1):n1,(n2+1):(2*n2)]=
    sigma1*sigma3*M(dist_mat,nu1,a1)+sigma2*sigma4*M(dist_mat,nu2,a2)
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*M(dist_mat,nu1,a1)+sigma4^2*M(dist_mat,nu2,a2)
  
  
  del_a1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_a1[1:n1,1:n2]=sigma1^2*M_a(dist_mat,nu1,a1)
  del_a1[(n1+1):(2*n1),1:n2]=del_a1[(1):n1,(n2+1):(2*n2)]=sigma1*sigma3*M_a(dist_mat,nu1,a1)
  del_a1[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma3^2*M_a(dist_mat,nu1,a1)
  
  del_a2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_a2[1:n1,1:n2]=sigma2^2*M_a(dist_mat,nu2,a2)
  del_a2[(n1+1):(2*n1),1:n2]=del_a2[(1):n1,(n2+1):(2*n2)]=sigma2*sigma4*M_a(dist_mat,nu2,a2)
  del_a2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma4^2*M_a(dist_mat,nu2,a2)
  
  
  inv_K=inv(K)
  
  #sigma1,sigma2,sigma3,sigma4,nu1=1,nu2=1,tau,a1,a2
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma1),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma2),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma3%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma3),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma4%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma4),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_tau%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_tau),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_a1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_a1),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_a2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_a2)
  )
  
  return(grad)
  
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
    #print(params)
    #print(c(f(params)))
  }
  #print("iter= ")
  #print(iter)
  return(params)
}





NU=rbind(c(0,1),c(1,0),c(2,2),c(5,5),c(5,2),c(30,30))
data("meuse")
meuse_df <- as.data.frame(meuse)
x=meuse_df$x
y=meuse_df$y

X=cbind(x,y)
X=scale(X)


Y=log(cbind(meuse_df$cadmium,meuse_df$zinc))
#Y=Y-mean(Y)


train_inds=list(sample(nrow(X),20))
for( i in 2:20){
  train_ind=train_inds[[i-1]]
  train_inds=append(train_inds, list(c(train_ind,sample((1:155)[-train_ind],5))))
}

a1=l1=.1;a2=l2=.1
sigma1=sigma2=sigma3=sigma4=1
tau=.1
#validity_check(rho,nu1,nu2,nu12,a1,a2,a12)

RMSEs_matern=matrix(nrow=nrow(NU),ncol=length(train_inds))
LogS_matern=matrix(nrow=nrow(NU),ncol=length(train_inds))
RMSEs_standard=LogS_standard=numeric(length(train_inds))

for (j in 1:length(train_inds)){
  train_ind=train_inds[[j]]
  test_ind=-train_ind
  
  Xtr=as.matrix(X[train_ind,])
  Ytr=Y[train_ind,]
  Xte=X[-train_ind,]
  Yte=Y[-train_ind,]
  
 
  mean_Y=mean(Ytr)
  Ytr=Ytr-mean_Y
  initial_par=c(sigma1,sigma2,sigma3,sigma4,tau,l1,l2)
  
  opt_par = Adam(
    function(x) {
      log_likelihood(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],x[6], x[7])
    },
    function(x) {
      grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],x[6], x[7])
    },
    initial_par,
    max_iter = 100,lr=.001
  )
  #opt_par=initial_par
  Ktetr_standard=cov_mat(Xte,Xtr,opt_par[1],
                         opt_par[2],
                         opt_par[3],
                         opt_par[4],opt_par[6],
                         opt_par[7]
  )
  
  Ktetr_trtr_inv_standard=Ktetr_standard%*%inv(cov_mat(Xtr,Xtr,opt_par[1],
                                                       opt_par[2],
                                                       opt_par[3],
                                                       opt_par[4],opt_par[6],
                                                       opt_par[7]
  )
  
  +diag(opt_par[5]^2,nrow=2*nrow(Xtr)))
  
  posterior_mean_standard=mean_Y+Ktetr_trtr_inv_standard%*%c(Ytr[,1],Ytr[,2])
  
  
  
  (RMSEs_standard[j]=mean(apply((cbind(posterior_mean_standard[1:(0.5*length(posterior_mean_standard))],
                                       posterior_mean_standard[(0.5*length(posterior_mean_standard)+1):
                                                                 (length(posterior_mean_standard))])
                                 -Yte)^2,1,sum))^.5)
  posterior_cov_standard= cov_mat(Xte,Xte,opt_par[1],
                                  opt_par[2],
                                  opt_par[3],
                                  opt_par[4],opt_par[6],
                                  opt_par[7])-Ktetr_trtr_inv_standard%*%t(Ktetr_standard)
  
  
  (LogS_standard[j]=(t(c(Yte[,1],Yte[,2])-posterior_mean_standard)%*%solve(posterior_cov_standard+
                                                                             diag(opt_par[5]^2,2*nrow(Xte)))%*%
                       (c(Yte[,1],Yte[,2])-posterior_mean_standard)
                     +determinant(posterior_cov_standard+diag(opt_par[5]^2,2*nrow(Xte)))$modulus[1]+
                       log(2*pi)*nrow(Xte))/nrow(Xte))
  
  for (i in 1:nrow(NU)){
    nu1=NU[i,1];nu2=NU[i,2]
    
    initial_par=c(sigma1,sigma2,sigma3,sigma4,tau,a1,a2)
    
    opt_par=Adam(function(x){log_likelihood_matern(Ytr,Xtr,x[1],x[2],x[3],x[4],nu1,nu2,
                                                          x[5],x[6],x[7])},
                        function(x){grad_matern(Xtr,Ytr,x[1],x[2],x[3],x[4],
                                                nu1=nu1,nu2=nu2,
                                                x[5],x[6],x[7])},
                        initial_par,
                        lr=.001,max_iter = 100
    )
    opt_tau=opt_par[5]
    opt_par=c(opt_par[1:4],nu1,nu2,opt_par[6:7])
    
    Ktetr_matern=matern_cov_mat(Xte,Xtr,opt_par[1],
                                opt_par[2],
                                opt_par[3],
                                opt_par[4],
                                opt_par[5],
                                opt_par[6],opt_par[7],opt_par[8])
    
    Ktetr_trtr_inv_matern=Ktetr_matern%*%inv(matern_cov_mat(Xtr,Xtr,opt_par[1],
                                                            opt_par[2],
                                                            opt_par[3],
                                                            opt_par[4],
                                                            opt_par[5],
                                                            opt_par[6],opt_par[7],
                                                            opt_par[8])+
                                               diag(opt_tau^2,nrow=2*nrow(Xtr)))
    
    posterior_mean_matern=mean_Y+Ktetr_trtr_inv_matern%*%c(Ytr[,1],Ytr[,2])
    
    
    
    (RMSEs_matern[i,j]=mean(apply((cbind(posterior_mean_matern[1:(0.5*length(posterior_mean_matern))],
                                         posterior_mean_matern[(0.5*length(posterior_mean_matern)+1):
                                                                 (length(posterior_mean_matern))])
                                   -Yte)^2,1,sum))^.5)
    posterior_cov_matern= matern_cov_mat(Xte,Xte,opt_par[1],
                                         opt_par[2],
                                         opt_par[3],
                                         opt_par[4],
                                         opt_par[5],
                                         opt_par[6],opt_par[7],
                                         opt_par[8])-
      Ktetr_trtr_inv_matern%*%t(Ktetr_matern)
    
    
    (LogS_matern[i,j]=(t(c(Yte[,1],Yte[,2])-posterior_mean_matern)%*%solve(posterior_cov_matern+
                                                                             diag(opt_tau^2,2*nrow(Xte)))%*%
                         (c(Yte[,1],Yte[,2])-posterior_mean_matern)
                       +determinant(posterior_cov_matern+diag(opt_tau^2,2*nrow(Xte)))$modulus[1]+
                         log(2*pi)*nrow(Xte))/nrow(Xte))
    
  }
}

res=list(RMSEs_standard, LogS_standard,RMSEs_matern,LogS_matern)

save(list = "res", file = paste0("Spatial_stats_project_LMC", id, ".rda"))
