id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

library(sp)
library(pracma)
cov_mat=function(x1,x2,l1,sigma1,l2,sigma2,nu1=0,nu2=0){
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=distmat(t(x1),x2)
  }else{
    dist_mat=distmat(x1,x2)
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

log_likelihood_grad=function(ytr,xtr,l1,sigma1,l2,sigma2,sigma_obs,
                             equivariant=FALSE,fundamental=FALSE,nu1=0,nu2=nu1){
  if(length(ytr)>2){
    if(ncol(ytr)==2){
      ytr=c(ytr[,1],ytr[,2])
    }
  }
  ntr=ifelse(length(xtr)==2,1,nrow(xtr))
  
  if(equivariant){
    if(fundamental){
      Ktr=eq_cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,fundamental = 1,nu1=nu1,nu2=nu2)+
        sigma_obs*diag(nrow=2*ntr)
    }else{
      Ktr=eq_cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,nu1=nu1,nu2=nu2)+
        sigma_obs*diag(nrow=2*ntr)
    }
  }else{
    Ktr=cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,nu1=nu1,nu2=nu2)+
      sigma_obs*diag(nrow=2*ntr)
  }
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-ntr*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}

grad_standard=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs_2){
  
  K=cov_mat(xtr,xtr,l1_2^.5,sigma1_2^.5,l2_2^.5,sigma2_2^.5)+
    diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  #K=cov_mat(xtr,xtr,opt_par[1],opt_par[2],opt_par[3],opt_par[4])+
  # diag(opt_par[5],nrow=2*nrow(xtr))
  
  n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
  
  if(n1==1){
    dist_mat=distmat(t(xtr),xtr)
  }else{
    dist_mat=distmat(xtr,xtr)
  }
  del_cov1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov1[1:n1,1:n2]=sigma1_2*exp(-.5*dist_mat^2/(l1_2))*dist_mat^2/(2*l1_2^2)
  del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov2[1:n1,1:n2]=0
  del_cov2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2_2*exp(-.5*dist_mat^2/(l2_2))*dist_mat^2/(2*l2_2^2)
  
  del_cov3=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov3[1:n1,1:n2]=exp(-.5*dist_mat^2/(l1_2))
  del_cov3[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_cov4=matrix(0,nrow=2*n1,ncol=2*n2)
  del_cov4[1:n1,1:n2]=0
  del_cov4[(n1+1):(2*n1),(n2+1):(2*n2)]=exp(-.5*dist_mat^2/(l2_2))
  
  
  inv_K=solve(K)
  
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_cov1),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov3%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_cov3),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_cov2),
         
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov4%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_cov4),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K))
  
  return(grad)
  
}
norm=function(x){sum(x^2)^.5}
M=function(h,n,a){
  k=0:n
  exp(-a*norm(h))*sum(factorial(n+k)/factorial(2*n)*
                        (factorial(n)/(factorial(k)*factorial(n-k)))*
                        (2*a*norm(h))^(n-k))
}
M=function(H,n,a){
  k=0:n
  l=length(H)
  ncH=ncol(H)
  matrix(sapply(H,function(h){exp(-a*h)*sum(factorial(n+k)/factorial(2*n)*
                                              (factorial(n)/(factorial(k)*factorial(n-k)))*
                                              (2*a*h)^(n-k))}),ncol=ifelse(l>1,ncH,1),
         byrow=1)
  
}


matern_cov_mat=function(x1,x2,sigma1,sigma2,rho,nu1,nu2,nu12,a1,a2,a12){
  if(validity_check(rho,nu1,nu2,nu12,a1,a2,a12)==0){print("invalid parameters")}
  
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
  cov[1:n1,1:n2]=sigma1^2*M(dist_mat,nu1,a1)
  cov[(n1+1):(2*n1),1:n2]=cov[(1):n1,(n2+1):(2*n2)]=rho*sigma1*sigma2*
    M(dist_mat,nu12,a12)
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*M(dist_mat,nu2,a2)
  return(cov)
  
}




log_likelihood_matern=function(ytr,xtr,sigma1,sigma2,rho,nu1,nu2,nu12,tau,a1,a2,a12){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  Ktr=matern_cov_mat(xtr,xtr,sigma1,sigma2,rho,nu1,nu2,nu12,a1,a2,a12)+
    tau^2*diag(nrow = 2*nrow(xtr))
  
  
  ll=-0.5*t(ytr)%*%solve(Ktr)%*%ytr-0.5*determinant(Ktr,logarithm = 1)$modulus-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  #print(c(ll))
  return(ll)
  
  
}


grad_matern=function(xtr,ytr,sigma1,sigma2,rho,nu1=1,nu2=1,nu12=2,tau,a1,a2,a12){
  x1=x2=xtr
  K=matern_cov_mat(xtr,xtr,sigma1,sigma2,rho,nu1,nu2,nu12,a1,a2,a12)+
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
  del_sigma1[(n1+1):(2*n1),1:n2]=del_sigma1[(1):n1,(n2+1):(2*n2)]=rho*sigma2*
    M(dist_mat,nu12,a12)
  del_sigma1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_sigma2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_sigma2[1:n1,1:n2]=0
  del_sigma2[(n1+1):(2*n1),1:n2]=del_sigma2[(1):n1,(n2+1):(2*n2)]=rho*sigma1*
    M(dist_mat,nu12,a12)
  del_sigma2[(n1+1):(2*n1),(n2+1):(2*n2)]=2*sigma2*M(dist_mat,nu2,a2)
  
  del_rho=matrix(0,nrow=2*n1,ncol=2*n2)
  del_rho[1:n1,1:n2]=0
  del_rho[(n1+1):(2*n1),1:n2]=del_rho[(1):n1,(n2+1):(2*n2)]=sigma2*sigma1*
    M(dist_mat,nu12,a12)
  del_rho[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  del_tau=2*tau*diag(nrow=2*nrow(xtr))
  
  del_a1=matrix(0,nrow=2*n1,ncol=2*n2)
  del_a1[1:n1,1:n2]=sigma1^2*M_a(dist_mat,nu1,a1)
  del_a1[(n1+1):(2*n1),1:n2]=del_a1[(1):n1,(n2+1):(2*n2)]=0
  del_a1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  del_a2=matrix(0,nrow=2*n1,ncol=2*n2)
  del_a2[1:n1,1:n2]=0
  del_a2[(n1+1):(2*n1),1:n2]=del_a2[(1):n1,(n2+1):(2*n2)]=0
  del_a2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*M_a(dist_mat,nu2,a2)
  
  del_a12=matrix(0,nrow=2*n1,ncol=2*n2)
  del_a12[1:n1,1:n2]=0
  del_a12[(n1+1):(2*n1),1:n2]=del_a12[(1):n1,(n2+1):(2*n2)]=
    rho*sigma2*sigma1*M_a(dist_mat,nu12,a12)
  del_a12[(n1+1):(2*n1),(n2+1):(2*n2)]=0
  
  
  inv_K=inv(K)
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma1),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_sigma2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_sigma2),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_rho%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_rho),
             
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_tau%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_tau),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_a1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_a1),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_a2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_a2),
             -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_a12%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_a12))
  
  return(grad)
  
}

validity_check=function(rho,nu1,nu2,nu12,a1,a2,a12){
  rho_bound=a1^(nu1+.5)*a2^(nu2+.5)/(a12^(2*nu12+1))*gamma(nu12+.5)/
    (gamma(nu1+.5)^.5*gamma(nu2+.5)^.5)*(
      exp(1)*(a12^2-.5*(a1^2+a2^2))/(nu12-.5*(nu1+nu2))
    )^(nu12-.5*(nu1+nu2))
  ifelse((nu12>=0.5*(nu1+nu2))&(a12^2>=0.5*(a1^2+a2^2))&
           (abs(rho)<=rho_bound),1,0)
  
  
}
rho_bound=function(nu1,nu2,nu12,a1,a2,a12){
  a1^(nu1+.5)*a2^(nu2+.5)/(a12^(2*nu12+1))*gamma(nu12+.5)/
    (gamma(nu1+.5)^.5*gamma(nu2+.5)^.5)*(
      exp(1)*(a12^2-.5*(a1^2+a2^2))/(nu12-.5*(nu1+nu2))
    )^(nu12-.5*(nu1+nu2))
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


Adam_matern = function(f, grad_f, params_init, lr = lr, beta1 = 0.9, beta2 = 0.999,
                       epsilon = 1e-8, max_iter = max_iter, 
                       tol = 1e-10,alpha=FALSE,alpha2=FALSE) {
  
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
    a1=params_new[(length(params_new)-2)]
    a2=params_new[(length(params_new)-1)]
    a12=params_new[(length(params_new))]
    rho=params_new[3]
    
    a1_old=params[(length(params)-2)]
    a2_old=params[(length(params)-1)]
    a12_old=params[(length(params))]
    rho_old=params[3]
    
    a1_inter=a1;a2_inter=a2;a12_inter=a12;rho_inter=rho
    j=0
    while(validity_check(rho_inter,nu1,nu2,nu12,a1_inter,a2_inter,a12_inter)==0){
      a1_inter=a1;a2_inter=a2;a12_inter=a12;rho_inter=rho
      
      a1_inter=ifelse(rbinom(1,1,.5)==1,a1_inter,a1_old)
      a2_inter=ifelse(rbinom(1,1,.5)==1,a2_inter,a2_old)
      a12_inter=ifelse(rbinom(1,1,.5)==1,a12_inter,a12_old)
      rho_inter=ifelse(rbinom(1,1,.5)==1,rho_inter,rho_old)
      # print(validity_check(rho_inter,nu1,nu2,nu12,a1_inter,a2_inter,a12_inter))
#      max_iter=max_iter+1
      j=j+1
      if(j==10000){
        break
        print(j)
        
      }
      
      
    }
    if(j==10000){
      
      params_new[(length(params_new)-2)]=params[(length(params)-2)]
      params_new[(length(params_new)-1)]=params[(length(params)-1)]
      params_new[(length(params_new))]=params[(length(params))]
      params_new[3]=params[3]
      break
    }else{
      params_new[(length(params_new)-2)]=a1_inter
      params_new[(length(params_new)-1)]=a2_inter
      params_new[(length(params_new))]=a12_inter
      params_new[3]=rho_inter
    }
    
    
    params=params_new
    
    a=(validity_check(params[3],nu1,nu2,nu12,params[5],params[6],params[7]))
    if(a==0){print(a)}
    
    #print(params)
    #print(c(f(params)))
  }
  #print("iter= ")
  #print(iter)
  return(params)
}



NU=rbind(c(0,0,1),c(1,1,2),c(2,2,4),c(5,5,10),c(10,10,20),c(30,30,60))
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

a1=.1;a2=.1;a12=.2
sigma1=sigma2=1
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
  
  tau=.1
  mean_Y=mean(Ytr)
  Ytr=Ytr-mean_Y
  initial_par=c(1/a1,sigma1,1/a2,sigma2,tau)
  
  opt_par = Adam(
    function(x) {
      log_likelihood_grad(
        Ytr, Xtr, x[1], x[2], x[3], x[4], x[5])
    },
    function(x) {
      grad_standard(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5])
    },
    initial_par[1:5],
    max_iter = 50,lr=.001
  )^.5
  #opt_par=initial_par
  Ktetr_standard=cov_mat(Xte,Xtr,opt_par[1],
                         opt_par[2],
                         opt_par[3],
                         opt_par[4]
  )
  
  Ktetr_trtr_inv_standard=Ktetr_standard%*%inv(cov_mat(Xtr,Xtr,opt_par[1],
                                                       opt_par[2],
                                                       opt_par[3],
                                                       opt_par[4],
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
                                  opt_par[4])-Ktetr_trtr_inv_standard%*%t(Ktetr_standard)
  
  
  (LogS_standard[j]=(t(c(Yte[,1],Yte[,2])-posterior_mean_standard)%*%solve(posterior_cov_standard+
                                                                          diag(opt_par[5]^2,2*nrow(Xte)))%*%
                    (c(Yte[,1],Yte[,2])-posterior_mean_standard)
                  +determinant(posterior_cov_standard+diag(opt_par[5]^2,2*nrow(Xte)))$modulus[1]+
                    log(2*pi)*nrow(Xte))/nrow(Xte))
  
  for (i in 1:nrow(NU)){
    nu1=NU[i,1];nu2=NU[i,2];nu12=NU[i,3]
    rho=min(cor(Ytr)[2],rho_bound(nu1,nu2,nu12,a1,a2,a12))
    print(validity_check(rho,nu1,nu2,nu12,a1,a2,a12))
    
    initial_par=c(sigma1,sigma2,rho,tau,a1,a2,a12)
    
    opt_par=Adam_matern(function(x){log_likelihood_matern(Ytr,Xtr,x[1],x[2],x[3],nu1,nu2,nu12,
                                                          x[4],x[5],x[6],x[7])},
                        function(x){grad_matern(Xtr,Ytr,x[1],x[2],x[3],
                                                nu1=nu1,nu2=nu2,nu12=nu12,
                                                x[4],x[5],x[6],x[7])},
                        initial_par,
                        lr=.001,max_iter = 50
    )
    opt_tau=opt_par[4]
    opt_par=c(opt_par[1:3],nu1,nu2,nu12,opt_par[5:7])
    
    Ktetr_matern=matern_cov_mat(Xte,Xtr,opt_par[1],
                                opt_par[2],
                                opt_par[3],
                                opt_par[4],
                                opt_par[5],
                                opt_par[6],opt_par[7],
                                opt_par[8],
                                opt_par[9])
    
    Ktetr_trtr_inv_matern=Ktetr_matern%*%inv(matern_cov_mat(Xtr,Xtr,opt_par[1],
                                                            opt_par[2],
                                                            opt_par[3],
                                                            opt_par[4],
                                                            opt_par[5],
                                                            opt_par[6],opt_par[7],
                                                            opt_par[8],
                                                            opt_par[9])+
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
                                         opt_par[8],
                                         opt_par[9])-
      Ktetr_trtr_inv_matern%*%t(Ktetr_matern)
    
    
    (LogS_matern[i,j]=(t(c(Yte[,1],Yte[,2])-posterior_mean_matern)%*%solve(posterior_cov_matern+
                                                                        diag(opt_tau^2,2*nrow(Xte)))%*%
                     (c(Yte[,1],Yte[,2])-posterior_mean_matern)
                   +determinant(posterior_cov_matern+diag(opt_tau^2,2*nrow(Xte)))$modulus[1]+
                     log(2*pi)*nrow(Xte))/nrow(Xte))
    
  }
}

res=list(RMSEs_standard, LogS_standard,RMSEs_matern,LogS_matern)

save(list = "res", file = paste0("Spatial_stats_project", id, ".rda"))
