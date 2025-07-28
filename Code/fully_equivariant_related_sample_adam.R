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

jit=1e-8
fund_cov_mat=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  norm=function(x){sum(x^2)^.5}
  if(norm(x1[4:6]-x1[1:3])>=norm(x1[7:9]-x1[1:3])){
    x_H_bar_1=x1[4:6]-x1[1:3]
    x_H_bar_2=x1[7:9]-x1[1:3]
    r1=norm(x_H_bar_1)
    r2=norm(x_H_bar_2)
  }else{
    x_H_bar_1=x1[7:9]-x1[1:3]
    x_H_bar_2=x1[4:6]-x1[1:3]
    r1=norm(x_H_bar_1)
    r2=norm(x_H_bar_2)
  }
  
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
  
  A1=c(r1,(psi_3%*%x_H_tilde_2)[1],(psi_3%*%x_H_tilde_2)[2])
  
  
  if(norm(x2[4:6]-x2[1:3])>=norm(x2[7:9]-x2[1:3])){
    x_H2_bar_1=x2[4:6]-x2[1:3]
    x_H2_bar_2=x2[7:9]-x2[1:3]
    r1_2=norm(x_H2_bar_1)
    r2_2=norm(x_H2_bar_2)
  }else{
    x_H2_bar_1=x2[7:9]-x2[1:3]
    x_H2_bar_2=x2[4:6]-x2[1:3]
    r1_2=norm(x_H2_bar_1)
    r2_2=norm(x_H2_bar_2)
  }
  
  theta1_2=atan2(x_H2_bar_1[1],x_H2_bar_1[2])
  x_H2_tilde_1=c(0,sin(theta1_2)*x_H2_bar_1[1]+cos(theta1_2)*x_H2_bar_1[2],x_H2_bar_1[3])
  phi1_2=-atan2(x_H2_tilde_1[3],x_H2_tilde_1[2])
  
  psi_1_2=matrix(c(cos(theta1_2),-sin(theta1_2),0,sin(theta1_2),cos(theta1_2),0,0,0,1),
                 nrow=3,byrow = 1)
  psi_2_2=matrix(c(1,0,0,0,cos(phi1_2),-sin(phi1_2),0,sin(phi1_2),cos(phi1_2)),
                 nrow=3,byrow = 1)
  
  x_H2_tilde_2=psi_2_2%*%psi_1_2%*%x_H2_bar_2
  
  phi2_2=-atan2(x_H2_tilde_2[3],x_H2_tilde_2[1])
  psi_3_2=matrix(c(cos(phi2_2),0,-sin(phi2_2),0,1,0,sin(phi2_2),0,cos(phi2_2)),
                 nrow=3,byrow = 1)
  
  A2=c(r1_2,(psi_3_2%*%x_H2_tilde_2)[1],(psi_3_2%*%x_H2_tilde_2)[2])
  
  t(psi_3%*%psi_2%*%psi_1)%*%(exp(diag(c(-(norm(A1-A2))^2/(2*l1^2),
                                         -(norm(A1-A2))^2/(2*l2^2),
                                         -(norm(A1-A2))^2/(2*l3^2))))*
                                diag(c(sigma1^2,sigma2^2,sigma3^2)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
  
  
}

cross_fund_cov_mat=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  norm=function(x){sum(x^2)^.5}
  cov= matrix(0, ncol =3 * ifelse(length(as.matrix(x2)) == 9, 1,nrow(x2)),
              nrow = 3 * ifelse(length(as.matrix(x1)) == 9, 
                                1, nrow(x1)))
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      if(length(as.matrix(x2)) == 9){
        if(length(as.matrix(x1)) == 9){
          cov1= fund_cov_mat(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3)
          
        }else{
          cov1= fund_cov_mat(x1[i,],x2,l1,sigma1,l2,sigma2,l3,sigma3)
          
        }
      }else if(length(as.matrix(x1)) == 9){
        cov1= fund_cov_mat(x1,x2[j,],l1,sigma1,l2,sigma2,l3,sigma3)
        
      }else{
        cov1= fund_cov_mat(x1[i,],x2[j,],l1,sigma1,l2,sigma2,l3,sigma3)
        
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


ind_train=sample(331,30)
ind_test=-ind_train

Xtr=as.matrix(original_struct[ind_train,X_names])
Ytr=as.matrix(original_struct[ind_train,Y_names])

Xte=as.matrix(original_struct[ind_test,X_names])
Yte=as.matrix(original_struct[ind_test,Y_names])
for(dat in list(rot_0,rot_1,rot_2,rot_3,rot_4,rot_5,rot_6,rot_7,rot_8,rot_9)){
  Xtr=rbind(Xtr,as.matrix(dat[ind_train,X_names]))
  Ytr=rbind(Ytr,as.matrix(dat[ind_train,Y_names]))
  
  Xte=rbind(Xte,as.matrix(dat[ind_test,X_names]))
  Yte=rbind(Yte,as.matrix(dat[ind_test,Y_names]))
}

log_likelihood_fund=function(ytr,xtr,l1,sigma1,l2,sigma2,l3,sigma3){
  if(length(ytr)>3){
    if(ncol(ytr)==3){
      ytr=c(ytr[,1],ytr[,2],ytr[,3])
    }
  }
  
  
  ntr=ifelse(length(xtr)==9,1,nrow(xtr))
  Ktr=cross_fund_cov_mat(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*ntr)
  
  
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
  
  K=cross_fund_cov_mat(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*nrow(xtr))
  
  K_l1 =  K_l2 = K_l3 = K_sigma1 = K_sigma2 = K_sigma3 =
    matrix(0, ncol =3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)),
           nrow = 3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)))
  
  inv_K=solve(K)
  
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      x1=xtr[i,]
      x2=xtr[j,]
      if(norm(x1[4:6]-x1[1:3])>=norm(x1[7:9]-x1[1:3])){
        x_H_bar_1=x1[4:6]-x1[1:3]
        x_H_bar_2=x1[7:9]-x1[1:3]
        r1=norm(x_H_bar_1)
        r2=norm(x_H_bar_2)
      }else{
        x_H_bar_1=x1[7:9]-x1[1:3]
        x_H_bar_2=x1[4:6]-x1[1:3]
        r1=norm(x_H_bar_1)
        r2=norm(x_H_bar_2)
      }
      
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
      
      A1=c(r1,(psi_3%*%x_H_tilde_2)[1],(psi_3%*%x_H_tilde_2)[2])
      
      
      if(norm(x2[4:6]-x2[1:3])>=norm(x2[7:9]-x2[1:3])){
        x_H2_bar_1=x2[4:6]-x2[1:3]
        x_H2_bar_2=x2[7:9]-x2[1:3]
        r1_2=norm(x_H2_bar_1)
        r2_2=norm(x_H2_bar_2)
      }else{
        x_H2_bar_1=x2[7:9]-x2[1:3]
        x_H2_bar_2=x2[4:6]-x2[1:3]
        r1_2=norm(x_H2_bar_1)
        r2_2=norm(x_H2_bar_2)
      }
      
      theta1_2=atan2(x_H2_bar_1[1],x_H2_bar_1[2])
      x_H2_tilde_1=c(0,sin(theta1_2)*x_H2_bar_1[1]+cos(theta1_2)*x_H2_bar_1[2],x_H2_bar_1[3])
      phi1_2=-atan2(x_H2_tilde_1[3],x_H2_tilde_1[2])
      
      psi_1_2=matrix(c(cos(theta1_2),-sin(theta1_2),0,sin(theta1_2),cos(theta1_2),0,0,0,1),
                     nrow=3,byrow = 1)
      psi_2_2=matrix(c(1,0,0,0,cos(phi1_2),-sin(phi1_2),0,sin(phi1_2),cos(phi1_2)),
                     nrow=3,byrow = 1)
      
      x_H2_tilde_2=psi_2_2%*%psi_1_2%*%x_H2_bar_2
      
      phi2_2=-atan2(x_H2_tilde_2[3],x_H2_tilde_2[1])
      psi_3_2=matrix(c(cos(phi2_2),0,-sin(phi2_2),0,1,0,sin(phi2_2),0,cos(phi2_2)),
                     nrow=3,byrow = 1)
      
      A2=c(r1_2,(psi_3_2%*%x_H2_tilde_2)[1],(psi_3_2%*%x_H2_tilde_2)[2])
      
      
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%diag(c(
        sigma1^2*exp(-(norm(A1-A2))^2/(2*l1^2))/(l1^3)*(norm(A1-A2))^2,0,0))%*%
        psi_3_2%*%psi_2_2%*%psi_1_2
      
      cov2=t(psi_3%*%psi_2%*%psi_1)%*%diag(c(0,
                                             sigma2^2*exp(-(norm(A1-A2))^2/(2*l2^2))/(l2^3)*(norm(A1-A2))^2,0))%*%
        psi_3_2%*%psi_2_2%*%psi_1_2
      cov3=t(psi_3%*%psi_2%*%psi_1)%*%diag(c(0,0,
                                             sigma3^2*exp(-(norm(A1-A2))^2/(2*l3^2))/(l3^3)*(norm(A1-A2))^2))%*%
        psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      
      cov4= t(psi_3%*%psi_2%*%psi_1)%*%
        diag(c(2*sigma1*exp(-(norm(A1-A2))^2/(2*l1^2)),0,0))%*%
        psi_3_2%*%psi_2_2%*%psi_1_2
      cov5= t(psi_3%*%psi_2%*%psi_1)%*%
        diag(c(0,2*sigma2*exp(-(norm(A1-A2))^2/(2*l2^2)),0))%*%
        psi_3_2%*%psi_2_2%*%psi_1_2
      cov6= t(psi_3%*%psi_2%*%psi_1)%*%
        diag(c(0,0,2*sigma3*exp(-(norm(A1-A2))^2/(2*l3^2))))%*%
        psi_3_2%*%psi_2_2%*%psi_1_2
      
      x1=x2=xtr
      
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
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l1%*%inv_K%*%c(
    ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l1),
    -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma1%*%inv_K%*%c(
      ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_sigma1),
    -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l2%*%inv_K%*%(
      c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_l2),
    
    -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma2%*%inv_K%*%(
      c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_sigma2),
    -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l3%*%inv_K%*%c(
      ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l3),
    -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma3%*%inv_K%*%(
      c(ytr[,1],ytr[,2],ytr[,3]))+Trace(inv_K%*%K_sigma3)
  )
  
  return(grad)
  
  
}



special_rotation_angles=rbind(rep(2*pi,3),
                              c(2.0855206959894113, -5.855942700212618, 1.7635657246972212),
                              c(0.3806374921400987, -2.493801425248193, 5.4444235829280405),
                              c(2.111937675196926, -1.0284754309126032, 4.835224398435446),
                              c(3.567474808732965, -2.551963613501586, 0.4896806445127213),
                              c(4.550038933593765, -2.694427117348039, 4.257027634446941),
                              c(1.8939230198156405, -1.1832644670865402, 3.6014066278768766),
                              c(4.831141257674098, -0.44935010319804064, 5.947625845273707),
                              c(3.8091501738418696, -0.9957447632831032, 5.845416470491046),
                              c(2.3493017469783752, -6.132935370986087, 5.983612965600155),
                              c(3.1422456290913, -3.730433851277036, 5.776461351230881))
coordinate_transform=read.table("coordinate_transform.txt")

cov_mat3d_2=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  x1=as.numeric(x1);x2=as.numeric(x2)
  dat=list(original_struct,rot_0,rot_1,rot_2,
              rot_3,rot_4,rot_5,
              rot_6,rot_7,rot_8,rot_9)
  for(i in 1:11){
    if(sum(x1%in%(as.matrix(dat[[i]])))==9){
      ang=special_rotation_angles[i,]
      phi1=ang[1];phi2=ang[2];phi3=ang[3]
      
      rho_x1=matrix(c(1,0,0,0,cos(phi1),-sin(phi1),0,sin(phi1),cos(phi1)),
                 nrow=3,byrow = 1)%*%matrix(c(cos(phi2),0,-sin(phi2),0,1,0,
                                              sin(phi2),0,cos(phi2)),
                                            nrow=3,byrow = 1)%*%
        matrix(c(cos(phi3),-sin(phi3),0,sin(phi3),cos(phi3),0,0,0,1),
               nrow=3,byrow = 1)
     
      
      special_rep_x1=rho_x1%*%as.numeric(coordinate_transform[which(abs(as.numeric(
        apply(as.matrix(dat[[i]])[,4:12]-x1,1,sum)))<=1e-14),-1])
    
    }
    if(sum(x2%in%(as.matrix(dat[[i]])))==9){
      
      ang=special_rotation_angles[i,]
      phi1=ang[1];phi2=ang[2];phi3=ang[3]
      
      rho_x2=matrix(c(1,0,0,0,cos(phi1),-sin(phi1),0,sin(phi1),cos(phi1)),
                    nrow=3,byrow = 1)%*%matrix(c(cos(phi2),0,-sin(phi2),0,1,0,
                                                 sin(phi2),0,cos(phi2)),
                                               nrow=3,byrow = 1)%*%
        matrix(c(cos(phi3),-sin(phi3),0,sin(phi3),cos(phi3),0,0,0,1),
               nrow=3,byrow = 1)
      
      special_rep_x2=rho_x2%*%as.numeric(coordinate_transform[which(abs(as.numeric(
        apply(as.matrix(dat[[i]])[,4:12]-x2,1,sum)))<=1e-14),-1])
      
    }
  }
  norm=function(x){sum(x^2)^.5}
  x1=special_rep_x1;x2=special_rep_x2
  return(exp(diag(c(-(norm(x1)-norm(x2))^2/(2*l1^2),
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


log_likelihood_standard2=function(ytr,xtr,l1,sigma1,l2,sigma2,l3,sigma3){
  if(length(ytr)>3){
    if(ncol(ytr)==3){
      ytr=c(ytr[,1],ytr[,2],ytr[,3])
    }
  }
  
  
  ntr=ifelse(length(xtr)==9,1,nrow(xtr))
  Ktr=cross_cov_mat3d_2(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*ntr)
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-ntr*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}


jit=1e-8
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
      
      x1=as.numeric(x1[i,]);x2=as.numeric(x2[j,])
      dat=list(original_struct,rot_0,rot_1,rot_2,
               rot_3,rot_4,rot_5,
               rot_6,rot_7,rot_8,rot_9)
      for(l in 1:11){
        if(sum(x1%in%(as.matrix(dat[[l]])))==9){
          ang=special_rotation_angles[l,]
          phi1=ang[1];phi2=ang[2];phi3=ang[3]
          
          rho_x1=matrix(c(1,0,0,0,cos(phi1),-sin(phi1),0,sin(phi1),cos(phi1)),
                        nrow=3,byrow = 1)%*%matrix(c(cos(phi2),0,-sin(phi2),0,1,0,
                                                     sin(phi2),0,cos(phi2)),
                                                   nrow=3,byrow = 1)%*%
            matrix(c(cos(phi3),-sin(phi3),0,sin(phi3),cos(phi3),0,0,0,1),
                   nrow=3,byrow = 1)
          
          
          special_rep_x1=rho_x1%*%as.numeric(coordinate_transform[which(abs(as.numeric(
            apply(as.matrix(dat[[l]])[,4:12]-x1,1,sum)))<=1e-15),-1])
          
        }
        if(sum(x2%in%(as.matrix(dat[[l]])))==9){
          
          ang=special_rotation_angles[l,]
          phi1=ang[1];phi2=ang[2];phi3=ang[3]
          
          rho_x2=matrix(c(1,0,0,0,cos(phi1),-sin(phi1),0,sin(phi1),cos(phi1)),
                        nrow=3,byrow = 1)%*%matrix(c(cos(phi2),0,-sin(phi2),0,1,0,
                                                     sin(phi2),0,cos(phi2)),
                                                   nrow=3,byrow = 1)%*%
            matrix(c(cos(phi3),-sin(phi3),0,sin(phi3),cos(phi3),0,0,0,1),
                   nrow=3,byrow = 1)
          
          special_rep_x2=rho_x2%*%as.numeric(coordinate_transform[which(abs(as.numeric(
            apply(as.matrix(dat[[l]])[,4:12]-x2,1,sum)))<=1e-15),-1])
          
        }
      }
      
      x1=special_rep_x1;x2=special_rep_x2
      
      cov=exp(diag(c(-(norm(x1)-norm(x2))^2/(2*l1^2),
                     -(norm(x1)-norm(x2))^2/(2*l2^2),
                     -(norm(x1)-norm(x2))^2/(2*l3^2))))*
        diag(c(sigma1^2,sigma2^2,sigma3^2))
      
      cov1= diag(c(exp(-(norm(x1)-norm(x2))^2/(2*l1^2)),0,0))*
        diag(c(sigma1^2,sigma2^2,sigma3^2))*(norm(x1)-norm(x2))^2/(l1^3)
      
      cov2= diag(c(0,exp(-(norm(x1)-norm(x2))^2/(2*l2^2)),0))*
        diag(c(sigma1^2,sigma2^2,sigma3^2))*(norm(x1)-norm(x2))^2/(l2^3)
      
      
      cov3= diag(c(0,0,exp(-(norm(x1)-norm(x2))^2/(2*l3^2))))*
        diag(c(sigma1^2,sigma2^2,sigma3^2))*(norm(x1)-norm(x2))^2/(l3^3)
      
      cov4= diag(c(exp(-(norm(x1)-norm(x2))^2/(2*l1^2)),0,0))*2*sigma1
      
      
      cov5= diag(c(0,exp(-(norm(x1)-norm(x2))^2/(2*l2^2)),0))*2*sigma2
      
      cov6= diag(c(0,0,exp(-(norm(x1)-norm(x2))^2/(2*l3^2))))*2*sigma3
      
      x1=x2=xtr
      
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


opt_par=Adam(function(x) {
  log_likelihood_standard2(
    Ytr,Xtr, x[1], x[2], x[3], x[4], x[5],x[6])
},function(x) {
  grad_standard2(Xtr, Ytr, x[1], x[2], x[3], x[4], x[5],
                 x[6])
},rep(1,6),lr=.1,max_iter=20)


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


(rmse2=mean(apply((cbind(posterior_mean2[1:(length(posterior_mean2)/3)],
                         posterior_mean2[(length(posterior_mean2)/3+1):
                                           (2*length(posterior_mean2)/3)],
                         posterior_mean2[(2*length(posterior_mean2)/3+1):
                                           (length(posterior_mean2))])-Yte)^2,1,sum)
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







res=rbind(c(rmse2),c(LogS2))

save(list = "res", file = paste0("Special_transform_related_sample_adam", id, ".rda"))




