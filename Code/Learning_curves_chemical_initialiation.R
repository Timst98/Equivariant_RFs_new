id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

library(pracma)
library(doParallel)
library(cubature)
library(foreach)
library(progress)

Export=c('inv')
Packages=c()

initial_par=c(runif(2,0.001,5))
notcluster=is.na(id)
samplesize=50
data_set_size=1261-410
train_ind=sample(data_set_size,10)
ntrain=10

#lr=0.1;max_iter=1000
lr=0.001;max_iter=1000
jit=1e-8
#lr=0.1;max_iter=1;ntrain=1;train_ind=sample(data_set_size,10)
#lr=0.1;max_iter=100;ntrain=2;test_ind=sample(851,10)
#lr=0.1;max_iter=100

if(max_iter==1000){
  ntrain=6;x_values=c(20,30,40,50,80,110)
}


X_names=c("X.X_O","X.Y_O","X.Z_O","X.X_H1","X.Y_H1","X.Z_H1",
          "X.X_H2","X.Y_H2","X.Z_H2")
Y_names=c("Y.X_DIPOLE_MP2","Y.Y_DIPOLE_MP2","Y.Z_DIPOLE_MP2")

data <- readLines("results_mp2/original_struct.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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

data <- readLines("results_mp2/original_struct.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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

data <- readLines("results_mp2/rot_0.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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



data <- readLines("results_mp2/rot_1.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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


data <- readLines("results_mp2/rot_2.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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



data <- readLines("results_mp2/rot_3.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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

data <- readLines("results_mp2/rot_4.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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


data <- readLines("results_mp2/rot_5.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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



data <- readLines("results_mp2/rot_6.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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



data <- readLines("results_mp2/rot_7.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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


data <- readLines("results_mp2/rot_8.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
    loc3=loc3[!is.na(loc3)]
    
    X_H2 <- append(X_H2, loc3[1])
    Y_H2 <- append(Y_H2, loc3[2])
    Z_H2 <- append(Z_H2, loc3[3])
  }
}

# Combine the vectors into matrices
Y <- cbind(X_DIPOLE_MP2, Y_DIPOLE_MP2, Z_DIPOLE_MP2)
X <- cbind(X_O, Y_O, Z_O, X_H1, Y_H1, Z_H1, X_H2, Y_H2, Z_H2)

rot_8=data.frame(Y=Y,X=X)


data <- readLines("results_mp2/rot_9.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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



data <- readLines("results_mp2/rot_3.xyz")

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
    
    loc1=as.numeric(unlist(strsplit(data[i+12], " ")))
    loc1=loc1[!is.na(loc1)]
    X_O <- append(X_O, loc1[1])
    Y_O <- append(Y_O, loc1[2])
    Z_O <- append(Z_O, loc1[3])
    
    loc2= as.numeric(unlist(strsplit(data[i+13], " ")))
    loc2=loc2[!is.na(loc2)]
    X_H1 <- append(X_H1,loc2[1])
    Y_H1 <- append(Y_H1, loc2[2])
    Z_H1 <- append(Z_H1, loc2[3])
    
    loc3=as.numeric(unlist(strsplit(data[i+14], " ")))
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



######## REMOVE DUPLICATES

projection_fund=function(x1){
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
  return(A1)
}
X=as.matrix(original_struct[X_names])
Y=as.matrix(original_struct[Y_names])

A=t(apply(X,1,projection_fund))
#A=matrix(rnorm(3*3641,mean(A),sd(A)),ncol=3)


a=diag(1,nrow(A),nrow(A))
for(i in 1:nrow(A)){
  for(j in 1:nrow(A)){
    if(i!=j){a[i,j]=sum((A[i,]-A[j,])^2)^.5}
  }
}

### DUPLICATES CHECK
threshold=1e-10
duplicates=c()
for(i in 1:1261){
  ind=which(a[i,]<threshold)
  duplicates=c(duplicates,ind[ind>i])
}
####################################################################



fund_cov_mat_scalar=function(x1,x2,l,sigma){
  norm=function(x){sum(x^2)^.5}
  if(norm(x1[4:6]-x1[1:3])>norm(x1[7:9]-x1[1:3])){
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
  
  
  if(norm(x2[4:6]-x2[1:3])>norm(x2[7:9]-x2[1:3])){
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
  
  sigma^2*exp(-(norm(A1-A2))^2/(2*l^2))*t(psi_3%*%psi_2%*%psi_1)%*%psi_3_2%*%psi_2_2%*%psi_1_2
  
  
}

cross_fund_cov_mat_scalar=function(x1,x2,l,sigma){
  norm=function(x){sum(x^2)^.5}
  cov= matrix(0, ncol =3 * ifelse(length(as.matrix(x2)) == 9, 1,nrow(x2)),
              nrow = 3 * ifelse(length(as.matrix(x1)) == 9, 
                                1, nrow(x1)))
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      if(length(as.matrix(x2)) == 9){
        if(length(as.matrix(x1)) == 9){
          cov1= fund_cov_mat_scalar(x1,x2,l,sigma)
          
        }else{
          cov1= fund_cov_mat_scalar(x1[i,],x2,l,sigma)
          
        }
      }else if(length(as.matrix(x1)) == 9){
        cov1= fund_cov_mat_scalar(x1,x2[j,],l,sigma)
        
      }else{
        cov1= fund_cov_mat_scalar(x1[i,],x2[j,],l,sigma)
        
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

log_likelihood_fund_scalar=function(ytr,xtr,l,sigma){
  if(length(ytr)>3){
    if(ncol(ytr)==3){
      ytr=c(ytr[,1],ytr[,2],ytr[,3])
    }
  }
  
  
  ntr=ifelse(length(xtr)==9,1,nrow(xtr))
  Ktr=cross_fund_cov_mat_scalar(xtr,xtr,l,sigma)+
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

cov_mat3d_scalar=function(x1,x2,l,sigma){
  norm=function(x){sum(x^2)^.5}
  
  diag(c(sigma^2*exp(-norm(x1-x2)^2/(2*l^2)),
         sigma^2*exp(-norm(x1-x2)^2/(2*l^2)),
         sigma^2*exp(-norm(x1-x2)^2/(2*l^2))))
}



cross_cov_mat3d_scalar=function(x1,x2,l,sigma){
  cov= matrix(0, ncol =3 * ifelse(length(as.matrix(x2)) == 9, 1,nrow(x2)),
              nrow = 3 * ifelse(length(as.matrix(x1)) == 9, 1, nrow(x1)))
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      
      cov1= cov_mat3d_scalar(x1[i,],x2[j,],l,sigma)
      
      
      
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





grad_fund_scalar=function(xtr,ytr,l,sigma){
  Trace=function(x){sum(diag(x))}
  x1=x2=xtr
  norm=function(x){sum(x^2)^.5}
  
  K=cross_fund_cov_mat_scalar(xtr,xtr,l,sigma)+
    diag(jit,nrow=3*nrow(xtr))
  
  K_l=K_sigma=
    matrix(0, ncol =3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)),
           nrow = 3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)))
  
  inv_K=inv(K)
  
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      x1=xtr[i,]
      x2=xtr[j,]
      if(norm(x1[4:6]-x1[1:3])>norm(x1[7:9]-x1[1:3])){
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
      
      
      if(norm(x2[4:6]-x2[1:3])>norm(x2[7:9]-x2[1:3])){
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
      
      
      
      cov1=sigma^2*exp(-(norm(A1-A2))^2/(2*l^2))/(l^3)*(norm(A1-A2))^2*t(psi_3%*%psi_2%*%psi_1)%*%
        psi_3_2%*%psi_2_2%*%psi_1_2
      
      cov2=2*sigma*exp(-(norm(A1-A2))^2/(2*l^2))*t(psi_3%*%psi_2%*%psi_1)%*%
        psi_3_2%*%psi_2_2%*%psi_1_2
      
      x1=x2=xtr
      
      K_l[i, j] = cov1[1,1]
      K_l[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,2]
      
      K_l[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,3]
      
      K_l[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[2,1]
      K_l[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[3,1]
      
      K_l[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,2]
      K_l[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,3]
      
      K_l[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,2]
      
      K_l[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,3]
      
      
      
      K_sigma[i, j] = cov2[1,1]
      K_sigma[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
              ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,2]
      
      K_sigma[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
              2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,3]
      
      K_sigma[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[2,1]
      K_sigma[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[3,1]
      
      K_sigma[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,2]
      K_sigma[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,3]
      
      K_sigma[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
              ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,2]
      
      K_sigma[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
              2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,3]
      
    }
  }
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l%*%inv_K%*%c(
    ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l),
    -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma%*%inv_K%*%c(
      ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_sigma)
  )
  
  return(grad)
  
  
}


grad_standard_scalar=function(xtr,ytr,l,sigma){
  Trace=function(x){sum(diag(x))}
  x1=x2=xtr
  
  norm=function(x){sum(x^2)^.5}
  
  K=cross_cov_mat3d_scalar(xtr,xtr,l,sigma)+
    diag(jit,nrow=3*nrow(xtr))
  
  
  K_l = K_sigma =
    matrix(0, ncol =3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)),
           nrow = 3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)))
  
  inv_K=solve(K)
  
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      x1=x1[i,];x2=x2[j,]
      cov1= cov_mat3d_scalar(x1,x2,l,sigma)*(norm(x1-x2))^2/(l^3)
      
      cov2= cov_mat3d_scalar(x1,x2,l,1)*sigma*2
      
      x1=x2=xtr
      
      K_l[i, j] = cov1[1,1]
      K_l[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,2]
      
      K_l[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,3]
      
      K_l[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[2,1]
      K_l[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov1[3,1]
      
      K_l[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,2]
      K_l[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[1,3]
      
      K_l[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[3,2]
      
      K_l[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
          2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov1[2,3]
      
      
      
      K_sigma[i, j] = cov2[1,1]
      K_sigma[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
              ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,2]
      
      K_sigma[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
              2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,3]
      
      K_sigma[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[2,1]
      K_sigma[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i, j] = cov2[3,1]
      
      K_sigma[i,ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,2]
      K_sigma[i,2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[1,3]
      
      K_sigma[2*ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
              ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[3,2]
      
      K_sigma[ifelse(length(as.matrix(x1))==9,1,nrow(x1)) + i,
              2*ifelse(length(as.matrix(x2))==9,1,nrow(x2)) + j] = cov2[2,3]
      
    }
  }
  
  grad=0.5*c(-t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_l%*%inv_K%*%c(
    ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_l),
    -t(c(ytr[,1],ytr[,2],ytr[,3]))%*%inv_K%*%K_sigma%*%inv_K%*%c(
      ytr[,1],ytr[,2],ytr[,3])+Trace(inv_K%*%K_sigma)
  )
  
  return(grad)
  
  
  
}


log_likelihood_standard_scalar=function(ytr,xtr,l,sigma){
  if(length(ytr)>3){
    if(ncol(ytr)==3){
      ytr=c(ytr[,1],ytr[,2],ytr[,3])
    }
  }
  
  
  ntr=ifelse(length(xtr)==9,1,nrow(xtr))
  Ktr=cross_cov_mat3d(xtr,xtr,l,sigma)+
    diag(jit,nrow=3*ntr)
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*determinant(Ktr)$modulus[1]-ntr*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}


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

log_likelihood_fund=function(ytr,xtr,l1,sigma1,l2,sigma2,l3,sigma3){
  if(length(ytr)>3){
    if(ncol(ytr)==3){
      ytr=c(ytr[,1],ytr[,2],ytr[,3])
    }
  }
  
  
  ntr=ifelse(length(xtr)==9,1,nrow(xtr))
  Ktr=cross_fund_cov_mat(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*ntr)
  
  
  ll=-0.5*t(ytr)%*%inv(Ktr)%*%ytr-0.5*log(det(Ktr))-ntr*log(2*pi)
  #print(c(ll))
  ll=ifelse(is.nan(ll),10^6,-ll)
  #ll=ifelse(ll>0,ll,1e6)
  
  return(ll)
  
  
}

cov_mat3d=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  norm=function(x){sum(x^2)^.5}
  
  diag(c(sigma1^2*exp(-norm(x1-x2)^2/(2*l1^2)),
         sigma2^2*exp(-norm(x1-x2)^2/(2*l2^2)),
         sigma3^2*exp(-norm(x1-x2)^2/(2*l3^2))))
}

cross_cov_mat3d=function(x1,x2,l1,sigma1,l2,sigma2,l3,sigma3){
  cov= matrix(0, ncol =3 * ifelse(length(as.matrix(x2)) == 9, 1,nrow(x2)),
              nrow = 3 * ifelse(length(as.matrix(x1)) == 9, 1, nrow(x1)))
  
  for(i in 1:(length(as.matrix(x1))/9)){
    for (j in 1:(length(as.matrix(x2))/9)){
      
      cov1= cov_mat3d(x1[i,],x2[j,],l1,sigma1,l2,sigma2,l3,sigma3)
      
      
      
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

grad_fund=function(xtr,ytr,l1,sigma1,l2,sigma2,l3,sigma3){
  Trace=function(x){sum(diag(x))}
  x1=x2=xtr
  norm=function(x){sum(x^2)^.5}
  
  K=cross_fund_cov_mat(xtr,xtr,l1,sigma1,l2,sigma2,l3,sigma3)+
    diag(jit,nrow=3*nrow(xtr))
  
  K_l1 =  K_l2 = K_l3 = K_sigma1 = K_sigma2 = K_sigma3 =
    matrix(0, ncol =3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)),
           nrow = 3 * ifelse(length(as.matrix(xtr)) == 9, 1, nrow(xtr)))
  
  inv_K=inv(K)
  
  
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
      x1=x1[i,];x2=x2[j,]
      cov1= cov_mat3d(x1,x2,l1,sigma1,l2,0 ,l3,0)*(norm(x1-x2))^2/(l1^3)
      
      cov2= cov_mat3d(x1,x2,l1,0,l2,sigma2 ,l3,0)*(norm(x1-x2))^2/(l2^3)
      
      cov3= cov_mat3d(x1,x2,l1,0,l2,0 ,l3,sigma3)*(norm(x1-x2))^2/(l3^3)
      
      cov4= cov_mat3d(x1,x2,l1,1,l2,0,l3,0)*2*sigma1
      
      
      cov5= cov_mat3d(x1,x2,l1,0,l2,1,l3,0)*2*sigma2
      
      cov6= cov_mat3d(x1,x2,l1,0,l2,0,l3,1)*2*sigma3
      
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




dat=list(original_struct,rot_0,rot_1,rot_2,rot_3,rot_4,rot_5,rot_6,rot_7,rot_8,rot_9)


rot_ind=sample(10,1261,replace = TRUE)
swap_ind=sample(2,1261,replace = TRUE)


X=(t(sapply(1:1261,function(i){
  if(swap_ind[i]==1){
    as.matrix(dat[[rot_ind[i]]][X_names])[i,]
  }else{
    as.matrix(dat[[rot_ind[i]]][X_names])[i,c(1:3,7:9,4:6)]
  }
}))+t(sapply(1:1261,function(i){
  a=runif(3,-2,2)
  rep(a,3)
})))[-duplicates,]


Y=t(sapply(1:1261,function(i){
  as.matrix(dat[[rot_ind[i]]][Y_names])[i,]
}))[-duplicates,]





if(max_iter==1000){
  train_size_stepsize=c(10,10,10,10,30,30)
}else{
  train_size_stepsize=rep(10,10)
}

numcores=30
cl=makeCluster(numcores)
registerDoParallel(cl)

Res=foreach(i=1:numcores,.combine='rbind',.export = Export,.packages=Packages)%dopar%{
  res=numeric(8)
  train_ind=sample(data_set_size,samplesize)
  test_ind=(1:data_set_size)[-train_ind]
  
  Xtr=X[train_ind,]
  Xte=X[test_ind,]
  Ytr=Y[train_ind,]
  Yte=Y[test_ind,]
  
  
 
  opt_par=Adam(function(x) {
    log_likelihood_standard_scalar(
      Ytr,Xtr,x[1], x[2])
  },function(x) {
    grad_standard_scalar(Xtr, Ytr, x[1], x[2])
  },initial_par,lr=lr,max_iter=max_iter)
  
  
  res[1:2]=opt_par
  Ktetr2=cross_cov_mat3d_scalar(Xte,Xtr,opt_par[1],
                                opt_par[2])
  
  Ktetr_trtr_inv2=Ktetr2%*%solve(cross_cov_mat3d_scalar(Xtr,Xtr,opt_par[1],
                                                        opt_par[2])+diag(jit,nrow=3*nrow(Xtr)))
  
  
  
  posterior_mean=posterior_mean2=Ktetr_trtr_inv2%*%c(Ytr[,1],Ytr[,2],Ytr[,3])
  
  
  (rmse2=mean(apply((cbind(posterior_mean2[1:(length(posterior_mean)/3)],
                                      posterior_mean2[(length(posterior_mean)/3+1):
                                                        (2*length(posterior_mean)/3)],
                                      posterior_mean2[(2*length(posterior_mean)/3+1):
                                                        (length(posterior_mean))])-Yte)^2,1,sum)
  )^.5)
  
  
  posterior_cov2= cross_cov_mat3d_scalar(Xte,Xte,opt_par[1],
                                         opt_par[2])-
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
  
  
  res[3]=log(rmse2)
  res[4]=LogS2
  
  
  opt_par=Adam(function(x) {
    log_likelihood_fund_scalar(
      Ytr,Xtr, x[1], x[2])
  },function(x) {grad_fund_scalar(Xtr, Ytr, x[1], x[2])},
  initial_par,lr=lr,max_iter = max_iter)
  
  res[5:6]=opt_par
  
  Ktetr_eq=cross_fund_cov_mat_scalar(Xte,Xtr,opt_par[1],
                                     opt_par[2])
  
  Ktetr_trtr_inv_eq=Ktetr_eq%*%inv(cross_fund_cov_mat_scalar(Xtr,Xtr,opt_par[1],
                                                             opt_par[2])+
                                     diag(jit,nrow=3*nrow(Xtr)))
  
  
  
  posterior_mean_eq=Ktetr_trtr_inv_eq%*%c(Ytr[,1],Ytr[,2],Ytr[,3])
  
  
  (rmse_eq=mean(apply((cbind(posterior_mean_eq[1:(length(posterior_mean_eq)/3)],
                                        posterior_mean_eq[(length(posterior_mean_eq)/3+1):
                                                            (2*length(posterior_mean_eq)/3)],
                                        posterior_mean_eq[(2*length(posterior_mean_eq)/3+1):
                                                            (length(posterior_mean_eq))])-Yte)^2,1,sum)
  )^.5)
  
  
  posterior_cov_eq= cross_fund_cov_mat_scalar(Xte,Xte,opt_par[1],
                                              opt_par[2])-
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
  
  
  
  res[7]=log(rmse_eq)
  res[8]=LogS_eq
  
  as.numeric(res)
  
}
stopCluster(cl)

res=c(initial_par,apply(Res,2,mean))
save(list = "res", file = paste0("learning_curves_chemical_initialization_", id, ".rda"))

