id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
library(pracma)
library(cubature)
library(progress)
notcluster=is.na(id)
fund_training=0

#lr=0.1;max_iter=1000
lr=0.001;max_iter=1000
jit=1e-8
Nsplit=10


if(max_iter==1000){
  ntrain=6
  train_size_stepsize=c(10,10,10,10,30,30)
}else{
  train_size_stepsize=rep(10,10)
  ntrain=10
}
#max_iter=1


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



Adam=function(f, grad_f, params_init, lr = lr, beta1 = 0.9, beta2 = 0.999,
              epsilon = 1e-8, max_iter = max_iter, tol = 1e-10) {
  
  params = params_init
  
  m = numeric(length(params))
  v = numeric(length(params))
  
  
  iter = 0
  
  if (notcluster) {
    pb = progress_bar$new(
      format = "  Progress [:bar] :percent in :elapsed",
      total = max_iter,
      clear = FALSE,
      width = 60
    )
  }
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
    if (notcluster) {
      pb$tick()
    }
  }
  #print("iter= ")
  #print(iter)
  return(params)
}




fund_cov_mat=function(x1,x2,sigma1,sigma2,sigma3,l1,l2,l3){
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

cross_fund_cov_mat=function(x1,x2,sigma1,sigma2,sigma3,l1,l2,l3){
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
cross_cov_mat3d_2=function(x1,x2,sigma1,sigma2,sigma3,l1,l2,l3){
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



cov_mat3d_2=function(x1,x2,sigma1,sigma2,sigma3,l1,l2,l3){
  norm=function(x){sum(x^2)^.5}
  (exp(diag(c(-(norm(x1-x2))^2/(2*l1^2),
              -(norm(x1-x2))^2/(2*l2^2),
              -(norm(x1-x2))^2/(2*l3^2))))*
      diag(c(sigma1^2,sigma2^2,sigma3^2)))
}


validation_ll_3d=function(yval,ytr,xval,xtr,sigma1,sigma2,sigma3,l1,l2,l3){
  nval=nrow(xval)
  ntr=nrow(xtr)
  norm=function(x){sum(x^2)^.5}
  K_val_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3 * nval)
  
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      x1=xval;x2=xtr
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,
                          norm((x1[i,]-x2[j,])/l2)^2,
                          norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(sigma1^2,sigma2^2,sigma3^2))
      
      K_val_tr[i, j] = cov1[1,1]
      K_val_tr[nval + i,
               ntr + j] = cov1[2,2]
      
      K_val_tr[2*nval + i,
               2*ntr + j] = cov1[3,3]
      
      K_val_tr[nval + i, j] = cov1[2,1]
      K_val_tr[2*nval+i, j] = cov1[3,1]
      
      K_val_tr[i,ntr + j] = cov1[1,2]
      K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      K_val_tr[2*nval + i,
               ntr + j] = cov1[3,2]
      
      K_val_tr[nval + i,
               2*ntr+ j] = cov1[2,3]
      
      
    }
  }
  
  
  K_val= matrix(0, ncol =3 * nval,
                nrow = 3 * nval)
  x1=xval;x2=xval
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,
                          norm((x1[i,]-x2[j,])/l2)^2,
                          norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(sigma1^2,sigma2^2,sigma3^2))
      
      K_val[i, j] = cov1[1,1]
      K_val[nval + i,
            nval + j] = cov1[2,2]
      
      K_val[2*nval + i,
            2*nval + j] = cov1[3,3]
      
      K_val[nval + i, j] = cov1[2,1]
      K_val[2*nval + i, j] = cov1[3,1]
      
      K_val[i,nval + j] = cov1[1,2]
      K_val[i,2*nval + j] = cov1[1,3]
      
      K_val[2*nval + i,
            nval + j] = cov1[3,2]
      
      K_val[nval + i,
            2*nval+ j] = cov1[2,3]
      
      
    }
  }
  
  K_tr= matrix(0, ncol =3 * ntr,
               nrow = 3 * ntr)
  
  x1=xtr;x2=xtr
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,
                          norm((x1[i,]-x2[j,])/l2)^2,
                          norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(sigma1^2,sigma2^2,sigma3^2))
      
      K_tr[i, j] = cov1[1,1]
      K_tr[ntr + i,
           ntr + j] = cov1[2,2]
      
      K_tr[2*ntr + i,
           2*ntr + j] = cov1[3,3]
      
      K_tr[ntr + i, j] = cov1[2,1]
      K_tr[2*ntr + i, j] = cov1[3,1]
      
      K_tr[i,ntr + j] = cov1[1,2]
      K_tr[i,2*ntr + j] = cov1[1,3]
      
      K_tr[2*ntr + i,
           ntr + j] = cov1[3,2]
      
      K_tr[ntr + i,
           2*ntr+ j] = cov1[2,3]
      
      
    }
  }
  
  #  K_tr=K_tr+diag(tau^2,nrow=2*ntr)
  
  ll=norm(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))^2+
    determinant(K_val-K_val_tr%*% solve(K_tr)%*%t(K_val_tr))$modulus
  #ll=determinant(K_val-K_val_tr%*% solve(K_tr)%*%t(K_val_tr))$modulus
  return(ll)
  
  
}


validation_ll_fund_3d=function(yval,ytr,xval,xtr,sigma1,sigma2,sigma3,l1,l2,l3){
  nval=nrow(xval)
  ntr=nrow(xtr)
  norm=function(x){sum(x^2)^.5}
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
    return(list(A1=A1,psi_1=psi_1,psi_2=psi_2,psi_3=psi_3))
  }
  
  
  
  K_val_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*diag(c(sigma1^2,sigma2^2,sigma3^2)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      K_val_tr[i, j] = cov1[1,1]
      K_val_tr[nval + i,
               ntr + j] = cov1[2,2]
      
      K_val_tr[2*nval + i,
               2*ntr + j] = cov1[3,3]
      
      K_val_tr[nval + i, j] = cov1[2,1]
      K_val_tr[2*nval+i, j] = cov1[3,1]
      
      K_val_tr[i,ntr + j] = cov1[1,2]
      K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      K_val_tr[2*nval + i,
               ntr + j] = cov1[3,2]
      
      K_val_tr[nval + i,
               2*ntr+ j] = cov1[2,3]
      
      
    }
  }
  
  
  K_val= matrix(0, ncol =3 * nval,
                nrow = 3 * nval)
  x1=xval;x2=xval
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(sigma1^2,sigma2^2,sigma3^2)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      K_val[i, j] = cov1[1,1]
      K_val[nval + i,
            nval + j] = cov1[2,2]
      
      K_val[2*nval + i,
            2*nval + j] = cov1[3,3]
      
      K_val[nval + i, j] = cov1[2,1]
      K_val[2*nval + i, j] = cov1[3,1]
      
      K_val[i,nval + j] = cov1[1,2]
      K_val[i,2*nval + j] = cov1[1,3]
      
      K_val[2*nval + i,
            nval + j] = cov1[3,2]
      
      K_val[nval + i,
            2*nval+ j] = cov1[2,3]
      
      
    }
  }
  
  K_tr= matrix(0, ncol =3 * ntr,
               nrow = 3 * ntr)
  
  
  x1=xtr;x2=xtr
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(sigma1^2,sigma2^2,sigma3^2)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      K_tr[i, j] = cov1[1,1]
      K_tr[ntr + i,
           ntr + j] = cov1[2,2]
      
      K_tr[2*ntr + i,
           2*ntr + j] = cov1[3,3]
      
      K_tr[ntr + i, j] = cov1[2,1]
      K_tr[2*ntr + i, j] = cov1[3,1]
      
      K_tr[i,ntr + j] = cov1[1,2]
      K_tr[i,2*ntr + j] = cov1[1,3]
      
      K_tr[2*ntr + i,
           ntr + j] = cov1[3,2]
      
      K_tr[ntr + i,
           2*ntr+ j] = cov1[2,3]
      
      
    }
  }
  
  K_tr=K_tr+diag(jit,nrow=3*ntr)
  K_val=K_val+diag(jit,nrow=3*nval)
  
  ll=norm(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))^2+
    determinant(K_val-K_val_tr%*% solve(K_tr)%*%t(K_val_tr))$modulus
  #ll=determinant(K_val-K_val_tr%*% solve(K_tr)%*%t(K_val_tr))$modulus
  return(ll)
  
  
}

grad_val_3d=function(yval,ytr,xval,xtr,sigma1,sigma2,sigma3,l1,l2,l3){
  nval=nrow(xval)
  ntr=nrow(xtr)
  norm=function(x){sum(x^2)^.5}
  
  K_val_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,
                          norm((x1[i,]-x2[j,])/l2)^2,
                          norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(sigma1^2,sigma2^2,sigma3^2))
      
      K_val_tr[i, j] = cov1[1,1]
      K_val_tr[nval + i,
               ntr + j] = cov1[2,2]
      
      K_val_tr[2*nval + i,
               2*ntr + j] = cov1[3,3]
      
      K_val_tr[nval + i, j] = cov1[2,1]
      K_val_tr[2*nval + i, j] = cov1[3,1]
      
      K_val_tr[i,ntr + j] = cov1[1,2]
      K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      K_val_tr[2*nval + i,
               ntr + j] = cov1[3,2]
      
      K_val_tr[nval + i,
               2*ntr+ j] = cov1[2,3]
      
      
    }
  }
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,0,0)
      ))*
        diag(c(sigma1^2*norm(xval[i,]-xtr[j,])^2/(l1^3),
               0,0))
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval+i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  K_val= matrix(0, ncol =3 * nval,
                nrow = 3 * nval)
  x1=xval;x2=xval
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,
                          norm((x1[i,]-x2[j,])/l2)^2,
                          norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(sigma1^2,sigma2^2,sigma3^2))
      
      K_val[i, j] = cov1[1,1]
      K_val[nval + i,
            nval + j] = cov1[2,2]
      
      K_val[2*nval + i,
            2*nval + j] = cov1[3,3]
      
      K_val[nval + i, j] = cov1[2,1]
      K_val[2*nval + i, j] = cov1[3,1]
      
      K_val[i,nval + j] = cov1[1,2]
      K_val[i,2*nval + j] = cov1[1,3]
      
      K_val[2*nval + i,
            nval + j] = cov1[3,2]
      
      K_val[nval + i,
            2*nval+ j] = cov1[2,3]
      
      
    }
  }
  
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,0,0)
      ))*
        diag(c(sigma1^2*norm(x1[i,]-x2[j,])^2/(l1^3),
               0,0))
      
      
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  K_tr= matrix(0, ncol =3 * ntr,
               nrow = 3 * ntr)
  
  x1=xtr;x2=xtr
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,
                          norm((x1[i,]-x2[j,])/l2)^2,
                          norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(sigma1^2,sigma2^2,sigma3^2))
      
      K_tr[i, j] = cov1[1,1]
      K_tr[ntr + i,
           ntr + j] = cov1[2,2]
      
      K_tr[2*ntr + i,
           2*ntr + j] = cov1[3,3]
      
      K_tr[ntr + i, j] = cov1[2,1]
      K_tr[2*ntr + i, j] = cov1[3,1]
      
      K_tr[i,ntr + j] = cov1[1,2]
      K_tr[i,2*ntr + j] = cov1[1,3]
      
      K_tr[2*ntr + i,
           ntr + j] = cov1[3,2]
      
      K_tr[ntr + i,
           2*ntr+ j] = cov1[2,3]
      
      
    }
  }
  
  
  
  
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,0,0)
      ))*
        diag(c(sigma1^2*norm(x1[i,]-x2[j,])^2/(l1^3),
               0,0))
      
      
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  trace=function(x){sum(diag(x))}
  
  grad_l1=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(0,norm((x1[i,]-x2[j,])/l2)^2,0)
      ))*
        diag(c(0,sigma2^2*norm(x1[i,]-x2[j,])^2/(l2^3),
               0))
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval+i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xval;x2=xval
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      cov1=exp(-.5*diag(c(0,norm((x1[i,]-x2[j,])/l2)^2,0)
      ))*
        diag(c(0,sigma2^2*norm(x1[i,]-x2[j,])^2/(l2^3),
               0))
      
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xtr;x2=xtr
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(0,norm((x1[i,]-x2[j,])/l2)^2,0)
      ))*diag(c(0,sigma2^2*norm(x1[i,]-x2[j,])^2/(l2^3),
                0))
      
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  
  
  grad_l2=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(0,0,norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(0,0,sigma3^2*norm(x1[i,]-x2[j,])^2/(l3^3)
        ))
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval + i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xval;x2=xval
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      cov1=exp(-.5*diag(c(0,0,norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(0,0,sigma3^2*norm(x1[i,]-x2[j,])^2/(l3^3)
        ))
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xtr;x2=xtr
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(0,0,norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(0,0,sigma3^2*norm(x1[i,]-x2[j,])^2/(l3^3)
        )) 
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  
  
  grad_l3=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,0,0)
      ))*
        diag(c(sigma1*2,0,0
        ))
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval + i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xval;x2=xval
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,0,0)
      ))*
        diag(c(sigma1*2,0,0
        ))
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xtr;x2=xtr
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(norm((x1[i,]-x2[j,])/l1)^2,0,0)
      ))*
        diag(c(sigma1*2,0,0
        ))
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  
  
  grad_sigma1=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(0,norm((x1[i,]-x2[j,])/l2)^2,0)
      ))*
        diag(c(0,sigma2*2,0
        ))
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval + i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xval;x2=xval
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      cov1=exp(-.5*diag(c(0,norm((x1[i,]-x2[j,])/l2)^2,0)
      ))*
        diag(c(0,sigma2*2,0
        ))
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xtr;x2=xtr
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(0,norm((x1[i,]-x2[j,])/l2)^2,0)
      ))*
        diag(c(0,sigma2*2,0
        ))
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  
  
  grad_sigma2=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(0,0,norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(0,0,sigma3*2
        ))
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval + i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xval;x2=xval
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      cov1=exp(-.5*diag(c(0,0,norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(0,0,sigma3*2
        ))
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xtr;x2=xtr
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      cov1=exp(-.5*diag(c(0,0,norm((x1[i,]-x2[j,])/l3)^2)
      ))*
        diag(c(0,0,sigma3*2
        ))
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  
  
  grad_sigma3=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  
  
  
  return(c(grad_sigma1,
           grad_sigma2,
           grad_sigma3,
           grad_l1,
           grad_l2,
           grad_l3
  ))
  
  
}







grad_val_fund_3d=function(yval,ytr,xval,xtr,sigma1,sigma2,sigma3,l1,l2,l3){
  nval=nrow(xval)
  ntr=nrow(xtr)
  norm=function(x){sum(x^2)^.5}
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
    return(list(A1=A1,psi_1=psi_1,psi_2=psi_2,psi_3=psi_3))
  }
  
  K_val_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(sigma1^2,sigma2^2,sigma3^2)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      K_val_tr[i, j] = cov1[1,1]
      K_val_tr[nval + i,
               ntr + j] = cov1[2,2]
      
      K_val_tr[2*nval + i,
               2*ntr + j] = cov1[3,3]
      
      K_val_tr[nval + i, j] = cov1[2,1]
      K_val_tr[2*nval + i, j] = cov1[3,1]
      
      K_val_tr[i,ntr + j] = cov1[1,2]
      K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      K_val_tr[2*nval + i,
               ntr + j] = cov1[3,2]
      
      K_val_tr[nval + i,
               2*ntr+ j] = cov1[2,3]
      
      
    }
  }
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*diag(c(sigma1^2*norm((A1-A2))^2/(l1^3),0,0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval+i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  K_val= matrix(0, ncol =3 * nval,
                nrow = 3 * nval)
  x1=xval;x2=xval
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*diag(c(sigma1^2,sigma2^2,sigma3^2)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      K_val[i, j] = cov1[1,1]
      K_val[nval + i,
            nval + j] = cov1[2,2]
      
      K_val[2*nval + i,
            2*nval + j] = cov1[3,3]
      
      K_val[nval + i, j] = cov1[2,1]
      K_val[2*nval + i, j] = cov1[3,1]
      
      K_val[i,nval + j] = cov1[1,2]
      K_val[i,2*nval + j] = cov1[1,3]
      
      K_val[2*nval + i,
            nval + j] = cov1[3,2]
      
      K_val[nval + i,
            2*nval+ j] = cov1[2,3]
      
      
    }
  }
  
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*diag(c(sigma1^2*norm((A1-A2))^2/(l1^3),0,0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  K_tr= matrix(0, ncol =3 * ntr,
               nrow = 3 * ntr)
  x1=xtr;x2=xtr
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(sigma1^2,sigma2^2,sigma3^2)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      K_tr[i, j] = cov1[1,1]
      K_tr[ntr + i,
           ntr + j] = cov1[2,2]
      
      K_tr[2*ntr + i,
           2*ntr + j] = cov1[3,3]
      
      K_tr[ntr + i, j] = cov1[2,1]
      K_tr[2*ntr + i, j] = cov1[3,1]
      
      K_tr[i,ntr + j] = cov1[1,2]
      K_tr[i,2*ntr + j] = cov1[1,3]
      
      K_tr[2*ntr + i,
           ntr + j] = cov1[3,2]
      
      K_tr[ntr + i,
           2*ntr+ j] = cov1[2,3]
      
      
    }
  }
  
  K_tr=K_tr+diag(jit,nrow=3*ntr)
  K_val=K_val+diag(jit,nrow=3*nval)
  
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*diag(c(sigma1^2*norm((A1-A2))^2/(l1^3),0,0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  trace=function(x){sum(diag(x))}
  
  grad_l1=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*diag(c(0,sigma2^2*norm((A1-A2))^2/(l2^3),0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval+i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xval;x2=xval
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,sigma2^2*norm((A1-A2))^2/(l2^3),0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xtr;x2=xtr
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,sigma2^2*norm((A1-A2))^2/(l2^3),0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  
  
  grad_l2=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,0,sigma3^2*norm((A1-A2))^2/(l3^3))))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval + i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xval;x2=xval
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,0,sigma3^2*norm((A1-A2))^2/(l3^3))))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xtr;x2=xtr
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,0,sigma3^2*norm((A1-A2))^2/(l3^3))))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  
  
  grad_l3=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(sigma1*2,0,0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval + i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xval;x2=xval
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(sigma1*2,0,0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xtr;x2=xtr
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(sigma1*2,0,0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  
  
  grad_sigma1=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,sigma2*2,0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval + i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xval;x2=xval
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,sigma2*2,0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xtr;x2=xtr
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,sigma2*2,0)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  
  
  grad_sigma2=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  
  
  del_K_val_tr= matrix(0, ncol =3 * ntr,
                       nrow = 3 * nval)
  x1=xval;x2=xtr
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,0,sigma3*2)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      
      del_K_val_tr[i, j] = cov1[1,1]
      del_K_val_tr[nval + i,
                   ntr + j] = cov1[2,2]
      
      del_K_val_tr[2*nval + i,
                   2*ntr + j] = cov1[3,3]
      
      del_K_val_tr[nval + i, j] = cov1[2,1]
      del_K_val_tr[2*nval + i, j] = cov1[3,1]
      
      del_K_val_tr[i,ntr + j] = cov1[1,2]
      del_K_val_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_val_tr[2*nval + i,
                   ntr + j] = cov1[3,2]
      
      del_K_val_tr[nval + i,
                   2*ntr+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xval;x2=xval
  
  del_K_val= matrix(0, ncol =3 * nval,
                    nrow = 3 * nval)
  for(i in 1:nval){
    #if(i%%100==0){print(i)}
    
    for (j in 1:nval){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,0,sigma3*2)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      
      del_K_val[i, j] = cov1[1,1]
      del_K_val[nval + i,
                nval + j] = cov1[2,2]
      
      del_K_val[2*nval + i,
                2*nval + j] = cov1[3,3]
      
      del_K_val[nval + i, j] = cov1[2,1]
      del_K_val[2*nval + i, j] = cov1[3,1]
      
      del_K_val[i,nval + j] = cov1[1,2]
      del_K_val[i,2*nval + j] = cov1[1,3]
      
      del_K_val[2*nval + i,
                nval + j] = cov1[3,2]
      
      del_K_val[nval + i,
                2*nval+ j] = cov1[2,3]
      
      
      
    }
  }
  
  
  
  x1=xtr;x2=xtr
  del_K_tr= matrix(0, ncol =3 * ntr,
                   nrow = 3* ntr)
  
  for(i in 1:ntr){
    #if(i%%100==0){print(i)}
    
    for (j in 1:ntr){
      
      p1=projection_fund(x1[i,])
      p2=projection_fund(x2[j,])
      
      psi_3=p1$psi_3; psi_2=p1$psi_2;psi_1=p1$psi_1;A1=p1$A1
      psi_3_2=p2$psi_3; psi_2_2=p2$psi_2;psi_1_2=p2$psi_1;A2=p2$A1
      
      cov1=t(psi_3%*%psi_2%*%psi_1)%*%(exp(-.5*diag(c(norm((A1-A2)/l1)^2,
                                                      norm((A1-A2)/l2)^2,
                                                      norm((A1-A2)/l3)^2)
      ))*
        diag(c(0,0,sigma3*2)))%*%psi_3_2%*%psi_2_2%*%psi_1_2
      
      
      
      del_K_tr[i, j] = cov1[1,1]
      del_K_tr[ntr + i,
               ntr + j] = cov1[2,2]
      
      del_K_tr[2*ntr + i,
               2*ntr + j] = cov1[3,3]
      
      del_K_tr[ntr + i, j] = cov1[2,1]
      del_K_tr[2*ntr + i, j] = cov1[3,1]
      
      del_K_tr[i,ntr + j] = cov1[1,2]
      del_K_tr[i,2*ntr + j] = cov1[1,3]
      
      del_K_tr[2*ntr + i,
               ntr + j] = cov1[3,2]
      
      del_K_tr[ntr + i,
               2*ntr+ j] = cov1[2,3]
      
      
      
      
      
      
    }
  }
  
  
  
  grad_sigma3=-2*t(c(yval[,1],yval[,2],yval[,3])-K_val_tr%*%solve(K_tr)%*%(c(ytr[,1],ytr[,2],ytr[,3])))%*%
    ((del_K_val_tr%*%solve(K_tr)-K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr))%*%c(ytr[,1],ytr[,2],ytr[,3]))+
    trace(solve(K_val-K_val_tr%*%
                  solve(K_tr)%*%t(K_val_tr))%*%(
                    del_K_val-del_K_val_tr%*%solve(K_tr)%*%t(K_val_tr)+
                      K_val_tr%*%solve(K_tr)%*%del_K_tr%*%solve(K_tr)%*%t(K_val_tr)-
                      K_val_tr%*%solve(K_tr)%*%t(del_K_val_tr)
                  ))
  
  
  
  
  
  
  return(c(grad_sigma1,
           grad_sigma2,
           grad_sigma3,
           grad_l1,
           grad_l2,
           grad_l3
  ))
  
  
}


agg_validation_ll_3d=function(ytr,xtr,sigma1,sigma2,sigma3,l1,l2,l3,ratio=.25,nsplit=10,seed=1){
  set.seed(seed)
  ntr=nrow(xtr)
  
  nval=floor(ratio*ntr)
  
  inds=sapply(1:nsplit,function(x){sample(ntr,nval)})
  
  mean(apply(inds,2,function(x){
    validation_ll_3d(ytr[x,],ytr[-x,],xtr[x,],xtr[-x,],sigma1,sigma2,sigma3,l1,l2,l3)
  }))
  
}



agg_validation_ll_fund_3d=function(ytr,xtr,sigma1,sigma2,sigma3,l1,l2,l3,ratio=.25,nsplit=10,seed=1){
  set.seed(seed)
  ntr=nrow(xtr)
  
  nval=floor(ratio*ntr)
  
  inds=sapply(1:nsplit,function(x){sample(ntr,nval)})
  
  mean(apply(inds,2,function(x){
    validation_ll_fund_3d(ytr[x,],ytr[-x,],xtr[x,],xtr[-x,],sigma1,sigma2,sigma3,l1,l2,l3)
  }))
  
}

grad_agg_validation_ll_3d=function(ytr,xtr,sigma1,sigma2,sigma3,l1,l2,l3,
                                   ratio=.25,nsplit=10,seed=1){
  set.seed(seed)
  ntr=nrow(xtr)
  
  nval=floor(ratio*ntr)
  
  inds=sapply(1:nsplit,function(x){sample(ntr,nval)})
  
  apply((apply(inds,2,function(x){
    grad_val_3d(ytr[x,],ytr[-x,],xtr[x,],xtr[-x,],sigma1,sigma2,sigma3,l1,l2,l3)
  })),1,mean)
  
}


grad_agg_validation_ll_fund_3d=function(ytr,xtr,sigma1,sigma2,sigma3,l1,l2,l3,
                                        ratio=.25,nsplit=10,seed=1){
  set.seed(seed)
  ntr=nrow(xtr)
  
  nval=floor(ratio*ntr)
  
  inds=sapply(1:nsplit,function(x){sample(ntr,nval)})
  
  apply((apply(inds,2,function(x){
    grad_val_fund_3d(ytr[x,],ytr[-x,],xtr[x,],xtr[-x,],sigma1,sigma2,sigma3,l1,l2,l3)
  })),1,mean)
  
}


if(fund_training){
  RMSES=LogS=matrix(nrow=ntrain,ncol=2);Opt_par=list()
  colnames(RMSES)<-c("standard","fundamental")
  colnames(LogS)<-c("standard","fundamental")
  
  
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
  
  data_set_size=nrow(X)
  train_ind=sample(data_set_size,10)
  tic()
  for( i in 1:ntrain){
    
    train_ind=c(train_ind,sample((1:data_set_size)[-train_ind],10))
    test_ind=(1:data_set_size)[-train_ind]
    
    Xtr=X[train_ind,]
    Xte=X[test_ind,]
    Ytr=Y[train_ind,]
    Yte=Y[test_ind,]
    seed=sample(10000,1)
    
    OPT_PAR=matrix(nrow=6,ncol=ncol(RMSES))
    
    opt_par=Adam(function(x) {
      agg_validation_ll_3d(
        Ytr,Xtr, x[1], x[2], x[3], x[4], x[5],x[6],nsplit=Nsplit,seed=seed)
    },function(x) {grad_agg_validation_ll_3d(Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],
                                             x[6],nsplit = Nsplit,seed=seed)},rep(1,6),lr=lr,
    max_iter = max_iter)
    
    OPT_PAR[,1]=opt_par
    
    Ktetr=cross_cov_mat3d_2(Xte,Xtr,opt_par[1],
                            opt_par[2],
                            opt_par[3],
                            opt_par[4],
                            opt_par[5],
                            opt_par[6])
    
    Ktetr_trtr_inv=Ktetr%*%inv(cross_cov_mat3d_2(Xtr,Xtr,opt_par[1],
                                                 opt_par[2],
                                                 opt_par[3],
                                                 opt_par[4],
                                                 opt_par[5],
                                                 opt_par[6])+
                                 diag(jit,nrow=3*nrow(Xtr)))
    
    
    
    posterior_mean=Ktetr_trtr_inv%*%c(Ytr[,1],Ytr[,2],Ytr[,3])
    
    
    (RMSES[i,1]=mean(apply((cbind(posterior_mean[1:(length(posterior_mean)/3)],
                                  posterior_mean[(length(posterior_mean)/3+1):
                                                   (2*length(posterior_mean)/3)],
                                  posterior_mean[(2*length(posterior_mean)/3+1):
                                                   (length(posterior_mean))])-Yte)^2,1,sum)
    )^.5)
    
    
    posterior_cov= cross_cov_mat3d_2(Xte,Xte,opt_par[1],
                                     opt_par[2],
                                     opt_par[3],
                                     opt_par[4],opt_par[5],
                                     opt_par[6])-
      Ktetr_trtr_inv%*%t(Ktetr)
    
    
    (LogS[i,1]=(t(c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean)%*%
                  solve(posterior_cov
                        +diag(1e-6,3*nrow(Xte))
                  )%*%
                  (c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean)
                +determinant(posterior_cov
                             +diag(1e-6,3*nrow(Xte))
                )$modulus[1]+
                  log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    
    
    opt_par=Adam(function(x) {
      agg_validation_ll_fund_3d(
        Ytr,Xtr, x[1], x[2], x[3], x[4], x[5],x[6],nsplit=Nsplit,seed=seed)
    },function(x) {
      grad_agg_validation_ll_fund_3d(Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],
                                     x[6],nsplit=Nsplit,seed=seed)
    },rep(1,6),lr=lr,max_iter=max_iter)
    
    OPT_PAR[,2]=opt_par
    
    Ktetr2=cross_fund_cov_mat(Xte,Xtr,opt_par[1],
                              opt_par[2],
                              opt_par[3],
                              opt_par[4],
                              opt_par[5],
                              opt_par[6])
    
    Ktetr_trtr_inv2=Ktetr2%*%solve(cross_fund_cov_mat(Xtr,Xtr,opt_par[1],
                                                      opt_par[2],
                                                      opt_par[3],
                                                      opt_par[4],
                                                      opt_par[5],
                                                      opt_par[6])+diag(jit,nrow=3*nrow(Xtr)))
    
    
    
    posterior_mean=posterior_mean2=Ktetr_trtr_inv2%*%c(Ytr[,1],Ytr[,2],Ytr[,3])
    
    
    (RMSES[i,2]=mean(apply((cbind(posterior_mean2[1:(length(posterior_mean)/3)],
                                  posterior_mean2[(length(posterior_mean)/3+1):
                                                    (2*length(posterior_mean)/3)],
                                  posterior_mean2[(2*length(posterior_mean)/3+1):
                                                    (length(posterior_mean))])-Yte)^2,1,sum)
    )^.5)
    
    
    posterior_cov2= cross_fund_cov_mat(Xte,Xte,opt_par[1],
                                       opt_par[2],
                                       opt_par[3],
                                       opt_par[4],opt_par[5],
                                       opt_par[6])-
      Ktetr_trtr_inv2%*%t(Ktetr2)
    
    
    (LogS[i,2]=(t(c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean2)%*%
                  solve(posterior_cov2+diag(1e-6,3*nrow(Xte))
                  )%*%
                  (c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean2)+
                  determinant(posterior_cov2+diag(1e-6,3*nrow(Xte))
                  )$modulus[1]+
                  log(2*pi)*nrow(Xte))/nrow(Xte))
    
    
    
    Opt_par[[i]]=OPT_PAR
  }
  toc()
  res=list(RMSES,LogS,OPT_PAR)
  save(list = "res", file = paste0("learning_curves_chemical_VAL_LL_", id, ".rda"))
}else{
  

RMSES=LogS=numeric(ntrain);Opt_par=matrix(ncol=ntrain,nrow=6)


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

data_set_size=nrow(X)
train_ind=sample(data_set_size,10)
tic()
print(c(lr,max_iter,ntrain,nrow(X),Nsplit))
for( i in 1:ntrain){
  
  train_ind=c(train_ind,sample((1:data_set_size)[-train_ind],10))
  test_ind=(1:data_set_size)[-train_ind]
  
  Xtr=X[train_ind,]
  Xte=X[test_ind,]
  Ytr=Y[train_ind,]
  Yte=Y[test_ind,]
  seed=sample(10000,1)
  
  
  
  opt_par=Adam(function(x) {
    agg_validation_ll_3d(
      Ytr,Xtr, x[1], x[2], x[3], x[4], x[5],x[6],nsplit=Nsplit,seed=seed)
  },function(x) {grad_agg_validation_ll_3d(Ytr, Xtr, x[1], x[2], x[3], x[4], x[5],
                                           x[6],nsplit = Nsplit,seed=seed)},rep(1,6),lr=lr,
  max_iter = max_iter)
  Opt_par[,i]=opt_par
  
  
  Ktetr=cross_cov_mat3d_2(Xte,Xtr,opt_par[1],
                          opt_par[2],
                          opt_par[3],
                          opt_par[4],
                          opt_par[5],
                          opt_par[6])
  
  Ktetr_trtr_inv=Ktetr%*%inv(cross_cov_mat3d_2(Xtr,Xtr,opt_par[1],
                                               opt_par[2],
                                               opt_par[3],
                                               opt_par[4],
                                               opt_par[5],
                                               opt_par[6])+
                               diag(jit,nrow=3*nrow(Xtr)))
  
  
  
  posterior_mean=Ktetr_trtr_inv%*%c(Ytr[,1],Ytr[,2],Ytr[,3])
  
  
  (RMSES[i]=mean(apply((cbind(posterior_mean[1:(length(posterior_mean)/3)],
                                posterior_mean[(length(posterior_mean)/3+1):
                                                 (2*length(posterior_mean)/3)],
                                posterior_mean[(2*length(posterior_mean)/3+1):
                                                 (length(posterior_mean))])-Yte)^2,1,sum)
  )^.5)
  
  
  posterior_cov= cross_cov_mat3d_2(Xte,Xte,opt_par[1],
                                   opt_par[2],
                                   opt_par[3],
                                   opt_par[4],opt_par[5],
                                   opt_par[6])-
    Ktetr_trtr_inv%*%t(Ktetr)
  
  
  (LogS[i]=(t(c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean)%*%
                solve(posterior_cov
                      +diag(1e-6,3*nrow(Xte))
                )%*%
                (c(Yte[,1],Yte[,2],Yte[,3])-posterior_mean)
              +determinant(posterior_cov
                           +diag(1e-6,3*nrow(Xte))
              )$modulus[1]+
                log(2*pi)*nrow(Xte))/nrow(Xte))
  
  
print(i)  
}
toc()
res=list(RMSES,LogS,Opt_par)
save(list = "res", file = paste0("learning_curves_chemical_VAL_LL_", id, ".rda"))
}
