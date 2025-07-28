l=(1:10000)
train_inds=1:20
RMSEs_standard=LogS_standard=matrix(nrow=length(l),ncol=length(train_inds))
RMSEs_matern1=LogS_matern1=matrix(nrow=length(l),ncol=length(train_inds))
RMSEs_matern2=LogS_matern2=matrix(nrow=length(l),ncol=length(train_inds))
RMSEs_matern3=LogS_matern3=matrix(nrow=length(l),ncol=length(train_inds))
RMSEs_matern4=LogS_matern4=RMSEs_matern5=LogS_matern5=
  RMSEs_matern6=LogS_matern6=matrix(nrow=length(l),ncol=length(train_inds))

for (id in l) {
  file_name <- paste0(
    "Spatial_stats_project", id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    RMSEs_standard[id,]=as.numeric(res[[1]])
    LogS_standard[id,]=as.numeric(res[[2]])
    RMSEs_matern1[id,]=as.numeric(res[[3]][1,])
    LogS_matern1[id,]=as.numeric(res[[4]][1,])
    RMSEs_matern2[id,]=as.numeric(res[[3]][2,])
    LogS_matern2[id,]=as.numeric(res[[4]][2,])
    
    RMSEs_matern3[id,]=as.numeric(res[[3]][3,])
    LogS_matern3[id,]=as.numeric(res[[4]][3,])
    
    RMSEs_matern4[id,]=as.numeric(res[[3]][4,])
    LogS_matern4[id,]=as.numeric(res[[4]][4,])
    
    RMSEs_matern5[id,]=as.numeric(res[[3]][5,])
    LogS_matern5[id,]=as.numeric(res[[4]][5,])
    
    RMSEs_matern6[id,]=as.numeric(res[[3]][6,])
    LogS_matern6[id,]=as.numeric(res[[4]][6,])
    
    
  } else {
    cat("File not found:", file_name, "\n")
  }
}
results <- list(RMSEs_standard,LogS_standard,
                RMSEs_matern1,LogS_matern1,
                RMSEs_matern2,LogS_matern2,
                RMSEs_matern3,LogS_matern3,
                RMSEs_matern4,LogS_matern4,RMSEs_matern5,LogS_matern5,
                  RMSEs_matern6,LogS_matern6)
save(list = "results", file =
       "Spatial_stats_project.rda")



for(i in 1:155){
file_name <- paste0(
    "SS_LOOCV", id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    RMSEs_standard[id,]=as.numeric(res[[1]])
    LogS_standard[id,]=as.numeric(res[[2]])
    RMSEs_matern1[id,]=as.numeric(res[[3]][1,])
    LogS_matern1[id,]=as.numeric(res[[4]][1,])
    RMSEs_matern2[id,]=as.numeric(res[[3]][2,])
    LogS_matern2[id,]=as.numeric(res[[4]][2,])

    RMSEs_matern3[id,]=as.numeric(res[[3]][3,])
    LogS_matern3[id,]=as.numeric(res[[4]][3,])

    RMSEs_matern4[id,]=as.numeric(res[[3]][4,])
    LogS_matern4[id,]=as.numeric(res[[4]][4,])

    RMSEs_matern5[id,]=as.numeric(res[[3]][5,])
    LogS_matern5[id,]=as.numeric(res[[4]][5,])

    RMSEs_matern6[id,]=as.numeric(res[[3]][6,])
    LogS_matern6[id,]=as.numeric(res[[4]][6,])


  } else {
    cat("File not found:", file_name, "\n")
  }
}
results <- list(RMSEs_standard,LogS_standard,
                RMSEs_matern1,LogS_matern1,
                RMSEs_matern2,LogS_matern2,
                RMSEs_matern3,LogS_matern3,
                RMSEs_matern4,LogS_matern4,RMSEs_matern5,LogS_matern5,
                  RMSEs_matern6,LogS_matern6)
save(list = "results", file =
       "Spatial_stats_project.rda")


}
