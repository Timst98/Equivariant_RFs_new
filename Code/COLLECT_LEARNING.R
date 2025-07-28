
l=1:400
results_standard=results_eq=results_ST <- matrix(nrow=length(l),ncol=3)
non_existent=c()
for (id in l) {
  file_name <- paste0(
    "Learning_curves_originial_struct_lr_minus1"
    #"Learning_curves_originial_struct_lr_minus2"
    #"Learning_curves_originial_struct_lr_minus3"
    , id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    results_standard[id,] <- res[[1]]
    results_eq[id,] <- res[[2]]
    results_ST[id,] <- res[[3]]
    
    
  } else {
    non_existent=c(non_existent,id)
    cat("File not found:", file_name, "\n")
  }
}
results <- list(results_standard,results_eq,results_ST,non_existent)
save(list = "results", file =
       "Learning_curves_original_struct_01.rda"
     #"Scores_fully_equivariant_related_sample.rda"
)

results_standard=results_eq=results_ST <- matrix(nrow=length(l),ncol=3)
non_existent=c()
for (id in l) {
  file_name <- paste0(
    #"Learning_curves_originial_struct_lr_minus1"
    "Learning_curves_originial_struct_lr_minus2"
    #"Learning_curves_originial_struct_lr_minus3"
    , id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    results_standard[id,] <- res[[1]]
    results_eq[id,] <- res[[2]]
    results_ST[id,] <- res[[3]]
    
    
  } else {
    non_existent=c(non_existent,id)
    cat("File not found:", file_name, "\n")
  }
}
results <- list(results_standard,results_eq,results_ST,non_existent)
save(list = "results", file =
       "Learning_curves_original_struct_001.rda"
     #"Scores_fully_equivariant_related_sample.rda"
)
results_standard=results_eq=results_ST <- matrix(nrow=length(l),ncol=3)
non_existent=c()

for (id in l) {
  file_name <- paste0(
    #"Learning_curves_originial_struct_lr_minus1"
    #"Learning_curves_originial_struct_lr_minus2"
    "Learning_curves_originial_struct_lr_minus3"
    , id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    results_standard[id,] <- res[[1]]
    results_eq[id,] <- res[[2]]
    results_ST[id,] <- res[[3]]
    
    
  } else {
    non_existent=c(non_existent,id)
    cat("File not found:", file_name, "\n")
  }
}
results <- list(results_standard,results_eq,results_ST,non_existent)
save(list = "results", file =
       "Learning_curves_original_struct_0001.rda"
     #"Scores_fully_equivariant_related_sample.rda"
)



l=1:400
results_standard=results_eq=results_ST <- matrix(nrow=length(l),ncol=3)
non_existent=c()
for (id in l) {
  file_name <- paste0(
    "Learning_curves_related_lr_minus1"
    #"Learning_curves_originial_struct_lr_minus2"
    #"Learning_curves_originial_struct_lr_minus3"
    , id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    results_standard[id,] <- res[[1]]
    results_eq[id,] <- res[[2]]
    results_ST[id,] <- res[[3]]
    
    
  } else {
    non_existent=c(non_existent,id)
    cat("File not found:", file_name, "\n")
  }
}
results <- list(results_standard,results_eq,results_ST,non_existent)
save(list = "results", file =
       "Learning_curves_related_01.rda"
     #"Scores_fully_equivariant_related_sample.rda"
)

results_standard=results_eq=results_ST <- matrix(nrow=length(l),ncol=3)
non_existent=c()
for (id in l) {
  file_name <- paste0(
    #"Learning_curves_originial_struct_lr_minus1"
    "Learning_curves_related_lr_minus2"
    #"Learning_curves_originial_struct_lr_minus3"
    , id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    results_standard[id,] <- res[[1]]
    results_eq[id,] <- res[[2]]
    results_ST[id,] <- res[[3]]
    
    
  } else {
    non_existent=c(non_existent,id)
    cat("File not found:", file_name, "\n")
  }
}
results <- list(results_standard,results_eq,results_ST,non_existent)
save(list = "results", file =
       "Learning_curves_related_001.rda"
     #"Scores_fully_equivariant_related_sample.rda"
)
results_standard=results_eq=results_ST <- matrix(nrow=length(l),ncol=3)

non_existent=c()
for (id in l) {
  file_name <- paste0(
    #"Learning_curves_originial_struct_lr_minus1"
    #"Learning_curves_originial_struct_lr_minus2"
    "Learning_curves_related_lr_minus3"
    , id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    results_standard[id,] <- res[[1]]
    results_eq[id,] <- res[[2]]
    results_ST[id,] <- res[[3]]
    
    
  } else {
    non_existent=c(non_existent,id)
    cat("File not found:", file_name, "\n")
  }
}
results <- list(results_standard,results_eq,results_ST, non_existent)
save(list = "results", file =
       "Learning_curves_related_0001.rda"
     #"Scores_fully_equivariant_related_sample.rda"
)



l=1:400
results_standard=results_eq=results_ST <- matrix(nrow=length(l),ncol=3)
    non_existent=c()
for (id in l) {
  file_name <- paste0(
    "Learning_curves_random_lr_minus1"
    #"Learning_curves_originial_struct_lr_minus2"
    #"Learning_curves_originial_struct_lr_minus3"
    , id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    results_standard[id,] <- res[[1]]
    results_eq[id,] <- res[[2]]
    results_ST[id,] <- res[[3]]
    
    
  } else {
    non_existent=c(non_existent,id)
    cat("File not found:", file_name, "\n")
  }
}
results <- list(results_standard,results_eq,results_ST, non_existent)
save(list = "results", file =
       "Learning_curves_random_01.rda"
     #"Scores_fully_equivariant_related_sample.rda"
)

results_standard=results_eq=results_ST <- matrix(nrow=length(l),ncol=3)
    non_existent=c()
for (id in l) {
  file_name <- paste0(
    #"Learning_curves_originial_struct_lr_minus1"
    "Learning_curves_random_lr_minus2"
    #"Learning_curves_originial_struct_lr_minus3"
    , id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    results_standard[id,] <- res[[1]]
    results_eq[id,] <- res[[2]]
    results_ST[id,] <- res[[3]]
    
    
  } else {
    non_existent=c(non_existent,id)
    cat("File not found:", file_name, "\n")
  }
}
results <- list(results_standard,results_eq,results_ST,    non_existent)
save(list = "results", file =
       "Learning_curves_random_001.rda"
     #"Scores_fully_equivariant_related_sample.rda"
)
results_standard=results_eq=results_ST <- matrix(nrow=length(l),ncol=3)

    non_existent=c()
for (id in l) {
  file_name <- paste0(
    #"Learning_curves_originial_struct_lr_minus1"
    #"Learning_curves_originial_struct_lr_minus2"
    "Learning_curves_random_minus3"
    , id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    results_standard[id,] <- res[[1]]
    results_eq[id,] <- res[[2]]
    results_ST[id,] <- res[[3]]
    
    
  } else {
    non_existent=c(non_existent,id)
    cat("File not found:", file_name, "\n")
  }
}
results <- list(results_standard,results_eq,results_ST,    non_existent)
save(list = "results", file =
       "Learning_curves_random_0001.rda"
     #"Scores_fully_equivariant_related_sample.rda"
)

