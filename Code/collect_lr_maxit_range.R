l <- (1:2420)
#l=c(1:5,1000)
results <- vector("list", length(l))

for (id in l) {
  file_name <- paste0("rmses_LR_MAXITER_range_1_2_", id, ".rda")
  
  if (file.exists(file_name)) {
    load(file_name)
    results[[id]] <- res
  } else {
    cat("File not found:", file_name, "\n")
  }
}

results <- do.call(rbind, results)
save(list = "results", file = "results_lr_maxiter_range_1_2_.rda")
