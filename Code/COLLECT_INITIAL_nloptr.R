
l=1:249
results <- vector("list", length(l))
for (id in l) {
 load(paste0("rmses_DECOMP_inital_nloptr", id, ".rda"))
 results[[id]] <- res
}
results <- do.call(rbind, results)
save(list = "results", file = "results_initial_nloptr.rda")
