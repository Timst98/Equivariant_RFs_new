

l=(1:100)

results <- vector("list", length(l))
for (id in l) {
file_name <- paste0(
#"learning_curves_NMF_no_train_translations_long_" # Take from l=9!!
"learning_curves_NMF_"
, id, ".rda")
  if (file.exists(file_name)) {
    load(file_name)
    results[[id]] <- res
  } else {
    cat("File not found:", file_name, "\n")
  }
}
results <- do.call(rbind, results)
save(list = "results", file =

"learning_curves_NMF.rda"
)
