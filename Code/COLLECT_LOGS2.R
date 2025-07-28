
l=(1:6000)

results <- vector("list", length(l))
for (id in l) {
file_name <- paste0(
"Scores_ST_equivariant_related_sample_adam"
#"Fully_equivariant_related_sample"
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
 "Scores_ST_related_sample.rda"
#"Scores_fully_equivariant_related_sample.rda"
)


results <- vector("list", length(l))
for (id in l) {
file_name <- paste0(
"Scores_ST_equivariant_related_sample_adam"
#"Fully_equivariant_related_sample"
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
 "Scores_ST_related_sample_adam.rda"
#"Scores_fully_equivariant_related_sample.rda"
)


results <- vector("list", length(l))
for (id in l) {
file_name <- paste0(
"Scores_ST_equivariant_random_sample"
#"Fully_equivariant_related_sample"
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
 "Scores_ST_random_sample.rda"
#"Scores_fully_equivariant_related_sample.rda"
)


results <- vector("list", length(l))
for (id in l) {
file_name <- paste0(
"Scores_ST_equivariant_random_sample_adam"
#"Fully_equivariant_related_sample"
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
 "Scores_ST_random_sample_adam.rda"
#"Scores_fully_equivariant_related_sample.rda"
)

