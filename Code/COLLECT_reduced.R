RES=list("Scores_adam_standard_helm_reduced.rda","Scores_adam_mix_reduced.rda",
         "Scores_adam_mix2_reduced.rda","Scores_adam_range_reduced.rda",
         "Scores_adam_range2_reduced.rda","Scores_adam_range_flexible_reduced.rda",

"Scores_adam_standard_helm_GULF.rda","Scores_adam_mix_GULF.rda",
              "Scores_adam_mix2_GULF.rda","Scores_adam_range_gulf.rda",
              "Scores_adam_range2_gulf.rda","Scores_adam_range_flexible_gulf.rda",
              "Scores_adam_flexible_one_center_gulf.rda",
              "Scores_adam_range_one_center_norange_gulf.rda",
              "Scores_adam_flexible_one_center_gulf.rda")

results <- vector("list", length(RES))
for (id in RES) {
  file_name <- paste0(id)
  if (file.exists(file_name)) {
    load(file_name)
    results[[id]] <- res
  } else {
    cat("File not found:", file_name, "\n")
  }
}
results <- do.call(rbind, results)
save(list = "results", file = "Adam_Scores_reduced.rda")
