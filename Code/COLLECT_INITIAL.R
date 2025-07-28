
l=(1:8000)
results <- 0
i=0
for (id in l) {
file_name <- paste0("Differences_", id, ".rda")
  if (file.exists(file_name)) {
    tryCatch({load(file_name)
         i=i+1
                         
         results=results+ res
         },error=function(x){print("error")}
         
) } else {
    cat("File not found:", file_name, "\n")
  }
}
results <- results/i
save(list = "results", file = "Differences.rda")
