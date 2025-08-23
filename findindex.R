# find_identical_rows_with_submatrix_include_self
findindex <- function(Pz, w) {
  n <- nrow(Pz)
  q <- ncol(Pz)
  
  if (w < 1 || w > q) {
    stop("w must be between 1 and the number of columns of Pz")
  }
  
  # Generate all column combinations
  comb_list <- combn(q, w, simplify = FALSE)
  num_combinations <- length(comb_list)
  
  # Initialize result list
  result_list <- vector("list", n)
  for (i in 1:n) {
    result_list[[i]] <- matrix(nrow = 0, ncol = 2)
    colnames(result_list[[i]]) <- c("row_index", "submatrix_index")
  }
  
  # Process each column combination
  for (k in 1:num_combinations) {
    cols <- comb_list[[k]]
    sub_mat <- Pz[, cols, drop = FALSE]
    
    # Convert rows to strings for comparison
    row_strings <- apply(sub_mat, 1, function(x) paste(x, collapse = "\r"))
    
    # For each row, find identical rows (including self)
    for (i in 1:n) {
      identical_rows <- which(row_strings == row_strings[i])
      
      if (length(identical_rows) > 0) {
        # Add results to matrix
        new_entries <- cbind(identical_rows, rep(k, length(identical_rows)))
        result_list[[i]] <- rbind(result_list[[i]], new_entries)
      }
    }
  }
  
  return(result_list)
}

# Create example design matrix with user-specified number of rows
set.seed(123)
n_rows <- 81  # User can change this value to specify the number of rows
Pz <- matrix(sample(c("A", "B", "C"), n_rows * 3, replace = TRUE), nrow = n_rows, ncol = 3)
colnames(Pz) <- c("Var1", "Var2", "Var3")
rownames(Pz) <- paste0("Row", 1:n_rows)

print("Design matrix:")
print(Pz)

# Using the function (including self)
w <- 2
result_include_self <- findindex(Pz, w)

# Print results
print("Results (including self):")
for (i in 1:length(result_include_self)) {
  cat(paste0("\nIdentical row information for row ", i, ":\n"))
  if (nrow(result_include_self[[i]]) > 0) {
    print(result_include_self[[i]])
  } else {
    cat("  No identical rows\n")
  }
}

# Verify correctness of results
cat("\nVerification (including self):\n")
# All column combinations
comb_list <- combn(3, 2, simplify = FALSE)
for (k in 1:length(comb_list)) {
  cols <- comb_list[[k]]
  sub_mat <- Pz[, cols, drop = FALSE]
  cat(paste0("Submatrix ", k, " (Columns", paste(cols, collapse = ","), "):\n"))
  print(sub_mat)
  
  row_strings <- apply(sub_mat, 1, function(x) paste(x, collapse = "\r"))
  for (i in 1:n_rows) {
    identical_rows <- which(row_strings == row_strings[i])
    cat(paste0("  Identical rows for row ", i, ": ", paste(identical_rows, collapse = ","), "\n"))
  }
  cat("\n")
}