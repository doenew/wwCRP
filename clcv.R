clcv <- function(matrix_data, w) {
  # Check input validity
  if (!is.matrix(matrix_data) && !is.data.frame(matrix_data)) {
    stop("Input must be a matrix or data frame")
  }
  
  if (w > ncol(matrix_data)) {
    stop("w cannot be greater than the number of columns of the matrix")
  }
  
  # Get all possible combinations of w columns
  col_combinations <- combn(ncol(matrix_data), w, simplify = FALSE)
  
  # Calculate the number of distinct level combinations for each column combination
  result_vector <- sapply(col_combinations, function(cols) {
    # Extract the selected columns
    selected_cols <- matrix_data[, cols, drop = FALSE]
    
    # Convert each row to a string representing the level combination
    combinations <- apply(selected_cols, 1, function(row) {
      paste(row, collapse = "-")
    })
    
    # Calculate the number of distinct combinations
    length(unique(combinations))
  })
  
  # Add names to the vector (using column names)
  if (!is.null(colnames(matrix_data))) {
    names(result_vector) <- sapply(col_combinations, function(cols) {
      paste(colnames(matrix_data)[cols], collapse = "-")
    })
  }
  
  return(result_vector)
}
