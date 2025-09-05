#' Classify Rows of a Matrix by Level Combinations of Variables
#'
#' @param matrix_data Input matrix or data frame
#' @param cols Optional parameter, specifies column indices or names for classification. If NULL, all columns are used
#' @return A numeric vector indicating the category of each row
#'
#' @examples
#' # Create example matrix
#' mat <- matrix(c(1,1,2,2,1,2,1,2), ncol=2)
#' # Classify
#' category_vector <- classify_rows(mat)
#' print(category_vector)
classify_rows <- function(matrix_data, cols = NULL) {
  # Check input
  if (!is.matrix(matrix_data) && !is.data.frame(matrix_data)) {
    stop("Input must be a matrix or data frame")
  }
  
  # If columns are specified, extract them
  if (!is.null(cols)) {
    if (is.character(cols)) {
      # Select by column names
      if (!all(cols %in% colnames(matrix_data))) {
        stop("Specified column names do not exist")
      }
      matrix_data <- matrix_data[, cols, drop = FALSE]
    } else if (is.numeric(cols)) {
      # Select by column indices
      if (any(cols > ncol(matrix_data)) || any(cols < 1)) {
        stop("Column indices out of range")
      }
      matrix_data <- matrix_data[, cols, drop = FALSE]
    } else {
      stop("cols parameter must be a character vector (column names) or numeric vector (column indices)")
    }
  }
  
  # Convert each row to a string representation of level combinations
  row_combinations <- apply(matrix_data, 1, function(row) {
    paste(row, collapse = "-")
  })
  
  # Create factor and convert to numeric categories
  factor_combinations <- factor(row_combinations)
  category_vector <- as.numeric(factor_combinations)
  
  # Add category labels as attributes
  category_labels <- levels(factor_combinations)
  attr(category_vector, "category_labels") <- category_labels
  
  return(category_vector)
}
