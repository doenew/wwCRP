clcv <- function(matrix_data, w) {
  # 检查输入有效性
  if (!is.matrix(matrix_data) && !is.data.frame(matrix_data)) {
    stop("输入必须是矩阵或数据框")
  }
  
  if (w > ncol(matrix_data)) {
    stop("w不能大于矩阵的列数")
  }
  
  # 获取所有可能的w列组合
  col_combinations <- combn(ncol(matrix_data), w, simplify = FALSE)
  
  # 计算每种列组合的非重复水平组合数
  result_vector <- sapply(col_combinations, function(cols) {
    # 提取选定的列
    selected_cols <- matrix_data[, cols, drop = FALSE]
    
    # 将每行转换为一个字符串表示的水平组合
    combinations <- apply(selected_cols, 1, function(row) {
      paste(row, collapse = "-")
    })
    
    # 计算非重复组合数
    length(unique(combinations))
  })
  
  # 为向量添加名称（使用列名）
  if (!is.null(colnames(matrix_data))) {
    names(result_vector) <- sapply(col_combinations, function(cols) {
      paste(colnames(matrix_data)[cols], collapse = "-")
    })
  }
  
  return(result_vector)
}