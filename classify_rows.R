#' 将矩阵的行按变量水平组合分类
#'
#' @param matrix_data 输入矩阵或数据框
#' @param cols 可选参数，指定用于分类的列索引或列名。如果为NULL，则使用所有列
#' @return 一个数值向量，表示每行所属的类别
#'
#' @examples
#' # 创建示例矩阵
#' mat <- matrix(c(1,1,2,2,1,2,1,2), ncol=2)
#' # 分类
#' category_vector <- classify_rows(mat)
#' print(category_vector)
classify_rows <- function(matrix_data, cols = NULL) {
  # 检查输入
  if (!is.matrix(matrix_data) && !is.data.frame(matrix_data)) {
    stop("输入必须是矩阵或数据框")
  }
  
  # 如果指定了列，则提取这些列
  if (!is.null(cols)) {
    if (is.character(cols)) {
      # 按列名选择
      if (!all(cols %in% colnames(matrix_data))) {
        stop("指定的列名不存在")
      }
      matrix_data <- matrix_data[, cols, drop = FALSE]
    } else if (is.numeric(cols)) {
      # 按列索引选择
      if (any(cols > ncol(matrix_data)) || any(cols < 1)) {
        stop("列索引超出范围")
      }
      matrix_data <- matrix_data[, cols, drop = FALSE]
    } else {
      stop("cols参数必须是字符向量(列名)或数值向量(列索引)")
    }
  }
  
  # 将每行转换为一个字符串表示的水平组合
  row_combinations <- apply(matrix_data, 1, function(row) {
    paste(row, collapse = "-")
  })
  
  # 创建因子并转换为数值类别
  factor_combinations <- factor(row_combinations)
  category_vector <- as.numeric(factor_combinations)
  
  # 添加类别标签作为属性
  category_labels <- levels(factor_combinations)
  attr(category_vector, "category_labels") <- category_labels
  
  return(category_vector)
}