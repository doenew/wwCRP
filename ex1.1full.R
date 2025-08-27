# 
setwd("F:\\Marginally coupled mixture desgins")
source("wwCRP.R")
source("Gibbssamplingformixtureswithlinearconstraints.R")
source("uniformcriteriaformixtures.R")
source("classify_rows.R")
source("ITM.R")

library(DoE.base)
library(vcd)
library(MaxPro)
# 3 factors 2 levels  
base_design <- fac.design(
  nfactors = 3,
  nlevels = c(2,2,2),
  replications = 8,
  randomize = FALSE
)
Pz0=data.matrix(base_design)
Pz = Pz0[,-4]

p=3
a=c(0,0,0)
b=c(1,1,1)
# trainning sample for new methods
ts <- Gibbs_l(a,b,n = 100000)


unifullD105=c()
unifullD205=c()
unifullD305=c()
unifullD108=c()
unifullD208=c()
unifullD308=c()
unifullMPD=c()
unifullFFFD=c()

library(foreach)
library(doParallel)

# 定义unicri函数
unicri <- function(design, region){
  mud <- sqrt(sum((colMeans(design) - colMeans(region))^2))
  sdd <- sqrt(sum((apply(design, 2, sd) - apply(region, 2, sd))^2))
  tdesign <- t(design)
  drd <- apply(region, 1, function(x){
    min((colSums((x - tdesign) * (x - tdesign))))
  })
  RMSD <- sqrt(mean(drd))
  MaD <- max(sqrt(drd))
  MiD <- min(dist(design))
  
  return(data.frame(mu = mud, sd = sdd, RMSD = RMSD, MaD = MaD, MiD = MiD))
}

# 设置基础路径
base_path <- "F:/Marginally coupled mixture desgins/ex1.1/"

# 设置并行计算
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# 导出必要的变量和函数到所有节点
clusterExport(cl, c("unicri", "ts", "base_path"))

# 并行计算
result_list <- foreach(kk = 1:50, .combine = 'c', .packages = c()) %dopar% {
  # 使用完整路径读取文件
  fullD105 <- read.csv(paste0(base_path, "1wCDM0.5", kk, ".csv"))
  fullD205 <- read.csv(paste0(base_path, "2wCDM0.5", kk, ".csv"))
  fullD305 <- read.csv(paste0(base_path, "3wCDM0.5", kk, ".csv"))
  fullD108 <- read.csv(paste0(base_path, "1wCDM0.8", kk, ".csv"))
  fullD208 <- read.csv(paste0(base_path, "2wCDM0.8", kk, ".csv"))
  fullD308 <- read.csv(paste0(base_path, "3wCDM0.8", kk, ".csv"))
  fullMPD  <- read.csv(paste0(base_path, "MPD", kk, ".csv"))
  fullFFFD <- read.csv(paste0(base_path, "FFFD", kk, ".csv"))
  
  # 计算unicri并返回列表
  list(
    unifullD105 = unicri(fullD105, ts$sample),
    unifullD205 = unicri(fullD205, ts$sample),
    unifullD305 = unicri(fullD305, ts$sample),
    unifullD108 = unicri(fullD108, ts$sample),
    unifullD208 = unicri(fullD208, ts$sample),
    unifullD308 = unicri(fullD308, ts$sample),
    unifullMPD  = unicri(fullMPD, ts$sample),
    unifullFFFD = unicri(fullFFFD[,1:3], ts$sample)
  )
}

# 停止并行集群
stopCluster(cl)

# 修复extract_and_combine函数
extract_and_combine <- function(result_list, index) {
  # 提取每个结果的指定索引
  items <- result_list[seq(index, length(result_list), by = 8)]
  # 合并所有数据框
  do.call(rbind, items)
}

# 提取并合并结果
unifullD105 <- extract_and_combine(result_list, 1)
unifullD205 <- extract_and_combine(result_list, 2)
unifullD305 <- extract_and_combine(result_list, 3)
unifullD108 <- extract_and_combine(result_list, 4)
unifullD208 <- extract_and_combine(result_list, 5)
unifullD308 <- extract_and_combine(result_list, 6)
unifullMPD  <- extract_and_combine(result_list, 7)
unifullFFFD <- extract_and_combine(result_list, 8)

for(j in 1:5){
  results64=data.frame(D105=unifullD105[,j],
                       D205= unifullD205[,j],
                       D305= unifullD205[,j],
                       D108=unifullD108[,j],
                       D208 =unifullD208[,j],
                       D308 =unifullD308[,j],
                       MPD = unifullMPD[,j],
                       FFFD=unifullFFFD[,j])
  boxplot( results64)
  write.csv(results64,paste0(base_path, "re64", j, ".csv"))
}

