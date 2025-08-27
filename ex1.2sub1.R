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
  nfactors = 2,
  nlevels = c(2,5),
  replications = 5,
  randomize = FALSE
)
Pz0=data.matrix(base_design)
Pz = Pz0[,-3]

p=3
a=c(0.1,0.2,0.1)
b=c(0.7,0.8,0.6)
# trainning sample for new methods
ts <- Gibbs_l(a,b,n = 100000)


unisub1D105=c()
unisub1D205=c()
unisub1D305=c()
unisub1D108=c()
unisub1D208=c()
unisub1D308=c()
unisub1MPD=c()
unisub1FFFD=c()

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
base_path <- "F:/Marginally coupled mixture desgins/ex1.2/"

# 设置并行计算
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# 导出必要的变量和函数到所有节点
clusterExport(cl, c("unicri", "ts", "base_path"))

# 并行计算
result_list <- foreach(kk = 1:50, .combine = 'c', .packages = c()) %dopar% {
  # 使用完整路径读取文件
  sub1D105 <- read.csv(paste0(base_path, "1wCDM0.5", kk, ".csv"))
  sub1D205 <- read.csv(paste0(base_path, "2wCDM0.5", kk, ".csv"))
  sub1D108 <- read.csv(paste0(base_path, "1wCDM0.8", kk, ".csv"))
  sub1D208 <- read.csv(paste0(base_path, "2wCDM0.8", kk, ".csv"))
  sub1MPD  <- read.csv(paste0(base_path, "MPD", kk, ".csv"))
  sub1FFFD <- read.csv(paste0(base_path, "FFFD", kk, ".csv"))
  
  # 计算unicri并返回列表
  list(
    unisub1D105 = unicri(sub1D105[Pz[,1]==1,],ts$sample)+
    unicri(sub1D105[Pz[,1]==2,], ts$sample)+
    unicri(sub1D105[Pz[,2]==1,], ts$sample)+
    unicri(sub1D105[Pz[,2]==2,], ts$sample)+
    unicri(sub1D105[Pz[,2]==3,], ts$sample)+
    unicri(sub1D105[Pz[,2]==4,], ts$sample)+
    unicri(sub1D105[Pz[,2]==5,], ts$sample),
    unisub1D205= unicri(sub1D205[Pz[,1]==1,],ts$sample)+
    unicri(sub1D205[Pz[,1]==2,], ts$sample)+
    unicri(sub1D205[Pz[,2]==1,], ts$sample)+
    unicri(sub1D205[Pz[,2]==2,], ts$sample)+
    unicri(sub1D205[Pz[,2]==3,], ts$sample)+
    unicri(sub1D205[Pz[,2]==4,], ts$sample)+
    unicri(sub1D205[Pz[,2]==5,], ts$sample),
    unisub1D108 = unicri(sub1D108[Pz[,1]==1,],ts$sample)+
      unicri(sub1D108[Pz[,1]==2,], ts$sample)+
    unicri(sub1D108[Pz[,2]==1,], ts$sample)+
    unicri(sub1D108[Pz[,2]==2,], ts$sample)+
    unicri(sub1D108[Pz[,2]==3,], ts$sample)+
    unicri(sub1D108[Pz[,2]==4,], ts$sample)+
    unicri(sub1D108[Pz[,2]==5,], ts$sample),
    unisub1D208 = unicri(sub1D208[Pz[,1]==1,],ts$sample)+
    unicri(sub1D208[Pz[,1]==2,], ts$sample)+
    unicri(sub1D208[Pz[,2]==1,], ts$sample)+
    unicri(sub1D208[Pz[,2]==2,], ts$sample)+
    unicri(sub1D208[Pz[,2]==3,], ts$sample)+
    unicri(sub1D208[Pz[,2]==4,], ts$sample)+
    unicri(sub1D208[Pz[,2]==5,], ts$sample),
    unisub1MPD  = unicri(sub1MPD[Pz[,1]==1,],ts$sample)+
    unicri(sub1MPD[Pz[,1]==2,], ts$sample)+
    unicri(sub1MPD[Pz[,2]==1,], ts$sample)+
    unicri(sub1MPD[Pz[,2]==2,], ts$sample)+
    unicri(sub1MPD[Pz[,2]==3,], ts$sample)+
    unicri(sub1MPD[Pz[,2]==4,], ts$sample)+
    unicri(sub1MPD[Pz[,2]==5,], ts$sample),
    unisub1FFFD = unicri(sub1FFFD[sub1FFFD[,4]=="L1", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,4]=="L2", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,5]=="L1", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,5]=="L2", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,5]=="L3", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,5]=="L4", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,5]=="L5", 1:3], ts$sample)
  )
}

# 停止并行集群
stopCluster(cl)

# 修复extract_and_combine函数
extract_and_combine <- function(result_list, index) {
  # 提取每个结果的指定索引
  items <- result_list[seq(index, length(result_list), by = 6)]
  # 合并所有数据框
  do.call(rbind, items)
}

# 提取并合并结果
unisub1D105 <- extract_and_combine(result_list, 1)
unisub1D205 <- extract_and_combine(result_list, 2)
unisub1D108 <- extract_and_combine(result_list, 3)
unisub1D208 <- extract_and_combine(result_list, 4)
unisub1MPD  <- extract_and_combine(result_list, 5)
unisub1FFFD <- extract_and_combine(result_list, 6)

for(j in 1:5){
  results50=data.frame(D105=unisub1D105[,j],
                       D205= unisub1D205[,j],
                       D108=unisub1D108[,j],
                       D208 =unisub1D208[,j],
                       MPD = unisub1MPD[,j],
                       FFFD=unisub1FFFD[,j])
  boxplot( results50)
  write.csv(results50,paste0(base_path, "sub1re50", j, ".csv"))
}

### n=100
# 设置基础路径
base_path <- "F:/Marginally coupled mixture desgins/ex10/"

unisub1D105=c()
unisub1D205=c()
unisub1D305=c()
unisub1D108=c()
unisub1D208=c()
unisub1D308=c()
unisub1MPD=c()
unisub1FFFD=c()
# 设置并行计算
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# 导出必要的变量和函数到所有节点
clusterExport(cl, c("unicri", "ts", "base_path"))

# 并行计算
result_list <- foreach(kk = 1:50, .combine = 'c', .packages = c()) %dopar% {
  # 使用完整路径读取文件
  sub1D105 <- read.csv(paste0(base_path, "1wCDM0.5", kk, ".csv"))
  sub1D205 <- read.csv(paste0(base_path, "2wCDM0.5", kk, ".csv"))
  sub1D108 <- read.csv(paste0(base_path, "1wCDM0.8", kk, ".csv"))
  sub1D208 <- read.csv(paste0(base_path, "2wCDM0.8", kk, ".csv"))
  sub1MPD  <- read.csv(paste0(base_path, "MPD", kk, ".csv"))
  sub1FFFD <- read.csv(paste0(base_path, "FFFD", kk, ".csv"))
  
  # 计算unicri并返回列表
  list(
    unisub1D105 = unicri(sub1D105[Pz[,1]==1,],ts$sample)+
      unicri(sub1D105[Pz[,1]==2,], ts$sample)+
      unicri(sub1D105[Pz[,2]==1,], ts$sample)+
      unicri(sub1D105[Pz[,2]==2,], ts$sample)+
      unicri(sub1D105[Pz[,2]==3,], ts$sample)+
      unicri(sub1D105[Pz[,2]==4,], ts$sample)+
      unicri(sub1D105[Pz[,2]==5,], ts$sample),
    unisub1D205= unicri(sub1D205[Pz[,1]==1,],ts$sample)+
      unicri(sub1D205[Pz[,1]==2,], ts$sample)+
      unicri(sub1D205[Pz[,2]==1,], ts$sample)+
      unicri(sub1D205[Pz[,2]==2,], ts$sample)+
      unicri(sub1D205[Pz[,2]==3,], ts$sample)+
      unicri(sub1D205[Pz[,2]==4,], ts$sample)+
      unicri(sub1D205[Pz[,2]==5,], ts$sample),
    unisub1D108 = unicri(sub1D108[Pz[,1]==1,],ts$sample)+
      unicri(sub1D108[Pz[,1]==2,], ts$sample)+
      unicri(sub1D108[Pz[,2]==1,], ts$sample)+
      unicri(sub1D108[Pz[,2]==2,], ts$sample)+
      unicri(sub1D108[Pz[,2]==3,], ts$sample)+
      unicri(sub1D108[Pz[,2]==4,], ts$sample)+
      unicri(sub1D108[Pz[,2]==5,], ts$sample),
    unisub1D208 = unicri(sub1D208[Pz[,1]==1,],ts$sample)+
      unicri(sub1D208[Pz[,1]==2,], ts$sample)+
      unicri(sub1D208[Pz[,2]==1,], ts$sample)+
      unicri(sub1D208[Pz[,2]==2,], ts$sample)+
      unicri(sub1D208[Pz[,2]==3,], ts$sample)+
      unicri(sub1D208[Pz[,2]==4,], ts$sample)+
      unicri(sub1D208[Pz[,2]==5,], ts$sample),
    unisub1MPD  = unicri(sub1MPD[Pz[,1]==1,],ts$sample)+
      unicri(sub1MPD[Pz[,1]==2,], ts$sample)+
      unicri(sub1MPD[Pz[,2]==1,], ts$sample)+
      unicri(sub1MPD[Pz[,2]==2,], ts$sample)+
      unicri(sub1MPD[Pz[,2]==3,], ts$sample)+
      unicri(sub1MPD[Pz[,2]==4,], ts$sample)+
      unicri(sub1MPD[Pz[,2]==5,], ts$sample),
    unisub1FFFD = unicri(sub1FFFD[sub1FFFD[,4]=="L1", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,4]=="L2", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,5]=="L1", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,5]=="L2", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,5]=="L3", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,5]=="L4", 1:3], ts$sample)+
      unicri(sub1FFFD[sub1FFFD[,5]=="L5", 1:3], ts$sample)
  )
}

# 停止并行集群
stopCluster(cl)

# 修复extract_and_combine函数
extract_and_combine <- function(result_list, index) {
  # 提取每个结果的指定索引
  items <- result_list[seq(index, length(result_list), by = 6)]
  # 合并所有数据框
  do.call(rbind, items)
}

# 提取并合并结果
unisub1D105 <- extract_and_combine(result_list, 1)
unisub1D205 <- extract_and_combine(result_list, 2)
unisub1D108 <- extract_and_combine(result_list, 3)
unisub1D208 <- extract_and_combine(result_list, 4)
unisub1MPD  <- extract_and_combine(result_list, 5)
unisub1FFFD <- extract_and_combine(result_list, 6)

for(j in 1:5){
  results100=data.frame(D105=unisub1D105[,j],
                        D205= unisub1D205[,j],
                        D108=unisub1D108[,j],
                        D208 =unisub1D208[,j],
                        MPD = unisub1MPD[,j],
                        FFFD=unisub1FFFD[,j])
  boxplot( results100)
#  write.csv(results100,paste0(base_path, "sub1re100", j, ".csv"))
}




