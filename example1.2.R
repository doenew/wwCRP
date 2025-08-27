setwd("F:/Marginally coupled mixture desgins")
source("wwCRP.R")
source("Gibbssamplingformixtureswithlinearconstraints.R")
source("uniformcriteriaformixtures.R")
source("classify_rows.R")
source("InverseTransformMethod.R")
#library(randtoolbox)
library(DoE.base)
library(vcd)
library(MaxPro)
# 2 factors 2,5 levels  
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
bp=0.5 # balance parameter

for(kk in 1:50){
 # D3=wwCRPParallel(Pz,w = 3,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D2=wwCRPParallel(Pz,w = 2,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D1=wwCRPParallel(Pz,w = 1,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
 # mypath3 <- paste0("F:\\Marginally coupled mixture desgins\\ex1.2\\3wCDM0.5",kk,".csv")
#  write.csv(D3,mypath3,row.names = F)
  mypath2 <- paste0("F:\\Marginally coupled mixture desgins\\ex1.2\\2wCDM0.5",kk,".csv")
  write.csv(D2,mypath2,row.names = F)
  mypath1 <- paste0("F:\\Marginally coupled mixture desgins\\ex1.2\\1wCDM0.5",kk,".csv")
  write.csv(D1,mypath1,row.names = F)
  # MPD
  rand_design_part1=apply(matrix(rep(seq(from=0,to=1,length=nrow(Pz)),p-1),ncol=p-1),2,sample)
  InitialDesign=cbind(rand_design_part1,Pz)
  obj=MaxProQQ(InitialDesign, p_nom=3) 
  MPD=ITM(obj$Design[,1:(p-1)],a,b)
  mypathMPD <- paste0("F:\\Marginally coupled mixture desgins\\ex1.2\\MPD",kk,".csv")
  write.csv(MPD,mypathMPD,row.names = F)
}

bp=0.8

for(kk in 1:50){
#  D3=wwCRPParallel(Pz,w = 3,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D2=wwCRPParallel(Pz,w = 2,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D1=wwCRPParallel(Pz,w = 1,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
 # mypath3 <- paste0("F:\\Marginally coupled mixture desgins\\ex1.2\\3wCDM0.8",kk,".csv")
 # write.csv(D3,mypath3,row.names = F)
  mypath2 <- paste0("F:\\Marginally coupled mixture desgins\\ex1.2\\2wCDM0.8",kk,".csv")
  write.csv(D2,mypath2,row.names = F)
  mypath1 <- paste0("F:\\Marginally coupled mixture desgins\\ExampleFFFD\\1wCDM0.8",kk,".csv")
  write.csv(D1,mypath1,row.names = F)
}


base_design <- fac.design(
  nfactors = 2,
  nlevels = c(2,5),
  replications = 10,
  randomize = FALSE
)
Pz0=data.matrix(base_design)
Pz = Pz0[,-3]

p=3
a=c(0.1,0.2,0.1)
b=c(0.7,0.8,0.6)
# trainning sample for new methods
ts <- Gibbs_l(a,b,n = 100000)
bp=0.5 # balance parameter

for(kk in 1:50){
  # D3=wwCRPParallel(Pz,w = 3,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D2=wwCRPParallel(Pz,w = 2,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D1=wwCRPParallel(Pz,w = 1,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  # mypath3 <- paste0("F:\\Marginally coupled mixture desgins\\ex10\\3wCDM0.5",kk,".csv")
  #  write.csv(D3,mypath3,row.names = F)
  mypath2 <- paste0("F:\\Marginally coupled mixture desgins\\ex10\\2wCDM0.5",kk,".csv")
  write.csv(D2,mypath2,row.names = F)
  mypath1 <- paste0("F:\\Marginally coupled mixture desgins\\ex10\\1wCDM0.5",kk,".csv")
  write.csv(D1,mypath1,row.names = F)
  # MPD
  rand_design_part1=apply(matrix(rep(seq(from=0,to=1,length=nrow(Pz)),p-1),ncol=p-1),2,sample)
  InitialDesign=cbind(rand_design_part1,Pz)
  obj=MaxProQQ(InitialDesign, p_nom=3) 
  MPD=ITM(obj$Design[,1:(p-1)],a,b)
  mypathMPD <- paste0("F:\\Marginally coupled mixture desgins\\ex10\\MPD",kk,".csv")
  write.csv(MPD,mypathMPD,row.names = F)
}

bp=0.8

for(kk in 1:50){
  #  D3=wwCRPParallel(Pz,w = 3,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D2=wwCRPParallel(Pz,w = 2,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D1=wwCRPParallel(Pz,w = 1,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  # mypath3 <- paste0("F:\\Marginally coupled mixture desgins\\ex1.2\\3wCDM0.8",kk,".csv")
  # write.csv(D3,mypath3,row.names = F)
  mypath2 <- paste0("F:\\Marginally coupled mixture desgins\\ex10\\2wCDM0.8",kk,".csv")
  write.csv(D2,mypath2,row.names = F)
  mypath1 <- paste0("F:\\Marginally coupled mixture desgins\\ex10\\1wCDM0.8",kk,".csv")
  write.csv(D1,mypath1,row.names = F)
}
