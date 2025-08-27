setwd("F:/Marginally coupled mixture desgins")
source("wwCRP.R")
source("Gibbssamplingformixtureswithlinearconstraints.R")
source("uniformcriteriaformixtures.R")
source("classify_rows.R")
source("ITM.R")
#library(randtoolbox)
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
a=rep(0,p)
b=rep(1,p)
# trainning sample for new methods
ts <- Gibbs_l(a,b,n = 100000)
bp=0.5 # balance parameter

for(kk in 1:50){
  D3=wwCRPParallel(Pz,w = 3,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D2=wwCRPParallel(Pz,w = 2,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D1=wwCRPParallel(Pz,w = 1,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  mypath3 <- paste0("F:\\Marginally coupled mixture desgins\\ExampleFFFD\\3wCDM0.5",kk,".csv")
  write.csv(D3,mypath3,row.names = F)
  mypath2 <- paste0("F:\\Marginally coupled mixture desgins\\ExampleFFFD\\2wCDM0.5",kk,".csv")
  write.csv(D2,mypath2,row.names = F)
  mypath1 <- paste0("F:\\Marginally coupled mixture desgins\\ExampleFFFD\\1wCDM0.5",kk,".csv")
  write.csv(D1,mypath1,row.names = F)
  # MPD
  rand_design_part1=apply(matrix(rep(seq(from=0,to=1,length=nrow(Pz)),p-1),ncol=p-1),2,sample)
  InitialDesign=cbind(rand_design_part1,Pz)
  obj=MaxProQQ(InitialDesign, p_nom=3) 
  MPD=ITM(obj$Design[,1:(p-1)])
  mypathMPD <- paste0("F:\\Marginally coupled mixture desgins\\ExampleFFFD\\MPD",kk,".csv")
  write.csv(MPD,mypathMPD,row.names = F)
}

bp=0.8

for(kk in 1:50){
  D3=wwCRPParallel(Pz,w = 3,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D2=wwCRPParallel(Pz,w = 2,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  D1=wwCRPParallel(Pz,w = 1,ts=ts$sample,lambda =bp,avetheta = F,T = 200,ifparallel = T)
  mypath3 <- paste0("F:\\Marginally coupled mixture desgins\\ExampleFFFD\\3wCDM0.8",kk,".csv")
  write.csv(D3,mypath3,row.names = F)
  mypath2 <- paste0("F:\\Marginally coupled mixture desgins\\ExampleFFFD\\2wCDM0.8",kk,".csv")
  write.csv(D2,mypath2,row.names = F)
  mypath1 <- paste0("F:\\Marginally coupled mixture desgins\\ExampleFFFD\\1wCDM0.8",kk,".csv")
  write.csv(D1,mypath1,row.names = F)
}

D3=wwCRPParallel(Pz,w = 3,ts=ts$sample,lambda =bp,avetheta = F,T = 300,ifparallel = T)
ternaryplot(D3,col=data.matrix(classify_rows(Pz)),pch=data.matrix(classify_rows(Pz))+10,cex=1)
ternaryplot(D3,col=data.matrix(classify_rows(Pz[,1:2])),pch=data.matrix(classify_rows(Pz[,1:2]))+13,cex=1)
ternaryplot(D3,col=Pz[,1],pch=Pz[,1]+13,cex=1)

D2=wwCRPParallel(Pz,w = 2,ts=ts$sample,lambda =bp,avetheta = F,T = 300,ifparallel = T)
ternaryplot(D2,col=data.matrix(classify_rows(Pz[,1:2])),pch=data.matrix(classify_rows(Pz[,1:2]))+13,cex=1)
ternaryplot(D2,col=data.matrix(classify_rows(Pz[,c(1,3)])),pch=data.matrix(classify_rows(Pz[,c(1,3)]))+13,cex=1)
ternaryplot(D2,col=data.matrix(classify_rows(Pz[,2:3])),pch=data.matrix(classify_rows(Pz[,2:3]))+13,cex=1)

D1=wwCRPParallel(Pz,w =1,ts=ts$sample,lambda =bp,avetheta = F,T = 300,ifparallel = T)
for(i in 1:3)
{
  ternaryplot(D1,col=Pz[,i],pch=Pz[,i]+15,cex=1)
}

write.csv(cbind(D3,Pz),"F:\\Marginally coupled mixture desgins\\ExampleFFFD\\w3.csv",row.names=FALSE)

#MaxPro design
rand_design_part1=apply(matrix(rep(seq(from=0,to=1,length=nrow(Pz)),p-1),ncol=p-1),2,sample)
InitialDesign=cbind(rand_design_part1,Pz)
obj=MaxProQQ(InitialDesign, p_nom=3,nstarts = 10) 
MPD=ITM(obj$Design[,1:(p-1)])
ternaryplot(MPD,col=data.matrix(classify_rows(Pz)),pch=data.matrix(classify_rows(Pz))+13,cex=1)
ternaryplot(MPD,col=data.matrix(classify_rows(Pz[,1:2])),pch=data.matrix(classify_rows(Pz[,1:2]))+13,cex=1)
ternaryplot(MPD,col=Pz[,1],pch=Pz[,1]+13,cex=1)

write.csv(cbind(MPD,Pz),"F:\\Marginally coupled mixture desgins\\ExampleFFFD\\MPD.csv",row.names=FALSE)

# random number
RN=ITM(matrix(runif(nrow(Pz)*(p-1)),ncol=2))
ternaryplot(RN,col=data.matrix(classify_rows(Pz)),pch=data.matrix(classify_rows(Pz))+13,cex=1)
ternaryplot(RN,col=data.matrix(classify_rows(Pz[,1:2])),pch=data.matrix(classify_rows(Pz[,1:2]))+13,cex=1)
ternaryplot(RN,col=Pz[,1],pch=Pz[,1]+13,cex=1)

do.call(rbind,apply(matrix(1:8,ncol=1),1,function(i){unicri(D3[classify_rows(Pz)==i,],region = ts$sample)}))
do.call(rbind,apply(matrix(1:8,ncol=1),1,function(i){unicri(MPD[classify_rows(Pz)==i,],region = ts$sample)}))
do.call(rbind,apply(matrix(1:8,ncol=1),1,function(i){unicri(RN[classify_rows(Pz)==i,],region = ts$sample)}))

rbind(unicri(D3,region = ts$sample),
unicri(MPD,region = ts$sample),
unicri(RN,region = ts$sample))