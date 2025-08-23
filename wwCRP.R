# $w$-way coupled representative points based on Hybrid energy distance criterion
# $$E^{\boldsymbol{\lambda}^z}\left(F,F_{\mathcal{P}_n}\right)$$

# Pz: Design of $q$ qualitative variables, q\leq 2, otherwise use the function for Sliced representative points directly 
# ts: training sample y_1,...y_N
## w=1 +++++ Marginally coupled representative points (default)
## w=2 ----- Doubly coupled representative points
## w=q ***** Sliced representative points
# lambda: balanced parameters
# T: maximum number of iterations
# ifparallel: TRUE means to parallel update points in MM iterations

# setwd("/correct/directory/") # If your source file reports an error, please use setwd()
source("clcv.R") 
source("findindex.R") 

wwCRPParallel=function(Pz,ts,w=1,lambda=1/2,avetheta=TRUE,T=200,ifparallel=TRUE)
{ # Confirm the number of slices and the number of points in each slice
  n=nrow(Pz)   # the total number of representative points 
  q=ncol(Pz)   # the number of qualitative variables
  if(q==1){stop("please use slicedHEDParallel fucntion 
                When there is only one qualitative variable")}
  # Among each of the w variables selection, rows that have same level combination as row i 
  Indexlist <- findindex(Pz,w) 
  lcv <- clcv(Pz,w)
  if(avetheta){
    theta <- 1/choose(q,w)           #K^w
  }else{theta<-lcv/sum(lcv)}
  # Initial points
  N <- nrow(ts) # the number of points in training sample
  p <- ncol(ts) # the number of quantitative variables
  Pn <- ts[sample(1:N,n,replace = FALSE),]
  
  if(ifparallel==TRUE){
    if (!require("foreach")) install.packages("foreach")
    if (!require("doParallel")) install.packages("doParallel")
    cl <- makeCluster(detectCores()-1) # Adjust the number of threads according to the actual situation
    registerDoParallel(cl)
    cat(paste0("Parallel computing has been enabled, using ", getDoParWorkers(), " worker threads\n"))
  }
  trts=t(ts)  #The transposition of training samples is helpful for subsequent calculations
  
  #MM algorithm iteration
  Pnnew <- Pn
  for(t in 1:T){
    if(ifparallel==TRUE){
      Pnnew_list <- foreach(i = 1:n, .combine = rbind) %dopar% {
        xdy <- Pn[i,]-trts      # Result are saved by row
        normxy <- sqrt(colSums(xdy*xdy))
        normxy[normxy==0] <- 1  # Handle zero norms
        xdx <- t(Pn[i,]-t(Pn))
        normxx <- sqrt(rowSums(xdx*xdx))
        normxx[normxx==0] <- 1  # Handle zero norms
        xdx=xdx/normxx  # Normalization
        
        # within subdesigns
        group_indices <- split(Indexlist[[i]][, 1], Indexlist[[i]][, 2])
        submean <- t(sapply(group_indices, function(rows) {
          colMeans(xdx[rows, , drop = FALSE])
        }))
        
        
        # Update point
        update <- 1/mean(1/normxy)*(colMeans(ts/normxy)+(1-lambda)*colSums(theta*submean)+mean(lambda)*colMeans(xdx))
       
      }
      Pnnew <- as.matrix(Pnnew_list)  
    } else {
      
      for(i in 1:n){
        xdy <- Pn[i,]-trts      # Result are saved by row
        normxy <- sqrt(colSums(xdy*xdy))
        normxy[normxy==0] <- 1  # Handle zero norms
        xdx <- t(Pn[i,]-t(Pn))
        normxx <- sqrt(rowSums(xdx*xdx))
        normxx[normxx==0] <- 1  # Handle zero norms
        xdx=xdx/normxx  # Normalization
        # within subdesigns
        group_indices <- split(Indexlist[[i]][, 1], Indexlist[[i]][, 2])
        submean <- t(sapply(group_indices, function(rows) {
          colMeans(xdx[rows, , drop = FALSE])
        }))
        # Update point
        Pnnew[i,] <- 1/mean(1/normxy)*(colMeans(ts/normxy)+colMeans((1-lambda)*submean)+mean(lambda)*colMeans(xdx))          
      }
    }
    Pn <- Pnnew 
    cat(sprintf("Iteration %d/%d complete\n", t, T))
  }
  
  # Close parallel cluster
  if(ifparallel==TRUE){
    stopCluster(cl)
    registerDoSEQ()  # Restore to serial computing
  }
  
  return(Pn=Pn)
}


# Example
library(randtoolbox)
library(DoE.base)
# 3 factors 2 levels  
base_design <- fac.design(
  nfactors = 3,
  nlevels = c(2,2,2),
  replications = 15,
  randomize = FALSE
)
Pz0=data.matrix(base_design)
Pz = Pz0[,-4]

D3=wwCRPParallel(Pz,w = 3,ts=sobol(100000,2,scrambling = 3),lambda =0.8,avetheta = F,T = 300,ifparallel = T)
plot(D3,col=data.matrix(classify_rows(Pz)),pch=data.matrix(classify_rows(Pz))+10,cex=2,xlim=c(0,1),ylim=c(0,1))

D2=wwCRPParallel(Pz,w = 2,ts=sobol(100000,2,scrambling = 3),lambda =0.5,avetheta = F,T = 300,ifparallel = T)
plot(D2,col=data.matrix(classify_rows(Pz[,1:2])),pch=data.matrix(classify_rows(Pz[,1:2]))+13,cex=2,xlim=c(0,1),ylim=c(0,1))
plot(D2,col=data.matrix(classify_rows(Pz[,c(1,3)])),pch=data.matrix(classify_rows(Pz[,c(1,3)]))+13,cex=2,xlim=c(0,1),ylim=c(0,1))
plot(D2,col=data.matrix(classify_rows(Pz[,2:3])),pch=data.matrix(classify_rows(Pz[,2:3]))+13,cex=2,xlim=c(0,1),ylim=c(0,1))

D1=wwCRPParallel(Pz,w =1,ts=sobol(100000,2,scrambling = 3),lambda =0.5,avetheta = F,T = 300,ifparallel = T)
for(i in 1:3)
{
  plot(D1,col=Pz[,i],pch=Pz[,i]+15,cex=2,xlim=c(0,1),ylim=c(0,1))
}

for (i in 1:2) {
  plot(D[Pz[,i]==1,],pch=16,cex=2,xlim=c(0,1),ylim=c(0,1))
  points(D[Pz[,i]==2,],pch=17,col=2,cex=2,xlim=c(0,1),ylim=c(0,1))
}
for (i in 3) {
  plot(D[Pz[,i]==1,],pch=16,cex=2,xlim=c(0,1),ylim=c(0,1))
  points(D[Pz[,i]==2,],pch=17,col=2,cex=2,xlim=c(0,1),ylim=c(0,1))
  points(D[Pz[,i]==3,],pch=15,col=3,cex=2,xlim=c(0,1),ylim=c(0,1))
  points(D[Pz[,i]==4,],pch=18,col=4,cex=2,xlim=c(0,1),ylim=c(0,1))
  
  }
















