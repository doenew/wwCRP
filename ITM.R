ITM=function(A){
  q=ncol(A)+1
  n=nrow(A)
  s=1/((q-1):1)
  C1=t(t(A)^s)
  C2=cbind(1,C1)
  EMP=matrix(0,n,q)
  for(j in 2:q){
    EMP[,j]=apply(C2[,1:j],1,prod) }
  EMP[,1]=1
  AT=cbind((1-C1)*EMP[,1:(q-1)],EMP[,q])
  AT }