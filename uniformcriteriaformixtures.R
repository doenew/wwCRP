#mean(mu), standard variance(sd), distance(RMSD), maximum distance (MaD),minimum between design(MiD)
print("unicri:For the uniformity evaluation of the design in the mixture experimental region, 
      the first parameter is the design to be evaluated, and the second parameter is a uniform 
      random sample with a sufficiently large sample size in the experimental region.")
unicri=function(design,region){
  mud=sqrt(sum((colMeans(design)-colMeans(region))^2))
  sdd=sqrt(sum((apply(design,2,sd)-apply(region,2,sd))^2))
  tdesign=t(design)
  drd=apply(region,1,function(x){
    min((colSums((x-tdesign)*(x-tdesign))))
  })
  RMSD=sqrt(mean(drd))
  MaD=max(sqrt(drd))
  MiD=min(dist(design))
  
  return(data.frame(mu=mud,sd=sdd,RMSD=RMSD,MaD=MaD,MiD=MiD))
}
