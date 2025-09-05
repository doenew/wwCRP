library(foreach)
library(doParallel)
library(DiceKriging)  # Make sure to load the DiceKriging package
library(DoE.base)

### y_3
model_basic <- function(x) {
  # 1. Base term (no z interaction)
  base <- 2.87*x[1] + 1.08*x[2] + 2.01*x[3] - 1.14*x[1]*x[2] - 1.00*x[1]*x[3] + 0.20*x[2]*x[3] + 3.18*x[1]*x[2]*x[3]
  
  # 2. Single z interaction terms (x[4], x[5], x[6] interacting with x respectively)
  z1term <- (0.49*x[1] + 0.18*x[2] + 0.25*x[3] - 0.81*x[1]*x[2] - 0.54*x[1]*x[3] - 0.14*x[2]*x[3] + 0.07*x[1]*x[2]*x[3]) * x[4]
  z2term <- (0.71*x[1] + 0.25*x[2] + 0.40*x[3] - 0.59*x[1]*x[2] - 0.01*x[1]*x[3] + 0.07*x[2]*x[3] - 1.41*x[1]*x[2]*x[3]) * x[5]
  z3term <- (-0.09*x[1] - 0.08*x[2] + 0.01*x[3] + 0.10*x[1]*x[2] - 0.03*x[1]*x[3] - 0.19*x[2]*x[3] + 0.11*x[1]*x[2]*x[3]) * x[6]
  
  # 3. Double z interaction terms (x[4]x[5], x[4]x[6], x[5]x[6] interacting with x respectively)
  z1z2term <- (0.07*x[1] - 0.03*x[2] + 0.00*x[3] - 0.06*x[1]*x[2] + 0.14*x[1]*x[3] + 0.23*x[2]*x[3] + 1.74*x[1]*x[2]*x[3]) * x[4]*x[5]
  z1z3term <- (-0.05*x[1] - 0.05*x[2] + 0.17*x[3] + 0.14*x[1]*x[2] - 0.27*x[1]*x[3] - 0.25*x[2]*x[3] - 0.71*x[1]*x[2]*x[3]) * x[4]*x[6]
  z2z3term <- (0.10*x[1] - 0.03*x[2] - 0.05*x[3] - 0.19*x[1]*x[2] - 0.43*x[1]*x[3] + 0.12*x[2]*x[3] + 1.77*x[1]*x[2]*x[3]) * x[5]*x[6]
  
  # 4. Triple z interaction term (x[4]x[5]x[6] interacting with x)
  z1z2z3term <- (0.04*x[1] - 0.04*x[2] - 0.04*x[3] - 0.09*x[1]*x[2] - 0.12*x[1]*x[3] + 0.27*x[2]*x[3] - 1.33*x[1]*x[2]*x[3]) * x[4]*x[5]*x[6]
  
  
  base + z1term + z2term + z3term + z1z2term + z1z3term + z2z3term + z1z2z3term
}

# # y_2
# model_basic <- function(x) {
#   # 1. Base term (no z interaction)
#   base <- 2.87*x[1] + 1.08*x[2] + 2.01*x[3] - 1.14*x[1]*x[2] - 1.00*x[1]*x[3] + 0.20*x[2]*x[3] + 3.18*x[1]*x[2]*x[3]
#   
#   # 2. Single z interaction terms (x[4], x[5], x[6] interacting with x respectively)
#   z1term <- (0.49*x[1] + 0.18*x[2] + 0.25*x[3] - 0.81*x[1]*x[2] - 0.54*x[1]*x[3] - 0.14*x[2]*x[3] + 0.07*x[1]*x[2]*x[3]) * x[4]
#   z2term <- (0.71*x[1] + 0.25*x[2] + 0.40*x[3] - 0.59*x[1]*x[2] - 0.01*x[1]*x[3] + 0.07*x[2]*x[3] - 1.41*x[1]*x[2]*x[3]) * x[5]
#   z3term <- (-0.09*x[1] - 0.08*x[2] + 0.01*x[3] + 0.10*x[1]*x[2] - 0.03*x[1]*x[3] - 0.19*x[2]*x[3] + 0.11*x[1]*x[2]*x[3]) * x[6]
#   
#   # 3. Double z interaction terms (x[4]x[5], x[4]x[6], x[5]x[6] interacting with x respectively)
#   z1z2term <- (0.07*x[1] - 0.03*x[2] + 0.00*x[3] - 0.06*x[1]*x[2] + 0.14*x[1]*x[3] + 0.23*x[2]*x[3] + 1.74*x[1]*x[2]*x[3]) * x[4]*x[5]
#   z1z3term <- (-0.05*x[1] - 0.05*x[2] + 0.17*x[3] + 0.14*x[1]*x[2] - 0.27*x[1]*x[3] - 0.25*x[2]*x[3] - 0.71*x[1]*x[2]*x[3]) * x[4]*x[6]
#   z2z3term <- (0.10*x[1] - 0.03*x[2] - 0.05*x[3] - 0.19*x[1]*x[2] - 0.43*x[1]*x[3] + 0.12*x[2]*x[3] + 1.77*x[1]*x[2]*x[3]) * x[5]*x[6]
#   
#   # 4. Total (error term ε needs to be passed externally or assumed 0 for testing)
#   base + z1term + z2term +  z3term + z1z2term + z1z3term + z2z3term
# }
#  # y_1
# model_basic <- function(x) {
#   # Base term
#   base <- 2.87*x[1] + 1.08*x[2] + 2.01*x[3] - 1.14*x[1]*x[2] - 1.00*x[1]*x[3] + 0.20*x[2]*x[3] + 3.18*x[1]*x[2]*x[3]
#   # x[4] interaction term
#   z1term <- (0.49*x[1] + 0.18*x[2] + 0.25*x[3] - 0.81*x[1]*x[2] - 0.54*x[1]*x[3] - 0.14*x[2]*x[3] + 0.07*x[1]*x[2]*x[3]) * x[4]
#   # x[5] interaction term
#   z2term <- (0.71*x[1] + 0.25*x[2] + 0.40*x[3] - 0.59*x[1]*x[2] - 0.01*x[1]*x[3] + 0.07*x[2]*x[3] - 1.41*x[1]*x[2]*x[3]) * x[5]
#   # x[6] interaction term
#   z3term <- (-0.09*x[1] - 0.08*x[2] + 0.01*x[3] + 0.10*x[1]*x[2] - 0.03*x[1]*x[3] - 0.19*x[2]*x[3] + 0.11*x[1]*x[2]*x[3]) * x[6]
#   # Total
#   base + z1term + z2term + z3term
# }
# 
# # y_0
# model_basic <- function(x) {
#   2.87*x[1] + 1.08*x[2] + 2.01*x[3] - 1.14*x[1]*x[2] - 1.00*x[1]*x[3] + 0.20*x[2]*x[3] + 3.18*x[1]*x[2]*x[3]
# }


# 3 factors 2 levels  
base_design <- fac.design(
  nfactors = 3,
  nlevels = c(2,2,2),
  replications =1,
  randomize = FALSE
)
Pz0=data.matrix(base_design)
Pz = Pz0[,-4]*2-3   # Change the levels to -1 and 1

# Set base path
base_path <- "J:/Marginally coupled mixture desgins/ex2small/"

# Set up parallel computing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Export necessary variables and functions to all nodes
clusterExport(cl, c("base_path", "Pz", "Gibbs_l", "model_basic", "km", "predict.km"))

# Parallel computation
results <- foreach(kk = 1:50, .combine = 'rbind', .packages = c("DiceKriging")) %dopar% {
  # Generate test data
  testdata0 <- Gibbs_l(rep(0,3), rep(1,3), n = 18000)
  testdata <- cbind(testdata0$sample[10001:18000,], Pz[rep(1:8, times=1000),])
  Ytestdata <- apply(testdata, 1, model_basic)
  colnames(testdata) <- c("V1", "V2", "V3", "A", "B", "C")
  
  #################################################
  timesn=2 # 2,4,8,12 for n=16,32,64,96############
  #################################################
  
  # Read and process data
  fullD105 <- cbind(read.csv(paste0(base_path, "1wCDM0.5", kk, ".csv")), Pz[rep(1:8, times=timesn),])
  fullD205 <- cbind(read.csv(paste0(base_path, "2wCDM0.5", kk, ".csv")), Pz[rep(1:8, times=timesn),])
  fullD305 <- cbind(read.csv(paste0(base_path, "3wCDM0.8", kk, ".csv")), Pz[rep(1:8, times=timesn),])
  fullD108 <- cbind(read.csv(paste0(base_path, "1wCDM0.8", kk, ".csv")), Pz[rep(1:8, times=timesn),])
  fullD208 <- cbind(read.csv(paste0(base_path, "2wCDM0.8", kk, ".csv")), Pz[rep(1:8, times=timesn),])
  fullD308 <- cbind(read.csv(paste0(base_path, "3wCDM0.8", kk, ".csv")), Pz[rep(1:8, times=timesn),])
  fullMPD  <- cbind(read.csv(paste0(base_path, "MPD", kk, ".csv")), Pz[rep(1:8, times=timesn),])
  
  FFFD0 <- read.csv(paste0(base_path, "FFFD", kk, ".csv"))
  num_mat <- matrix(NA, nrow = nrow(FFFD0[,4:6]), ncol = ncol(FFFD0[,4:6]))
  num_mat[FFFD0[,4:6] == "L1"] <- -1
  num_mat[FFFD0[,4:6] == "L2"] <- 1
  fullFFFD <- cbind(FFFD0[,1:3], num_mat)
  colnames(fullFFFD) <- c("V1", "V2", "V3", "A", "B", "C")

  # Calculate MSPE
  calculate_mspe <- function(design_data, test_data, y_test) {
    y_design <- apply(design_data, 1, model_basic)
    design_model <- km(design = design_data, response = y_design)
    predictions <- predict.km(design_model, newdata = test_data, type = "SK")$mean
    mean((predictions - y_test)^2)
  }
  
  # Calculate MSPE for each design
#  fullD0MSPE <- calculate_mspe(fullD0, testdata, Ytestdata)
  fullD105MSPE <- calculate_mspe(fullD105, testdata, Ytestdata)
  fullD205MSPE <- calculate_mspe(fullD205, testdata, Ytestdata)
  fullD305MSPE <- calculate_mspe(fullD305, testdata, Ytestdata)
  fullD108MSPE <- calculate_mspe(fullD108, testdata, Ytestdata)
  fullD208MSPE <- calculate_mspe(fullD208, testdata, Ytestdata)
  fullD308MSPE <- calculate_mspe(fullD308, testdata, Ytestdata)
  fullMPDMSPE  <- calculate_mspe(fullMPD, testdata, Ytestdata)
  fullFFFDMSPE <- calculate_mspe(fullFFFD, testdata, Ytestdata)
  
  # Return results
  c(fullD105MSPE, fullD205MSPE, fullD305MSPE, fullD108MSPE, 
    fullD208MSPE, fullD308MSPE, fullMPDMSPE, fullFFFDMSPE)
}

# Stop parallel cluster
stopCluster(cl)

# Extract results


allfullD105MSPE <- results[,1]
allfullD205MSPE <- results[, 2]
allfullD305MSPE <- results[,3]
allfullD108MSPE <- results[,4]
allfullD208MSPE <- results[, 5]
allfullD308MSPE <- results[, 6]
allfullMPDMSPE  <- results[, 7]
allfullFFFDMSPE <- results[, 8]

# Combine results and draw boxplot
allfullMSPE <- cbind(
  
  allfullD105MSPE,
  allfullD205MSPE,
  allfullD305MSPE,
  allfullD108MSPE,
  allfullD208MSPE,
  allfullD308MSPE,
  allfullMPDMSPE,
   allfullFFFDMSPE
)

boxplot(allfullMSPE)
colMeans(allfullMSPE)
write.csv(allfullMSPE,file =paste0("J:/Marginally coupled mixture desgins/ex2/small3.csv"),row.names = FALSE  )
