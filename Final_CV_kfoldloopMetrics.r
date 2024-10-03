library(caret)
library(Metrics)
library(mgcv)
library(sdmTMB)
library(sdmTMBextra)
library(ggplot2)
library(INLA)
library(tidyr)
library(dplyr)

###############################################################
##                 MULTI-STATS CODE
###############################################################
wd <- setwd("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Papers\\Drones and Prey distribution\\Prediction2024\\CV_final")

###____Datasets:
xdat.Full <- readRDS("xdat.Full_finalCV.RDS")
xdat.21_23 <- readRDS("xdat.21-23_finalCV.RDS")
xdat.Agg <- readRDS("Polygon_sum_finalCV.RDS")

xdat.Full$Year <- as.factor(xdat.Full$Year)

write.csv(xdat.Full,"xdat.Full.csv")
write.csv(xdat.21_23,"xdat.21-23.csv")
write.csv(xdat.Agg,"xdat.Agg.csv")


str(xdat.Full)
str(xdat.21_23)
str(xdat.Agg)

colnames(xdat.21_23)

####################################################
# GAM cross validation + metrics
###################################################
cval.gam=function(xdat, k, xformula){
  
  kdat=xdat[sample(1:nrow(xdat)), ] 
  folds <- cut(seq(1,nrow(xdat)),breaks=k,labels=FALSE)
  #MAPEkres=NULL
  R2db=NULL
  RMSEdb=NULL
  MAEdb=NULL
  likdb=NULL
  
  for(i in 1:k){
    start=Sys.time()
    cat(i)
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- kdat[testIndexes, ]
    trainData <- kdat[-testIndexes, ]
    
    mx=gam(formula = xformula,
           data=trainData,
           family = gaussian(link = "identity"))
    testData$pred=predict(mx, newdata=testData)
    # MAPE <- mape(testData$pred , testData$log_NASC )
     xlik = as.numeric(logLik(mx))
    xR2 = R2(testData$pred,testData$log_NASC)
    xRMSE = RMSE(testData$pred,testData$log_NASC) 
    xMAE = MAE(testData$pred,testData$log_NASC) 
    #  MAPEkres=c(MAPEkres, MAPE )
     likdb = c(likdb, xlik )
    R2db = c(R2db, xR2 )
    RMSEdb = c(RMSEdb, xRMSE )
    MAEdb = c(MAEdb, xMAE )
  }

  result.db <- as.data.frame(list( likdb, R2db, MAEdb,RMSEdb))
  colnames(result.db) <- c("ML", "R2", "MAE", "RMSE")
  result.db$rowK <- 1:length(result.db$R2)
  result.db <- pivot_longer(result.db, cols= (1:4), values_to = "value" , names_to = "metrics")
  result.db <-result.db %>% group_by(metrics) %>% dplyr::summarise(meanv = mean(value))
 # result.db$Model <- "GAM"
  result.db$metrics <- factor(result.db$metrics, levels = c("ML", "R2", "MAE", "RMSE"))
  return(result.db)
}


####################################################
# MAIN sdmTMB cross validation + metrics 
###################################################
cval.TMB.main=function(xdat, k, xformula){
  
  kdat=xdat[sample(1:nrow(xdat)), ]
  folds <- cut(seq(1,nrow(xdat)),breaks=k,labels=FALSE)
  # MAPEkres=NULL
  R2db=NULL
  RMSEdb=NULL
  MAEdb=NULL
  likdb=NULL
  
  for(i in 1:k){
    start=Sys.time()
    cat(i)
    #Segement your data by fold using the which() function 
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- kdat[testIndexes, ]
    trainData <- kdat[-testIndexes, ]
    #xmesh <- make_mesh(trainData, xy_cols = c("Xkm", "Ykm"), cutoff = xsize)
    bound.outer = 5
    #bound.outer = diff(range(dati2$X))/3
    max.edge=0.95
    
    ##supermesh
    mesh2 = inla.mesh.2d(loc=cbind(trainData$X, trainData$Y),
                         max.edge = c(1,20)*max.edge,
                         cutoff =3, 
                         offset = c(max.edge, bound.outer))
    
    xmesh <- make_mesh(trainData, c("X", "Y"), mesh = mesh2)
    
    mx <- sdmTMB(
      formula = xformula, 
      data = trainData,
      mesh = xmesh,
      family = student(link = "identity", df=4),
      spatiotemporal = "AR1", 
      spatial = "on",
      time = "Month")
    #xpreds=predict(mx, newdata=testData)
    
    xpreds=predict(mx, newdata=testData, nsim=99)
    testData$pred=(apply(xpreds, 1, median))
    # MAPE <- mape(testData$pred , testData$log_NASC )
    xlik = as.numeric(mx[["model"]][["objective"]][1]) 
    xR2 = R2(testData$pred,testData$log_NASC)
    xRMSE = RMSE(testData$pred,testData$log_NASC)
    xMAE = MAE(testData$pred,testData$log_NASC)
    # MAPEkres=c(MAPEkres, MAPE )
    R2db = c(R2db, xR2 )
    RMSEdb = c(RMSEdb, xRMSE )
    MAEdb = c(MAEdb, xMAE )
    likdb = c(likdb, xlik )
  }
  #  return(mean(kres))
  result.db <- as.data.frame(list( likdb, R2db, MAEdb,RMSEdb))
  colnames(result.db) <- c( "ML", "R2", "MAE", "RMSE")
  result.db$rowK <- 1:length(result.db$R2)
  result.db <- pivot_longer(result.db, cols= (1:4), values_to = "value" , names_to = "metrics")
  result.db <-result.db %>% group_by(metrics) %>% dplyr::summarise(meanv = mean(value))
 # result.db$Model <- "sdmGLMM"
  result.db$metrics <- factor(result.db$metrics, levels = c("ML", "R2", "MAE", "RMSE"))
  return(result.db)
}



####################################################
# 2021-2023 sdmTMB cross validation + metrics 
###################################################
cval.TMB.21_23=function(xdat, k, xformula){
  
  kdat=xdat[sample(1:nrow(xdat)), ] 
  folds <- cut(seq(1,nrow(xdat)),breaks=k,labels=FALSE)
  # MAPEkres=NULL
  R2db=NULL
  RMSEdb=NULL
  MAEdb=NULL
  likdb=NULL
  
  for(i in 1:k){
    start=Sys.time()
    cat(i)
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- kdat[testIndexes, ]
    trainData <- kdat[-testIndexes, ]
    #xmesh <- make_mesh(trainData, xy_cols = c("Xkm", "Ykm"), cutoff = xsize)
    bound.outer = 5
    #bound.outer = diff(range(dati2$X))/3
    max.edge=0.95
    
    ##supermesh
    mesh2 = inla.mesh.2d(loc=cbind(trainData$X, trainData$Y),
                         max.edge = c(1,20)*max.edge,
                         cutoff =3, 
                         offset = c(max.edge, bound.outer))
    
    xmesh <- make_mesh(trainData, c("X", "Y"), mesh = mesh2)
    
    mx <- sdmTMB(
      formula = xformula, 
      data = trainData,
      mesh = xmesh,
      family = student(link = "identity", df=4),
      spatiotemporal = "AR1", 
      spatial = "on",
      time = "Month")
    #xpreds=predict(mx, newdata=testData)
    
    xpreds=predict(mx, newdata=testData, nsim=99)
    testData$pred=(apply(xpreds, 1, median))
    # MAPE <- mape(testData$pred , testData$log_NASC )
    xlik = as.numeric(mx[["model"]][["objective"]][1])
    xR2 = R2(testData$pred,testData$log_NASC)
    xRMSE = RMSE(testData$pred,testData$log_NASC)
    xMAE = MAE(testData$pred,testData$log_NASC)
    # MAPEkres=c(MAPEkres, MAPE )
    R2db = c(R2db, xR2 )
    RMSEdb = c(RMSEdb, xRMSE )
    MAEdb = c(MAEdb, xMAE )
    likdb = c(likdb, xlik )
  }
  #  return(mean(kres))
  result.db <- as.data.frame(list( likdb, R2db, MAEdb,RMSEdb))
  colnames(result.db) <- c( "ML", "R2", "MAE", "RMSE")
  result.db$rowK <- 1:length(result.db$R2)
  result.db <- pivot_longer(result.db, cols= (1:4), values_to = "value" , names_to = "metrics")
  result.db <-result.db %>% group_by(metrics) %>% dplyr::summarise(meanv = mean(value))
  # result.db$Model <- "sdmGLMM"
  result.db$metrics <- factor(result.db$metrics, levels = c("ML","R2", "MAE", "RMSE"))
  return(result.db)
}


####################################################
# Aggregated sdmTMB cross validation + metrics
###################################################
cval.TMB.agg=function(xdat, k, xformula){
  
  kdat=xdat[sample(1:nrow(xdat)), ] 
  folds <- cut(seq(1,nrow(xdat)),breaks=k,labels=FALSE)
  R2db=NULL
  RMSEdb=NULL
  MAEdb=NULL
  likdb=NULL
  
  for(i in 1:k){
    start=Sys.time()
    cat(i)
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- kdat[testIndexes, ]
    trainData <- kdat[-testIndexes, ]
    #xmesh <- make_mesh(trainData, xy_cols = c("Xkm", "Ykm"), cutoff = xsize)
    bound.outer = 3
    max.edge=0.95
    
    mesh2 = inla.mesh.2d(loc=cbind(trainData$X, trainData$Y),
                         max.edge = c(1,20)*max.edge,
                         cutoff =10, 
                         offset = c(max.edge, bound.outer))
    
    
    xmesh <- make_mesh(trainData, c("X", "Y"), mesh = mesh2)
    
    #_____sdmTMB model_______________________
    mx <- sdmTMB(
      formula = xformula, ###
      data = trainData,
      mesh = xmesh,
      family = student(link = "identity", df=4),
      spatiotemporal = "IID", 
      spatial = "on",
      time = "Month")

    
        
    xpreds=predict(mx, newdata=testData, nsim=99)
    testData$pred=(apply(xpreds, 1, median))
    # MAPE <- mape(testData$pred , testData$log_NASC )
    xlik = as.numeric(mx[["model"]][["objective"]][1])
    xR2 = R2(testData$pred,testData$m_log_NASC)
    xRMSE = RMSE(testData$pred,testData$m_log_NASC)
    xMAE = MAE(testData$pred,testData$m_log_NASC)
    # MAPEkres=c(MAPEkres, MAPE )
    R2db = c(R2db, xR2 )
    RMSEdb = c(RMSEdb, xRMSE )
    MAEdb = c(MAEdb, xMAE )
    likdb = c(likdb, xlik )
  }
  #  return(mean(kres))
  result.db <- as.data.frame(list( likdb, R2db, MAEdb,RMSEdb))
  colnames(result.db) <- c( "ML", "R2", "MAE", "RMSE")
  result.db$rowK <- 1:length(result.db$R2)
  result.db <- pivot_longer(result.db, cols= (1:4), values_to = "value" , names_to = "metrics")
  result.db <-result.db %>% group_by(metrics) %>% dplyr::summarise(meanv = mean(value))
  # result.db$Model <- "sdmGLMM"
  result.db$metrics <- factor(result.db$metrics, levels = c("ML","R2", "MAE", "RMSE"))
  return(result.db)
}



##___Model formulas___

# GAM "no" environment: 
Model1= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(Year, bs="re") + s(Month, bs="re")

# GAM full
Model2= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm  + Mean_current_speed + s(Year, bs="re") + s(Month, bs="re")

# sdmTMB "no" environment: 
Model3= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth) + (1|Year) 

## sdmTMB full 
Model4= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

## USV 2021-2023 sdmTMB full 
Model5= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(SB_CTTemp) + s(SO_surf) + s(SB_ECO3_CHL) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

## CMSI 2021-2023 sdmTMB full 
Model6= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

##____________Aggregation_______________ 
Model7= m_log_NASC ~  s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

##____________Model for comparing to aggregated 
Model8= log_NASC ~  s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

##____________Model for actual prediction on CMSI grid:
Model9= log_NASC ~  s(Hour,bs = 'cc') + s(Depth)+ Mean_thetao_SST + SO_surf + Chl_surf + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)



###___________Store_________________
#Model 1
gam.store1=NULL
for(x in 1:10){
  temp <- cval.gam(xdat = xdat.Full, k=8, xformula=Model1)
  gam.store1=rbind(gam.store1, temp)
}
gam.store1$Model <- "Model 1"

#Model 2
gam.store2=NULL
for(x in 1:10){
  temp <- cval.gam(xdat = xdat.Full, k=8, xformula=Model2)
  gam.store2=rbind(gam.store2, temp)
}
gam.store2$Model <- "Model 2"

#Model 3
tmb.store3 <- NULL
for (x in 1) {
  success <- FALSE
  while (!success) {
    tryCatch({
      temp <- cval.TMB.main(xdat=xdat.Full, k=8, xformula = Model3)
      tmb.store3= rbind(tmb.store3, temp)
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred. Retrying...\n")
      Sys.sleep(5)  
    })
  }
}

tmb.store3$Model <- "Model 3"

#Model 4
tmb.store4 <- NULL
for (x in 1) {
  success <- FALSE
  while (!success) {
    tryCatch({
      temp <- cval.TMB.main(xdat=xdat.Full, k=8, xformula = Model4)
      tmb.store4= rbind(tmb.store4, temp)
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred. Retrying...\n")
      Sys.sleep(5)  
    })
  }
}

tmb.store4$Model <- "Model 4"

#Model 5: USV 21-23
tmb.store5 <- NULL
for (x in 1) {
  success <- FALSE
  while (!success) {
    tryCatch({
      temp <- cval.TMB.21_23(xdat=xdat.21_23, k=8, xformula = Model5)
      tmb.store5= rbind(tmb.store5, temp)
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred. Retrying...\n")
      Sys.sleep(5) 
    })
  }
}

tmb.store5$Model <- "Model 5"

#Model 6: CMSI 21-23
tmb.store6 <- NULL
for (x in 1) {
  success <- FALSE
  while (!success) {
    tryCatch({
      temp <- cval.TMB.21_23(xdat=xdat.21_23, k=8, xformula = Model6)
      tmb.store6= rbind(tmb.store6, temp)
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred. Retrying...\n")
      Sys.sleep(5) 
    })
  }
}

tmb.store6$Model <- "Model 6"



#Model 7 (prev 8): Aggregation
tmb.store7 <- NULL
for (x in 1:5) {
  success <- FALSE
  while (!success) {
    tryCatch({
      temp <- cval.TMB.agg(xdat=xdat.Agg, k=4, xformula = Model7)
      tmb.store7= rbind(tmb.store7, temp)
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred. Retrying...\n")
      Sys.sleep(1) 
    })
  }
}

tmb.store7$Model <- "Model 7"

saveRDS(tmb.store7,"tmb_storeAggMod7.RDS")

# #Model 9: Pred_GRID
# tmb.store9 <- NULL
# for (x in 1:1) {
#   success <- FALSE
#   while (!success) {
#     tryCatch({
#       temp <- cval.TMB.main(xdat=xdat.Full, k=8, xformula = Model9)
#       tmb.store9= rbind(tmb.store9, temp)
#       success <- TRUE
#     }, error = function(e) {
#       cat("Error occurred. Retrying...\n")
#       Sys.sleep(5)  # Wait for 5 seconds before retrying
#     })
#   }
# }
# 
# tmb.store9$Model <- "Model 9"

#Model 10: Aggregation
tmb.store10 <- NULL
for (x in 1:5) {
  success <- FALSE
  while (!success) {
    tryCatch({
      temp <- cval.TMB.main(xdat=xdat.Full, k=4, xformula = Model10)
      tmb.store10= rbind(tmb.store10, temp)
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred. Retrying...\n")
      Sys.sleep(1)  
    })
  }
}

tmb.store10$Model <- "Model 10"



##############################################################################
##############################################################################
#### boxplot by metrics
#boxplt.db  <- rbind(gam.store1,gam.store2,tmb.store3,tmb.store4,tmb.store5,tmb.store6,tmb.store7,tmb.store8,tmb.store9, tmb.store10)

##Extract model 10, rename to model 8 and refit to boxplt.db :) Remove Model 7 and ad the new metrics for aggregated data as model 7.

##Rename boxplot 10 to 8
boxplt.db10 <- boxplt.db[boxplt.db$Model==c("Model 10"),]
boxplt.db10$Model <- "Model 8"

# Ad new aggdata model 7
boxplt.db <- rbind(boxplt.db, boxplt.db10)
boxplt.db <- boxplt.db[!boxplt.db$Model==c("Model 7"),]
boxplt.db  <- rbind(boxplt.db,tmb_store7)

# boxplt.db$LL <- boxplt.db$ML
boxplt.db_x <- boxplt.db[boxplt.db$metrics=="ML",]
boxplt.db_x$metrics <- "LL"
boxplt.db <- boxplt.db[!boxplt.db$metrics=="ML",]
boxplt.db  <- rbind(boxplt.db,boxplt.db_x)
tail(boxplt.db)

boxplt.db8 <- boxplt.db[boxplt.db$Model==c("Model 8") & boxplt.db$metrics=="LL",]
boxplt.db8$meanv <- -boxplt.db8$meanv
boxplt.dbLL <- boxplt.db[boxplt.db$metrics=="LL",]
boxplt.dbLL <- boxplt.dbLL[!boxplt.dbLL$Model==c("Model 8"),]
boxplt.dbLL <- rbind(boxplt.dbLL, boxplt.db8)
boxplt.db <- boxplt.db[!boxplt.db$metrics=="LL",]
boxplt.db  <- rbind(boxplt.db,boxplt.dbLL)



saveRDS(boxplt.db, "boxplt.db.RDS")

# boxplot.db <- readRDS("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Papers\\Drones and Prey distribution\\Prediction2024\\CV_final\\boxplt.db.10iter.RDS")

boxplt.db <- readRDS("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Papers\\Drones and Prey distribution\\Prediction2024\\CV_final\\boxplt.db.RDS")

library(dplyr)
mean.metrics <- boxplt.db %>% group_by(Model,metrics ) %>% dplyr::summarise(meanv = mean(meanv))
mean.metrics

library(readr)
# write_csv(boxplot.db, "boxplot_db10iter.csv")
write_csv(boxplot.db, "boxplot_db5iter.csv")
# write_csv(mean.metrics, "mean_metrics10iter.csv")
write_csv(mean.metrics, "mean_metrics5iter.csv")

library(ggplot2)
jpeg("boxplot_metrics5iter.jpeg",width = 150, height = 150, units = "mm", res = 600)
ggplot(boxplt.db) + geom_boxplot(aes(Model,meanv, color= Model)) + facet_wrap(~metrics, scales = "free") + theme_bw()+ ylab("mean value of k-fold")+   theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7))
dev.off()






