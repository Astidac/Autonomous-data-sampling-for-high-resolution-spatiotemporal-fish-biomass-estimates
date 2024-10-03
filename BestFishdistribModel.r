
### Main script for manusscript November 2023
cval=function(xdat, k, formula){
  
  kdat=xdat[sample(1:nrow(xdat)), ] 
  folds <- cut(seq(1,nrow(xdat)),breaks=k,labels=FALSE)
  kres=NULL
  
  for(i in 1:k){
    start=Sys.time()
    cat(i)
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- kdat[testIndexes, ]
    trainData <- kdat[-testIndexes, ]
    
    mx=gam(formula = formula,
           data=trainData,
           family=gaussian())
    testData$pred=predict(mx, newdata=testData)
    xres=cor(testData$log_NASC,testData$pred)^2
    kres=c(kres, xres)
    
  }
  return(kres)
}

cval.TMB=function(xdat, k, xformula){
  
  kdat=xdat[sample(1:nrow(xdat)), ] 
  folds <- cut(seq(1,nrow(xdat)),breaks=k,labels=FALSE)
  kres=NULL
  
  for(i in 1:k){
    start=Sys.time()
    cat(i)
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- kdat[testIndexes, ]
    trainData <- kdat[-testIndexes, ]
    #xmesh <- make_mesh(trainData, xy_cols = c("Xkm", "Ykm"), cutoff = xsize)
    #xmesh <- make_mesh(trainData, xy_cols = c("X", "Y"), cutoff = xsize) 

    bound.outer = 5
    #bound.outer = diff(range(dati2$X))/3
    max.edge=0.95
    
    mesh2 = inla.mesh.2d(loc=cbind(xdat.model$X, xdat.model$Y),
                         max.edge = c(1,20)*max.edge,
                         cutoff =3, 
                         offset = c(max.edge, bound.outer))
    
    
    xmesh <- make_mesh(xdat.model, c("X", "Y"), mesh = mesh2)
    
    mx <- sdmTMB(
      formula = kformula, 
      data = xdat.model,
      mesh = xmesh,
      family = student(link = "identity", df=4),
      spatiotemporal = "AR1", 
      spatial = "on",
      time = "Month")
    
    xpreds=predict(mx, newdata=testData, nsim=99)
    testData$pred=(apply(xpreds, 1, median))
    xres=cor(testData$log_NASC,testData$pred)^2
    kres=c(kres, xres)
    
  }
  return(mean(kres))
}

#dati2 <- readRDS("C:/Users/adcn0002/OneDrive - Sveriges lantbruksuniversitet/Doktorand SLU/Papers/Drones and Prey distribution/Prediction2024/CV_final/Prediction/dati2_rightDepth.RDS")

wd <- setwd("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Papers\\Drones and Prey distribution\\Prediction2024\\")

wd

library(mgcv)
library(sdmTMB)
library(sdmTMBextra)
library(ggplot2)
library(mgcViz)
library(visreg)
library(INLA)

xdat.model=dati2[,c('X','Y','log_NASC',"Week", "Hour",'Depth','Year','Month', "Mean_thetao_SST", "Chl_surf", "Mean_NS_vo_Wclm", "Mean_EW_uo_Wclm","SO_surf","Mean_current_speed", "JDay")] 

# xdat.model=dati2[,c('X','Y','log_NASC',"Hour",'Depth','Year','Month', "Mean_thetao_SST", "Chl_surf", "Mean_NS_vo_Wclm", "Mean_EW_uo_Wclm","SO_surf","Mean_current_speed")] 

#saveRDS(xdat.model,"xdat.Full_finalCV.RDS")
#saveRDS(dati2,"dati2_rightDepth.RDS")

#xdat.model=datix[,c('X','Y','log_NASC',"Week", "Hour",'Depth_comb','Year','Month', "SB_CTTemp", "Chl_comb", "Mean_NS_vo_Wclm", "Mean_EW_uo_Wclm","SO_surf","Mean_current_speed")] 
#kformula= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb)+ s(SB_CTTemp) + s(Chl_comb) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year)
#kformula= log_NASC ~  s(Week)  + s(Hour) + s(Depth_comb)+ s(SB_CTTemp) + s(Chl_comb) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf) + s(Mean_current_speed) + Year
# gam.check(mgam1)

xformula6= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth) + s(Year, bs="re") + s(Month, bs="re")
xformula6= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm  + Mean_current_speed + s(Year, bs="re") + s(Month, bs="re") 

mgam1 <- gam(formula = xformula6,   data=xdat.model,  family=gaussian())
acf_residuals <- acf(residuals(mgam1))

resids <- residuals(mgam1) 
qqnorm(resids, main="Model 1");qqline(resids)

gam.store=NULL
for(x in 1){
  gam.store=c(gam.store,mean(cval(xdat = xdat.model, k=4, formula=xformula6)) )  
}

saveRDS(gam.store,"gam.store_wEnv.RDS")

visreg(mgam1)

# kformula= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth) + Mean_NS_vo_Wclm + Mean_EW_uo_Wclm + Mean_current_speed + (1|Year)
kformula= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth) +  (1|Year)
kformula= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)
# kformula= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ Mean_thetao_SST + SO_surf + Chl_surf + Mean_NS_vo_Wclm + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

##____________Model for comparing Agg to non-agg model est:
kformula= m_log_NASC ~  s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

kformula= log_NASC ~  s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

# # Pred-agg test SAME
# kformula= log_NASC ~ s(Depth) + Mean_thetao_SST + SO_surf + Chl_surf +
#   s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm + Mean_current_speed +
#   (1 | Year)


###____Mesh for sdmTMB's____________________________
bound.outer = 5
#bound.outer = diff(range(dati2$X))/3
max.edge=0.95
mesh2 = inla.mesh.2d(loc=cbind(xdat.model$X, xdat.model$Y),
                     max.edge = c(1,20)*max.edge,
                     cutoff =10, 
                     offset = c(max.edge, bound.outer))

xmesh <- make_mesh(xdat.model, c("X", "Y"), mesh = mesh2)

#xdat.model$Year <- as.factor(xdat.model$Year)
# xdat.model$Month <- as.numeric(xdat.model$Month)
# xdat.model$Month <- xdat.model$Month+3

mx <- sdmTMB(
  formula = kformula,
  data = xdat.model,
  mesh = xmesh,
  family = student(link = "identity", df=4),
  spatiotemporal = "AR1", 
  spatial = "on",
  time = "Month")

sanity(mx)
summary(mx)
tidy(mx, conf.int = TRUE)

log_likelihood <- logLik(mx)
print(log_likelihood)

predicted <- predict(mx, newdata = xdat.model)
observed <- xdat.model$m_log_NASC
residuals <- observed - predicted$est
rmse <- sqrt(mean(residuals^2))
print(rmse)

absolute_residuals <- abs(observed - predicted$est)
mae <- mean(absolute_residuals)
print(mae)

resids <- residuals(mx) 
qqnorm(resids, main="Model 9");qqline(resids)
acf_residuals <- acf(resids)

# rq_res <- residuals(mx)
# rq_res <- rq_res[is.finite(rq_res)] 

saveRDS(resids, "resids_Model3NoEnv.RDS")
saveRDS(mx, "mx_fullcomptoagg_smooth.RDS")
wd

tmb.store <- NULL
for (x in 1:5) {
  success <- FALSE
  while (!success) {
    tryCatch({
      tmb.store <- c(tmb.store, cval.TMB(xdat = xdat.model, k = 8, xformula = kformula, xsize = 5))
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred. Retrying...\n")
      Sys.sleep(1)
    })
  }
}

saveRDS(tmb.store,"tmb.store_noEnv.RDS")
mean(tmb.store_fullfinalfishmod)

#visreg::visreg(mx, scale = "response", ylim = c(1.8, 4), add = TRUE) #if it doesnt work, remove ad=TRUE
visreg::visreg(mx, scale = "response", ylim = c(1.8, 5), add = TRUE) 

# visreg::visreg(mx, xvar='Week', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Hour', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='SB_CTTemp', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='SO_surf', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Chl_comb', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Depth_comb', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Mean_NS_vo_Wclm', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Mean_EW_uo_Wclm', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Mean_current_speed', by='Year', scale='response', overlay = TRUE)

plot(mesh2)

#_________________Final prediction on area_________________________________
# readRDS("df.gridNew.RDS")
# newpred <- dplyr::select(df.gridNew,"X","Y", "Mean_NS_vo_Wclm","Year","Month","Mean_EW_uo_Wclm","SB_CTTemp","Depth_comb","Chl_comb","SO_surf","Mean_current_speed","Hour", "lat", "lon")
# colnames(newpred) <- c("X","Y","Mean_NS_vo_Wclm","Year","Month","Mean_EW_uo_Wclm","Mean_thetao_SST","Depth","Chl_surf","SO_surf","Mean_current_speed","Hour", "lat", "lon")
# head(newpred)
# saveRDS(newpred,"newpred.RDS")

# newpred <- readRDS("newpred.RDS")
# summary(polygon_sum)

##_________Final prediction
temp <- readRDS("polygon_sum_Meshsized.RDS")
temp$Hour <- 1
newpred <- dplyr::select(temp,"X","Y", "NS_vo_Wclm","Year","Month","EW_uo_Wclm","thetao_SST","max_depth_","Chl_surf","SO_surf","current_sp","Hour", "lat", "lon")
colnames(newpred) <- c("X","Y","Mean_NS_vo_Wclm","Year","Month","Mean_EW_uo_Wclm","Mean_thetao_SST","Depth","Chl_surf","SO_surf","Mean_current_speed","Hour", "lat", "lon")
summary(newpred)
saveRDS(newpred,"newpred2.RDS")

newpred <- readRDS("newpred2.RDS")

#xdat.Agg$Month <- as.factor(xdat.Agg$Month)
newpred$Month <- as.factor(newpred$Month)
#newpred$Month <- as.numeric(newpred$Month+3)
newpred$Month <- as.integer(newpred$Month)
newpred$Hour <- as.integer(newpred$Hour)
#newpred$Month <- as.factor(newpred$Month)
newpred$Year <- as.factor(newpred$Year)

# # # # # # # # # # # # # # # # # # #
p <- predict(mx, newdata = newpred) # full Grid pred
p2 <- predict(mx, newdata = xdat.model) #  full Self pred

p3 <- predict(mx, newdata = xdat.Agg) #  Agg on Self pred
p4 <- predict(mx, newdata = newpred) #  Agg on grid pred
# # # # # # # # # # # # # # # # # # #

saveRDS(p,"p_full_noweekhour_grid_smooth.RDS")
saveRDS(p2,"p__noweekhour_OnSelf_full_smooth.RDS")

saveRDS(p3,"p_agg_OnSelf_full_smooth.RDS")
saveRDS(p4,"p_agg_grid_smooth.RDS")

summary(p$est)
summary(p2$est)
summary(p3$est)
summary(p4$est)
summary(dati2$log_NASC)
summary(Polygon_sum$m_log_NASC)


wd <- setwd("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Papers\\Drones and Prey distribution\\Prediction2024\\CV_final\\Prediction\\")


###______raw data plotted on map for MS___________________
head(p)
raw_pl <- ggplot(bc_coast_proj)  +
  geom_point(data = dati2, aes(x = X * 1000, y = Y * 1000, col = log_NASC)) +
  scale_colour_viridis_c()+
  geom_sf() +
  xlim(600000, 700000 ) +
  ylim(6280000,6430000 ) +
  theme_light() +
  labs(fill = "Raw") +
  labs(x = "Longitude", y = "Latitude")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(Year ~ Month)+ ggtitle("Raw data by month")

raw_pl


###______Random effects plotted on map for MS___________________
p1$resids <- residuals(mx) 
p2$resids <- residuals(mx3)
p3$resids <- residuals(mx2)


##Model 9
est_resid1 <- ggplot(bc_coast_proj)  +
  geom_point(data = p1, aes(x = X * 1000, y = Y * 1000, col = resids)) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = mean(p2$resids))+
  geom_sf() +
  xlim(600000, 700000 ) +
  ylim(6280000,6430000 ) +
  theme_light() +
  labs(fill = "Predicted") +
  labs(x = "Longitude", y = "Latitude")+ 
  facet_grid(~ Year)+ ggtitle("Model 9")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

est_resid1

saveRDS(est_resid1, "est_resid1.RDS")

#Model 8
est_resid2 <- ggplot(bc_coast_proj)  +
  geom_point(data = p2, aes(x = X * 1000, y = Y * 1000, col = resids)) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = mean(p2$resids))+
  geom_sf() +
  xlim(600000, 700000 ) +
  ylim(6280000,6430000 ) +
  theme_light() +
  labs(fill = "Predicted") +
  labs(x = "Longitude", y = "Latitude")+ 
  facet_grid(~ Year)+ ggtitle("Model 8")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
est_resid2

#Model 7
est_resid3 <- ggplot(bc_coast_proj)  +
  geom_point(data = p3, aes(x = X * 1000, y = Y * 1000, col = resids)) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = mean(p3$resids))+
  geom_sf() +
  xlim(600000, 700000 ) +
  ylim(6280000,6430000 ) +
  theme_light() +
  labs(fill = "Predicted") +
  labs(x = "Longitude", y = "Latitude")+ 
  facet_grid(~ Year)+ ggtitle("Model 7")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

est_resid3
est_resid <- grid.arrange(est_resid3,est_resid2,est_resid1,  nrow = 3)

#jpeg("Resids_Model9.jpeg",width = 500, height = 550, units = "mm", res = 1000)
jpeg("Resids_Model7_9.jpeg",width = 500, height = 550, units = "mm", res = 1000)
print(est_resid1)
dev.off()

##____________prediction on same data to test distribution_____________
# pnew <- predict(mx, newdata = xdat.model) 
# pnew <- p4 
head(pnew)

histp4 <- ggplot() +
  geom_histogram(data = p4, aes(x = est, fill = "Variable 1"), alpha = 0.5, bins = 30) +
  geom_histogram(data = dati2, aes(x = log_NASC, fill = "Variable 2"), alpha = 0.5, bins = 30) +
  labs(x = "Value", y = "Density", fill = "Variable") +
  scale_fill_manual(values = c("Variable 1" = "blue", "Variable 2" = "red"), 
                    labels = c("Variable 1" = "Estimate", "Variable 2" = "Observation")) +
  theme_minimal()+ 
  ggtitle("Model 7")

histp <- ggplot() +
  geom_histogram(data = p, aes(x = est, fill = "Variable 1"), alpha = 0.5, bins = 30) +
  geom_histogram(data = dati2, aes(x = log_NASC, fill = "Variable 2"), alpha = 0.5, bins = 30) +
  labs(x = "Value", y = "Density", fill = "Variable") +
  scale_fill_manual(values = c("Variable 1" = "blue", "Variable 2" = "red"), 
                    labels = c("Variable 1" = "Estimate", "Variable 2" = "Observation")) +
  ggtitle("Model 8")+
  theme_minimal()


library(gridExtra)
grid.arrange(histp4, histp, nrow = 1)

xdat.model$estnew <- pnew$est
xdat.model$pred_discr <- xdat.model$log_NASC-xdat.model$estnew
hist(xdat.model$pred_discr, breaks=100)

# df.grid$pred <-  p$est
# df.gridNew$pred <-  p$est


##___________Match and test against agg pred___________:
p2$Month <- as.factor(p2$Month)
pcomp <- left_join(p2,p3, by=c("Month", "Year", "X", "Y"))
pcomp$disc_est <- pcomp$est.x-pcomp$est.y

hist(pcomp$disc_est, breaks=120, main="Discrepancy estimates")
saveRDS(pcomp, "pcomp_fullVSagg_db.RDS")

 ## Compare agg est with non-agg est
ggplot() +
  geom_histogram(data = pcomp, aes(x = est.x, fill = "Variable 1"), alpha = 0.5, bins = 30) +
  geom_histogram(data = pcomp, aes(x = est.y, fill = "Variable 2"), alpha = 0.5, bins = 30) +
  geom_histogram(data = dati2, aes(x = log_NASC, fill = "Variable 3"), alpha = 0.5, bins = 30) +
  # geom_smooth(data = dati2, aes(x = log_NASC, y = ..density..), method = "lm", color = "black", size = 1, formula = y ~ x) +  # Add smoothed curve
   labs(x = "log NASC", y = "Density", fill = "Variable") +
  scale_fill_manual(values = c("Variable 1" = "blue", "Variable 2" = "red"), 
                    labels = c("Variable 1" = "Estimate full res.", "Variable 2" = "Estimate aggregated")) +
  theme_minimal() #+
#facet_wrap(Month~Year)


##______________     regression between obs and est________
ggplot(data = pcomp, aes(x = est.y, y = est.x)) +
  geom_point(color = "black", shape=1) +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE, aes(color = "Model fit"), linewidth=1) +
  geom_abline(intercept = 0, slope = 1, color = "red", aes(color = "1-1 Diagonal"),linewidth=1) +
  scale_color_manual(values = c("blue" = "blue", "red" = "red"),
                     name = "Lines",
                     labels = c("Model fit", "1-1 Diagonal")) +
  
  ylim(-2.5,10)+
  xlim(-1,8)+

  labs(x = "Full resolution estimate", y = "Aggregated estimate", title = "Prediction estimates") +
  theme_minimal()


##__________________ Model est versus observation_______________________________
ggplot() +
  geom_histogram(data = p, aes(x = est, fill = "Variable 1"), alpha = 0.5, bins = 30) +
  geom_histogram(data = xdat.model, aes(x = log_NASC, fill = "Variable 2"), alpha = 0.5, bins = 30) +
  labs(x = "log NASC", y = "Density", fill = "Variable") +
  scale_fill_manual(values = c("Variable 1" = "blue", "Variable 2" = "red"), 
                    labels = c("Variable 1" = "Estimate", "Variable 2" = "Observation")) +
  theme_minimal() #+
  #facet_wrap(Month~Year)


library(ggplot2)

ggplot(data = p2, aes(x = est, y = log_NASC)) +
  geom_point(color = "black", shape=1) +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE, aes(color = "Model fit"), linewidth=1) +
  geom_abline(intercept = 0, slope = 1, color = "red", aes(color = "1-1 Diagonal"),linewidth=1) +
  scale_color_manual(values = c("blue" = "blue", "red" = "red"),
                     name = "Lines",
                     labels = c("Model fit", "1-1 Diagonal")) +
  
  ylim(-2.5,10)+
  xlim(-1,8)+
  labs(x = "Estimated log NASC", y = "Observed log NASC", title = "Prediction Model 10") +
  theme_minimal()

##_________Prediction plotted on field
# max_lat <- 57.8#max(datix$Lat_M)
# min_lat <- 56.8#min(datix$Lat_M)
# max_lon <- 18.0#max(datix$Lon_M)
# min_lon <- 16.95#min(datix$Lon_M)
# 
# library(dplyr)
# 
# filtered_data <- p %>%
#   filter(lat >= min_lat & lat <= max_lat & lon >= min_lon & lon <= max_lon)
# 
# filtered_data <- filtered_data %>%
#  filter(Depth<(-20))

#plot(filtered_data$X, filtered_data$Y)
p4$Month <- factor(p4$Month, levels = c(4, 5, 6, 7), labels = c("April", "May", "June", "July"))

est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = p4, aes(x = X * 1000, y = Y * 1000, col = exp(est))) +
  scale_color_viridis_c(
    trans = "sqrt",
    na.value = "red", limits = c(0, quantile(exp(p4$est), 0.995))
  ) +
  theme_light() +
  labs(fill = "Predicted", col = "NASC (nmi^2 m^2)") +
  labs(x = "Longitude", y = "Latitude") +
  facet_grid(Year ~ Month) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")  

# Set plot limits
est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)

est_pl

# save final comparison plot
jpeg("Est_AllYears_Agg_smooth.jpeg",width = 500, height = 320, units = "mm", res = 1000)
print(est_pl)
dev.off()


####______________CMSI vs USV environment: models________________:
head(dati2)
xdat.model=dati2[,c('X','Y','log_NASC', "Week", "Hour",'Depth','Year','Month', "SB_CTTemp", "SB_CTCond",  "SB_ECO3_CHL", "Mean_thetao_SST",  "Mean_NS_vo_Wclm", "Mean_EW_uo_Wclm","Mean_current_speed", "SO_surf", "Chl_surf")] 
xdat.model <- xdat.model[!xdat.model$Year %in% c(2019, 2020), ]

xdat.model <- xdat.model[!is.na(xdat.model$SB_ECO3_CHL), ]
xdat.model <- xdat.model[!is.na(xdat.model$SB_CTTemp), ]
 
xdat.model$Year <- as.factor(xdat.model$Year)
xdat.model$Year <- droplevels(xdat.model$Year)

#saveRDS(xdat.model,"xdat.21-23_finalCV.RDS")

##Sailbuoy 
##Old mod:
# formi=  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb, bs='bs')+ s(SB_CTTemp, bs = 'bs') + s(SB_CTCond) + s(SB_ECO3_CHL) +  s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + Mean_current_speed + (1|Year) #s(SB_ECO3_PC)+ removed due to 5000+ NA's; s(SB_ECO3_Beta) 

#USV
formi=  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(SB_CTTemp) + s(SO_surf) + s(SB_ECO3_CHL) +  s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm + Mean_current_speed + (1|Year)

#Copernicus
formi= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

mgam1 <- gam(formula = formi,   data=xdat.model,  family=gaussian())
summary(mgam1)

# gam.store=NULL
# for(x in 1:2){
#   gam.store=c(gam.store,mean(cval(xdat = xdat.model, k=8, formula=xformula6)) )  
# }

###____Mesh____________________________
bound.outer = 5
max.edge=0.95
mesh3 = inla.mesh.2d(loc=cbind(xdat.model$X, xdat.model$Y),
                     max.edge = c(1,20)*max.edge,
                     cutoff =3, 
                     offset = c(max.edge, bound.outer))

columns_to_convert <- c(
  "SB_CTTemp", "SB_CTCond", "SB_ECO3_CHL",
  "Mean_NS_vo_Wclm", "Mean_EW_uo_Wclm",  "Mean_current_speed", "SO_surf", "Chl_surf"
)

xdat.model[columns_to_convert] <- lapply(xdat.model[columns_to_convert], as.numeric)
tmesh <- make_mesh(xdat.model, c("X", "Y"), mesh = mesh3)


#Model
modi <- sdmTMB(
  formula = formi,
  data = xdat.model,
  mesh = tmesh,
  family = student(link = "identity", df=4),
  spatiotemporal = "ar1",
  spatial = "on",
  time = "Month",
  control = sdmTMBcontrol(newton_loops = 1))


sanity(modi)
summary(modi)

resids <- residuals(modi) 
qqnorm(resids, main="Model 7");qqline(resids)

##Save USV
saveRDS(modi,"modi_modelUSVenvir.RDS")
saveRDS(resids,"resids_modelUSVenvir.RDS")

#Save CMSI
saveRDS(modi,"modi_modelCMSIenvir.RDS")
saveRDS(resids,"resids_modelCMSIenvir.RDS")


tmb.store <- NULL
for (x in 1) {
  success <- FALSE
  while (!success) {
    tryCatch({
      tmb.store <- c(tmb.store, cval.TMB(xdat = xdat.model, k = 4, xformula = kformula, xsize = 5))
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred. Retrying...\n")
      Sys.sleep(5)  
    })
  }
}
#saveRDS(tmb.storeSB21,"tmb.storeSB21.RDS")
#saveRDS(tmb.storeCP21,"tmb.storeCP21.RDS")


jpeg("boxplot_R2_TMB_SB21env.jpeg",width = 150, height = 150, units = "mm", res = 600)
boxplot(tmb.store, tmb.storeNoenv,names=c("sdmTMB","sdmTMB SB Environment"))
dev.off()

t.test(tmb.store, tmb.storeSB21)


###_____________________Predd model that match aggregated model__________
#Model 10: Aggregation
tmb.store10 <- NULL
for (x in 1) {
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


#####________________LOG NASC HIST FOR ECHOGRAM FIG __________
sub2020 <- dati2[dati2$Year==2021,]
sub20201 <- sub2020[sub2020$Month==5,]
plot(sub20201$Lon_M, sub20201$Lat_M)

write.csv(sub20201,"sub2021May.csv")

## ___ Distribution plot with Upper and lower extreme, etc.______
summary_stats <- summary(dati2$log_NASC)
lower_extreme <- min(dati2$log_NASC)
first_quartile <- quantile(dati2$log_NASC, probs = 0.25)
mean_value <- mean(dati2$log_NASC)
third_quartile <- quantile(dati2$log_NASC, probs = 0.75)
upper_extreme <- max(dati2$log_NASC)

p <- ggplot(dati2, aes(x = log_NASC)) +
  geom_histogram(binwidth = 0.1, fill = "white", color = "black") +
  xlab(print(expression("Log NASC"~(m^2~nmi^{-2})))) +
  ylab("Frequency")+
  ylim(0,1850)+
  theme_light() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

p <- p +
  geom_vline(xintercept = lower_extreme, color = "orange", linetype = "dashed") +
  geom_vline(xintercept = first_quartile, color = "orange", linetype = "dashed") +
  geom_vline(xintercept = mean_value, color = "orange", linetype = "dashed") +
  geom_vline(xintercept = third_quartile, color = "orange", linetype = "dashed") +
  geom_vline(xintercept = upper_extreme, color = "orange", linetype = "dashed") +
  
  annotate("text", x = lower_extreme, y = 1850, label = "Lower extreme", vjust = 0, hjust = 1, color = "orange", angle = 90) +
  annotate("text", x = first_quartile, y = 1850, label = "1. Quart", vjust = 0, hjust = 1, color = "orange", angle = 90) +
  annotate("text", x = mean_value, y = 1850, label = "Mean", vjust = 0, hjust = 1, color = "orange", angle = 90) +
  annotate("text", x = third_quartile, y = 1850, label = "3. Quart", vjust = 0, hjust = 1, color = "orange", angle = 90) +
  annotate("text", x = upper_extreme, y = 1800, label = "Upper extreme", vjust = 0, hjust = 1, color = "orange", angle = 90)

print(p)
dati21 <- dati2[dati2$Year=="2021",]

##________Chl________________
gg <- ggplot(dati21, aes(x = JDay)) +
  geom_point(aes(y = SB_ECO3_CHL), size = 1.2, type=3, color = "#4B8BBE") +
  geom_point(aes(y = Mean_Chl_Daily),size = 1.2, color = "#77AC30") +
  #geom_smooth(aes(y = SB_ECO3_CHL, linetype = "SB_ECO3_CHL"),  se = FALSE, color = "#4B8BBE", linewidth = 1) +
  #geom_smooth(aes(y = Mean_Chl_Daily, linetype = "Mean_Chl_Daily"), se = FALSE, color = "#77AC30", linewidth = 1) +
  labs( x = "Day of the year", y = "Chlorophyll (??g/l)") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_manual(name = "Chlorophyll Variables",
                     values = c("SB_ECO3_CHL" = "#4B8BBE", "Mean_Chl_Daily" = "#77AC30"),
                     labels = c("SB_ECO3_CHL" = "USV", "Mean_Chl_Daily" = "CMSI")) +
  scale_linetype_manual(name = "Chlorophyll Variables",
                        values = c("SB_ECO3_CHL" = "solid", "Mean_Chl_Daily" = "solid"),
                        labels = c("SB_ECO3_CHL" = "USV", "Mean_Chl_Daily" = "CMSI"))
gg
gg3 <- gg
ggsave("Discrep_SB_CMSI_ChlP.png", gg, width = 6, height = 3, units = "in")

### salinity_______________
dati21 <- dati21[dati21$SB_CTCond>1,]

gg <- ggplot(dati21, aes(x = JDay)) +
  geom_point(aes(y = SB_CTCond), size = 1, type=3, color = "#4B8BBE") +
  geom_point(aes(y = Mean_SS_Daily, col="red"), size = 1, type=3, color = "#77AC30") +
  # geom_line(aes(y = SB_CTCond, color = "SB_CTCond"), linewidth = 1) +
  # geom_line(aes(y = Mean_SS_Daily, color = "Mean_SS_Daily"), linewidth = 1) +
  #geom_smooth(aes(y = SB_CTCond, linetype = "SB_CTCond"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
  #geom_smooth(aes(y = Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
  labs( x = "Day of the year", y = "Salinity (mmhos/cm)") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_manual(name = "Salinity Variables",
                     values = c("SB_CTCond" = "#4B8BBE", "Mean_SS_Daily" = "#77AC30"),
                     labels = c("SB_CTCond" = "USV", "Mean_SS_Daily" = "CMSI")) +
  scale_linetype_manual(name = "Salinity Variables",
                        values = c("SB_CTCond" = "solid", "Mean_SS_Daily" = "solid"),
                        labels = c("SB_CTCond" = "USV", "Mean_SS_Daily" = "CMSI"))

gg
gg2 <- gg
ggsave("Discrep_SB_CMSI_SSP.png", gg, width = 6, height = 3, units = "in")

summary(dati21$SB_CTCond)

###____ST______________
gg <- ggplot(dati21, aes(x = JDay)) +
  geom_point(aes(y = SB_CTTemp), size = 1, type=3, color = "#4B8BBE") +
  geom_point(aes(y = Mean_ST_Daily, col="red" ), size = 1, type=3, color = "#77AC30") +
  #geom_smooth(aes(y = SB_CTTemp, linetype = "SB_CTTemp"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
  #geom_smooth(aes(y = Mean_ST_Daily, linetype = "Mean_ST_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
  labs( x = "Day of the year", y = "Temperature (\u00B0C)") +
  theme_minimal() +
  # theme(legend.position = "right") +
  scale_color_manual(name = "Temperature variables",
                     values = c("SB_CTTemp" = "#4B8BBE", "Mean_ST_Daily" = "#77AC30"),
                     labels = c("SB_CTTemp" = "USV", "Mean_ST_Daily" = "CMSI")) +
  scale_linetype_manual(name = "Temperature variables",
                        values = c("SB_CTTemp" = "solid", "Mean_ST_Daily" = "solid"),
                        labels = c("SB_CTTemp" = "USV", "Mean_ST_Daily" = "CMSI"))

gg
gg1 <- gg
ggsave("Discrep_SB_CMSI_STP1.png", gg, width = 6, height = 3, units = "in")


library(gridExtra)
gg4 <- grid.arrange(gg1, gg2, gg3, nrow = 3, heights = c(1, 1, 1))
#ggsave("Discrep_SB_CMSI_ALL.png", gg4, width = 6, height = 3, units = "in") 


###________________Salinity year by year from sailor...___________________________________________
dati3 <- dati2[dati2$Year==c(2021,2022,2023),] 
library(ggplot2)

gg <- ggplot(dati3, aes(x = JDay)) +
  geom_smooth(aes(y = SB_CTCond, linetype = "SB_CTCond"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
  geom_smooth(aes(y = Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
  labs(title = "Discrepancy salinity", x = "Day of the year", y = "Salinity (mmhos/cm)") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_manual(name = "Salinity Variables",
                     values = c("SB_CTCond" = "#4B8BBE", "Mean_SS_Daily" = "#77AC30"),
                     labels = c("SB_CTCond" = "USV", "Mean_SS_Daily" = "CMSI")) +
  scale_linetype_manual(name = "Salinity Variables",
                        values = c("SB_CTCond" = "solid", "Mean_SS_Daily" = "solid"),
                        labels = c("SB_CTCond" = "USV", "Mean_SS_Daily" = "CMSI"))+
  facet_wrap(~Year, scales = "free_y", ncol = 1) 

gg
ggsave("Sailor_SS_year.png", gg, width = 6, height = 3, units = "in")



##_________________Discrepancy environmental variables_____________________
db <- dati2
db2 <- db[!db$Year==c(2019, 2020),]

db$SST_disc <- db$SB_CTTemp - db$Mean_thetao_SST
db$SS_disc <- db$SB_CTCond - db$SO_surf
db2$Chl_disc <- db2$SB_ECO3_CHL - db2$Chl_surf

plot(db2$JDay,db2$Chl_disc)
plot(db$JDay,db$SS_disc)
plot(db$JDay,db$SST_disc)


#________Salinity__________
gg5 <- ggplot(db, aes(x = JDay, y = SS_disc, group = Year, color = factor(Year))) +
  #geom_smooth(aes(y = dati2$SB_CTCond, linetype = "SB_CTCond"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
  # geom_smooth(aes(y = dati3$Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
  geom_point(position = position_jitter(width = 0.2), size = 0.5, show.legend = FALSE) +
  #geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
  labs(title = "", x = "JDay", y = "Salinity (PSU)") +
  theme_minimal() +
  #ylim(5,13)+
  #theme(legend.position = "top") +
  scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30","#2AA07A", "#E7AC30"))

gg5
ggsave("Discr_SS_allyears.png", gg, width = 6, height = 3, units = "in")



##_________Temperature
gg4 <- ggplot(db, aes(x = JDay, y = SST_disc, group = Year, color = factor(Year))) +
  #geom_smooth(aes(y = dati2$SB_CTTemp, linetype = "SB_CTTemp"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
  # geom_smooth(aes(y = dati3$Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
  geom_point(position = position_jitter(width = 0.2), size = 0.5, show.legend = FALSE) +
  #geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
  labs(title = "", x = "JDay", y = "Temperature ('C)") +
  theme_minimal() +
  #ylim(0,13)+
  #theme(legend.position = "top") +
  scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30","#2AA07A", "#E7AC30"))  

gg4
ggsave("Discr_ST_allyears.png", gg, width = 6, height = 3, units = "in")


##_________Chlorophyll_________________
gg6 <- ggplot(db2, aes(x = JDay, y = db2$Chl_disc, group = Year, color = factor(Year))) +
  #geom_smooth(aes(y = dati2$SB_CTTemp, linetype = "SB_CTTemp"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
  # geom_smooth(aes(y = dati3$Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
  geom_point(position = position_jitter(width = 0.2), size = 0.5, show.legend = FALSE) +
  #geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
  labs(title = "", x = "JDay", y = "Chlorophyll-a") +
  theme_minimal() +
  #theme(legend.position = "top") +
  scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30","#2AA07A", "#E7AC30"))  

gg6
ggsave("Discr_Chl_allyears.png", gg, width = 6, height = 3, units = "in")


library(gridExtra)
left_column <- arrangeGrob(gg1, gg2, gg3, nrow = 3)
right_column <- arrangeGrob(gg4, gg5, gg6, nrow = 3)
ggall <- grid.arrange(left_column, right_column, nrow = 1)
print(ggall)

ggsave("Discrepancy_fig.png", ggall, width = 6, height = 6, units = "in")
dev.off()


##_________station BY38 KARLSVDJ__________
sharkweb_data <- read.delim("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Data\\Main_data\\Environment\\sharkweb_data.txt")
library(dplyr)
smhi <- sharkweb_data
smhi <- smhi[smhi$Sample.maximum.depth<7,]
smhi <- smhi %>%
  mutate(Value = gsub(",", ".", Value))

ggplot(smhi, aes(x = Month, y = Value, group = Year, color = factor(Year))) +
  # geom_smooth(aes(y = dati3$SB_CTCond, linetype = "SB_CTCond"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
  # geom_smooth(aes(y = dati3$Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
  #geom_point(position = position_jitter(width = 0.2), size = 2) +
  geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
  labs(title = "BY38 KARLSVDJ, SMHI Salinity (PSU)", x = "Month", y = "Salinity (PSU)") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30")) 

smhi2 <- smhi


# ##__________________________USV Envir variables all years______________
# #________Salinity__________
# dati2 <- dati2[!is.na(dati2$SB_CTCond),]
# gg <- ggplot(dati2, aes(x = JDay, y = SB_CTCond, group = Year, color = factor(Year))) +
#   #geom_smooth(aes(y = dati2$SB_CTCond, linetype = "SB_CTCond"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
#   # geom_smooth(aes(y = dati3$Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
#   geom_point(position = position_jitter(width = 0.2), size = 0.5) +
#   #geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
#   labs(title = "USV sea surface salinity", x = "JDay", y = "Salinity (PSU)") +
#   theme_minimal() +
#   ylim(5,13)+
#   #theme(legend.position = "top") +
#   scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30","#2AA07A", "#E7AC30"))  # Adjust colors as needed
# 
# gg
# ggsave("Sailor_SS_allyears.png", gg, width = 6, height = 3, units = "in")

# ##_________Temperature
# gg <- ggplot(dati2, aes(x = JDay, y = SB_CTTemp, group = Year, color = factor(Year))) +
#   #geom_smooth(aes(y = dati2$SB_CTTemp, linetype = "SB_CTTemp"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
#   # geom_smooth(aes(y = dati3$Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
#   geom_point(position = position_jitter(width = 0.2), size = 0.5) +
#   #geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
#   labs(title = "USV sea surface temperature", x = "JDay", y = "Temperature ('C)") +
#   theme_minimal() +
#   #ylim(0,13)+
#   #theme(legend.position = "top") +
#   scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30","#2AA07A", "#E7AC30"))  # Adjust colors as needed
# 
# gg
# ggsave("Sailor_ST_allyears.png", gg, width = 6, height = 3, units = "in")


# ##_________Chlorophyll_________________
# gg <- ggplot(dati2, aes(x = JDay, y = SB_ECO3_CHL, group = Year, color = factor(Year))) +
#   #geom_smooth(aes(y = dati2$SB_CTTemp, linetype = "SB_CTTemp"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
#   # geom_smooth(aes(y = dati3$Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
#   geom_point(position = position_jitter(width = 0.2), size = 0.5) +
#   geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
#   labs(title = "USV sea surface Chlorophyll", x = "JDay", y = "Chlorophyll-a") +
#   theme_minimal() +
#   #theme(legend.position = "top") +
#   scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30","#2AA07A", "#E7AC30"))  # Adjust colors as needed
# 
# gg
# ggsave("Sailor_Chl_allyears_line.png", gg, width = 6, height = 3, units = "in")


# ##__________________________CMSI Envir variables all years______________
# #________Salinity__________
# dati2 <- dati2[!is.na(dati2$SB_CTCond),]
# 
# gg <- ggplot(dati2, aes(x = JDay, y = SO_surf, group = Year, color = factor(Year))) +
#   #geom_smooth(aes(y = dati2$SB_CTCond, linetype = "SB_CTCond"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
#   # geom_smooth(aes(y = dati3$Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
#   geom_point(position = position_jitter(width = 0.2), size = 0.5) +
#   #geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
#   labs(title = "CMSI sea surface salinity", x = "JDay", y = "Salinity (PSU)") +
#   theme_minimal() +
#   ylim(5,13)+
#   #theme(legend.position = "top") +
#   scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30","#2AA07A", "#E7AC30"))  # Adjust colors as needed
# 
# gg
# ggsave("CMSI_SS_allyears.png", gg, width = 6, height = 3, units = "in")
# 
# 
# 
# ##_________Temperature
# gg <- ggplot(dati2, aes(x = JDay, y = Mean_thetao_SST, group = Year, color = factor(Year))) +
#   #geom_smooth(aes(y = dati2$SB_CTTemp, linetype = "SB_CTTemp"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
#   # geom_smooth(aes(y = dati3$Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
#   geom_point(position = position_jitter(width = 0.2), size = 0.5) +
#   #geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
#   labs(title = "CMSI sea surface temperature", x = "JDay", y = "Temperature ('C)") +
#   theme_minimal() +
#   #ylim(0,13)+
#   #theme(legend.position = "top") +
#   scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30","#2AA07A", "#E7AC30"))  # Adjust colors as needed
# 
# gg
# ggsave("CMSI_ST_allyears.png", gg, width = 6, height = 3, units = "in")
# 
# 
# ##_________Chlorophyll_________________
# gg <- ggplot(dati2, aes(x = JDay, y = Chl_surf, group = Year, color = factor(Year))) +
#   #geom_smooth(aes(y = dati2$SB_CTTemp, linetype = "SB_CTTemp"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
#   # geom_smooth(aes(y = dati3$Mean_SS_Daily, linetype = "Mean_SS_Daily"), method = "loess", se = FALSE, color = "#77AC30", size = 1) +
#   geom_point(position = position_jitter(width = 0.2), size = 0.5) +
#   #geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
#   labs(title = "CMSI sea surface Chlorophyll", x = "JDay", y = "Chlorophyll-a") +
#   theme_minimal() +
#   #theme(legend.position = "top") +
#   scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30","#2AA07A", "#E7AC30"))  # Adjust colors as needed
# 
# gg
# 
# ggsave("CMSI_Chl_allyears.png", gg, width = 6, height = 3, units = "in")


## CALCULATION % percent of trip was in potential acoustic deadzone
#### How much of data is shallower than 100m at daytime.
filtered_data <- dati2[as.numeric(dati2$Hour) > 3 & as.numeric(dati2$Hour) < 21 & dati2$Depth > -80, ]
summary(filtered_data)
percentage <- nrow(filtered_data) / nrow(dati2) * 100
print(percentage)


##___________New shapefile for predictions:________________
library(sf)
shapefile <- st_read("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Papers\\Drones and Prey distribution\\Prediction2024\\polygon_sum_sh_ExportFeature2.shp") #old file: polygon_sum_sh_Smooth.shp

polygon_geometry <- st_geometry(shapefile)
plot(polygon_geometry)

library(foreign)
data <- read.dbf("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Papers\\Drones and Prey distribution\\Prediction2024\\polygon_sum_cl_E_ExportTable.dbf")
saveRDS(data,"polygon_sum_Meshsized.RDS")

###_______________________FIXED MAP FOR PREDICTIONS, RAW AND RESIDUALS_________________________:
#p$Month <- factor(p$Month, levels = c(4, 5, 6, 7), labels = c("April", "May", "June", "July"))
p$Month <- factor(p$Month, levels = c(4, 5, 6, 7), labels = c("April", "May", "June", "July"))

est_pl <- ggplot() +
  # geom_sf(data = polygon_geometry) +  
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = p, aes(x = X * 1000, y = Y * 1000, col = exp(est))) +
  scale_color_viridis_c(
    trans = "sqrt",
    na.value = "red", limits = c(0, quantile(exp(p$est), 0.995))
  ) +
  theme_light() +
  labs(fill = "Predicted", col = "NASC (nmi^2 m^2)") +  
  labs(x = "Longitude", y = "Latitude") +
  facet_grid(Year ~ Month) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")  

est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)
est_pl


##__________Depth map________________________________________
est_pl <- ggplot() +
  # geom_sf(data = polygon_geometry) +  
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = p2, aes(x = X * 1000, y = Y * 1000, col = Depth), size=3) +
  theme_light() +
  labs(fill = "Depth", col = "Depth (m)") +  
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")  

est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)
est_pl
ggsave("Depth_by_points.png", est_pl, width = 3, height = 6, units = "in")


##____NS current_____________
est_pl <- ggplot() +
  # geom_sf(data = polygon_geometry) +  
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = p2, aes(x = X * 1000, y = Y * 1000, col = Mean_NS_vo_Wclm), size=2.5) +
  theme_light() +
  labs(fill = "North-South Current", col = "North-South Current (ms)") +  
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") + 
  facet_grid(Year ~ Month)   
  
est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)

est_pl
ggsave("NS_Current.png", est_pl, width = 3, height = 6, units = "in")


##____EW current_____________
est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = p2, aes(x = X * 1000, y = Y * 1000, col = Mean_EW_uo_Wclm), size=2.5) +
  theme_light() +
  labs(fill = "East-West Current", col = "East-West Current (ms)") +  
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") + 
  facet_grid(Year ~ Month)   

est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)

est_pl
ggsave("EW_Current.png", est_pl, width = 3, height = 6, units = "in")


##____Current speed_____________
est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = p2, aes(x = X * 1000, y = Y * 1000, col = Mean_current_speed), size=2.5) +
  theme_light() +
  labs(fill = "Current speed", col = "Current speed (ms)") +  
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") + 
  facet_grid(Year ~ Month)  

est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)

est_pl
ggsave("Current_speed.png", est_pl, width = 3, height = 6, units = "in")


##____Temperature_____________
est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = p2, aes(x = X * 1000, y = Y * 1000, col = Mean_thetao_SST), size=2.5) +
  theme_light() +
  labs(fill = "Seasurface temperature", col = "Mean_thetao_SST") +  
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") + 
  facet_grid(Year ~ Month)+ 
  scale_x_continuous(breaks = function(x) seq(from = floor(min(X * 1000)), to = ceiling(max(X * 1000)), by = 50000))

est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)

#ggsave("SST_grid.png", est_pl, width = 8, height = 6, units = "in")
jpeg("SST_grid.png.jpeg",width = 150, height = 150, units = "mm", res = 600)
est_pl
dev.off()

#____Salinity_____________
est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = p2, aes(x = X * 1000, y = Y * 1000, col = SO_surf), size=2.5) +
  theme_light() +
  labs(fill = "Seasurface Salinity", col = "Seasurface Salinity") +  
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") + 
  facet_grid(Year ~ Month)   

est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)

est_pl
ggsave("So_grid.png", est_pl, width = 3, height = 6, units = "in")

##____Chlorophyll_____________
est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = p2, aes(x = X * 1000, y = Y * 1000, col = Chl_surf), size=2.5) +
  theme_light() +
  labs(fill = "Seasurface Chlorophyll-a", col = "Seasurface Chlorophyll-a") +  
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") + 
  facet_grid(Year ~ Month)

est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)

est_pl
ggsave("Chl_grid.png", est_pl, width = 3, height = 6, units = "in")

##_______________________Hotspot plot RAW____________________
quantile_value <- quantile(dati2$log_NASC, probs = 0.75, na.rm = TRUE) 
datip <- dati2[dati2$log_NASC > quantile_value, ]
print(datip)

est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = dati2, aes(x = X * 1000, y = Y * 1000), color = "black") +  
  geom_point(data = datip, aes(x = X * 1000, y = Y * 1000), color = "red") +   
  scale_color_viridis_c() +
  theme_light() +
  labs(fill = "NASC", col = "NASC (nmi^2 m^2)") +
  labs(x = "Longitude", y = "Latitude") +
  facet_grid(Year ~ Month) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")  

est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)

jpeg("observed_upper_Quant0.75.jpeg",width = 500, height = 320, units = "mm", res = 1000)
print(est_pl)
dev.off()

###______________Get that mesh
ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = p, aes(x = X * 1000, y = Y * 1000, color = "darkblue"), alpha = 0.5)+
  scale_color_identity() +  
  theme_void()

dati22 <- dati2[dati2$Year=="2022",]
dati22 <- dati22[dati22$Month==5,]

write.csv(dati22,"dati22_arcgis.csv")
