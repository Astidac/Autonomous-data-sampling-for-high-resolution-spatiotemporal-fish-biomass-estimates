


#____________________________________________________________________________________________________________
#
#     Statistical analyses for:
#
#    "Autonomous data sampling for high-resolution spatiotemporal fish biomass estimates"
#     Astrid A. Carlsen , Michele Casini, Francesco Masnadi, Olof Olsson, Aron Hejdst??m, Jonas Hentati-Sundberg 
#     Accepted for publication 08.10.2024, Ecological Informatics
#____________________________________________________________________________________________________________
### Obs, if I in any code reffer to "dati2" this is the working name for xdat.full, and can be exchanged for this!


### Cross validation functions:
# GAM (Model 1 and 2)
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


#CV function for Models 3-9: SdmTMB models ### OBS distribution family must be changed for Model 7 from AR1 to IID 
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
    bound.outer = 5
    max.edge=0.95
    
    mesh2 = inla.mesh.2d(loc=cbind(xdat.model$X, xdat.model$Y),
                         max.edge = c(1,20)*max.edge,
                         cutoff =3, 
                         offset = c(max.edge, bound.outer))
    
    
    xmesh <- make_mesh(xdat.model, c("X", "Y"), mesh = mesh2) #mesh relates to dataset and must be updated betwen model 7 and the others
    
    mx <- sdmTMB(
      formula = kformula, 
      data = xdat.model,
      mesh = xmesh,
      family = student(link = "identity", df=4),
      spatiotemporal = "AR1", ### IID for model 7
      spatial = "on",
      time = "Month")
    
    xpreds=predict(mx, newdata=testData, nsim=99)
    testData$pred=(apply(xpreds, 1, median))
    xres=cor(testData$log_NASC,testData$pred)^2
    kres=c(kres, xres)
    
  }
  return(mean(kres))
}

# set Working directory
wd <- setwd("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Papers\\Drones and Prey distribution\\Prediction2024\\")

# get data
xdat.full <- readRDS("xdat.full.RDS") # full resolution dataset
xdat.21_23 <- readRDS("xdat.21_23.RDS") #xdat.full for 2021-2023
xdat.agg <- readRDS("xdat.agg.RDS") #2x2km aggregated data /month

#library
library(mgcv)
library(sdmTMB)
library(sdmTMBextra)
library(ggplot2)
library(mgcViz)
library(visreg)
library(INLA)

#working data
#For models 1-4 and 8-9.
xdat.model=xdat.full[,c('X','Y','log_NASC',"Week", "Hour",'Depth','Year','Month', "Mean_thetao_SST", "Chl_surf", "Mean_NS_vo_Wclm", "Mean_EW_uo_Wclm","SO_surf","Mean_current_speed", "JDay")] 


#if necessary
#xdat.model$Year <- as.factor(xdat.model$Year)
# xdat.model$Month <- as.numeric(xdat.model$Month)
# xdat.model$Month <- xdat.model$Month+3

##____________________________Modelling:___________________________________
# For each of the models, run through the code below to get the statistical values and figures of interest. 
# Most of the code refers to the sdmTMB ('GLMM') models, reflecting the content of the paper


### Model 1: Formula for baseline GAM without evironmental variables
xformula6= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth) + s(Year, bs="re") + s(Month, bs="re")

### Model 2: Formula for GAM with evironmental variables
xformula6= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm  + Mean_current_speed + s(Year, bs="re") + s(Month, bs="re") 

## Model
mgam1 <- gam(formula = xformula6,   data=xdat.model,  family=gaussian())
##Temporal autocorrelation
acf_residuals <- acf(residuals(mgam1))

#residuals
resids <- residuals(mgam1) 
qqnorm(resids, main="Model 1");qqline(resids)

# Cross validation 
gam.store=NULL
for(x in 1){
  gam.store=c(gam.store,mean(cval(xdat = xdat.model, k=4, formula=xformula6)) )  
}

#saveRDS(gam.store,"gam.store_wEnv.RDS")
#visreg
visreg(mgam1)

# Model 3: Baseline GLMM without evironmental variables
kformula= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth) +  (1|Year)

# Model 4: GLMM with evironmental variables
kformula= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

### Model 5 and 6 are separated below

# Model 7:
kformula= m_log_NASC ~  s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

# Model 8:
kformula= log_NASC ~ + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)

# Model 9: 
kformula= log_NASC ~  s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)




###____For sdmTMB's____________________________

## Cross validation
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

# saveRDS(tmb.store,"tmb.store_noEnv.RDS")
mean(tmb.store_fullfinalfishmod) #get mean CV value 

# INLA mesh
bound.outer = 5
max.edge=0.95
mesh2 = inla.mesh.2d(loc=cbind(xdat.model$X, xdat.model$Y),
                     max.edge = c(1,20)*max.edge,
                     cutoff =10, 
                     offset = c(max.edge, bound.outer))

xmesh <- make_mesh(xdat.model, c("X", "Y"), mesh = mesh2)

# Visual inspection
#plot(mesh2)


## Run model for summary 
mx <- sdmTMB(
  formula = kformula,
  data = xdat.model,
  mesh = xmesh,
  family = student(link = "identity", df=4),
  spatiotemporal = "AR1", 
  spatial = "on",
  time = "Month")

## Test sanity and get summary of model
sanity(mx)
summary(mx)
tidy(mx, conf.int = TRUE)
log_likelihood <- logLik(mx) # Calculate LL
print(log_likelihood) 

## Predict on same data to evaluate fit
predicted <- predict(mx, newdata = xdat.model)
observed <- xdat.model$m_log_NASC
residuals <- observed - predicted$est
rmse <- sqrt(mean(residuals^2)) # Calculate RMSE
print(rmse)

absolute_residuals <- abs(observed - predicted$est) # Calculate absolute residuals
mae <- mean(absolute_residuals)
print(mae)

resids <- residuals(mx) 
qqnorm(resids, main="Model 9");qqline(resids)
acf_residuals <- acf(resids) #Testing temporal autocorrelation

# saveRDS(resids, "resids_Model3NoEnv.RDS")
# saveRDS(mx, "mx_fullcomptoagg_smooth.RDS")

#Visreg figure
visreg::visreg(mx, scale = "response", ylim = c(1.8, 5), add = TRUE) 

# Additional variable responses 
# visreg::visreg(mx, xvar='Week', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Hour', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='SB_CTTemp', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='SO_surf', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Chl_comb', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Depth_comb', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Mean_NS_vo_Wclm', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Mean_EW_uo_Wclm', by='Year', scale='response', overlay = TRUE)
# visreg::visreg(mx, xvar='Mean_current_speed', by='Year', scale='response', overlay = TRUE)



#_________________Final prediction on area_________________________________
# readRDS("df.gridNew.RDS")
# newpred <- dplyr::select(df.gridNew,"X","Y", "Mean_NS_vo_Wclm","Year","Month","Mean_EW_uo_Wclm","SB_CTTemp","Depth_comb","Chl_comb","SO_surf","Mean_current_speed","Hour", "lat", "lon")
# colnames(newpred) <- c("X","Y","Mean_NS_vo_Wclm","Year","Month","Mean_EW_uo_Wclm","Mean_thetao_SST","Depth","Chl_surf","SO_surf","Mean_current_speed","Hour", "lat", "lon")
# head(newpred)
# saveRDS(newpred,"newpred.RDS")

# newpred <- readRDS("newpred.RDS")
# summary(polygon_sum)



##_________Final prediction________________________________
## get dataframe with environmental variables 
temp <- readRDS("polygon_sum_Meshsized.RDS")
temp$Hour <- 1 #choose hour for prediction. As time is in CEST we chose 01, representing the darkest hour
newpred <- dplyr::select(temp,"X","Y", "NS_vo_Wclm","Year","Month","EW_uo_Wclm","thetao_SST","max_depth_","Chl_surf","SO_surf","current_sp","Hour", "lat", "lon")
colnames(newpred) <- c("X","Y","Mean_NS_vo_Wclm","Year","Month","Mean_EW_uo_Wclm","Mean_thetao_SST","Depth","Chl_surf","SO_surf","Mean_current_speed","Hour", "lat", "lon")

#saveRDS(newpred,"newpred2.RDS")

### check compatibility of dataset variables
### If there is any deviation between classes (e.g. if xdat.model$Month is a factor while newpred$Month is numeric, R will crashupon trying to predict)
str(xdat.model)
str(newpred)



#______________________________________________________
#
##        Prediction on spatial grid
#
#______________________________________________________

# # # # # # # # # # # # # # # # # # #
p <- predict(mx, newdata = newpred) # full Grid pred
p2 <- predict(mx, newdata = xdat.model) #  full Self pred

p3 <- predict(mx, newdata = xdat.Agg) #  Agg on Self pred
p4 <- predict(mx, newdata = newpred) #  Agg on grid pred
# # # # # # # # # # # # # # # # # # #

# saveRDS(p,"p_full_noweekhour_grid_smooth.RDS")
# saveRDS(p2,"p__noweekhour_OnSelf_full_smooth.RDS")
# saveRDS(p3,"p_agg_OnSelf_full_smooth.RDS")
# saveRDS(p4,"p_agg_grid_smooth.RDS")

# Get summaries of estimates
summary(p$est)
summary(p2$est)
summary(p3$est)
summary(p4$est)

# Observation summaries
summary(xdat.full$log_NASC)
summary(Polygon_sum$m_log_NASC)





### Compare estimates from the two resolutions (model 7 nd 8) predicted on grid:
##Make sure variable classes are compatible
#p$Month <- as.factor(p$Month)
#p4$Month <- as.factor(p4$Month)

library(dplyr)
# join the two predictions for comparison
PComp <- left_join(
  p, p4,
  by = c("X", "Y", "Year", "Month")
) %>%
  select(X, Y, Year, Month, est.x, est.y)  # Select both est.x (from p) and est.y (from p4)

# Calculte discrepancy
PComp$est_discrep <- PComp$est.y-PComp$est.x
hist(PComp$est_discrep, breaks=100)

#Correlate estimates model 
cor(PComp$est.y,PComp$est.x) #[1] 0.8027277



##plot estimate fits against each other for correlation
ggplot(data = PComp, aes(x = est.x, y = est.y)) +
  geom_point(color = "black", shape=1) +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE, aes(color = "Model fit"), linewidth=1) +
  geom_abline(intercept = 0, slope = 1, color = "red", aes(color = "1-1 Diagonal"),linewidth=1) +
  scale_color_manual(values = c("blue" = "blue", "red" = "red"),
                     name = "Lines",
                     labels = c("Model fit", "1-1 Diagonal")) +
  
  ylim(-2.5,10)+
  xlim(-1,8)+
  labs(x = "Estimate Model 8", y = "Estimate Model 7", title = "Correltion model 7 and 8") +
  theme_minimal()



###______Random effects plotted on map for MS___________________
#Rename/rerun with models with individual names
# I'm sorry, I cannot remember why I did it in the oposite order but surely there must have been a reason?
p1$resids <- residuals(mod9)  #Model 9
p2$resids <- residuals(mod8) # Model 8
p3$resids <- residuals(mod7) #Model 7


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

jpeg("Resids_Model7_9.jpeg",width = 500, height = 550, units = "mm", res = 1000)
print(est_resid1)
dev.off()


##____________prediction on same data to test distribution_____________
histp4 <- ggplot() +
  geom_histogram(data = p4, aes(x = est, fill = "Variable 1"), alpha = 0.5, bins = 30) +
  geom_histogram(data = xdat.full, aes(x = log_NASC, fill = "Variable 2"), alpha = 0.5, bins = 30) +
  labs(x = "Value", y = "Density", fill = "Variable") +
  scale_fill_manual(values = c("Variable 1" = "blue", "Variable 2" = "red"), 
                    labels = c("Variable 1" = "Estimate", "Variable 2" = "Observation")) +
  theme_minimal()+ 
  ggtitle("Model 7")

histp <- ggplot() +
  geom_histogram(data = p, aes(x = est, fill = "Variable 1"), alpha = 0.5, bins = 30) +
  geom_histogram(data = xdat.full, aes(x = log_NASC, fill = "Variable 2"), alpha = 0.5, bins = 30) +
  labs(x = "Value", y = "Density", fill = "Variable") +
  scale_fill_manual(values = c("Variable 1" = "blue", "Variable 2" = "red"), 
                    labels = c("Variable 1" = "Estimate", "Variable 2" = "Observation")) +
  ggtitle("Model 8")+
  theme_minimal()


## Plots
library(gridExtra)
grid.arrange(histp4, histp, nrow = 1)


#xdat.model$estnew <- pnew$est
xdat.model$pred_discr <- xdat.model$log_NASC-xdat.model$est #estnew
hist(xdat.model$pred_discr, breaks=100)

# df.grid$pred <-  p$est
# df.gridNew$pred <-  p$est


##___________Match and test against agg pred___________:
p2$Month <- as.factor(p2$Month)
pcomp <- left_join(p2,p3, by=c("Month", "Year", "X", "Y")) #Model 7 and 8 Back-predicted on model data 
pcomp$disc_est <- pcomp$est.x-pcomp$est.y

hist(pcomp$disc_est, breaks=120, main="Discrepancy estimates")
saveRDS(pcomp, "pcomp_fullVSagg_db.RDS")

## Compare Model 7 est with model 8 est
ggplot() +
  geom_histogram(data = pcomp, aes(x = est.x, fill = "Variable 1"), alpha = 0.5, bins = 30) +
  geom_histogram(data = pcomp, aes(x = est.y, fill = "Variable 2"), alpha = 0.5, bins = 30) +
  geom_histogram(data = xdat.full, aes(x = log_NASC, fill = "Variable 3"), alpha = 0.5, bins = 30) +
  labs(x = "log NASC", y = "Density", fill = "Variable") +
  scale_fill_manual(values = c("Variable 1" = "blue", "Variable 2" = "red"), 
                    labels = c("Variable 1" = "Estimate full res.", "Variable 2" = "Estimate aggregated")) +
  theme_minimal() #+
#facet_wrap(Month~Year)



##______________regression between Model 7 est with model 8 est________
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
  theme_minimal() #+facet_wrap(Month~Year)


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

##_________Prediction plotted on grid/field
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

### ____Model 7 predicted on field:____________________ 
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


##_____________________________________________________________
##          CMSI vs USV environment:      
##
###          Model 5 and Model 6
#______________________________________________________________

#rerun data
#You can either use xdat.full with this code or xdat.21_23 where these few rows already have been run
xdat.model=xdat.full[,c('X','Y','log_NASC', "Week", "Hour",'Depth','Year','Month', "SB_CTTemp", "SB_CTCond",  "SB_ECO3_CHL", "Mean_thetao_SST",  "Mean_NS_vo_Wclm", "Mean_EW_uo_Wclm","Mean_current_speed", "SO_surf", "Chl_surf")] 

xdat.model <- xdat.model[!xdat.model$Year %in% c(2019, 2020), ]
xdat.model <- xdat.model[!is.na(xdat.model$SB_ECO3_CHL), ]
xdat.model <- xdat.model[!is.na(xdat.model$SB_CTTemp), ]

xdat.model$Year <- as.factor(xdat.model$Year)
xdat.model$Year <- droplevels(xdat.model$Year)

#saveRDS(xdat.model,"xdat.21_23.RDS")


# Same procedure as above: choose one model to run, and then the following code
# Model 5 formula USV
formi=  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(SB_CTTemp) + s(SO_surf) + s(SB_ECO3_CHL) +  s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm + Mean_current_speed + (1|Year)

# Model 6 formula Copernicus
formi= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth)+ s(Mean_thetao_SST) + s(SO_surf) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm +  Mean_current_speed + (1|Year)


# Cross validation (using function on top of the page)
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


## Model summary:
##  
###____Mesh____________________________
bound.outer = 5
max.edge=0.95
mesh3 = inla.mesh.2d(loc=cbind(xdat.model$X, xdat.model$Y),
                     max.edge = c(1,20)*max.edge,
                     cutoff =3, 
                     offset = c(max.edge, bound.outer))

##Make sure all columns have the necessary classes
columns_to_convert <- c(
  "SB_CTTemp", "SB_CTCond", "SB_ECO3_CHL",
  "Mean_NS_vo_Wclm", "Mean_EW_uo_Wclm",  "Mean_current_speed", "SO_surf", "Chl_surf"
)

xdat.model[columns_to_convert] <- lapply(xdat.model[columns_to_convert], as.numeric)
tmesh <- make_mesh(xdat.model, c("X", "Y"), mesh = mesh3)


#Model (5/6)
modi <- sdmTMB(
  formula = formi,
  data = xdat.model,
  mesh = tmesh,
  family = student(link = "identity", df=4),
  spatiotemporal = "ar1",
  spatial = "on",
  time = "Month",
  control = sdmTMBcontrol(newton_loops = 1))

# sanity and summary
sanity(modi)
summary(modi)

resids <- residuals(modi) 
qqnorm(resids, main="Model 7");qqline(resids) ##update name for model...!

# ##Save USV
# saveRDS(modi,"modi_modelUSVenvir.RDS")
# saveRDS(resids,"resids_modelUSVenvir.RDS")
# 
# #Save CMSI
# saveRDS(modi,"modi_modelCMSIenvir.RDS")
# saveRDS(resids,"resids_modelCMSIenvir.RDS")


## overview outcomes from crossvalidations
jpeg("boxplot_R2_TMB_SB21env.jpeg",width = 150, height = 150, units = "mm", res = 600)
boxplot(tmb.store, tmb.storeNoenv,names=c("sdmTMB","sdmTMB SB Environment"))
dev.off()


# Student t-test difference between the two model performances
t.test(tmb.store, tmb.storeSB21) #for this, just rename the tmb-store of the USV (tmb.storeSB21; SB for sailbouy) before you run the CMSI data model 


# ###_____________________Model 10__________
# #Model 10: Aggregation
# tmb.store10 <- NULL
# for (x in 1) {
#   success <- FALSE
#   while (!success) {
#     tryCatch({
#       temp <- cval.TMB.main(xdat=xdat.Full, k=4, xformula = Model10)
#       tmb.store10= rbind(tmb.store10, temp)
#       success <- TRUE
#     }, error = function(e) {
#       cat("Error occurred. Retrying...\n")
#       Sys.sleep(1)  
#     })
#   }
# }
# 
# tmb.store10$Model <- "Model 10" 
# 

#####________________Fig 4 histogram: LOG NASC HIST FOR ECHOGRAM FIG __________

## ___quartiles:
summary_stats <- summary(xdat.full$log_NASC)
lower_extreme <- min(xdat.full$log_NASC)
first_quartile <- quantile(xdat.full$log_NASC, probs = 0.25)
mean_value <- mean(xdat.full$log_NASC)
third_quartile <- quantile(xdat.full$log_NASC, probs = 0.75)
upper_extreme <- max(xdat.full$log_NASC)

p <- ggplot(xdat.full, aes(x = log_NASC)) +
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



## ______________Figure 5:___________________________________
### Environmental variable Figures 
xdat.full1 <- xdat.full[xdat.full$Year=="2021",]


###____Figure 5a: Seasurface Temperature______________
gg <- ggplot(xdat.full1, aes(x = JDay)) +
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


###________________Figure 5b: Salinity year by year from sailor...___________________________________________
xdat.full1 <- xdat.full1[xdat.full1$SB_CTCond>1,]

gg <- ggplot(xdat.full1, aes(x = JDay)) +
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


##________Figure 5c: Chlorophyll________________

gg <- ggplot(xdat.full1, aes(x = JDay)) +
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


##_________________Discrepancy environmental variables_____________________
db2 <- xdat.full[!xdat.full$Year==c(2019, 2020),] #select out 2021-2023

#calculate discrepancies
db$SST_disc <- db$SB_CTTemp - db$Mean_thetao_SST
db$SS_disc <- db$SB_CTCond - db$SO_surf
db2$Chl_disc <- db2$SB_ECO3_CHL - db2$Chl_surf

##_________Fig 5d: Temperature
gg4 <- ggplot(db, aes(x = JDay, y = SST_disc, group = Year, color = factor(Year))) +
  #geom_smooth(aes(y = xdat.full$SB_CTTemp, linetype = "SB_CTTemp"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
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


#________Fig 5e: Salinity__________
gg5 <- ggplot(db, aes(x = JDay, y = SS_disc, group = Year, color = factor(Year))) +
  #geom_smooth(aes(y = xdat.full$SB_CTCond, linetype = "SB_CTCond"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
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


##_________Fig 5f: Chlorophyll_________________
gg6 <- ggplot(db2, aes(x = JDay, y = db2$Chl_disc, group = Year, color = factor(Year))) +
  #geom_smooth(aes(y = xdat.full$SB_CTTemp, linetype = "SB_CTTemp"), method = "loess", se = FALSE, color = "#4B8BBE", size = 1) +
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

# ggsave("Discrepancy_fig.png", ggall, width = 6, height = 6, units = "in")
dev.off()



### Comparison values from SMHI
##_________Figure A5.1: station BY38 KARLSVDJ__________
sharkweb_data <- read.delim("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Data\\Main_data\\Environment\\sharkweb_data.txt")
library(dplyr)
smhi <- sharkweb_data
smhi <- smhi[smhi$Sample.maximum.depth<7,]
smhi <- smhi %>%
  mutate(Value = gsub(",", ".", Value))

ggplot(smhi, aes(x = Month, y = Value, group = Year, color = factor(Year))) +
  geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Year)) +
  labs(title = "BY38 KARLSVDJ, SMHI Salinity (PSU)", x = "Month", y = "Salinity (PSU)") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("#4B8BBE", "#FFA07A", "#77AC30")) 

smhi2 <- smhi




## CALCULATION % percent of trip that was sampled in potential acoustic deadzone
#### How much of data is shallower than 100m at daytime (from 03-21Hours).
#Values can be varied to test different depth thresholds or time periods

filtered_data <- xdat.full[as.numeric(xdat.full$Hour) > 3 & as.numeric(xdat.full$Hour) < 21 & xdat.full$Depth > -80, ]
summary(filtered_data)
percentage <- nrow(filtered_data) / nrow(xdat.full) * 100
print(percentage)



##___________Get shapefiles for predictions:________________
# library(sf)
# shapefile <- st_read("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Papers\\Drones and Prey distribution\\Prediction2024\\polygon_sum_sh_ExportFeature2.shp") #old file: polygon_sum_sh_Smooth.shp
# 
# polygon_geometry <- st_geometry(shapefile)
# plot(polygon_geometry)
# 
# library(foreign)
# data <- read.dbf("C:\\Users\\adcn0002\\OneDrive - Sveriges lantbruksuniversitet\\Doktorand SLU\\Papers\\Drones and Prey distribution\\Prediction2024\\polygon_sum_cl_E_ExportTable.dbf")
# saveRDS(data,"polygon_sum_Meshsized.RDS")



# ###_______________________FIXED MAP FOR PREDICTIONS, RAW AND RESIDUALS_________________________:
# p$Month <- factor(p$Month, levels = c(4, 5, 6, 7), labels = c("April", "May", "June", "July"))
# 
# est_pl <- ggplot() +
#   # geom_sf(data = polygon_geometry) +  
#   geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
#   geom_point(data = p, aes(x = X * 1000, y = Y * 1000, col = exp(est))) +
#   scale_color_viridis_c(
#     trans = "sqrt",
#     na.value = "red", limits = c(0, quantile(exp(p$est), 0.995))
#   ) +
#   theme_light() +
#   labs(fill = "Predicted", col = "NASC (nmi^2 m^2)") +  
#   labs(x = "Longitude", y = "Latitude") +
#   facet_grid(Year ~ Month) +  
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "bottom")  
# 
# est_pl <- est_pl +
#   xlim(613939, 699950) +
#   ylim(6288500, 6413060)
# est_pl


##__________Depth map________________________________________
est_pl <- ggplot() +
  # geom_sf(data = polygon_geometry) +  
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = newpred, aes(x = X * 1000, y = Y * 1000, col = Depth), size=3) +
  theme_light() +
  labs(fill = "Depth", col = "Depth (m)") +  
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")  

est_pl <- est_pl +
  xlim(613939, 699950) +
  ylim(6288500, 6413060)
est_pl
#ggsave("Depth_by_points.png", est_pl, width = 3, height = 6, units = "in")


##____NS current_____________
est_pl <- ggplot() +
  # geom_sf(data = polygon_geometry) +  
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = newpred, aes(x = X * 1000, y = Y * 1000, col = Mean_NS_vo_Wclm), size=2.5) +
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
#ggsave("NS_Current.png", est_pl, width = 3, height = 6, units = "in")


##____EW current_____________
est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = newpred, aes(x = X * 1000, y = Y * 1000, col = Mean_EW_uo_Wclm), size=2.5) +
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
#ggsave("EW_Current.png", est_pl, width = 3, height = 6, units = "in")


##____Current speed_____________
est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = newpred, aes(x = X * 1000, y = Y * 1000, col = Mean_current_speed), size=2.5) +
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
#ggsave("Current_speed.png", est_pl, width = 3, height = 6, units = "in")


##____Temperature_____________
est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = newpred, aes(x = X * 1000, y = Y * 1000, col = Mean_thetao_SST), size=2.5) +
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
#jpeg("SST_grid.png.jpeg",width = 150, height = 150, units = "mm", res = 600)
est_pl
dev.off()

#____Salinity_____________
est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = newpred, aes(x = X * 1000, y = Y * 1000, col = SO_surf), size=2.5) +
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
#ggsave("So_grid.png", est_pl, width = 3, height = 6, units = "in")

##____Chlorophyll_____________
est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = newpred, aes(x = X * 1000, y = Y * 1000, col = Chl_surf), size=2.5) +
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
#ggsave("Chl_grid.png", est_pl, width = 3, height = 6, units = "in")


###_____________Figure Appendix A14:_____________________________
##_______________________Hotspot plot RAW____________________
quantile_value <- quantile(xdat.full$log_NASC, probs = 0.75, na.rm = TRUE) 
datip <- xdat.full[xdat.full$log_NASC > quantile_value, ]
print(datip)

est_pl <- ggplot() +
  geom_sf(data = bc_coast_proj, fill = "lightgrey") +  
  geom_point(data = xdat.full, aes(x = X * 1000, y = Y * 1000), color = "black") +  
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

# jpeg("observed_upper_Quant0.75.jpeg",width = 500, height = 320, units = "mm", res = 1000)
print(est_pl)
dev.off()

