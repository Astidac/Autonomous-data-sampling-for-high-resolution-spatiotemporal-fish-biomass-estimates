

#____________________________________________________________________________________________________________
#
#     Sample of model selections for:
#
#    "Autonomous data sampling for high-resolution spatiotemporal fish biomass estimates"
#     Astrid A. Carlsen , Michele Casini, Francesco Masnadi, Olof Olsson, Aron Hejdst??m, Jonas Hentati-Sundberg 
#     Accepted for publication 08.10.2024, Ecological Informatics
#____________________________________________________________________________________________________________

library(mgcv)
library(sdmTMB)

##___________________Cross validation:__________________________________________
#___________Results: 2) Models for backscatter__________________________________
## Here is a representative selection of the variable based model simplification process 
## stepwise analyses of several rounds.

## Minimum, maximum and smoothed maximum compositions are compared below, but all possible combinations of smooths were tested before the final composition in the paper were deemed the best 


#______ Models including only Copernicus environment:
xdat.model=xdat.full[,c('X','Y','log_NASC',"Distgroup", "Week", "Hour",'Depth_comb','Year','Date','Month', "Mean_thetao_SST", "Chl_surf", "Mean_NS_vo_Wclm", "Mean_EW_uo_Wclm","SO_surf","Mean_current_speed","JDay", "SB_CTTemp", "Chl_comb")] 

xdat.model <- xdat.model[!is.na(xdat.model$SB_CTTemp),]

#With Environmental variables from Copernicus
kformula= log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb)+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year)

xdat.model$Year <- as.factor(xdat.model$Year)

#____Mesh:
bound.outer = 5
max.edge=0.95
mesh2 = inla.mesh.2d(loc=cbind(xdat.model$X, xdat.model$Y),
                     max.edge = c(1,20)*max.edge,
                     cutoff =3, 
                     offset = c(max.edge, bound.outer))

xmesh <- make_mesh(xdat.model, c("X", "Y"),mesh= mesh2)


#____Baseline model with intercept only
set.seed(1)

cv_intercept <- sdmTMB_cv(
  log_NASC ~ 1, 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1))


## No smooths except cyclic and random error
set.seed(1) 
cv_depth <- sdmTMB_cv(
  log_NASC ~  Week  + s(Hour,bs = 'cc') + Depth_comb+ Mean_thetao_SST + Chl_surf + Mean_NS_vo_Wclm + Mean_EW_uo_Wclm + SO_surf  + Mean_current_speed + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1))

##All smooths s(bs='tp')
set.seed(1) 
cv_depthX <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb)+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)

head(xdat.full)
## Best_comb


## Final model structure
set.seed(1) 
cv_depthX1 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb)+ s(SB_CTTemp) + s(Chl_comb) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm + s(SO_surf)  + Mean_current_speed + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)

# model selection was based on the sum_loglik though all likelyhoods extracted were investigated for stability across folds
cv_intercept$sum_loglik
cv_depth$sum_loglik 
cv_depthX$sum_loglik
cv_depthX1$sum_loglik # ML criterion at convergence: 59790.179



#__________________________Variables_________________________________________________________________________________
## Models were stepwise tested for variables to include
##No current speed
set.seed(1) 
cv_depth2 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb)+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + (1|Year),
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)

##No Salinity____
set.seed(1) 
cv_depth3 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb)+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)


#No EW currents
set.seed(1) 
cv_depth4 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb)+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) +  s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)

summary(xdat.model)
#No NS currents
set.seed(1) 
cv_depth5 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb)+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)

#No Chl surface
set.seed(1) 
cv_depth6 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb)+ s(Mean_thetao_SST)  + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year),
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)

##No SSTemperature
set.seed(1) 
cv_depth7 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)



## No Depth
set.seed(1) 
cv_depth8 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') +  s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)



## No Week
set.seed(1) 
cv_depth9 <- sdmTMB_cv(
  log_NASC ~   s(Hour,bs = 'cc') + s(Depth_comb)+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)


## No Hour
set.seed(1) 
cv_depth10 <- sdmTMB_cv(
  log_NASC ~  s(Week)  +  s(Depth_comb)+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)

# List for investigation, highest value (in this case closest to zero) wins
cv_intercept$sum_loglik
cv_depth$sum_loglik 
cv_depthX$sum_loglik 
cv_depth2$sum_loglik
cv_depth3$sum_loglik
cv_depth4$sum_loglik
cv_depth5$sum_loglik

cv_depth6$sum_loglik
cv_depth7$sum_loglik
cv_depth8$sum_loglik
cv_depth9$sum_loglik
cv_depth10$sum_loglik




# Summary of sum_loglik vectors
sum_loglik_list <- list(
  cv_intercept$sum_loglik,
  cv_depth$sum_loglik,
  cv_depth2$sum_loglik,
  cv_depth3$sum_loglik,
  cv_depth4$sum_loglik,
  cv_depth5$sum_loglik,
  cv_depth6$sum_loglik,
  cv_depth7$sum_loglik,
  cv_depth8$sum_loglik,
  cv_depth9$sum_loglik,
  cv_depth10$sum_loglik
)
max_index <- which.max(sapply(sum_loglik_list, max))
vector_with_highest <- sum_loglik_list[[max_index]]
print(vector_with_highest) #cv_depth= Full model, wins with -67706.96
# with close ties to cv_depth2 "no current speed". Clear margin to the next models.


##_____________________________CV for smoothers (CMSI based models 1-4, 6-9)____________________
### stepwise selection of smooths:
###  Obs: During modelling, fine-tuning of smooths with different kappa values and smooth types were tested, but the simplest form was deemed the most reasonable

#week linear
set.seed(1) 
cv_depth2 <- sdmTMB_cv(
  log_NASC ~  Week  + s(Hour,bs = 'cc') + s(Depth_comb)+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + (1|Year),
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)

##depth linear
set.seed(1) 
cv_depth3 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + Depth_comb+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + (1|Year),
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)


#Temp linear
set.seed(1) 
cv_depth4 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + s(Depth_comb)+ Mean_thetao_SST + s(Chl_surf) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)


#...... You get the hint. we ended with the best smooth composition being the ones found in MS Table 2:

set.seed(1) 
cv_depth10 <- sdmTMB_cv(
  log_NASC ~  s(Week)  + s(Hour,bs = 'cc') + poly(Depth_comb,2)+ s(Mean_thetao_SST) + s(Chl_surf) + s(Mean_NS_vo_Wclm) + Mean_EW_uo_Wclm + s(SO_surf)  + current_speed + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)



##____________Rerun the same analyses for model 5 with USV-sampled environmental variables
## Model 5:
##  here model 5, with USV-measured environmental variables


set.seed(1) 
cv_depth <- sdmTMB_cv(
  log_NASC ~  s(Hour,bs = 'cc') + s(Week) +s(Depth_comb)+ s(SB_CTTemp) + s(Chl_comb) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
  data = xdat.model,  
  family = student(link = "identity", df=4),
  mesh= xmesh,
  k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
)

cv_depth$sum_loglik

hist(xdat.full$log_NASC, breaks = 50, xlab="Log NASC", main="")





# 
# _____________Repeat of step-wise selection of variables for USV-data based model (Model 5)__________________________
# ##No current speed
# set.seed(1) 
# cv_depth2 <- sdmTMB_cv(
#   log_NASC ~  s(Hour,bs = 'cc') + s(Week) +s(Depth_comb)+ s(SB_CTTemp) + s(Chl_comb) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  +  (1|Year),
#   data = xdat.model,  
#   family = student(link = "identity", df=4),
#  mesh= xmesh,
#   k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
# )
# 
# ##No Salinity
# set.seed(1) 
# cv_depth3 <- sdmTMB_cv(
#   log_NASC ~  s(Hour,bs = 'cc') + s(Week) +s(Depth_comb)+ s(SB_CTTemp) + s(Chl_comb) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(Mean_current_speed) +(1|Year),
#   data = xdat.model,  
#   family = student(link = "identity", df=4),
#  mesh= xmesh,
#   k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
# )
# 
# #No EW currents
# set.seed(1) 
# cv_depth4 <- sdmTMB_cv(
#   log_NASC ~  s(Hour,bs = 'cc') + s(Week) +s(Depth_comb)+ s(SB_CTTemp) + s(Chl_comb) + s(Mean_NS_vo_Wclm) +  s(SO_surf)  + s(Mean_current_speed) + (1|Year),
#   data = xdat.model,  
#   family = student(link = "identity", df=4),
#  mesh= xmesh,
#   k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
# )
# 
# 
# #No NS currents
# set.seed(1) 
# cv_depth5 <- sdmTMB_cv(
#   log_NASC ~  s(Hour,bs = 'cc') + s(Week) +s(Depth_comb)+ s(SB_CTTemp) + s(Chl_comb) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year) 
#   data = xdat.model,  
#   family = student(link = "identity", df=4),
#  mesh= xmesh,
#   k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
# )
# 
# #No Chl surface
# set.seed(1) 
# cv_depth6 <- sdmTMB_cv(
#   log_NASC ~  s(Hour,bs = 'cc') + s(Week) +s(Depth_comb)+ s(SB_CTTemp) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year)
#   data = xdat.model,  
#   family = student(link = "identity", df=4),
#  mesh= xmesh,
#   k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
# )
# 
# ##No SSTemperature
# set.seed(1) 
# cv_depth7 <- sdmTMB_cv(
#   log_NASC ~  s(Hour,bs = 'cc') + s(Week) +s(Depth_comb) + s(Chl_comb) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
#   data = xdat.model,  
#   family = student(link = "identity", df=4),
#  mesh= xmesh,
#   k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
# )
# 
# 
# 
# ## No Depth
# set.seed(1) 
# cv_depth8 <- sdmTMB_cv(
#   log_NASC ~  s(Hour,bs = 'cc') + s(Week) + s(SB_CTTemp) + s(Chl_comb) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
#   data = xdat.model,  
#   family = student(link = "identity", df=4),
#   mesh= xmesh,
#   k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
# )
# 
# 
# 
# ## No Week
# set.seed(1) 
# cv_depth9 <- sdmTMB_cv(
#   log_NASC ~  s(Hour,bs = 'cc')  +s(Depth_comb)+  s(SB_CTTemp) + s(Chl_comb) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
#   data = xdat.model,  
#   family = student(link = "identity", df=4),
#   mesh= xmesh,
#   k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
# )
# 
# 
# ## No Hour
# set.seed(1) 
# cv_depth10 <- sdmTMB_cv(
#   log_NASC ~  s(Depth_comb)+  s(SB_CTTemp) + s(Chl_comb) + s(Mean_NS_vo_Wclm) + s(Mean_EW_uo_Wclm) + s(SO_surf)  + s(Mean_current_speed) + (1|Year), 
#   data = xdat.model,  
#   family = student(link = "identity", df=4),
#   mesh= xmesh,
#   k_folds = 8,   spatiotemporal = "ar1",   spatial = "on",   time = "Month",   control = sdmTMBcontrol(newton_loops = 1)
# )
# 
# 
# 
# 
# ##Runn CV's and compare major effects with best smoother
# cv_intercept$sum_loglik
# cv_depth$sum_loglik
# cv_depth2$sum_loglik
# cv_depth3$sum_loglik
# cv_depth4$sum_loglik
# cv_depth5$sum_loglik
# cv_depth6$sum_loglik
# cv_depth7$sum_loglik
# cv_depth8$sum_loglik
# cv_depth9$sum_loglik
# cv_depth10$sum_loglik




##_________________For Model 7:______________________________
##Obs need to make some small edits in the CV-function, for correlation structure== "IID" and rerun dataset xdat.model <- xdat.agg
## Also here we ran on the models above though without Hour and Week. While the best model was not exactly the same across the three datasets, we decided to go with the full model (e.g. bet from full-resolution dataset models above) with adjusted smoothers as the common-ground variable combination for the manuscript. 

bound.outer = 5
max.edge=0.95
mesh2 = inla.mesh.2d(loc=cbind(xdat.model$X, xdat.model$Y),
                     max.edge = c(1,20)*max.edge,
                     cutoff =3, 
                     offset = c(max.edge, bound.outer))
xmesh <- make_mesh(xdat.model, c("X", "Y"), mesh = mesh2)
xdat.model$Year <- as.factor(xdat.model$Year)


## testing model 7:
mx <- sdmTMB(
  formula = kformula, ###
  data = xdat.model,
  mesh = xmesh,
  family = student(link = "identity", df=4),
  #knots =knots ,
  spatiotemporal = "iid",
  spatial = "on",
  time = "Month",
  control = sdmTMBcontrol(newton_loops = 1))

sanity(mx)
summary(mx)
tidy(mx, conf.int = TRUE)

resids <- residuals(mx) 
qqnorm(resids);qqline(resids)



## ____Loop for retriel untill no error
tmb.store <- NULL
for (x in 1) {
  success <- FALSE
  while (!success) {
    tryCatch({
      tmb.store <- c(tmb.store, cval.TMB(xdat = xdat.model, k = 4, xformula = kformula, xsize = 5))
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred. Retrying...\n")
      Sys.sleep(5)  # Wait for 5 seconds before retrying
    })
  }
}

residuals <- residuals(mx)
acf(residuals)
qqnorm(residuals);qqline(residuals)



