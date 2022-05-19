############## Manuscript Diversity and Distribution ######################

# Title:
# Functional biogeography of coastal marine invertebrates along the Southeastern Pacific Coast
# 
# Short title:Functional biogeography of coastal marine invertebrates
# 
# Authors
# D. Herrera, S. A. Navarrete, F. Labra, S. Castillo, L. F. Opazo


#####################################################################################################



# We run the analyses on the different hypotheses to explain the latitudinal
# gradient in different diversity measures along the Eastern pacific coast
# 
# 

# OK, we want to look at the potential effects of the following hypotheses and causal variables:
#        Hypothesis     Variable    
# 1)Species -area	        Area (km2)
# 2)Species - energy	    SST (°C)
# 3)Species - tolerance	  O2 (mgm-3)
# 4)Species - energy	    NO3 (mgm-3)
# 5)Species - energy	    C ( g/m-3 d)
# 6)Species - energy	    Chl a (mgm-3)

#  We need to take into consideration the autocorrelation structure of these 5 variables
# To do so, we will use spdev package to carry out the modelling
# This is what Rivadeneira et al (2011) state:
# 
# The association between species richness, oceanographic and ecological variables was explored using both traditional ordinary least square (OLS) and a simultaneous autoregressive
# model (SARerr) (Dormann et al., 2007; Kissling & Carl, 2008). Analyses were carried out using the library spdev in the software
# R (R Development Core Team, 2010), following the procedures detailed by Kissling & Carl (2008). The best models were
# selected using the Akaike information criterion (DAIC < 2).

# We will try an alternative approach that should not be impacted by variable collinearity.
# Random Forest


##### Loading required libraries #####
library(openxlsx)
library(caret)
library(randomForest)
library(usdm)               # To calculate variance inflation factors (VIF)
library(corrplot)           # To generate correlogram plot
library(RColorBrewer)       # To allow acces to the ColorBrewer library for plots
library(ggplot2)            # To plot results
library(ncf)                # TO calculate spatial autocorrelograms




##### Loading the data #####

data<- read.csv('data_env.csv', na.string='n/a')

# We subset the data
 
Predictors<-as.data.frame(cbind((data$Area_GEBCO),data$SST,data$O2,data$Chl_mean,data$Carbon,data$Nitrate))
colnames(Predictors)<-c("Area","SST","O2","Chl","C","NO3")


colnames(data)

# We build the data frames
NSp<-data$Species.richness
FSpe<-data$FSpe
FRich<-data$FRic
FRed<-data$FRed
FEve<-data$FEve
FD_Rao<-data$RaoQ

NSp_data<-cbind(NSp,Predictors)
FSpe_data<-cbind(FSpe,Predictors)
FRich_data<-cbind(FRich,Predictors)
FRed_data<-cbind(FRed,Predictors)
FEve_data<-cbind(FEve,Predictors)
FD_Rao_data<-cbind(FD_Rao,Predictors)



# We start by generating the visualization of the pairwise correlation
# among the different variables
# 

M0<-cor(Predictors)


cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat0 <- cor.mtest(Predictors)


corrplot(M0, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat0, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)



## We first evaluate collinearity in the environmental predictors 
#  of the latitudinal diversity gradient (LDG) on species richness
# and functional diversity


#  Here we assess the VIF for the variables 
vif(Predictors)
tvif<-vifstep(Predictors,th=100) 
tvif  ## 'best variables'
goodvars=tvif@results$Variables
env.target=cbind(Predictors[names(Predictors)%in%goodvars])
# set1=cbind(env.target,spp=spp)



#  The degree of predictor multicollinearity among variables, measured as the 
# variance inflation factor (Dormann et al., 2013; Naimi et al., 2014), 
# was above the commonly used threshold of ten. This was particularly high for
# log10(area), SST, Oxygen and marginally so for chl-a. 
# Exclusion of the variable with highest VIF (O2) shows that VIF for SST is 13.43 
# indicating a problematic amount of collinearity for this  variable (Gareth et al (2013).
# In contrast, all remaining variables have acceptable levels of collinearity.
# Gareth, J., Daniela, W., Trevor, H., & Robert, T. (2013). An introduction 
# to statistical learning: with applications in R. Spinger.


# We Fit the Ranger Random Forest Model for the different response variables, 
# as follows:

###################
# 1) Species Richness
library(ranger) # To fit conditional random forest model
library(pdp)    # To calculate partial dependence

## Ranger Random forest model with all predictors included
set.seed(11) # to get our same results
rg.general_NSp<-ranger(NSp~., data=NSp_data,importance = "permutation")
rg.general_NSp
rg.general_NSp.imp<-importance_pvalues(rg.general_NSp,method="altmann",formula = NSp~., data=NSp_data,num.permutations=10000) # ...this will take a while, go get a cup of coffee
rg.general_NSp.imp<-rg.general_NSp.imp[order(rownames(rg.general_NSp.imp)),]
rg.general_NSp.imp

write.csv(rg.general_NSp.imp,"Conditional_importance_NSp.imp.csv")

# Diagnostic analyses of the model
residuals.rg.model_NSp<-(NSp_data$NSp-predict(rg.general_NSp,NSp_data)$predictions)/sd(NSp_data$NSp)
spline.model.general_NSp<-spline.correlog(x=data$Longitude,y=data$Latitude,z=residuals.rg.model_NSp,latlon=T,resamp=5000,quiet=T,save=T)

# Residual plot data frame
NSP_plot<-data.frame(Lat=data$Latitude,
                     Residuals=residuals.rg.model_NSp
)

write.csv(NSP_plot,"Residual_plot_NSp.imp.csv")

# Partial dependence plots
NSp_p1 <- partial(rg.general_NSp, pred.var = "Area", plot = F, rug = TRUE)
NSp_p2 <- partial(rg.general_NSp, pred.var = "C", plot = F, rug = TRUE)
NSp_p3 <- partial(rg.general_NSp, pred.var = "Chl", plot = F, rug = TRUE)
NSp_p4 <- partial(rg.general_NSp, pred.var = "NO3", plot = F, rug = TRUE)
NSp_p5 <- partial(rg.general_NSp, pred.var = "O2", plot = F, rug = TRUE)
NSp_p6 <- partial(rg.general_NSp, pred.var = "SST", plot = F, rug = TRUE)

partial_NSp <- cbind(NSp_p1,NSp_p2,NSp_p3,NSp_p4,NSp_p5,NSp_p6)

write.csv(partial_NSp,"Partial_dependence_data_NSp.imp.csv")


# Autocorrelogram plot data frames
NSP_corrplot<-data.frame(Pred_x=as.vector(unlist(spline.model.general_NSp$real$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_NSp$real$predicted$y)))

NSP_corrboot<-data.frame(Boot_x=as.vector(unlist(spline.model.general_NSp$boot$boot.summary$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_NSp$real$predicted$y)),
                         UCI=as.vector(unlist(spline.model.general_NSp$boot$boot.summary$predicted$y[2,])))

# plot(spline.model.general_NSp)
NSP_fig<-ggplot(data=NSP_corrplot,aes(x=Pred_x,y=Pred_y))+geom_line(size=1)+
  geom_ribbon(data=NSP_corrboot,aes(x=Boot_x,ymin=UCI,ymax=NSP_corrplot$Pred_y+abs(NSP_corrplot$Pred_y-UCI)),alpha=0.15)+
  xlab("Distance (km)")+ylab("Moran's I") + geom_hline(yintercept=0)+xlim(0,6200)+
  ylim(-1,1)+theme_classic()+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8))
NSP_fig

Moran_NSP= cbind(NSP_corrplot, NSP_corrboot)
write.csv(Moran_NSP,"Moran_NSP.csv")

###################
# 2) Functional Specialization
# library(ranger)
## Ranger Random forest model with all predictors included
set.seed(11) # to get our same results
rg.general_FSpe<-ranger(FSpe~., data=FSpe_data,importance = "permutation")
rg.general_FSpe
rg.general_FSpe.imp<-importance_pvalues(rg.general_FSpe,method="altmann",formula = FSpe~., data=FSpe_data,num.permutations=10000) # ...this will take a while, go get a cup of coffee
rg.general_FSpe.imp<-rg.general_FSpe.imp[order(rownames(rg.general_FSpe.imp)),]
rg.general_FSpe.imp

write.csv(rg.general_FSpe.imp,"Conditional_importance_FSpe.imp.csv")


# Diagnostic analyses of the model
residuals.rg.model_FSpe<-(FSpe_data$FSpe-predict(rg.general_FSpe,FSpe_data)$predictions)/sd(FSpe_data$FSpe)
spline.model.general_FSpe<-spline.correlog(x=data$Longitude,y=data$Latitude,z=residuals.rg.model_FSpe,latlon=T,resamp=5000,quiet=T,save=T)

# Residual plot data frame
FSpe_plot<-data.frame(Lat=data$Latitude,
                     Residuals=residuals.rg.model_FSpe
)

write.csv(FSpe_plot,"Residual_plot_FSpe.imp.csv")


# Partial dependence plots
FSpe_p1 <- partial(rg.general_FSpe, pred.var = "Area", plot = F, rug = TRUE)
FSpe_p2 <- partial(rg.general_FSpe, pred.var = "C", plot = F, rug = TRUE)
FSpe_p3 <- partial(rg.general_FSpe, pred.var = "Chl", plot = F, rug = TRUE)
FSpe_p4 <- partial(rg.general_FSpe, pred.var = "NO3", plot = F, rug = TRUE)
FSpe_p5 <- partial(rg.general_FSpe, pred.var = "O2", plot = F, rug = TRUE)
FSpe_p6 <- partial(rg.general_FSpe, pred.var = "SST", plot = F, rug = TRUE)

partial_FSpe <- cbind(FSpe_p1,FSpe_p2,FSpe_p3,FSpe_p4,FSpe_p5,FSpe_p6)

# We save these as CSV files
write.csv(partial_FSpe,"partial_FSpe.csv",row.names = FALSE)


# Autocorrelogram plot data frames
FSpe_corrplot<-data.frame(Pred_x=as.vector(unlist(spline.model.general_FSpe$real$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_FSpe$real$predicted$y))                         )

FSpe_corrboot<-data.frame(Boot_x=as.vector(unlist(spline.model.general_FSpe$boot$boot.summary$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_FSpe$real$predicted$y)),
                         UCI=as.vector(unlist(spline.model.general_FSpe$boot$boot.summary$predicted$y[2,])))

# plot(spline.model.general_FSpe)
FSpe_fig<-ggplot(data=FSpe_corrplot,aes(x=Pred_x,y=Pred_y))+geom_line(size=1)+
  geom_ribbon(data=FSpe_corrboot,aes(x=Boot_x,ymin=UCI,ymax=FSpe_corrplot$Pred_y+abs(FSpe_corrplot$Pred_y-UCI)),alpha=0.15)+
  xlab("Distance (km)")+ylab("Moran's I") + geom_hline(yintercept=0)+xlim(0,6200)+
  ylim(-1,1)+theme_classic()+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8))
FSpe_fig

Moran_FSpe_corrboot= cbind(FSpe_corrboot, FSpe_corrboot)
write.csv(Moran_FSpe_corrboot,"Moran_FSpe_corrboot.csv")


###################
# 3) Functional Richness
rg.general_FRich<-ranger(FRich~., data=FRich_data,importance = "permutation")
rg.general_FRich
rg.general_FRich.imp<-importance_pvalues(rg.general_FRich,method="altmann",formula = FRich~., data=FRich_data,num.permutations=10000) # ...this will take a while, go get a cup of coffee
rg.general_FRich.imp<-rg.general_FRich.imp[order(rownames(rg.general_FRich.imp)),]
rg.general_FRich.imp

write.csv(rg.general_FRich.imp,"Conditional_importance_FRich.imp.csv")


# Diagnostic analyses of the model
residuals.rg.model_FRich<-(FRich_data$FRich-predict(rg.general_FRich,FRich_data)$predictions)/sd(FRich_data$FRich)
spline.model.general_FRich<-spline.correlog(x=data$Longitude,y=data$Latitude,z=residuals.rg.model_FRich,latlon=T,resamp=5000,quiet=T,save=T)

# Residual plot data frame
FRich_plot<-data.frame(Lat=data$Latitude,
                     Residuals=residuals.rg.model_FRich
)

write.csv(FRich_plot,"Residual_plot_FRich.imp.csv")



# Partial dependence plots
FRich_p1 <- partial(rg.general_FRich, pred.var = "Area", plot = F, rug = TRUE)
FRich_p2 <- partial(rg.general_FRich, pred.var = "C", plot = F, rug = TRUE)
FRich_p3 <- partial(rg.general_FRich, pred.var = "Chl", plot = F, rug = TRUE)
FRich_p4 <- partial(rg.general_FRich, pred.var = "NO3", plot = F, rug = TRUE)
FRich_p5 <- partial(rg.general_FRich, pred.var = "O2", plot = F, rug = TRUE)
FRich_p6 <- partial(rg.general_FRich, pred.var = "SST", plot = F, rug = TRUE)

partial_FRich <- cbind(FRich_p1,FRich_p2,FRich_p3,FRich_p4,FRich_p5,FRich_p6)

# We save these as CSV files
write.csv(partial_FRich,"partial_FRich.csv")


# Autocorrelogram plot data frames
FRich_corrplot<-data.frame(Pred_x=as.vector(unlist(spline.model.general_FRich$real$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_FRich$real$predicted$y))                         )

FRich_corrboot<-data.frame(Boot_x=as.vector(unlist(spline.model.general_FRich$boot$boot.summary$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_FRich$real$predicted$y)),
                         UCI=as.vector(unlist(spline.model.general_FRich$boot$boot.summary$predicted$y[2,])))

# plot(spline.model.general_FRich)
FRich_fig<-ggplot(data=FRich_corrplot,aes(x=Pred_x,y=Pred_y))+geom_line(size=1)+
  geom_ribbon(data=FRich_corrboot,aes(x=Boot_x,ymin=UCI,ymax=FRich_corrplot$Pred_y+abs(FRich_corrplot$Pred_y-UCI)),alpha=0.15)+
  xlab("Distance (km)")+ylab("Moran's I") + geom_hline(yintercept=0)+xlim(0,6200)+
  ylim(-1,1)+theme_classic()+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8))
FRich_fig

Moran_FRich_corrboot= cbind(FRich_corrplot, FRich_corrboot)
write.csv(Moran_FRich_corrboot,"Moran_FRich_corrboot.csv")


###################
# 4) Functional group redundancy
rg.general_FRed<-ranger(FRed~., data=FRed_data,importance = "permutation")
rg.general_FRed
rg.general_FRed.imp<-importance_pvalues(rg.general_FRed,method="altmann",formula = FRed~., data=FRed_data,num.permutations=10000) # ...this will take a while, go get a cup of coffee
rg.general_FRed.imp<-rg.general_FRed.imp[order(rownames(rg.general_FRed.imp)),]
rg.general_FRed.imp

write.csv(rg.general_FRed.imp,"Conditional_importance_FRed.imp.csv")


# Diagnostic analyses of the model
residuals.rg.model_FRed<-(FRed_data$FRed-predict(rg.general_FRed,FRed_data)$predictions)/sd(FRed_data$FRed)
spline.model.general_FRed<-spline.correlog(x=data$Longitude,y=data$Latitude,z=residuals.rg.model_FRed,latlon=T,resamp=5000,quiet=T,save=T)

# Residual plot data frame
FRed_plot<-data.frame(Lat=data$Latitude,
                     Residuals=residuals.rg.model_FRed
)

write.csv(FRed_plot,"Residual_plot_FRed.imp.csv")

# Partial dependence plots
FRed_p1 <- partial(rg.general_FRed, pred.var = "Area", plot = F, rug = TRUE)
FRed_p2 <- partial(rg.general_FRed, pred.var = "C", plot = F, rug = TRUE)
FRed_p3 <- partial(rg.general_FRed, pred.var = "Chl", plot = F, rug = TRUE)
FRed_p4 <- partial(rg.general_FRed, pred.var = "NO3", plot = F, rug = TRUE)
FRed_p5 <- partial(rg.general_FRed, pred.var = "O2", plot = F, rug = TRUE)
FRed_p6 <- partial(rg.general_FRed, pred.var = "SST", plot = F, rug = TRUE)

partial_FRed <- cbind(FRed_p1,FRed_p2,FRed_p3,FRed_p4,FRed_p5,FRed_p6)

# We save these as CSV files
write.csv(partial_FRed,"partial_FRed.csv")


# Autocorrelogram plot data frames
FRed_corrplot<-data.frame(Pred_x=as.vector(unlist(spline.model.general_FRed$real$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_FRed$real$predicted$y))                         )

FRed_corrboot<-data.frame(Boot_x=as.vector(unlist(spline.model.general_FRed$boot$boot.summary$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_FRed$real$predicted$y)),
                         UCI=as.vector(unlist(spline.model.general_FRed$boot$boot.summary$predicted$y[2,])))

# plot(spline.model.general_FRed)
FRed_fig<-ggplot(data=FRed_corrplot,aes(x=Pred_x,y=Pred_y))+geom_line(size=1)+
  geom_ribbon(data=FRed_corrboot,aes(x=Boot_x,ymin=UCI,ymax=FRed_corrplot$Pred_y+abs(FRed_corrplot$Pred_y-UCI)),alpha=0.15)+
  xlab("Distance (km)")+ylab("Moran's I") + geom_hline(yintercept=0)+xlim(0,6200)+
  ylim(-1,1)+theme_classic()+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8))
FRed_fig


Moran_FRed_corrboot= cbind(FRed_corrplot, FRed_corrboot)
write.csv(Moran_FRed_corrboot,"Moran_FRed_corrboot.csv")


###################

# 5) Functional group evenness
rg.general_FEve<-ranger(FEve~., data=FEve_data,importance = "permutation")
rg.general_FEve
rg.general_FEve.imp<-importance_pvalues(rg.general_FEve,method="altmann",formula = FEve~., data=FEve_data,num.permutations=10000) # ...this will take a while, go get a cup of coffee
rg.general_FEve.imp<-rg.general_FEve.imp[order(rownames(rg.general_FEve.imp)),]
rg.general_FEve.imp

write.csv(rg.general_FEve.imp,"Conditional_importance_FEve.imp.csv")

# Diagnostic analyses of the model
residuals.rg.model_FEve<-(FEve_data$FEve-predict(rg.general_FEve,FEve_data)$predictions)/sd(FEve_data$FEve)
spline.model.general_FEve<-spline.correlog(x=data$Longitude,y=data$Latitude,z=residuals.rg.model_FEve,latlon=T,resamp=5000,quiet=T,save=T)

# Residual plot data frame
FEve_plot<-data.frame(Lat=data$Latitude,
                     Residuals=residuals.rg.model_FEve
)

write.csv(FEve_plot,"Residual_plot_FEve.imp.csv")


# Partial dependence plots
FEve_p1 <- partial(rg.general_FEve, pred.var = "Area", plot = F, rug = TRUE)
FEve_p2 <- partial(rg.general_FEve, pred.var = "C", plot = F, rug = TRUE)
FEve_p3 <- partial(rg.general_FEve, pred.var = "Chl", plot = F, rug = TRUE)
FEve_p4 <- partial(rg.general_FEve, pred.var = "NO3", plot = F, rug = TRUE)
FEve_p5 <- partial(rg.general_FEve, pred.var = "O2", plot = F, rug = TRUE)
FEve_p6 <- partial(rg.general_FEve, pred.var = "SST", plot = F, rug = TRUE)

partial_FEve <- cbind(FEve_p1,FEve_p2,FEve_p3,FEve_p4,FEve_p5,FEve_p6)

# We save these as CSV files
write.csv(partial_FEve,"partial_FEve.csv")


# Autocorrelogram plot data frames
FEve_corrplot<-data.frame(Pred_x=as.vector(unlist(spline.model.general_FEve$real$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_FEve$real$predicted$y))                         )

FEve_corrboot<-data.frame(Boot_x=as.vector(unlist(spline.model.general_FEve$boot$boot.summary$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_FEve$real$predicted$y)),
                         UCI=as.vector(unlist(spline.model.general_FEve$boot$boot.summary$predicted$y[2,])))

# plot(spline.model.general_FEve)
FEve_fig<-ggplot(data=FEve_corrplot,aes(x=Pred_x,y=Pred_y))+geom_line(size=1)+
  geom_ribbon(data=FEve_corrboot,aes(x=Boot_x,ymin=UCI,ymax=FEve_corrplot$Pred_y+abs(FEve_corrplot$Pred_y-UCI)),alpha=0.15)+
  xlab("Distance (km)")+ylab("Moran's I") + geom_hline(yintercept=0)+xlim(0,6200)+
  ylim(-1,1)+theme_classic()+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8))
FEve_fig


Moran_FEve_corrboot= cbind(FEve_corrplot, FEve_corrboot)
write.csv(Moran_FEve_corrboot,"Moran_FEve_corrboot.csv")


###################
# 6) Functional diversity (Rao)
rg.general_FD_Rao<-ranger(FD_Rao~., data=FD_Rao_data,importance = "permutation")
rg.general_FD_Rao
rg.general_FD_Rao.imp<-importance_pvalues(rg.general_FD_Rao,method="altmann",formula = FD_Rao~., data=FD_Rao_data,num.permutations=10000) # ...this will take a while, go get a cup of coffee
rg.general_FD_Rao.imp<-rg.general_FD_Rao.imp[order(rownames(rg.general_FD_Rao.imp)),]
rg.general_FD_Rao.imp

write.csv(rg.general_FD_Rao.imp,"Conditional_importance_FD_Rao_imp.csv")


# Diagnostic analyses of the model
residuals.rg.model_FD_Rao<-(FD_Rao_data$FD_Rao-predict(rg.general_FD_Rao,FD_Rao_data)$predictions)/sd(FD_Rao_data$FD_Rao)
spline.model.general_FD_Rao<-spline.correlog(x=data$Longitude,y=data$Latitude,z=residuals.rg.model_FD_Rao,latlon=T,resamp=5000,quiet=T,save=T)

# Residual plot data frame
FD_Rao_plot<-data.frame(Lat=data$Latitude,
                     Residuals=residuals.rg.model_FD_Rao
)


write.csv(FD_Rao_plot,"Residual_plot_Rao.imp.csv")

# Partial dependence plots
FD_Rao_p1 <- partial(rg.general_FD_Rao, pred.var = "Area", plot = F, rug = TRUE)
FD_Rao_p2 <- partial(rg.general_FD_Rao, pred.var = "C", plot = F, rug = TRUE)
FD_Rao_p3 <- partial(rg.general_FD_Rao, pred.var = "Chl", plot = F, rug = TRUE)
FD_Rao_p4 <- partial(rg.general_FD_Rao, pred.var = "NO3", plot = F, rug = TRUE)
FD_Rao_p5 <- partial(rg.general_FD_Rao, pred.var = "O2", plot = F, rug = TRUE)
FD_Rao_p6 <- partial(rg.general_FD_Rao, pred.var = "SST", plot = F, rug = TRUE)

partial_FD_Rao <- cbind(FD_Rao_p1,FD_Rao_p2,FD_Rao_p3,FD_Rao_p4,FD_Rao_p5,FD_Rao_p6)

# We save these as CSV files
write.csv(partial_FD_Rao,"partial_FD_Rao.csv")


# Autocorrelogram plot data frames
FD_Rao_corrplot<-data.frame(Pred_x=as.vector(unlist(spline.model.general_FD_Rao$real$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_FD_Rao$real$predicted$y))                         )

FD_Rao_corrboot<-data.frame(Boot_x=as.vector(unlist(spline.model.general_FD_Rao$boot$boot.summary$predicted$x)),
                         Pred_y=as.vector(unlist(spline.model.general_FD_Rao$real$predicted$y)),
                         UCI=as.vector(unlist(spline.model.general_FD_Rao$boot$boot.summary$predicted$y[2,])))

# plot(spline.model.general_FD_Rao)
FD_Rao_fig<-ggplot(data=FD_Rao_corrplot,aes(x=Pred_x,y=Pred_y))+geom_line(size=1)+
  geom_ribbon(data=FD_Rao_corrboot,aes(x=Boot_x,ymin=UCI,ymax=FD_Rao_corrplot$Pred_y+abs(FD_Rao_corrplot$Pred_y-UCI)),alpha=0.15)+
  xlab("Distance (km)")+ylab("Moran's I") + geom_hline(yintercept=0)+xlim(0,6200)+
  ylim(-1,1)+theme_classic()+theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8))
FD_Rao_fig

Moran_FD_Rao_corrboot= cbind(FD_Rao_corrplot, FD_Rao_corrboot)
write.csv(Moran_FD_Rao_corrboot,"Moran_FD_Rao_corrboot.csv")


######## We now want to to calculate the predicted values #############

predicted<-cbind(predict(rg.general_NSp,NSp_data)$predictions,
      predict(rg.general_FSpe,FSpe_data)$predictions,
      predict(rg.general_FRich,FRich_data)$predictions,
      predict(rg.general_FRed,FRed_data)$predictions,
      predict(rg.general_FEve,FEve_data)$predictions,
      predict(rg.general_FD_Rao,FD_Rao_data)$predictions)

predicted
write.csv(predicted,"predicted.csv")



# Let's add all the residuals into a single data frame:

Residuals<-cbind(NSP_plot,FSpe_plot[,2],FRich_plot[,2],FRed_plot[,2],FEve_plot[,2],FD_Rao_plot[,2])
Residuals

Correlograma<-cbind(NSP_corrplot,FSpe_corrplot,FRich_corrplot,FRed_corrplot,FEve_corrplot,FD_Rao_corrplot)

write.csv(Correlograma,"Correlogram.csv")



Corr_Boot <- cbind(NSP_corrboot$Boot_x,NSP_corrboot$UCI,
      NSP_corrplot$Pred_y+abs(NSP_corrplot$Pred_y-NSP_corrboot$UCI),
      FSpe_corrboot$UCI,
      FSpe_corrplot$Pred_y+abs(FSpe_corrplot$Pred_y-FSpe_corrboot$UCI),FRich_corrboot$UCI,
      FRich_corrplot$Pred_y+abs(FRich_corrplot$Pred_y-FRich_corrboot$UCI),
      FRed_corrboot$UCI,
      FRed_corrplot$Pred_y+abs(FRed_corrplot$Pred_y-FRed_corrboot$UCI),
      FEve_corrboot$UCI,
      FEve_corrplot$Pred_y+abs(FEve_corrplot$Pred_y-FEve_corrboot$UCI),
      FD_Rao_corrboot$UCI,
      FD_Rao_corrplot$Pred_y+abs(FD_Rao_corrplot$Pred_y-FD_Rao_corrboot$UCI))
write.csv(Corr_Boot,"Corr_Boot.csv")


