# Functional biogeography of coastal marine invertebrates along the Southeastern Pacific Coast

### Authors: D. Herrera, S. A. Navarrete, F. Labra, S. P. Castillo, L. F. Opazo

Contact: [Felipe Opazo-Mella](felipe.opazomella@gmail.com)

The complete routine is available in function/main.R it can be opened downloading the whole repo and using the available R enviroment.
All the codes and outputs -including null models results- are avalible in the [drive folder](https://drive.google.com/drive/folders/1BKJpW3I3InCsigit8RL1H7ZfRCpYx1lG?usp=sharing)

#### 0. Preamble
To load the necessary packages to run this analysis, we use the pacman package
```{r}
if (!require("pacman")) install.packages("pacman")

pacman::p_load(foreach,doParallel, vegan, SYNCSA, tidyr,tseries,strucchange, svMisc, dplyr, mFD)
```
We also define useful functions available in this repository:
```{r}
'%ni%'<- Negate('%in%')
source('functions/FD_df.R')
```
Finally, we create folder for the output files:
```{r}
# Make directory for output
dir.create('data_output')
```
#### 1. Load data
The dataset used in this pipeline and in the manuscript is available and it is loaded into the routine via
```{r}
df0=read.csv('data_input/latdata.csv', check.names = F, row.names = 1)

```

#### 2. Observed funcitonal diversity metrics
To compute the functional diversity metrics, we create a function called ```.fdmetrics``` depending upon the ```mFD``` package:

```{r}

.fdmetrics<- function(df = NULL, nPC, maxPcoa,features, nom.features=NULL,ord.features=NULL,quan.features=NULL ,weight = NULL, coord=NULL){
  
  # df = dataframe of relative abundances with dimensions time x ecocodes.
  # nPC = number of PC axes to include in the computation of functional diversity indices
  # maxPcoa = maximum number of axis used for PCoa
  # features = max value for each feature
  # nom.features = indices of the features that are nominal. Default is NULL
  # ord.features = indices of the features that are ordinal Default is NULL
  # quan.features = indices of the features that are quantitative Default is NULL
  # weight = matrix of abundances with dimensions time x ecocodes. If NULL, it is calculated from df.
  # coord = matrix of features/traits with dimensions ecocodes x features. If NULL, it is calculated from df.
  
  if(any(quan.features == nom.features)) stop('A feature cannot be nominal and quantitative')
  if(any(quan.features == ord.features)) stop('A feature cannot be ordinal and quantitative')
  if(any(nom.features == ord.features)) stop('A feature cannot be ordinal and nominal')
  

  m0<- df
  df<-(t(df)/rowSums(t(df)))
  df<-FD_df(as.data.frame(df0), features = features)

  if(is.null(weight) && is.null(coord)){

    weight <- pivot_wider(df[order(as.numeric(df$time)),c('ecocode', 'time', 'abundance')],
                          names_from = 'ecocode', values_from = 'abundance',id_cols = 'time' )
    
    weight<-data.frame(weight, check.names = FALSE)
    rownames(weight)<-weight$time
    weight<-weight[-1]
    weight <- as.matrix(weight)
    rownames(weight)= paste0('t.', rownames(weight))

    coord  <- subset(df,select = 1:(nchar(as.character(df$ecocode[1]))+1))
    coord<- unique(coord)
    rownames(coord)<- coord[,1]
    coord<-coord[,-1]
 
    for(i in quan.features){
      coord[,i] = as.numeric(coord[,i])
    }

    if(!is.null(nom.features)){
      for (n in nom.features) {
        coord[,n] = as.factor(coord[,n])
      }
      
    }
    
    if(!is.null(ord.features)){
      for (i in ord.features) {
        coord[,i] = factor(coord[,i],ordered = TRUE )
      }
    }
  }
  
  coord_cat = data.frame(trait_name= colnames(coord), trait_type=NA, trait_weight=1, fuzzy_name=NA)
  
  coord_cat[nom.features,'trait_type'] = 'N'
  coord_cat[ord.features,'trait_type'] = 'O'
  coord_cat[quan.features,'trait_type'] = 'Q'

  sp_dist = mFD::funct.dist(sp_tr         = coord,
                            tr_cat        = coord_cat,
                            metric        = "gower",
                            scale_euclid  = "scale_center",
                            ordinal_var   = "classic",
                            weight_type   = "equal",
                            stop_if_NA    = TRUE)

  fspaces_quality <- mFD::quality.fspaces(
    sp_dist             = sp_dist,
    maxdim_pcoa         = maxPcoa,
    deviation_weighting = 'absolute',
    fdist_scaling       = FALSE,
    fdendro             = 'average')

  sp_faxes_coord <- fspaces_quality$details_fspaces$sp_pc_coord

  alpha_fd_indices <- mFD::alpha.fd.multidim(
    
    sp_faxes_coord   = sp_faxes_coord[, 1:nPC],
    asb_sp_w         = weight,
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  fd_ind_values <- alpha_fd_indices$functional_diversity_indices
  pc_axes <<- sp_faxes_coord[, 1:nPC]
  
  obsRao<-rao.diversity(traits=coord, weight)

  fd_ind_values$fred <-(obsRao$Simpson-obsRao$FunRao)/obsRao$Simpson
  fd_ind_values$raoQ <-obsRao$FunRao
  fd_ind_values$simpson<- obsRao$Simpson
  rownames(fd_ind_values) =substring(rownames(fd_ind_values), 3, 10000L)
  fd_ind_values = data.frame(time= rownames(fd_ind_values), fd_ind_values)
  
  return(fd_ind_values)
}

#
```
We run the function on the empirical dataset by calling:

```{r}
obsFD<- .fdmetrics(df = df0 , #dataframe of relative abundances with dimensions time x ecocodes
                   nPC=4,   # nPC = number of PC axes to include in the computation of functional diversity indices
                   maxPcoa=10, # maxPcoa = maximum number of axis used for PCoa
                   features = c(4,6,6,3,2,3,3,4), # features = max value for each feature
                   nom.features = 2:8, #vector of indices of nominal features in features
                   ord.features = 1, #vector of indices of ordinal features in features
                   quan.features = NULL #vector of indices of quantitative features in features
                   
)

write.csv(obsFD, file='data_output/obsFD.csv') #save observed functional diversity metrics
write.csv(pc_axes, file='data_output/pc_axes.csv') #save values of PCA's axes selected by nPC
```

#### 3. Null models
To run null models we create a function that parallelize the process

```{r, tidy=TRUE, tidy.opts=list(width.cutoff=60)}
.nulls<- function(null0, rep.nulls, maxPcoa, nPC, features, nom.features, ord.features, quan.features, numCores){

  if(any(quan.features == nom.features)) stop('A feature cannot be nominal and quantitative')
  if(any(quan.features == ord.features)) stop('A feature cannot be ordinal and quantitative')
  if(any(nom.features == ord.features)) stop('A feature cannot be ordinal and nominal')
 
  nullsdf<- data.frame()
  registerDoParallel(numCores)  # use multicore, set to the number of our cores

  foreach (p=1:rep.nulls, .combine = rbind) %dopar% {
    
    source('functions/FD_df.R')
    pacman::p_load(foreach,doParallel, vegan, SYNCSA, tidyr,tseries,strucchange, svMisc, dplyr, mFD)    
    
    print(paste0('null iteration: ',p, ' out of ',rep.nulls ))
    
    m0<- null0[,,p]    
    m1<-((m0)/rowSums((m0)))
    datam<-FD_df(as.data.frame(m1), features)
    df<-datam
    weight <- pivot_wider(df[order(as.numeric(df$time)),c('ecocode', 'time', 'abundance')],
                          names_from = 'ecocode', values_from = 'abundance',id_cols = 'time' )
    weight<-data.frame(weight, check.names = FALSE)
    rownames(weight)<-weight$time
    weight<-weight[-1]
    weight <- as.matrix(weight)
    rownames(weight)= paste0('t.', rownames(weight))
   
    coord  <- subset(df,select = 1:(nchar(as.character(df$ecocode[1]))+1))
    coord<- unique(coord)
    rownames(coord)<- coord[,1]
    coord<-coord[,-1]

    for(i in quan.features){
      
      coord[,i] = as.numeric(coord[,i])
      
    }
    
    for (n in nom.features) {
      
      coord[,n] = as.factor(coord[,n])
      
    }
   
    if(!is.null(ord.features)){
      
      for (i in ord.features) {
        
        coord[,i] = factor(coord[,i],ordered = TRUE )
      }
      
    }
    coord_cat = data.frame(trait_name= colnames(coord), trait_type=NA, trait_weight=1, fuzzy_name=NA)
    coord_cat[nom.features,'trait_type'] = 'N'
    coord_cat[ord.features,'trait_type'] = 'O'
    coord_cat[quan.features,'trait_type'] = 'Q'

    
    
    sp_dist = mFD::funct.dist(sp_tr         = coord,
                              tr_cat        = coord_cat,
                              metric        = "gower",
                              scale_euclid  = "scale_center",
                              ordinal_var   = "classic",
                              weight_type   = "equal",
                              stop_if_NA    = TRUE)
    
    
    
    fspaces_quality <- mFD::quality.fspaces(sp_dist = sp_dist, 
                              maxdim_pcoa = maxPcoa,
                              deviation_weighting = 'absolute',
                              fdist_scaling = FALSE,
                              fdendro = 'average')

    sp_faxes_coord <- fspaces_quality$details_fspaces$sp_pc_coord
    
    alpha_fd_indices <- mFD::alpha.fd.multidim(
      sp_faxes_coord   = sp_faxes_coord[, 1:nPC],
      asb_sp_w         = weight,
      scaling          = TRUE,
      check_input      = TRUE,
      details_returned = TRUE)
    
    
    
    fd_ind_values <- alpha_fd_indices$functional_diversity_indices
    
    
    
    obsRao<-rao.diversity(traits=coord, weight)
   
    fd_ind_values$fred <-(obsRao$Simpson-obsRao$FunRao)/obsRao$Simpson
    fd_ind_values$raoQ <-obsRao$FunRao
    fd_ind_values$simpson<- obsRao$Simpson
    rownames(fd_ind_values) =substring(rownames(fd_ind_values), 3, 10000L)
    fd_ind_values$time = rownames(fd_ind_values)
    mdimDF <- fd_ind_values
    nullFunDiv<-data.frame(null=p,mdimDF)
    nullsdf<-rbind(nullsdf, nullFunDiv)
  }  
}
```

 We run the null models code by calling:
 
    ```{r}
    rep.nulls<-5000
    method.nulls<- 'r0_samp'  #for other methods see ?commsim
    nm<-nullmodel(df0,method.nulls)  ## Df of counts NOT proportional abundance
    null<-simulate(nm, nsim =rep.nulls)



      nulldf<- .nulls(null0 = null,
                rep.nulls,
                maxPcoa = 10, nPC=4,
                features =  c(4,6,6,3,2,3,3,4), # features = max value for each feature
                nom.features = 2:8, #vector of indices of nominal features in features
                ord.features = 1, #vector of indices of ordinal features in features
                quan.features = NULL, #vector of indices of quantitative features in features
                numCores= 20) #parallel::detectCores() detect your number of cores, NOT SUGGESTED TO RUN AT MAX NUMBER OF CORES

    write.csv(nulldf, file='data_output/allnull.csv')
    ```

#### 4. Outputs

Each time we ran the whole analysis (e.g., ```source(main.R)```), a file is created in the ```data_output``` folder. Additionally, we created a function that summariseds the data from the null models:
```{r}
.nullsummary<-function(truenull){
  
  nullmodelsummary<- truenull%>%
    group_by(time)%>%
    mutate(com=as.numeric(time))%>%
    summarise(msimpson= mean(simpson),
              mraoQ= mean(raoQ),
              mfred = mean(fred),
              mfric = mean(fric),
              mfeve= mean(feve),
              mfspe= mean(fspe),
              mfdis = mean(fdis),
              mfdiv= mean(fdiv),
              mfori = mean(fori),
              mTaxRich = mean(sp_richn),
              
              sd.simpson= sd(simpson, na.rm = TRUE),
              sd.raoQ= sd(raoQ, na.rm = TRUE),
              sd.fred= sd(fred, na.rm = TRUE),
              sd.fric= sd(fric, na.rm = TRUE),
              sd.feve= sd(feve, na.rm = TRUE),
              sd.fspe= sd(fspe, na.rm = TRUE),
              sd.fdis= sd(fdis, na.rm = TRUE),
              sd.fdiv= sd(fdiv, na.rm = TRUE),
              sd.fori= sd(fori, na.rm = TRUE),
              sd.TaxRich= sd(sp_richn, na.rm = TRUE),
              
              n.fric = n())%>%
                      mutate(se.simpson = sd.simpson / sqrt(n.fric),
                             se.raoQ = sd.raoQ / sqrt(n.fric),
                             se.fred = sd.fred / sqrt(n.fric),
                             se.fric = sd.fric / sqrt(n.fric),
                             se.feve = sd.feve / sqrt(n.fric),
                             se.fspe = sd.fspe / sqrt(n.fric),
                             se.fdis = sd.fdis / sqrt(n.fric),
                             se.fdiv = sd.fdiv / sqrt(n.fric),
                             se.fori = sd.fori / sqrt(n.fric),
                             se.TaxRich = sd.TaxRich / sqrt(n.fric),

                             lower.ci.simpson = msimpson - qt(1 - (0.05 / 2), n.fric - 1) * se.simpson,
                             upper.ci.simpson = msimpson + qt(1 - (0.05 / 2), n.fric - 1) * se.simpson,
                             lower.ci.raoQ = mraoQ - qt(1 - (0.05 / 2), n.fric - 1) * se.raoQ,
                             upper.ci.raoQ = mraoQ + qt(1 - (0.05 / 2), n.fric - 1) * se.raoQ,
                             lower.ci.fred = mfred - qt(1 - (0.05 / 2), n.fric - 1) * se.fred,
                             upper.ci.fred = mfred + qt(1 - (0.05 / 2), n.fric - 1) * se.fred,
                             lower.ci.fric = mfric - qt(1 - (0.05 / 2), n.fric - 1) * se.fric,
                             upper.ci.fric = mfric + qt(1 - (0.05 / 2), n.fric - 1) * se.fric,
                             lower.ci.feve = mfeve - qt(1 - (0.05 / 2), n.fric - 1) * se.feve,
                             upper.ci.feve = mfeve + qt(1 - (0.05 / 2), n.fric - 1) * se.feve,
                             lower.ci.fspe = mfspe - qt(1 - (0.05 / 2), n.fric - 1) * se.fspe,
                             upper.ci.fspe = mfspe + qt(1 - (0.05 / 2), n.fric - 1) * se.fspe,
                             lower.ci.fdis = mfdis - qt(1 - (0.05 / 2), n.fric - 1) * se.fdis,
                             upper.ci.fdis = mfdis + qt(1 - (0.05 / 2), n.fric - 1) * se.fdis,
                             lower.ci.fdiv = mfdiv - qt(1 - (0.05 / 2), n.fric - 1) * se.fdiv,
                             upper.ci.fdiv = mfdiv + qt(1 - (0.05 / 2), n.fric - 1) * se.fdiv,
                             lower.ci.fori = mfori - qt(1 - (0.05 / 2), n.fric - 1) * se.fori,
                             upper.ci.fori = mfori + qt(1 - (0.05 / 2), n.fric - 1) * se.fori,
                             lower.ci.TaxRich = mTaxRich - qt(1 - (0.05 / 2), n.fric - 1) * se.TaxRich,
                             upper.ci.TaxRich = mTaxRich + qt(1 - (0.05 / 2), n.fric - 1) * se.TaxRich)

  return(nullmodelsummary)
  
}

nullsummary<- .nullsummary(truenull=nulldf)
write.csv(nullsummary, file='data_output/nullsummary.csv') #Saves the summary

```
