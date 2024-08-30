

# ---------------------------------------------------------------------------------------# 


# Project: Coping with Collapse: Functional Robustness of Coral-Reef Fish Network to Simulated Cascade Extinction


# Script 2: organize data for network analysis and produce a map (supp info)
# (Analyzes based on "occupancy_data.RData")

#   RANDOM ORDER OF CORAL REMOVAL -- SCENARIO CREATED AFTER THE 1ST ROUND OF REVIEW
#   The other scripts (also numbered with '2' bring scenarios of degree- and vulnerability-based removal)



# ---------------------------------------------------------------------------------------# 

rm(list=ls())
# load basic data & functions
require(here)
source(here ("R", "functions.R"))
source(here ("R", "packages.R"))
load(file = here ("processed_data","occupancy_data.RData"))

# set seed
set.seed (1001)


# -----------------------------------------------


# count total N fish 
ns_fish <- length(unique(network_data$peixe))
# total N corals
ns_coral <- length(unique(network_data$coral))


# coral associated fish (following Luza et al. 2022, sci rep) 
network_data <- network_data[which(network_data$low.coral>0 & 
                                     network_data$estimate.turf<=0),]
# N coral associated fish
ns_fish_associated <- length(unique(network_data$peixe)) # ncol(m_web)


# PARTITE A ----------------------------------------------

# transform data into matrix (matrix C x Fca)
m_web <- cast (formula = coral~peixe,value = "pred",
               fun.aggregate = mean,na.rm=T,fill=0,
               data = network_data)
labels_web <- m_web[,1]; m_web<- m_web[,-1] # labels
# adjust
rownames(m_web)<- labels_web
# ordering
m_web<-sortweb(m_web,sort.order="dec")



# functional trait space loss --------------------
# load and organize trait data ----------

source ("R/functions.R")
trait_dataset<- organize_traits (fish_size)

# find coral associated fish & match species and trait data 

sp_analyzed_response <- do.call(rbind, extracted_data) # melt data
# coral associated fish (following Luza et al. 2020, sci rep)
sp_analyzed_response<-(sp_analyzed_response[which(sp_analyzed_response$low.coral>0 & 
                                                    sp_analyzed_response$estimate.turf<=0),])
## bind trait data
sp_analyzed_response <- cbind (sp_analyzed_response,
                               trait_dataset [match (sp_analyzed_response$peixe,  # bind traits
                                                     trait_dataset$scientificName), 
                                              "Body_size"])

# ordering species according to body size
sp_analyzed_response <- sp_analyzed_response[order(sp_analyzed_response$Body_size),]
sp_analyzed_response$age<-as.numeric(sp_analyzed_response$age)
# first name letter to up
sp_analyzed_response$peixe<-firstup(sp_analyzed_response$peixe)
sp_analyzed_response$coral<-firstup(sp_analyzed_response$coral)

# a complete functional space per species of coral
total <- do.call (rbind,extracted_data) # rbind extracted data (with coefficients and CI)


# select traits to analysis
sel_traits <- trait_dataset[which(trait_dataset$scientificName %in% total$peixe),
                            c("Aspect_ratio","Trophic_level","Size_group",
                              "TempPref_max","Depth_max",
                              "log_actual_size",
                              "scientificName")]


# select and organize the data (max value among adult and juvenile)
sel_traits <- (sel_traits %>%
                      group_by(scientificName) %>%
                      summarize (log_actual_size = max(log_actual_size),
                                 Aspect_ratio  = max(Aspect_ratio),
                                 Trophic_level = max(Trophic_level),
                                 Size_group = max(Size_group),
                                 TempPref_max = max(TempPref_max),
                                 Depth_max = max(Depth_max))
)

# the correlation between traits
(cor(sel_traits[,-which(colnames(sel_traits) == "scientificName")], 
     use = "complete.obs"))#  correlation is fine


# PARTITE B --------------------------------------------------------

# secondary extinctions
# organize network data at the site scale
names (extracted_data_occ)<-coral_species

# calculate correlations of occupancy probability
df_occ <- lapply (extracted_data_occ, function (i) 
  melt (cor (i,method = "pearson"))
)
df_occ<-do.call(rbind, df_occ)# melt

# define combinations of spp
df_occ$comb <- paste (df_occ$Var1,df_occ$Var2,sep="_")

# summarize correlations
require(dplyr)
df_occ <- df_occ %>%
  group_by(comb) %>%
  summarize ( cor = mean(value))

# adjust spp names
network_data_occ <- cbind (df_occ, 
                           sp1=do.call(rbind, strsplit (df_occ$comb, "_"))[,1],
                           sp2=do.call(rbind, strsplit (df_occ$comb, "_"))[,2]
)

# back to a wide matrix format (Fca x Fcor)
m_web_occ <- cast (formula = sp1 ~ sp2,
                   data=network_data_occ,
                   value = "cor",
                   fill=0)
labels_web_occ <- m_web_occ[,1]; m_web_occ<- m_web_occ[,-1] # labels
# adjust
rownames(m_web_occ)<- labels_web_occ
# coral associated fish in the row
m_web_occ <- m_web_occ[which(rownames (m_web_occ) %in% colnames (m_web)),
                       which(colnames (m_web_occ) %in% colnames (m_web) == F)]
# ordering
m_web_occ<- m_web_occ [match (colnames(m_web),
                              rownames(m_web_occ)),]
#rownames(m_web_occ) == colnames(m_web)
# correlations lower than 0.8 <-- 0
m_web_occ[m_web_occ<0.8] <- 0 # conservative enought to define association
# rm empty cols
m_web_occ <- m_web_occ [, which(colSums(m_web_occ)>0)]


# ----------------------------------


# select traits for spp with some association to corals and coral -associated fish
sel_traits_fish <- trait_dataset[which(trait_dataset$scientificName %in% unique(unlist(dimnames(m_web_occ)))),
                            c("Aspect_ratio","Trophic_level","Size_group",
                              "TempPref_max","Depth_max",
                              "log_actual_size",
                              "scientificName")]

# select and organize the data
sel_traits_fish <- (sel_traits_fish %>%
  group_by(scientificName) %>%
  summarize (log_actual_size = max(log_actual_size),
             Aspect_ratio  = max(Aspect_ratio),
             Trophic_level = max(Trophic_level),
             Size_group = max(Size_group),
             TempPref_max = max(TempPref_max),
             Depth_max = max(Depth_max))
)


# trait space comprising the 63 spp. -------------------------

# distance matrix (gower)
gower_matrix <- daisy (apply (sel_traits_fish[,-which(colnames(sel_traits_fish) == "scientificName")],2,scale), 
                       metric=("gower"),
                       type = list (ordratio = "Size_group"))

# principal coordinate analysis --------------------------
# Building the functional space based on a PCOA 
pco<-dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
 

##  space (Province trait space) --------------------------
all_pool <- cbind (pco$li[,1:2],
              ext = F,
              sp = sel_traits_fish$scientificName)

# convex hull
a_pool <- all_pool [chull(all_pool[,1:2], y = NULL),] # its convex hull

# use species with some relationships (n=63)
fish_related <- total [which(total$peixe %in% sel_traits_fish$scientificName),]


# random removal -------------------------------------

# create random partite A
ndat <- 1000
random_order <- replicate (ndat,shuffle(length(rownames(m_web))))
m_web_random <- lapply (seq(1,ncol(random_order)), function (i)

  m_web [random_order [,i],] # set the order

)

# calculate functional trait space loss based on RFS (Luza et al. 2022) for each random set --------------------
RFS_corals <- lapply (m_web_random, function (r) # for each random set
  
  lapply (seq (1, ns_coral), function (ncoral) {  # create a sequence of removal
  
  # choose corals to remove
  rm_corals <- rownames (r)[1:ncoral]
  
  # coral-associated fish ( based on regression coefficients )
  coral_associated <- fish_related[which(fish_related$low.coral > 0 & 
                                           fish_related$estimate.turf <= 0 & 
                                           fish_related$coral %in% rm_corals),] 
  
  # reduced space
  setB<-cbind(all_pool, ext1=ifelse(all_pool$sp %in% 
                                 unique(coral_associated$peixe),T,F))
  pk <-setB[which(setB$ext1==F),]
  f <- pk [chull(pk, y = NULL),] # hull
  
  # quantifying reduction in functional space
  chull.poly.complete <- Polygon(a_pool[,1:2], hole=F) # complete space
  chull.area.complete <- chull.poly.complete@area # polygon area
  
  # coral
  chull.poly.coral <- if (nrow (f)==0) {0} else {Polygon(f[,1:2], hole=F)} # reduced space
  chull.area.coral <- if (nrow (f)==0)  {0} else {chull.poly.coral@area}  # area
  
  # calculating reductions across all corals
  RFS<-data.frame (corals=(chull.area.complete-chull.area.coral)/chull.area.complete) # reduction in functional space area
  
  # results into a dataframe
  res <- list (remaining_associated_SR = 1-length(unique(coral_associated$peixe))/ns_fish_associated,
               remaining_SR_all = 1-length(unique(coral_associated$peixe))/sum (dim(m_web_occ)),
                RFS = 1-(RFS),
               coral.associated =  coral_associated$peixe[order(coral_associated$peixe)],
               corals.removed = rm_corals,
               space.coral = f
               )
    ; # return
    res
  }
  
))

# calculate loss for each random data set in PARTITE A, and prepare data for analyses ------------------------

RFS_corals_random <- lapply (RFS_corals, function( RFS_corals) {

    # total loss
    FD_loss_total <- lapply  (RFS_corals, function (i) i$RFS)
    FD_loss_total<- unlist(FD_loss_total)
    
    # loss SR coral associated
    CA_loss_total <- lapply  (RFS_corals, function (i) i$remaining_associated_SR)
    CA_loss_total<- do.call (rbind, CA_loss_total)
    
    # all  SR 
    ALL_loss_total <- lapply  (RFS_corals, function (i) i$remaining_SR_all)
    ALL_loss_total<- do.call (rbind, ALL_loss_total)
    
    # corals
    loss_corals <- lapply  (RFS_corals, function (i) 1-length(i$corals.removed)/ns_coral)
    loss_corals<- do.call (rbind, loss_corals)
    
    # analysis dataset
    # first step of no loss
    analysis_dataset<- (data.frame (remain_all=1,
                                   remain_fish_associated=1,
                                   remain_coral=1,
                                   remain_FD = 1))
    
    # bind removals
    analysis_dataset<-rbind (analysis_dataset,
           data.frame (remain_all=ALL_loss_total,
                remain_fish_associated=CA_loss_total,
                remain_coral=loss_corals,
                remain_FD = FD_loss_total))
    
    analysis_dataset

})


# calculate loss for each random data set in PARTITE B ---------------------------------------

RFS_corals_secondary_extinctions <-  lapply (m_web_random, function (r) 
  
  lapply (seq (1, ns_coral), function (ncoral) { 
  
  # choose corals to remove
  rm_corals <- rownames (r)[1:ncoral]
  
  # coral-associated fish
  coral_associated <- fish_related[which(fish_related$low.coral > 0 & 
                                           fish_related$estimate.turf <= 0 & 
                                           fish_related$coral %in% rm_corals),] 
  
  # associated fish plus those associated with them
  ass_cs <- m_web_occ[unique(coral_associated$peixe),]
  ass_cs <- ass_cs[,colSums(ass_cs)>0]
  all_to_remove <- c(unique(coral_associated$peixe),colnames(ass_cs))
  
  # reduced space
  setB<-cbind(all_pool, ext1=ifelse(all_pool$sp %in% 
                                 all_to_remove,T,F))
  pk <-setB[which(setB$ext1==F),]
  f <- pk [chull(pk, y = NULL),]
  
  # quantifying reduction in functional space
  # https://chitchatr.wordpress.com/2015/01/23/calculating-the-area-of-a-convex-hull/
  chull.poly.complete <- Polygon(a_pool[,1:2], hole=F) # complete space
  chull.area.complete <- chull.poly.complete@area # area
  
  # coral
  chull.poly.coral <- if (nrow (f)==0) {0} else {Polygon(f[,1:2], hole=F)} # reduced space
  chull.area.coral <- if (nrow (f)==0)  {0} else {chull.poly.coral@area}  # area
  
  # calculating reductions across all corals
  RFS<-data.frame (corals=(chull.area.complete-chull.area.coral)/chull.area.complete) # reduction in functional space area
  
  # results
  res <- list (all_to_remove = all_to_remove,
               richness_remaining = length(all_to_remove)/sum(dim(m_web_occ)),
               richness_remaining_associated = length(unique(coral_associated$peixe))/sum(dim(m_web_occ)),
               RFS = 1-RFS,
               coral.associated =  coral_associated$peixe[order(coral_associated$peixe)],
               corals.removed = rm_corals,
               space.coral = f
  )
  ; # return
  res
  }

  )
)

# prepare data for analyses (bind partite A and B) ------------------------------------------

RFS_corals_secondary_extinctions_random <- lapply (seq(1,ndat), function (i) {

      # total loss
    loss_total_secondary <- lapply  (RFS_corals_secondary_extinctions[[i]], function (i) i$RFS)
    loss_total_secondary<- do.call (rbind, loss_total_secondary)
    
    # richness
    loss_total_secondary_SR <- lapply  (RFS_corals_secondary_extinctions[[i]], function (i) i$richness_remaining)
    loss_total_secondary_SR<- do.call (rbind,loss_total_secondary_SR)
    
    # richness associated fish
    loss_total_secondary_SR_CA <- lapply  (RFS_corals_secondary_extinctions[[i]], function (i) i$richness_remaining_associated)
    loss_total_secondary_SR_CA<- do.call (rbind,loss_total_secondary_SR_CA)

    
    # bind indirect extinctions in the analysis dataset
    analysis_dataset<-cbind(RFS_corals_random[[i]], 
                            data.frame  (remain_RFS_secondary = c(1,unlist(loss_total_secondary)),
                                         remain_SR_secondary = c(1,1-unlist (loss_total_secondary_SR)),
                                         remain_SR_secondary_associated = c(1,1-unlist (loss_total_secondary_SR_CA)),
                                         ext.lower=1,
                                         no = seq(0,ns_coral)/ns_coral
                            )
    )
    
    ;
    
    analysis_dataset

})
# analysis_dataset <- RFS_corals_secondary_extinctions_random[[10]]


# -----------------------------------------------
#               Robustness analysis
# -----------------------------------------------

rob_analysis_random <- lapply (RFS_corals_secondary_extinctions_random, function (analysis_dataset) {
      
  
  
  tryCatch({
      
      # run function of hyperbolic curve for each random data set
      hyper_curve_secondary <- fit.hyperbolica (analysis_dataset[,c("no", "ext.lower", "remain_all")],
                                      plot.it = T)
      hyper_curve_associated <- fit.hyperbolica (analysis_dataset[,c("no", "ext.lower", "remain_SR_secondary")],
                                                plot.it = T)
      hyper_curve_associated_RFS <- fit.hyperbolica (analysis_dataset[,c("no", "ext.lower", "remain_FD")],
                                                     plot.it = T)
      hyper_curve_secondary_RFS <- fit.hyperbolica (analysis_dataset[,c("no", "ext.lower", "remain_RFS_secondary")],
                                                plot.it = T)
      
      # predictions
      analysis_dataset$pred_remain_SR_all <- hyper_curve_secondary$preds
      analysis_dataset$pred_remain_SR_secondary <- hyper_curve_associated$preds
      analysis_dataset$pred_remain_RFS_associated <- hyper_curve_associated_RFS$preds
      analysis_dataset$pred_remain_RFS_secondary <- hyper_curve_secondary_RFS$preds
      
      # robustness
      # sr
      SR1 <- (integrate(splinefun(hyper_curve_secondary_RFS$x, 
                            analysis_dataset$pred_remain_SR_all), 
                            min(hyper_curve_secondary_RFS$x), 
                            max(hyper_curve_secondary_RFS$x), 
                            subdivisions = max(100L, 
                                               length(hyper_curve_secondary_RFS$x))))
      
      # SR secondary
      SR2<-(integrate(splinefun(hyper_curve_secondary_RFS$x, 
                          analysis_dataset$pred_remain_SR_secondary), 
                 min(hyper_curve_secondary_RFS$x), 
                 max(hyper_curve_secondary_RFS$x), 
                 subdivisions = max(100L, 
                                               length(hyper_curve_secondary_RFS$x))))
      
      # FD coral associated
      FD1<-(integrate(splinefun(hyper_curve_secondary_RFS$x, 
                          analysis_dataset$pred_remain_RFS_associated), 
                 min(hyper_curve_secondary_RFS$x), 
                 max(hyper_curve_secondary_RFS$x), 
                 subdivisions = max(100L, 
                                                length(hyper_curve_secondary_RFS$x))))
      
      # fd otehr
      FD2<-(integrate(splinefun(hyper_curve_secondary_RFS$x, 
                           analysis_dataset$pred_remain_RFS_secondary), 
                 min(hyper_curve_secondary_RFS$x), 
                 max(hyper_curve_secondary_RFS$x), 
                 subdivisions = max(100L, 
                                    length(hyper_curve_secondary_RFS$x))))
  
      # output with data and robustness estimate
      output <- list(analysis_dataset = analysis_dataset,
                     robustness = cbind (SR1=SR1$value,SR2=SR2$value,FD1=FD1$value,FD2=FD2$value)   )
      output
      
  }, error = function (e) return (NULL))
  
      
})


# some numerical errors appeared
# there is nothing special about the data generating the error, so we removed them
rob_analysis_random <- rob_analysis_random [(unlist(lapply (rob_analysis_random,is.null))!=T)]


# save random run
save (rob_analysis_random, file = here ("output", "random_scenario.RData"))

# averages
apply (
  do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) , 2, 
  mean,na.rm=T)
apply (
  do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) , 2, 
  quantile, 0.025)
apply (
  do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) , 2, 
  quantile, 0.975)


do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))  %>%
  melt() %>%
  ggplot (aes(x=Var2,y=value)) + 
  geom_violin()+
  geom_point (data= do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))  %>%
                melt() %>%
                group_by(Var2) %>%
                summarize(value=mean(value)),
              aes(x=Var2,y=(value)),position = position_dodge())


# clean work space
rm(list=ls())
