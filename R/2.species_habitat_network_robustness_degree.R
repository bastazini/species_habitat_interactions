# ---------------------------------------------# 


#       
#       SPECIES  - HABITAT INTERACTION NETWORKS
# (Analyzes based on "occupancy_data.RData")


#   REMOVAL OF CORALS IS BASED ON THEIR DEGREE - NUMBER OF CORAL-ASSOCIATED FISH


#   IN DESCENDING ORDER (HIGH DEGREE IS REMOVED FIRST)


#   The other scripts (also numbered with '2' bring scenarios of random- and vulnerability-based removal)



# ----------------------------------------------# 

# load basic data & functions
rm(list=ls())
require(here)
source(here ("R", "functions.R"))
source(here ("R", "packages.R"))
load(file = here ("processed_data","occupancy_data.RData"))

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

# ordering (degree) --------------------
m_web<-sortweb(m_web,sort.order="dec")


# functional trait space loss --------------------
# load and organize trait data ----------

source ("R/functions.R")
trait_dataset<- organize_traits (fish_size)
#unique(network_data$peixe) %in% unique(trait_dataset$scientificName)


#---------------------------------------
# 		LOSS OF TRAIT SPACE (RFS)
# ------------------------------------

# find coral associated fish
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
m_web_occ[m_web_occ<0.8] <- 0 # conservativo suficiente para definir associacao
# rm empty cols
m_web_occ <- m_web_occ [, which(colSums(m_web_occ)>0)]


# plot and save tripartite network ---------------------------------

pdf(file = here ("output","plot.pdf"),width=10,height=10)

plotweb2(data.matrix(m_web),
         data.matrix(m_web_occ),
         method = "normal",
         empty=T,
         col.interaction = "gray",
         ybig=1,
         labsize = 0.75,
         spacing=0.01,
         lab.space	=0.1,
         method2="normal",
         spacing2=0.01,
         empty2=T,
         col.interaction2 = "orange" ,
         col.pred2 = "orange",
         col.prey2 = "gray80"
)

dev.off()


# ----------------------------------

# select traits for spp with some association to corals and coral -associated fish
# select traits to analysis
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


# trait space comprising the 63 spp.
# distance matrix (gower)
gower_matrix <- daisy (apply (sel_traits_fish[,-which(colnames(sel_traits_fish) == "scientificName")],2,scale), 
                       metric=("gower"),
                       type = list (ordratio = "Size_group"))


# principal coordinate analysis
# Building the functional space based on a PCOA 
pco<-dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 


##  space
all_pool <- cbind (pco$li[,1:2],
              ext = F,
              sp = sel_traits_fish$scientificName)


# convex hull
a_pool <- all_pool [chull(all_pool[,1:2], y = NULL),] # its convex hull

# use species with some relationships (n=63)
fish_related <- total [which(total$peixe %in% sel_traits_fish$scientificName),]


# --------------------------------------
# removal based on coral degree


# calculate functional trait space loss
# based on RFS (Luza et al. 2022)
RFS_corals <- lapply (seq (1, ns_coral), function (ncoral) { 
  
  # choose corals to remove
  rm_corals <- rownames (m_web)[1:ncoral]
  
  # coral-associated fish
  coral_associated <- fish_related[which(fish_related$low.coral > 0 & 
                                         fish_related$estimate.turf <= 0 & 
                                         fish_related$coral %in% rm_corals
                                         ),] 
  
  # reduced space
  setB<-cbind(all_pool, ext1=ifelse(all_pool$sp %in% 
                                 unique(coral_associated$peixe),T,F))
  pk <-setB[which(setB$ext1==F),]
  f <- pk [chull(pk, y = NULL),] # hull
  
  # quantifying reduction in functional space
  # https://chitchatr.wordpress.com/2015/01/23/calculating-the-area-of-a-convex-hull/
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
})

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


# simulating secondary (indirect) extinctions -----------------

RFS_corals_secondary_extinctions <- lapply (seq (1, ns_coral), function (ncoral) { 
  
  # choose corals to remove
  rm_corals <- rownames (m_web)[1:ncoral]
  
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
})

# total loss
loss_total_secondary <- lapply  (RFS_corals_secondary_extinctions, function (i) i$RFS)
loss_total_secondary<- do.call (rbind, loss_total_secondary)

# richness
loss_total_secondary_SR <- lapply  (RFS_corals_secondary_extinctions, function (i) i$richness_remaining)
loss_total_secondary_SR<- do.call (rbind,loss_total_secondary_SR)

# richness associated fish
loss_total_secondary_SR_CA <- lapply  (RFS_corals_secondary_extinctions, function (i) i$richness_remaining_associated)
loss_total_secondary_SR_CA<- do.call (rbind,loss_total_secondary_SR_CA)

# bind indirect extinctions in the analysis dataset
analysis_dataset<-cbind(analysis_dataset, 
                        data.frame  (remain_RFS_secondary = c(1,unlist(loss_total_secondary)),
                                    remain_SR_secondary = c(1,1-unlist (loss_total_secondary_SR)),
                                    remain_SR_secondary_associated = c(1,1-unlist (loss_total_secondary_SR_CA)),
                                    ext.lower=1,
                                    no = seq(0,ns_coral)/ns_coral
                   )
        )


# run function of hyperbolic curve
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
(SR1 <- (integrate(splinefun(hyper_curve_secondary_RFS$x, 
                      analysis_dataset$pred_remain_SR_all), 
                      min(hyper_curve_secondary_RFS$x), 
                      max(hyper_curve_secondary_RFS$x), 
                      subdivisions = max(100L, 
                                         length(hyper_curve_secondary_RFS$x)))))
# SR secondary
(SR2 <- (integrate(splinefun(hyper_curve_secondary_RFS$x, 
                    analysis_dataset$pred_remain_SR_secondary), 
           min(hyper_curve_secondary_RFS$x), 
           max(hyper_curve_secondary_RFS$x), 
           subdivisions = max(100L, 
                                         length(hyper_curve_secondary_RFS$x)))))

# FD cora lassociated
(FD1 <- (integrate(splinefun(hyper_curve_secondary_RFS$x, 
                    analysis_dataset$pred_remain_RFS_associated), 
           min(hyper_curve_secondary_RFS$x), 
           max(hyper_curve_secondary_RFS$x), 
           subdivisions = max(100L, 
                                          length(hyper_curve_secondary_RFS$x)))))
# fd otehr
(FD2 <- (integrate(splinefun(hyper_curve_secondary_RFS$x, 
                     analysis_dataset$pred_remain_RFS_secondary), 
           min(hyper_curve_secondary_RFS$x), 
           max(hyper_curve_secondary_RFS$x), 
           subdivisions = max(100L, 
                              length(hyper_curve_secondary_RFS$x)))))

# analysis data frame
rob_analysis_degree <- data.frame (SR1=SR1$value,
                                   SR2=SR2$value,
                                   FD1=FD1$value,
                                   FD2=FD2$value) %>%
  melt() %>%
  mutate("scenario"="degree")


# -----------------------------

# save degree results


# save random run
save (analysis_dataset,rob_analysis_degree, file = here ("output", "degree_results.RData"))



# plot and save the robustness curve -------------------------------------------

png (here ("output", "ATC.png"),width = 20, height = 20,units="cm",res=300)

# plot
ggplot(analysis_dataset[order(analysis_dataset$remain_coral,decreasing=F),], 
       aes (x= (1-remain_coral), y = remain_fish_associated)) +
  
  # functional diversity
  # direct loss
  geom_point( aes(y=remain_FD),linetype=2,size=3,col="orange",fill="orange",shape=22,alpha=1) + 
  geom_line( aes(y=pred_remain_RFS_associated),linetype=1,size=1.2,col="orange") + 
  geom_ribbon( aes(ymax=pred_remain_RFS_associated),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.3,
               fill="orange") + 
  
  
  
  # indirect loss
  geom_point( aes(y=remain_RFS_secondary),linetype=2,size=3, col = "orange",shape=19,alpha=1) + 
  geom_line( aes(y=pred_remain_RFS_secondary),linetype=2,size=1.2, col = "orange") + 
  geom_ribbon( aes(ymax=pred_remain_RFS_secondary),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.4,
               fill="orange") + 
  
  
  # taxonomic diversity
  geom_ribbon( aes(ymax=pred_remain_SR_all),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.3,
               fill="black") + 
  geom_point( aes(y=remain_all),linetype=1,size=3,col="black",fill="black",shape=22,alpha=1) + 
  geom_line( aes(y=pred_remain_SR_all),linetype=1,size=1.2, col= "black") + 
  
  #geom_line( aes(y=remain_fish_associated*100),linetype=1,size=1.2,col="orange") +
  geom_point( aes(y=remain_SR_secondary),linetype=1,size=3,col="black",shape=19,alpha=1) +
  geom_line( aes(y=pred_remain_SR_secondary),linetype=2,size=1.2,col="black") +
  geom_ribbon( aes(ymax=pred_remain_SR_secondary),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.4,
               fill="black") + 
  
  
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Remaining reef fish species",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~.*1, name="Remaining functional diversity (trait space area)")
  ) + 
  xlab ("Proportion of corals eliminated") + 
  theme_bw() + 
  ggtitle ("")+
  theme (axis.title=element_text(size=14))


dev.off()

# create dir to receive files
save.image(here ("processed_data", "basic_data.RData"))

# clean work space
rm(list=ls())
