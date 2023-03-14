# ---------------------------------------------# 

#       SPECIES  - HABITAT INTERACTION NETWORKS

# ----------------------------------------------# 
# load basic data & functions
require(here)
source(here ("R", "functions.R"))
source(here ("R", "packages.R"))
load(file = here ("data","occupancy_data.RData"))

# coun total N fish
ns_fish <- length(unique(network_data$peixe))
# total N corals
ns_coral <- length(unique(network_data$coral))

# histogram of predictions
ggplot (network_data, aes (x=pred)) +
  geom_histogram(binwidth=0.1,fill="yellow", colour="black")+
  facet_wrap(~coral)

# coral associated fish (following Luza et al. 2022, sci rep)
network_data <- network_data[which(network_data$low.coral>0 & 
                                     network_data$estimate.turf<=0),]
# N coral associated fish
ns_fish_associated <- length(unique(network_data$peixe)) # ncol(m_web)

# transform data into matrix (matrix C x Fca)
m_web <- cast (formula = coral~peixe,value = "pred",
               fun.aggregate = mean,na.rm=T,fill=0,
               data = network_data)
labels_web <- m_web[,1]; m_web<- m_web[,-1] # labels
# adjust
rownames(m_web)<- labels_web
# ordering
m_web<-sortweb(m_web,sort.order="dec")

# secondary extinctions following bipartite
(ex <- second.extinct(m_web, participant="lower", method="random", nrep=10, 
                      details=F))
second.extinct(m_web, participant="lower", method="random", nrep=50, 
                      details=T)

bipartite::robustness(ex)
fit.hyperbolica(ex)



# ===============================================
# functional trait space loss
# organize trait data

source ("R/functions.R")
trait_dataset<- organize_traits (fish_size)

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
# the correlation between traits
(cor(sel_traits[,-which(colnames(sel_traits) == "scientificName")], 
     use = "complete.obs"))#  correlation is fine
# distance matrix (gower)
gower_matrix <- daisy (apply (sel_traits[,-which(colnames(sel_traits) == "scientificName")],2,scale), 
                       metric=("gower"),
                       type = list (ordratio = "Size_group"))
# principal coordinate analysis
# Building the functional space based on a PCOA 
pco<-dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 

#barplot(pco$eig) # barplot of eigenvalues for each axis 
(Inertia2<-(pco$eig[1]+pco$eig[2]+pco$eig[3]) /(sum(pco$eig))) # percentage of inertia explained by the two first axes
# estimate quality of f space
#quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
#                                       plot="quality_funct_space_I") # it will produce a plot (hosted in the root folder)
## only the frst axis
(Inertia.first <- (pco$eig[1]) /(sum(pco$eig)))
## only the frst axis
(Inertia.scnd <- (pco$eig[2]) /(sum(pco$eig)))
## only the frst axis
(Inertia.trd <- (pco$eig[3]) /(sum(pco$eig)))
Inertia.first+Inertia.scnd

## complete space
all <- cbind (pco$li[,1:2],
              ext = F,
              sp = sel_traits$scientificName)
a <- all [chull(all[,1:2], y = NULL),] # its convex hull

# match ordination and traits
total<- cbind (total,
               all[(match (total$peixe, all$sp)),c("A1", "A2")])
# coral associated space
c_assoc <-cbind(all, ext1=ifelse(all$sp %in% 
                                   unique( (sp_analyzed_response$peixe)),T,F))
c_assoc <-c_assoc[which(c_assoc$ext1==T),]
c_assoc_set <- c_assoc [chull(c_assoc, y = NULL),]

# functional space loss
RFS_corals <- lapply (seq (1, ns_coral), function (ncoral) { 
  
  # choose corals to remove
  rm_corals <- rownames (m_web)[1:ncoral]
  
  # coral-associated fish
  coral_associated <- total[which(total$low.coral > 0 & 
                                  total$estimate.turf <= 0 & 
                                  total$coral %in% rm_corals),] 
  # reduced space
  setB<-cbind(all, ext1=ifelse(all$sp %in% 
                                 unique(coral_associated$peixe),T,F))
  pk <-setB[which(setB$ext1==F),]
  f <- pk [chull(pk, y = NULL),] # hull
  
  # quantifying reduction in functional space
  # https://chitchatr.wordpress.com/2015/01/23/calculating-the-area-of-a-convex-hull/
  chull.poly.complete <- Polygon(a[,1:2], hole=F) # complete space
  chull.area.complete <- chull.poly.complete@area # polygon area
  
  # coral
  chull.poly.coral <- Polygon(f[,1:2], hole=F) # reduced space
  chull.area.coral <- chull.poly.coral@area # polygon area
  
  # calculating reductions across all corals
  RFS<-data.frame (corals=(chull.area.complete-chull.area.coral)/chull.area.complete) # reduction in functional space area
  
  # results into a dataframe
  res <- list (remaining_associated_SR = 1-length(unique(coral_associated$peixe))/ns_fish_associated,
               remaining_SR_all = 1-length(unique(coral_associated$peixe))/ns_fish,
                RFS = (RFS),
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
            remain_FD = 1-FD_loss_total))

# --------------------------------------------------------
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
# cora associated fish in the row
m_web_occ <- m_web_occ[which(rownames (m_web_occ) %in% colnames (m_web)),
                        which(colnames (m_web_occ) %in% colnames (m_web) == F)]
# ordering
m_web_occ<- m_web_occ [match (colnames(m_web),
                              rownames(m_web_occ)),]
#rownames(m_web_occ) == colnames(m_web)
# correlations lower than 0.8 <-- 0
m_web_occ[m_web_occ<0.8] <- 0
# rm empty cols
m_web_occ <- m_web_occ [, which(colSums(m_web_occ)>0)]


# plot tripartite network

pdf(file = here ("output","plot.pdf"),width=10,height=10)

plotweb2(data.matrix(m_web),
  data.matrix(m_web_occ),
  method = "normal",
  empty=T,
  labsize = 0.75,
  spacing=0.01,
  lab.space	=0.75,
  method2="normal",
  spacing2=0.01,
  empty2=T
)

dev.off()


# simulating indirect extinctions
RFS_corals_secondary_extinctions <- lapply (seq (1, ns_coral), function (ncoral) { 
  
  # choose corals to remove
  rm_corals <- rownames (m_web)[1:ncoral]
  # coral-associated fish
  coral_associated <- total[which(total$low.coral > 0 & 
                                    total$estimate.turf <= 0 & 
                                    total$coral %in% rm_corals),] 
  
  # associated fish plus those associated with them
  ass_cs <- m_web_occ[unique(coral_associated$peixe),]
  ass_cs <- ass_cs[,colSums(ass_cs)>0]
  all_to_remove <- c(unique(coral_associated$peixe),colnames(ass_cs))
  
  # reduced space
  setB<-cbind(all, ext1=ifelse(all$sp %in% 
                                 all_to_remove,T,F))
  pk <-setB[which(setB$ext1==F),]
  f <- pk [chull(pk, y = NULL),]
  
  # quantifying reduction in functional space
  # https://chitchatr.wordpress.com/2015/01/23/calculating-the-area-of-a-convex-hull/
  chull.poly.complete <- Polygon(a[,1:2], hole=F) # complete space
  chull.area.complete <- chull.poly.complete@area # area
  
  # coral
  chull.poly.coral <- Polygon(f[,1:2], hole=F) # reduced space
  chull.area.coral <- chull.poly.coral@area # area
  
  # calculating reductions across all corals
  RFS<-data.frame (corals=(chull.area.complete-chull.area.coral)/chull.area.complete) # reduction in functional space area
  
  # results
  res <- list (all_to_remove = all_to_remove,
               richness_remaining = length(all_to_remove)/ns_fish,
               richness_remaining_associated = length(unique(coral_associated$peixe))/ns_fish,
               RFS = RFS,
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
                        data.frame  (remain_RFS_secondary = c(1,1-unlist(loss_total_secondary)),
                                    remain_SR_secondary = c(1,1-unlist (loss_total_secondary_SR)),
                                    remain_SR_secondary_associated = c(1,1-unlist (loss_total_secondary_SR_CA)),
                                    ext.lower=1,
                                    no = seq(1,nrow (analysis_dataset))/nrow (analysis_dataset)
                   ))

# run function of hyperbolic curve
hyper_curve_secondary <- fit.hyperbolica (analysis_dataset[,c("no", "ext.lower", "remain_all")],
                                plot.it = F)
hyper_curve_associated <- fit.hyperbolica (analysis_dataset[,c("no", "ext.lower", "remain_SR_secondary")],
                                          plot.it = F)
hyper_curve_associated_RFS <- fit.hyperbolica (analysis_dataset[,c("no", "ext.lower", "remain_FD")],
                                               plot.it = F)
hyper_curve_secondary_RFS <- fit.hyperbolica (analysis_dataset[,c("no", "ext.lower", "remain_RFS_secondary")],
                                          plot.it = F)
# predictions
analysis_dataset$pred_remain_SR_all <- predict(hyper_curve_associated$model, 
                                               newdata = data.frame(x = seq(0, 
                                                                            1, 
                                                                            length=9)))
analysis_dataset$pred_remain_SR_secondary <- predict(hyper_curve_secondary$model, 
                                                     newdata = data.frame(x = seq(0, 
                                                                                  1, length=9)))
analysis_dataset$pred_remain_RFS_associated <- predict(hyper_curve_associated_RFS$model, 
                                                       newdata = data.frame(x = seq(0, 
                                                                                    1, length=9)))
analysis_dataset$pred_remain_RFS_secondary <- predict(hyper_curve_secondary_RFS$model, 
                                                      newdata = data.frame(x = seq(0, 
                                                                                  1, length=9)))

# robustness
# sr
(integrate(splinefun(hyper_curve_secondary_RFS$x, 
                      analysis_dataset$pred_remain_SR_all), 
                      min(hyper_curve_secondary_RFS$x), 
                      max(hyper_curve_secondary_RFS$x), 
                      subdivisions = max(100L, 
                                         length(hyper_curve_secondary_RFS$x))))
# SR secondary
(integrate(splinefun(hyper_curve_secondary_RFS$x, 
                                analysis_dataset$pred_remain_SR_secondary), 
           min(hyper_curve_secondary_RFS$x), 
           max(hyper_curve_secondary_RFS$x), 
           subdivisions = max(100L, 
                                         length(hyper_curve_secondary_RFS$x))))
# FD cora lassociated
(integrate(splinefun(hyper_curve_secondary_RFS$x, 
                                 analysis_dataset$pred_remain_RFS_associated), 
           min(hyper_curve_secondary_RFS$x), 
           max(hyper_curve_secondary_RFS$x), 
           subdivisions = max(100L, 
                                          length(hyper_curve_secondary_RFS$x))))
# fd otehr
(integrate(splinefun(hyper_curve_secondary_RFS$x, 
                     analysis_dataset$pred_remain_RFS_secondary), 
           min(hyper_curve_secondary_RFS$x), 
           max(hyper_curve_secondary_RFS$x), 
           subdivisions = max(100L, 
                              length(hyper_curve_secondary_RFS$x))))



# project the trait space loss
# plot
png (here ("output", "ATC.png"))
ggplot(analysis_dataset[order(analysis_dataset$remain_coral,decreasing=F),], 
       aes (x= (1-remain_coral), y = remain_fish_associated)) +
  
  #geom_line( aes(y=remain_fish)) + 
  geom_point( aes(y=remain_all),linetype=1,size=1.2) + 
  geom_line( aes(y=pred_remain_SR_all),linetype=1,size=1.2, col= "gray") + 
  
  #geom_line( aes(y=remain_fish_associated*100),linetype=1,size=1.2,col="orange") +
  geom_point( aes(y=remain_SR_secondary),linetype=1,size=1.2,col="red") +
  geom_line( aes(y=pred_remain_SR_secondary),linetype=1,size=1.2,col="pink") +
  
  # trait space
  geom_point( aes(y=remain_FD),linetype=2,size=1.2,col="black") + 
  geom_line( aes(y=pred_remain_RFS_associated),linetype=2,size=1.2,col="gray") + 
  
  geom_point( aes(y=remain_RFS_secondary),linetype=2,size=1.2, col = "red") + 
  geom_line( aes(y=pred_remain_RFS_secondary),linetype=2,size=1.2, col = "pink") + 
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Remaining reef fish species",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~.*1, name="Remaining functional diversity (trait space area)")
  ) + 
  xlab ("Proportion of corals eliminated") + 
  theme_bw() + 
  ggtitle ("Direct (black) and indirect (red) loss\nof fish SR (solid line) and FD (dashed line)")+
  theme (axis.title=element_text(size=14))

dev.off()



# --------------------------------------------------------------------------------------
# trait space illustrations
# space of indirect extinctions
all_to_remove <- lapply  (RFS_corals_secondary_extinctions, function (i) 
  i$all_to_remove)

c_direct_indirect <-cbind(all, ext1=ifelse(all$sp %in% 
                                             unique(unlist(all_to_remove)),T,F))
c_direct_indirect <-c_direct_indirect[which(c_direct_indirect$ext1==T),]
c_direct_indirect_set <- c_direct_indirect [chull(c_direct_indirect, y = NULL),]


# ----------------------- complete space
# kernel densities
require(ks)
# optimal bandwidth estimation
hpi_mi_d1 <- Hpi(x = all[,c("A1","A2")])

# kernel density estimation
est_mi_d1 <- kde(x = all[,c("A1","A2")], 
                 H = hpi_mi_d1, 
                 compute.cont = TRUE)  

# bandwidths for each point
den_mi_d1 <- list(est_mi_d1$eval.points[[1]], 
                  est_mi_d1$eval.points[[2]], 
                  est_mi_d1$estimate)
names(den_mi_d1) <- c("x", "y", "z")
dimnames(den_mi_d1$z) <- list(den_mi_d1$x, den_mi_d1$y)
dcc_mi_d1 <- melt(den_mi_d1$z)

# run kernel
cl_50_mi_d1 <- cl(df = den_mi_d1, prob = 0.50) # 0.5 probability kernel
cl_95_mi_d1 <- cl(df = den_mi_d1, prob = 0.95) # 0.95 probability kernel
cl_99_mi_d1 <- cl(df = den_mi_d1, prob = 0.99)# 0.99 probability kernel

## PCoA
# colour palette
col_pal <- colorRampPalette(c("red", "yellow", "white"))(100)

# plot PCoA
PCoA_plot_mi_d1 <- ggplot(dcc_mi_d1, aes(x = Var1, y = Var2)) + 
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal), 
                       limits = c(0,20)) +
  # points for species
  geom_point(data = all, 
             aes(x = A1,y=A2), size = 0.3, alpha = 0.5, colour = "grey20") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_mi_d1, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_mi_d1, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_mi_d1, colour = "grey70", size = 1) +
  #scale_fill_gradient(low = 'white', high = 'red') +
  coord_equal() +
  theme_classic()+
  xlim (c(-0.6,0.6))+
  ylim (c(-0.6,0.6))

# plot density
PCoA_plot_mi_d1 <- PCoA_plot_mi_d1 +geom_polygon(data=a, aes (A1,A2),
                              alpha=0.05,
                              fill="black",
                              colour = "gray92",
                              size=1,
                              linetype = 1) + # complete space
  geom_text_repel(data = a, aes (x=A1, y=A2, label=firstup(sp)),
                  size=3)


# ---------------- coral associated
# optimal bandwidth estimation for coral associated fish
hpi_mi_d_CA <- Hpi(x = all[-which(all$sp %in% unique( (sp_analyzed_response$peixe))),
                           c("A1","A2")])

# kernel density estimation
est_mi_CA <- kde(x = all[-which(all$sp %in% unique( (sp_analyzed_response$peixe))),
                         c("A1","A2")], 
                 H = hpi_mi_d_CA, 
                 compute.cont = TRUE)  

# bandwidths for each point
den_mi_CA <- list(est_mi_CA$eval.points[[1]], 
                  est_mi_CA$eval.points[[2]], 
                  est_mi_CA$estimate)
names(den_mi_CA) <- c("x", "y", "z")
dimnames(den_mi_CA$z) <- list(den_mi_CA$x, 
                              den_mi_CA$y)
dcc_mi_CA <- melt(den_mi_CA$z)

# run kernel
cl_50_mi_CA <- cl(df = den_mi_CA, prob = 0.50)# 0.5 probability kernel
cl_95_mi_CA <- cl(df = den_mi_CA, prob = 0.95)# 0.95 probability kernel
cl_99_mi_CA <- cl(df = den_mi_CA, prob = 0.99)# 0.99 probability kernel

## PCoA
# plot PCoA
PCoA_plot_mi_d2 <- ggplot(dcc_mi_CA, aes(x = Var1, y = Var2)) + 
  
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal), 
                       limits = c(0,20)) +
  
  # points for species
  geom_point(data = all[-which(all$sp %in% unique( (sp_analyzed_response$peixe))),], 
             aes(x = A1,y=A2), size = 0.3, alpha = 0.5, colour = "grey20") +
  
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_mi_d1, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_mi_d1, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_mi_d1, colour = "grey70", size = 1) +
  #scale_fill_gradient(low = 'white', high = 'red') +
  coord_equal() +
  theme_classic()+
  xlim (c(-0.6,0.6))+
  ylim (c(-0.6,0.6))

# plot density
PCoA_plot_mi_d2 <- PCoA_plot_mi_d2 +
  geom_polygon(data=c_assoc_set, aes (A1,A2),
               alpha=0.05,
               fill="black",
               colour = "gray92",
               size=1,
               linetype = 1)+
  geom_text_repel(data = c_assoc_set, aes (x=A1, y=A2, label=firstup(sp)),
                  size=3)
  
# ---------------------- elimination of all 
# optimal bandwidth estimation for coral associated fish
hpi_mi_d_all <- Hpi(x = all[-which(all$sp %in% unique(unlist(all_to_remove))),
                           c("A1","A2")])

# kernel density estimation
est_mi_all <- kde(x = all[-which(all$sp %in% unique(unlist(all_to_remove))),c("A1","A2")], 
                 H = hpi_mi_d_all, 
                 compute.cont = TRUE)  

# bandwidths for each point
den_mi_all <- list(est_mi_all$eval.points[[1]], 
                  est_mi_all$eval.points[[2]], 
                  est_mi_all$estimate)
names(den_mi_all) <- c("x", "y", "z")
dimnames(den_mi_all$z) <- list(den_mi_all$x, 
                              den_mi_all$y)
dcc_mi_all <- melt(den_mi_all$z)

# run kernel
cl_50_mi_all <- cl(df = den_mi_all, prob = 0.50)# 0.5 probability kernel
cl_95_mi_all <- cl(df = den_mi_all, prob = 0.95)# 0.95 probability kernel
cl_99_mi_all <- cl(df = den_mi_all, prob = 0.99)# 0.99 probability kernel

## PCoA
# plot PCoA
PCoA_plot_mi_d3 <- ggplot(dcc_mi_all, aes(x = Var1, y = Var2)) + 
  
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal), 
                       limits = c(0,20)) +
  
  # points for species
  geom_point(data = all[-which(all$sp %in% unique(network_data$peixe)),], 
             aes(x = A1,y=A2), size = 0.3, alpha = 0.5, colour = "grey20") +
  
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_mi_d1, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_mi_d1, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_mi_d1, colour = "grey70", size = 1) +
  #scale_fill_gradient(low = 'white', high = 'red') +
  coord_equal() +
  theme_classic()+
  xlim (c(-0.6,0.6))+
  ylim (c(-0.6,0.6))

# plot density
PCoA_plot_mi_d3 <- PCoA_plot_mi_d3 +
  geom_polygon(data=c_direct_indirect_set, aes (A1,A2),
               alpha=0.05,
               fill="black",
               colour = "gray92",
               size=1,
               linetype = 1)+
  geom_text_repel(data = c_direct_indirect_set, aes (x=A1, y=A2, label=firstup(sp)),
                  size=3)



# ----------------- plot of the difference
# difference between complete and coral associated
diff_kernels <- den_mi_CA$z - den_mi_d1$z # difference
diff_kernels_coord <- den_mi_d1
diff_kernels_coord$z <- diff_kernels

# melt to plot
dcc_mi_d1_test2 <- melt(diff_kernels_coord$z)
dcc_mi_d1_test2$value[dcc_mi_d1_test2$value>0] <-0 # only negative differences

# plot PCoA
PCoA_plot_mi_diff <- ggplot(dcc_mi_d1_test2, aes(x = Var1, y = Var2)) + 
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_colour_gradient2(
    low = ("#F46B2F"),
    mid = "white",
    high = ("#4268cb"),
    midpoint = 0,
    aesthetics = "fill",
    limits=c(-12,0)
  ) +
  
  # points for species
  geom_point(data = all, 
             aes(x = A1,y=A2),
             size = 0.3, alpha = 0.5, colour = "grey20") +
  coord_equal() +
  theme_classic() +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black",size =5),
        axis.title = element_text(colour = "black",size =7),
        legend.position = "top"
  ) +
  xlim (c(-0.6,0.6))+
  ylim (c(-0.6,0.6))
PCoA_plot_mi_diff

## correlations to project trait values into the ordination
correlations <- cor (data.frame(sel_traits[,-which(colnames(sel_traits) == "scientificName")],
                                pco$li[,1:3]),
                     use = "complete.obs")
correlations<-correlations [,c("A1","A2")]# interesting correlations

# clean plot  to receive correlations
clean_plot <- ggplot ()+ geom_point(data = all, 
                    aes(x = A1,y=A2),
                    size = 0.3, alpha = 0.5, colour = "grey20") +
  coord_equal() +
  theme_classic() +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black",size =5),
        axis.title = element_text(colour = "black",size =7),
        legend.position = "top"
  )+
  xlim (c(-0.6,0.6))+
  ylim (c(-0.6,0.6))

clean_plot <- clean_plot + geom_segment(aes(x = 0, y = 0, 
                 xend = correlations[1,1]*0.2, 
                 yend = correlations[1,2]*0.2),size = 1,
             color="black",
             arrow = arrow(length = unit(.35, "cm")))  + 
  ## annotate
  annotate(geom="text",x=correlations[1,1]*0.25,
           y=correlations[1,2]*0.24,label="Aspect ratio",
           color="black") +
  
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[2,1]*0.2, 
                   yend = correlations[2,2]*0.2),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm"))) + 
  annotate(geom="text",x=correlations[2,1]*0.25,
           y=correlations[2,2]*0.25,label="Trophic level",
           color="black") +
  
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[3,1]*0.2, 
                   yend = correlations[3,2]*0.2),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm"))) + 
  annotate(geom="text",x=correlations[3,1]*0.20,
           y=correlations[3,2]*0.25,label="Group size",
           color="black") +
  
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[4,1]*0.2, 
                   yend = correlations[4,2]*0.2),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm"))) + 
  annotate(geom="text",x=correlations[4,1]*0.25,
           y=correlations[4,2]*0.29,label="TºC max",
           color="black") + 
  
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[5,1]*0.2, 
                   yend = correlations[5,2]*0.2),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm"))) + 
  annotate(geom="text",x=correlations[5,1]*0.25,
           y=correlations[5,2]*0.23,label="Depth max",
           color="black") + 
  
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[6,1]*0.2, 
                   yend = correlations[6,2]*0.2),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm"))) + 
  annotate(geom="text",x=correlations[6,1]*0.18,
           y=correlations[6,2]*0.25,label="Body size",
           color="black") 


# -------------------- difference complete vs removal of all 

diff_kernels <- den_mi_all$z - den_mi_d1$z # difference
diff_kernels_coord <- den_mi_d1
diff_kernels_coord$z <- diff_kernels

# melt to plot
dcc_mi_d1_test2 <- melt(diff_kernels_coord$z)
dcc_mi_d1_test2$value[dcc_mi_d1_test2$value>0] <-0 # negative diff

# plot PCoA
PCoA_plot_mi_diff_2 <- ggplot(dcc_mi_d1_test2, aes(x = Var1, y = Var2)) + 
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_colour_gradient2(
    low = ("#F46B2F"),
    mid = "white",
    high = ("#4268cb"),
    midpoint = 0,
    aesthetics = "fill",
    limits=c(-12,0)
  ) +
  
  # points for species
  geom_point(data = all, 
             aes(x = A1,y=A2),
             size = 0.3, alpha = 0.5, colour = "grey20") +
  coord_equal() +
  theme_classic() +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black",size =5),
        axis.title = element_text(colour = "black",size =7),
        legend.position = "top"
  )+
  xlim (c(-0.6,0.6))+
  ylim (c(-0.6,0.6))

pdf (here ("output", "trait_space.pdf "),width=10,height=8)
grid.arrange(PCoA_plot_mi_d1+theme(legend.position = "top"),
             PCoA_plot_mi_d2+theme(legend.position = "top"),
             PCoA_plot_mi_d3+theme(legend.position = "top"),
             PCoA_plot_mi_diff,
             PCoA_plot_mi_diff_2,
            clean_plot, nrow=2,ncol=3)

dev.off()

# trait space occupancy
dcc_mi_d1$value
round (sum(dcc_mi_all$value>0)/sum(dcc_mi_d1$value>0),2)
round(sum(dcc_mi_CA$value>0)/sum(dcc_mi_d1$value>0),2)
