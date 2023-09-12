
# --------------------------------------------------------------------------------------


#             trait space representation/analyses


# --------------------------------------------------------------------------------------

load(here ("processed_data", "basic_data.RData"))

# create the complete trait space (all 113 species in the dataset)
# distance matrix (gower)
gower_matrix_complete <- daisy (apply (sel_traits[,-which(colnames(sel_traits) == "scientificName")],2,scale), 
                                metric=("gower"),
                                type = list (ordratio = "Size_group"))

# principal coordinate analysis
# Building the functional space based on a PCOA 
pco_complete<-dudi.pco(quasieuclid(gower_matrix_complete), 
                       scannf=F, 
                       nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 

#barplot(pco$eig) # barplot of eigenvalues for each axis 
(Inertia2<-(pco_complete$eig[1]+pco_complete$eig[2]+pco_complete$eig[3]) /(sum(pco_complete$eig))) # percentage of inertia explained by the two first axes

# estimate quality of f space
#quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
#                                       plot="quality_funct_space_I") # it will produce a plot (hosted in the root folder)
## only the frst axis
(Inertia.first <- (pco_complete$eig[1]) /(sum(pco_complete$eig)))
## only the frst axis
(Inertia.scnd <- (pco_complete$eig[2]) /(sum(pco_complete$eig)))
## only the frst axis
(Inertia.trd <- (pco_complete$eig[3]) /(sum(pco_complete$eig)))
Inertia.first+Inertia.scnd

## complete space
all_complete <- cbind (pco_complete$li[,1:2],
                       ext = F,
                       sp = sel_traits$scientificName)

# convex hull
a_complete <- all_complete [chull(all_complete[,1:2], y = NULL),] # its convex hull

# coral associated space
c_assoc <-cbind(all_complete, ext1=ifelse(all_complete$sp %in% 
                                            rownames(m_web_occ),T,F))
c_assoc <-c_assoc[which(c_assoc$ext1==T),] # the vector (choose species in the sets of coral associated)
c_assoc_set <- c_assoc [chull(c_assoc, y = NULL),] # the convex hull vertices

# coccurring fish
c_direct_indirect <-cbind(all_complete, ext1=ifelse(all_complete$sp %in% 
                                                      unlist(dimnames(m_web_occ)),T,F))
c_direct_indirect <-c_direct_indirect[which(c_direct_indirect$ext1==T),]# the vector (choose species in the sets of coccurrent spp)
c_direct_indirect_set <- c_direct_indirect [chull(c_direct_indirect, y = NULL),]

# not associated
not_c_assoc <-cbind(all_complete, ext1=ifelse(all_complete$sp %in% 
                                                unlist(dimnames(m_web_occ)) ==F,T,F))
not_c_assoc <-not_c_assoc[which(not_c_assoc$ext1==T),] # the vector (choose species in the sets of coral associated)
not_c_assoc_set <- not_c_assoc [chull(not_c_assoc, y = NULL),] # the convex hull vertices

# out not associated
out_not_c_assoc <-cbind(all_complete, ext1=ifelse(all_complete$sp %in% 
                                                unlist(dimnames(m_web_occ)) ==F,T,F))
out_not_c_assoc <-out_not_c_assoc[which(out_not_c_assoc$ext1==F),] # the vector (choose species in the sets of coral associated)
out_not_c_assoc_set <- out_not_c_assoc [chull(out_not_c_assoc, y = NULL),] # the convex hull vertices

# representation of each trait space
# coral associated space / complete space
Polygon(c_assoc_set[,1:2], hole=F)@area / Polygon(a_complete[,1:2], hole=F)@area 
# coocurring / complete
Polygon(c_direct_indirect_set[,1:2], hole=F)@area  / Polygon(a_complete[,1:2], hole=F)@area 
# not associated / complete
Polygon(not_c_assoc_set[,1:2], hole=F)@area  / Polygon(a_complete[,1:2], hole=F)@area 

# out of not associated / complete
Polygon(out_not_c_assoc_set[,1:2], hole=F)@area  / Polygon(a_complete[,1:2], hole=F)@area 



# ----------------------- complete space
# kernel densities
require(ks)
# optimal bandwidth estimation
hpi_mi_d1 <- Hpi(x = all_complete[,c("A1","A2")])

# kernel density estimation
est_mi_d1 <- kde(x = all_complete[,c("A1","A2")], 
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
  geom_point(data = all_complete, 
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

# plot convex hull

PCoA_plot_mi_d1 <- PCoA_plot_mi_d1 + geom_polygon(data=a_complete, aes (A1,A2),
                                                  alpha=0.05,
                                                  fill="black",
                                                  colour = "gray92",
                                                  size=1,
                                                  linetype = 1) + # complete space
  geom_text_repel(data = a_complete, aes (x=A1, y=A2, label=firstup(sp)),
                  size=3,col="black") + 
  theme(legend.position = "top")


# marginal plots 
PCoA_plot_mi_d1 <- ggMarginal(PCoA_plot_mi_d1, type="densigram",fill = "red",alpha=0.8)



# ---------------- coral associated
# optimal bandwidth estimation for coral associated fish
hpi_mi_d_CA <- Hpi(x = all_complete[-which(all_complete$sp %in% rownames(m_web_occ)),
                                    c("A1","A2")])

# kernel density estimation
est_mi_CA <- kde(x = all_complete[-which(all_complete$sp %in% rownames(m_web_occ)),
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
  geom_point(data = all_complete[-which(all_complete$sp %in% rownames(m_web_occ)),], 
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
                  size=3)+
  theme(legend.position = "top")

# plot density
PCoA_plot_mi_d2 <- ggMarginal(PCoA_plot_mi_d2, type="densigram",fill = "red",alpha=0.8)


# ---------------------- elimination of all 

# optimal bandwidth estimation for coral associated fish

hpi_mi_d_all <- Hpi(x = all_complete[-which(all_complete$sp %in% unlist(dimnames(m_web_occ))),
                                     c("A1","A2")])

# kernel density estimation
est_mi_all <- kde(x = all_complete[-which(all_complete$sp %in% unlist(dimnames(m_web_occ))),
                                   c("A1","A2")], 
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
  geom_point(data = all_complete[-which(all_complete$sp %in% unlist(dimnames(m_web_occ))),], 
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
                  size=3)+ 
  theme(legend.position = "top")

# plot density
PCoA_plot_mi_d3 <- ggMarginal(PCoA_plot_mi_d3, type="densigram",fill = "red",alpha=0.8)


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
  geom_point(data = all_complete[-which(all_complete$sp %in% rownames(m_web_occ)),], 
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
  ylim (c(-0.6,0.6))+
  theme(legend.position = "top")

# plot density
PCoA_plot_mi_diff <- ggMarginal(PCoA_plot_mi_diff, 
                                type="densigram",fill = "red",alpha=0.8)


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
  geom_point(data = all_complete[-which(all_complete$sp %in% unlist(dimnames(m_web_occ))),], 
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

# plot density
PCoA_plot_mi_diff_2 <- ggMarginal(PCoA_plot_mi_diff_2, 
                                  type="densigram",fill = "red",alpha=0.8)





## correlations to project trait values into the ordination
correlations <- cor (data.frame(sel_traits[,-which(colnames(sel_traits) == "scientificName")],
                                pco_complete$li[,1:3]),
                     use = "complete.obs")
correlations<-correlations [,c("A1","A2")]# interesting correlations

# clean plot  to receive correlations
clean_plot <- ggplot () +  +
  coord_equal() +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black",size =5),
        axis.title = element_text(colour = "black",size =7),
        legend.position = "top"
  )+
  xlim (c(-0.4,0.4))+
  ylim (c(-0.4,0.4))


# plot trait spaces
# plot density
clean_plot <- ggplot() +
  geom_polygon(data=a_complete, aes (A1,A2),
               alpha=0.5,
               fill="black",
               colour = "black",
               size=1,
               linetype = 2)+
  geom_polygon(data=not_c_assoc_set, aes (A1,A2),
               alpha=0.5,
               fill="white",
               colour = "white",
               size=0,
               linetype = 1)+
  geom_polygon(data=c_direct_indirect_set, aes (A1,A2),
               alpha=0.5,
               fill="green",
               colour = "gray92",
               size=0,
               linetype = 1)+
  geom_polygon(data=c_assoc_set, aes (A1,A2),
               alpha=0.5,
               fill="red",
               colour = "gray92",
               size=0,
               linetype = 1) +
  geom_polygon(data=out_not_c_assoc_set, aes (A1,A2),
               alpha=0.5,
               fill="pink",
               colour = "gray92",
               size=0,
               linetype = 1)+
  
  
  ## annotate
  annotate(geom="text",
           x=-0.15,
           y=0.25,
           label=paste (
             
             round (Polygon(not_c_assoc_set[,1:2], hole=F)@area  / Polygon(a_complete[,1:2], hole=F)@area,2)*100,
             "%", sep=""
           ),
           
           color="black")+
  annotate(geom="text",
           x=0.13,
           y=0.17,
           label=paste (
             
             round (Polygon(c_direct_indirect_set[,1:2], hole=F)@area  / Polygon(a_complete[,1:2], hole=F)@area,2)*100,
             "%", sep=""
           ),
           
           color="black")+
  annotate(geom="text",
           x=0.0,
           y=0.0,
           label=paste (
             
             round (Polygon(c_assoc_set[,1:2], hole=F)@area  / Polygon(a_complete[,1:2], hole=F)@area,2)*100,
             "%", sep=""
           ),
           
           color="black")  +
  geom_segment(aes(x = 0.2, y = 0.22, 
                   xend = 0.25, 
                   yend = 0.27),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm"))) + 
  annotate(geom="text",
           x=0.25,
           y=0.3,
           label=paste (
             
             100-round (Polygon(not_c_assoc_set[,1:2], hole=F)@area  / Polygon(a_complete[,1:2], hole=F)@area ,2)*100,
             "%", sep=""
           ))

# add points
clean_plot <- clean_plot + geom_point(data = all_complete, 
           aes(x = A1,y=A2),
           size = 2, alpha = 0.5, colour = "grey20")


clean_plot <- clean_plot + 
  
  
  # size
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[1,1]*0.5, 
                   yend = correlations[1,2]*0.5),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm")))  + 
  
  
  ## annotate
  annotate(geom="text",x=correlations[1,1]*0.5,
           y=correlations[1,2]*0.5,label="Body size",
           color="black") +
  
  # aspect ratio
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[2,1]*0.5, 
                   yend = correlations[2,2]*0.5),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm")))  + 
  
  
  ## annotate
  annotate(geom="text",x=correlations[2,1]*0.5,
           y=correlations[2,2]*0.5,label="Aspect ratio",
           color="black") +
  
  # trophic level
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[3,1]*0.5, 
                   yend = correlations[3,2]*0.5),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm"))) + 
  
  annotate(geom="text",x=correlations[3,1]*0.5,
           y=correlations[3,2]*0.5,
           label="Trophic level",
           color="black") +
  
  # group size
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[4,1]*0.3, 
                   yend = correlations[4,2]*0.3),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm"))) + 
  annotate(geom="text",x=correlations[4,1]*0.3,
           y=correlations[4,2]*0.3,label="Group size",
           color="black") +
  
  # temperature
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[5,1]*0.5, 
                   yend = correlations[5,2]*0.5),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm"))) + 
  annotate(geom="text",x=correlations[5,1]*0.50,
           y=correlations[5,2]*0.50,label="TºC max",
           color="black") + 
  
  # depth
  geom_segment(aes(x = 0, y = 0, 
                   xend = correlations[6,1]*0.5, 
                   yend = correlations[6,2]*0.5),size = 1,
               color="black",
               arrow = arrow(length = unit(.35, "cm"))) + 
  annotate(geom="text",x=correlations[6,1]*0.5,
           y=correlations[6,2]*0.5,label="Depth max",
           color="black")  +
  theme_classic() 


clean_plot <- ggMarginal(clean_plot, 
                        type="densigram",fill = "red",alpha=0.8)

# sasve


pdf (here ("output", "trait_space.pdf"),width=12,height=9,onefile = T)

grid.arrange(PCoA_plot_mi_d1,
             PCoA_plot_mi_d2,
             PCoA_plot_mi_d3,
             PCoA_plot_mi_diff,
             PCoA_plot_mi_diff_2,
             clean_plot,
             nrow=2,ncol=3)

dev.off()


# save the trait space (convex hulls)

pdf (here ("output", "trait_space_convex_hull.pdf"),
     width=5,height=5,onefile = T)

trait_spaces

dev.off()

