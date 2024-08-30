
# ------------------------------------------------------------------


# Project: Coping with Collapse: Functional Robustness of Coral-Reef Fish Network to Simulated Cascade Extinction


# Script 1: organize data for network analysis and produce a map (supp info)


# ------------------------------------------------------------------


rm(list=ls())

# load packages 
require(here); require(reshape); require(jagsUI); require(bipartite); require(igraph)
require(ggplot2);require(dplyr)

# load neeeded functions and additional packages
source(here ("R", "quality_funct_space_fromdist2.R"))
source(here ("R", "functions.R"))

# ----------------------------# 
# load fish and coral data
# ----------------------------# 


load (here ("data",
           "Data_fish_detection_LONGO_AUED.RData"))



# ------------------------------------------------- #
# site occupancy estimates and predictions
# ------------------------------------------------- #  

# load model output

load(here ("data",
           "samples_OCCcoral_PdepthTime_longo_RdmP.RData")) 


# chose standard deviations from the mean
chosen_sd_coral<- 4
chosen_sd_algae<- 0 # the average


## data to use in the plot
extracted_data <-   lapply (seq(1,length(coral_species)), function (cor)
  
  lapply (seq(1,length(samples_OCCcoral_PdepthTime_longo_RdmP[[cor]])) ,function (age) # across ages
    
    do.call(rbind,lapply(seq(1,length (fish_species[[cor]][[age]])), function (k)  # across fish species
      
      data.frame (
        
        # coral
        coral = coral_species[cor],
        
        # fish
        peixe = fish_species [[cor]][[age]][k],
        age = age,
        
        # regression coefficient and 90% credible interval
        # corals
        estimate.coral = mean (samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$sims.list$beta1 [,k]),
        low.coral = quantile (samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$sims.list$beta1 [,k], 0.05),
        high.coral = quantile (samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$sims.list$beta1 [,k],0.95),
        # turf
        estimate.turf = mean (samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$sims.list$beta2 [,k]),
        # estimate
        psi = samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$mean$n.occ [k]/36,#psi [,k]),
        p_exceedance = sum(samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$sims.list$psi[,,k]>0.5)/length(samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$sims.list$psi[,,k]), # 3000*36
        # predictions
        pred = plogis(samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$mean$intercept.psi[k]+
                        samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$mean$beta1[k]*chosen_sd_coral+
                        samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$mean$beta2[k]*chosen_sd_algae)
        
        
      ) # close df
    ) # close fish 
    ) # close fish
  ) # close age
) # close corals


# check how much these standard deviations means in terms of increasing coral cover
(round (apply(cob_corals,2,sd),4)*100)*chosen_sd_coral
sd(cob_algae$`calcareous turf`)*100

# melt these data within age
extracted_data <- lapply (extracted_data, function (coral)
  
  do.call(rbind, coral)
)

# adjusted wrong spp name
extracted_data <- lapply (extracted_data, function (i) {
  
  i$peixe[which(i$peixe == "Chaenopsis ocellata")] <- "Chaetodon ocellatus"
  
  i
}
)


# organize network data
network_data <- do.call(rbind,extracted_data)



# data for secondary/indirect extinctions 
## data to use in the plot
age=1
extracted_data_occ <-   lapply (seq(1,length(coral_species)), function (cor) 
  
  #lapply (seq(1,length(samples_OCCcoral_PdepthTime_longo_RdmP[[cor]])) ,function (age) # across ages
        
        do.call(cbind, lapply(seq(1,length (fish_species[[cor]][[age]])), function (k){  # across fish species
          
          # regression coefficient and 90% credible interval
          # corals
          psi <- data.frame(samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$mean$psi[,k])
          rownames (psi) <- cob_algae$sites_analysis
          colnames (psi) <- fish_species [[cor]][[age]][k]
          ;psi
          
        }
  
   ) # close fish
    ) # close do.call
  #) # close age
) # close corals


# adjusted wrong spp name
extracted_data_occ <- lapply (extracted_data_occ, function (i) {
  
  colnames(i)[which(colnames (i) == "Chaenopsis ocellata")] <- "Chaetodon ocellatus"
  
  i
}
)


#save
dir.create("processed_data")
save(cob_corals,coordenadas,network_data, extracted_data, fish_size,fish_species,coral_species,extracted_data_occ,
     file = here ("processed_data","occupancy_data.RData"))



# producing a map  -----------------------------------------------

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


# calculate values to map ------------------------------------------------------------

coordenadas <- data.frame (coordenadas,
            Corals = rowSums(cob_corals>0)/ncol(cob_corals),
            CAF = 
                    do.call(cbind, 
                            
                            lapply (extracted_data_occ, function (i)
                    
                                rowSums(i[,which(colnames(i) %in% rownames(m_web_occ))])/nrow(m_web_occ)
                    
                              )
                    ) %>%
                      
                      rowMeans() , 
                    COO = do.call(cbind, 
                                       
                                       lapply (extracted_data_occ, function (i)
                                         
                                         rowSums(i[,which(colnames(i) %in% colnames(m_web_occ))])/ncol(m_web_occ)
                                         
                                       )
                    ) %>%
              
              rowMeans()
)
       
# background

# mapa mundi
world <- ne_countries(scale = "medium", returnclass = "sf")

# crop mapa mundi
wm <- ggplot() + 
  geom_sf (data=world, size = 0.1, 
           fill= "gray90",colour="gray90") +
  coord_sf (xlim = c(-50, -30),  ylim = c(-30, 2), expand = FALSE) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightcyan",
                                        colour = "lightcyan"),
        axis.text.x = element_text(size=6),
        axis.ticks.x=element_line(size=1),
        axis.text.y = element_text(size=6),
        axis.ticks.y=element_line(size=1),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        title = element_blank()) 

# arrange
p1 <- wm + 
  
  geom_point(data = coordenadas,aes(x=decimalLongitude,
                                    y=decimalLatitude),
             size=4,alpha=0.3) + 
  geom_point(data=coordenadas[c(1,10,15,18,26,30,36),],
             aes(x=decimalLongitude,
                 y=decimalLatitude),
             size=4,col = "orange",stroke=2,shape=21)+
  labs(x="Longitude",y="Latitude")


# prop spp 
p2<-ggplot() +
  geom_bar(data=coordenadas[c(1,10,15,18,26,30,36),] %>%
             mutate("Site" = c("Abrolhos", "Costa dos Corais", "Espírito Santo",
                               "Ilha Bela", "Manuel Luis", "Parrachos", "Rocas")) %>%
             melt(id.vars = c("Group.1", "decimalLongitude","decimalLatitude", "Site")) %>%
             arrange(-decimalLatitude) ,
           aes (y=variable, x=value,fill=variable),
           stat = "identity")+
  scale_fill_manual(values= c("black", "orange", "gray"))+
  facet_wrap(~reorder(Site,-decimalLatitude),ncol=1)+
  theme_minimal()+
  xlim(c(0,1))+
  labs (y="",x="Proportion of the species")+
  theme(legend.position = "none",panel.grid =  element_blank(),
        strip.text = element_text(hjust = 0.02),
        axis.text.y =element_blank(),
        axis.title.x = element_text(size=9))


# arrange
dir.create("output")
pdf(here ("output","fig_map.pdf"),width = 7,height = 5)#,units = "cm")#,res=300)
grid.arrange(p1,
             p2,
             ncol=2,
             widths = c(0.75,0.25))

dev.off()

rm(list=ls())
