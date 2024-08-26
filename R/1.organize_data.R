
# ------------------------------------------------------------------


# script to organize data for network analysis and produce a map


# ------------------------------------------------------------------


rm(list=ls())

# load packages 
require(here); require(reshape); require(jagsUI); require(bipartite); require(igraph)
require(ggplot2);require(dplyr)

# load neeeded functions and additional packages
source(here ("\\.","Pos_Doc_Sinbiose","Review_within_the_group","coral_fish_project","R_comm_wide","packages.R"))
source(here ("\\.","Pos_Doc_Sinbiose","Review_within_the_group","coral_fish_project","R_comm_wide","functions.R"))
source(here ("\\.","Pos_Doc_Sinbiose","Review_within_the_group","coral_fish_project","R_comm_wide","quality_funct_space_fromdist2.R"))

# ----------------------------# 

# load fish and coral data

# ----------------------------# 

load (here("\\.","Pos_Doc_Sinbiose","Review_within_the_group",
           "coral_fish_project","output_comm_wide_R1",
           "Data_fish_detection_LONGO_AUED.RData"))

# load traits 
load(file = here ("data","occupancy_data.RData"))

# ------------------------------------------------- #
# site occupancy estimates and predictions
# ------------------------------------------------- #  

# load model output

load(here("\\.","Pos_Doc_Sinbiose","Review_within_the_group","coral_fish_project",
          "output_comm_wide_R1",
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
pdf(here ("output","fig_map.pdf"),width = 7,height = 5)#,units = "cm")#,res=300)
grid.arrange(p1,
             p2,
             ncol=2,
             widths = c(0.75,0.25))

dev.off()

# relating SR and coral cover

SR_all_fish <- do.call(rbind,lapply (extracted_data_occ, function (i) 
  (rowSums(i))))


# coral associated
trait_dataset<- organize_traits (fish_size)

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

# SR associated
associated <- lapply (extracted_data_occ, function (i) 
  
  i[,which(colnames(i) %in% unique(sp_analyzed_response$peixe))]

  )
SR_associated <- do.call(rbind,lapply (associated, function (i) 
  (rowSums(i))))

# FRic
comp_fish <- Reduce("+",extracted_data_occ)/8

all_traits <- sel_traits %>% 
  group_by(scientificName) %>%
  summarize (Aspect_ratio=mean(Aspect_ratio,na.rm=T),
             Trophic_level=mean(Trophic_level,na.rm=T),
             Size_group=mean(Size_group,na.rm=T),
             TempPref_max=mean(TempPref_max,na.rm=T),
             Depth_max=mean(Depth_max,na.rm=T),
             log_actual_size = mean(log_actual_size,na.rm=T ))


require(FD)
comp_fish<-(comp_fish[,which(colnames(comp_fish) %in% all_traits$scientificName)])
comp_fish<-comp_fish[,order(colnames(comp_fish))]
comp_fish <- ifelse (comp_fish>0.75,1,0) # discreticizing vals
comp_fish <- comp_fish[,which(colSums (comp_fish) >0)]
all_traits <- all_traits[which(all_traits$scientificName %in% colnames(comp_fish)),]
all_traits<-data.frame(all_traits)
rownames (all_traits)<- all_traits$scientificName; all_traits<-all_traits[,-1]




# distance matrix 
gower_matrix <- daisy (apply (all_traits,
                              2,scale), 
                       metric=("gower"),
                       type = list (ordratio = "Size_group")) 

dimnames(gower_matrix) <- rownames(all_traits)




FD_all <- dbFD (gower_matrix,
      comp_fish,
      corr="cailliez")

# coral ass
comp_fish_ass <- comp_fish[,which(colnames(comp_fish) %in% unique(sp_analyzed_response$peixe))]
comp_fish_ass <- comp_fish_ass[,which(colSums(comp_fish_ass)>0)]

all_traits_ass <- all_traits[which(rownames(all_traits) %in% colnames(comp_fish_ass)),]

FD_CA <- dbFD (all_traits_ass,
               comp_fish_ass ,
                corr="cailliez")


# df plots 
df_exploration <- data.frame (
  AllFish = colMeans(SR_all_fish),
  CoralAssoc = colMeans(SR_associated),
  #TS=FD_all$FRic,
  sum_coral_cover = rowSums(cob_corals))

df_exploration <- melt(df_exploration,id.vars = "sum_coral_cover")

ggplot(df_exploration, aes(x=sum_coral_cover,
                           y=value, 
                           group = variable,
                           color=variable))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_viridis_d(end=0.8)+theme_bw()+
  xlab ("Site-level total coral cover") + 
  ylab ("Expected species richness (sum of psi [site occ probability])")+
  annotate(geom="text",x=0.15,
           y=45,
           label="Beta=18.62, adjR2=0.16",
           color="purple4",alpha=0.5) +
  annotate(geom="text",
           x=0.15,
           y=30,
           label="Beta=35.89, adjR2=0.14",
           color="#2B7A0B",alpha=0.5) 


df2 <- data.frame (cob=rowSums(cob_corals),
                   fric=FD_all$FRic)

ggplot (df2,aes (x=cob,y=fric))+
  geom_point()+
  geom_smooth(method= "lm")

summary (lm (df2$fric~df2$cob))

with (df_exploration%>%
        filter(variable == "AllFish"),
      
      summary(lm(value~sum_coral_cover)))

with (df_exploration%>%
        filter(variable == "CoralAssoc"),
      
      summary(lm(value~sum_coral_cover)))


rm(list=ls())
