# script to organize data

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

load (here("\\.","Pos_Doc_Sinbiose","Review_within_the_group","coral_fish_project","output_comm_wide_R1","Data_fish_detection_LONGO_AUED.RData"))

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
dir.create("data")
save(network_data, extracted_data, fish_size,fish_species,coral_species,extracted_data_occ,
     file = here ("data","occupancy_data.RData"))


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
