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
extracted_data_occ <-   lapply (seq(1,length(coral_species)), function (cor) 
  
  # lapply (seq(1,length(samples_OCCcoral_PdepthTime_longo_RdmP[[cor]])) ,function (age) # across ages
  
  do.call(cbind, lapply(seq(1,length (fish_species[[cor]][[age]])), function (k){  # across fish species
    
    # regression coefficient and 90% credible interval
    # corals
    psi <- data.frame(samples_OCCcoral_PdepthTime_longo_RdmP[[cor]][[age]]$mean$psi[,k])
    rownames (psi) <- cob_algae$sites_analysis
    colnames (psi) <- fish_species [[cor]][[age]][k]
    ;psi
  }
  
  ) # close fish
  #  ) # close age
  )) # close corals

#save
dir.create("data")
save(network_data, extracted_data, fish_size,fish_species,coral_species,extracted_data_occ,
     file = here ("data","occupancy_data.RData"))
rm(list=ls())
