
# ==================================================
# organize trait dataset

organize_traits <- function (dataset) {
  
          trait_dataset <- dataset %>%
            
            dplyr::select(scientificName, 
                   family,
                   measurementValue,  # actually the measurement of size
                   Fish_age,
                   Body_size,
                   Aspect_ratio,Trophic_level, Size_group,
                   TempPref_max, Depth_max) %>%
            
            group_by(scientificName,Fish_age) %>% 
            
            summarise (family = unique(family),
                       actual_size=mean(as.numeric(measurementValue),na.rm=T),
                       Body_size = mean(Body_size,na.rm=T),
                       Aspect_ratio = mean(Aspect_ratio,na.rm=T),
                       Trophic_level = mean(Trophic_level,na.rm=T),
                       Size_group = unique(Size_group),
                       TempPref_max = mean(TempPref_max,na.rm=T),
                       Depth_max = mean(Depth_max,na.rm=T)) %>% # count the number of rows that a given answered appeared
            
            mutate (log_actual_size = log(actual_size),
                    log_Body_size = log(Body_size),
                    ordered_Size_group = ordered(Size_group, 
                                                 levels = c("sol","pair", "smallg",
                                                            "medg", "largeg"))) 
          
          
          # load trait data (to correct taxonomic inconsistency)
          #L.peixes$scientificName [which(L.peixes$scientificName == "Chaenopsis ocellata")] <- "Chaetodon ocellatus"
          
          traits <- read.csv (here ("\\.","Pos_Doc_Sinbiose","Review_within_the_group","coral_fish_project","data", "traits", "Atributos_especies_Atlantico_&_Pacifico_Oriental_2020_04_28.csv"),
                              sep = ";",h=T)
          
          # find and replace to correct
          # size group
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"), "Size_group"] <- traits[traits$Name == "Chaetodon ocellatus","Size_group"]
          # body size
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"), "Body_size"] <- as.numeric(traits[traits$Name == "Chaetodon ocellatus","Body_size"])
          # aspect ratio
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"), "Aspect_ratio"] <- as.numeric(gsub (",", ".", 
                                                                                                                          traits[traits$Name == "Chaetodon ocellatus","Aspect_ratio"]))
          # trophic level
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"), "Trophic_level"] <- as.numeric(gsub (",", ".", 
                                                                                                                           traits[traits$Name == "Chaetodon ocellatus","Trophic_level"]))
          # temp max
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"), "TempPref_max"] <- as.numeric(gsub (",", ".", 
                                                                                                                          traits[traits$Name == "Chaetodon ocellatus","TempPref_max"]))
          # depth
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"), "Depth_max"] <- as.numeric(gsub (",", ".", 
                                                                                                                       traits[traits$Name == "Chaetodon ocellatus","Depth_max"]))
          # actual size (from the dataset)
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"),"actual_size"]<- 14.5 # from Longo et al. 
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"),"log_actual_size"] <- log(14.5) # log actual size
          # log max tot body length
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"),"log_Body_size"] <- log(20)
          # schooling size
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"),"ordered_Size_group"] <- "pair"
          # chance the sp name itself
          trait_dataset[which(trait_dataset$scientificName == "Chaenopsis ocellata"),"scientificName"] <- "Chaetodon ocellatus"
          # and the family
          trait_dataset[which(trait_dataset$scientificName == "Chaetodon ocellatus"),"family"] <- "Chaetodontidae"
          
          
          trait_dataset[which(trait_dataset$scientificName == "Chaetodon ocellatus"),]
          
          # analyzed fish
          # number of analyzed species
          adult <- lapply (fish_species, function (i) i[[1]]) 
          adult<-unique(unlist(adult)) # adult
          adult[which(adult == "Chaenopsis ocellata")] <- "Chaetodon ocellatus"
          juvenile <- lapply (fish_species, function (i) i[[2]]) 
          juvenile<-unique(unlist(juvenile))#juvenile
          juvenile[which(juvenile == "Chaenopsis ocellata")] <- "Chaetodon ocellatus"
          table(juvenile %in% adult)# both
          
          # subset
          trait_dataset <- trait_dataset[which(trait_dataset$scientificName %in% unique(c(adult, juvenile))),]
          
          # ordering group size
          trait_dataset$Size_group <- sapply(trait_dataset$Size_group , function(x) {
            if (x=="sol") {1} 
            else if (x=="pair") {2} 
            else if (x=="smallg") {3} 
            else if (x=="medg") {4} 
            else if (x=="largeg") {5}}
          )

          
          return(trait_dataset)
}


###  kernel function from Cooke et al. 2019 (check it out here : https://github.com/03rcooke/hyper_pca/commit/2b7df79a30242d3d479e75382a8865df3f5a6f7d)

cl <- function(df, prob) {
  dx <- diff(df$x[1:2])
  dy <- diff(df$y[1:2])
  sz <- sort(df$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}


## toupper for species names

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}