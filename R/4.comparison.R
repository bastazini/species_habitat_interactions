#         -----------------


#     comparison of scenarios


#         -----------------
require(here)
source(here ("R", "packages.R"))

# load random
load(here ("output", "random_scenario.RData"))

# extract the analysis data set and add the run
analysis_dataset_random <- lapply (seq (1,length(rob_analysis_random)), function (i) {
  
  
  # extract the analysis data set and add the run
  df_analysis<-data.frame (rob_analysis_random[[i]]$analysis_dataset,
                           run = i)
  # reorder
  df_analysis <- df_analysis[order(df_analysis$remain_coral,decreasing=F),]
  df_analysis
  
}
)
analysis_dataset_random <- do.call(rbind,analysis_dataset_random) # melt results



# load degree
load(here ("output", "degree_results.RData"))

# load vulnerability results
load (file = here ("output", "vulnerability_results.RData"))

# degree robustness curves ------------------------------------------------------

main_plot <- ggplot(analysis_dataset[order(analysis_dataset$remain_coral,decreasing=F),], 
                    aes (x= (1-remain_coral), y = remain_fish_associated)) +
  
  # functional diversity
  # direct loss
  
  geom_line( aes(y=pred_remain_RFS_associated),linetype=1,size=1.2,col="orange") + 
  geom_ribbon( aes(ymax=pred_remain_RFS_associated),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.3,
               fill="orange") + 
  
  geom_point( aes(y=remain_FD),
              size=3,
              stroke=2,
              col="orange",
              fill="white",
              shape=23,
              alpha=1) +
  
  
  # indirect loss
  geom_line( aes(y=pred_remain_RFS_secondary),linetype=1,size=1.2, col = "orange3") + 
  geom_ribbon( aes(ymax=pred_remain_RFS_secondary),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.4,
               fill="orange3") + 
  
  geom_point( aes(y=remain_RFS_secondary),size=3,
              stroke=2,
              col="orange3",
              fill="white",
              shape=21,
              alpha=1) + 
  
  # taxonomic diversity
  geom_ribbon( aes(ymax=pred_remain_SR_all),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.3,
               fill="gray50") + 
  geom_line( aes(y=pred_remain_SR_all),
             linetype=1,size=1.2,col="gray50") + 
  geom_point( aes(y=remain_all),size=3,
              stroke=2,
              col="gray50",
              fill="white",
              shape=23,
              alpha=1) + 
  
  #geom_line( aes(y=remain_fish_associated*100),linetype=1,size=1.2,col="orange") +
  geom_line( aes(y=pred_remain_SR_secondary),linetype=1,size=1.2,col="black") +
  geom_ribbon( aes(ymax=pred_remain_SR_secondary),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.8,
               fill="black") + 
  geom_point( aes(y=remain_SR_secondary),
              size=3,
              stroke=2,
              col="black",
              fill="white",
              shape=21,
              alpha=1) +
  
  
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Remaining Fish Taxonomic Diversity",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~.*1, name="Remaining Fish Functional Diversity")
  ) + 
  xlab ("Proportion of corals eliminated") + 
  theme_bw() + 
  ggtitle ("(a) Degree centrality")+
  theme (axis.title=element_text(size=14),
         axis.text = element_text(size=12),
         plot.title = element_text(size=18))


# vulnerability curves
vuln_plot <- ggplot(analysis_dataset_vulnerability[order(analysis_dataset_vulnerability$remain_coral,decreasing=F),], 
                    aes (x= (1-remain_coral), y = remain_fish_associated)) +
  
  # functional diversity
  # direct loss
  
  geom_line( aes(y=pred_remain_RFS_associated),linetype=1,size=1.2,col="orange") + 
  geom_ribbon( aes(ymax=pred_remain_RFS_associated),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.3,
               fill="orange") + 
  
  geom_point( aes(y=remain_FD),
              size=3,
              stroke=2,
              col="orange",
              fill="white",
              shape=23,
              alpha=1) +
  
  
  # indirect loss
  geom_line( aes(y=pred_remain_RFS_secondary),linetype=1,size=1.2, col = "orange3") + 
  geom_ribbon( aes(ymax=pred_remain_RFS_secondary),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.4,
               fill="orange3") + 
  
  geom_point( aes(y=remain_RFS_secondary),size=3,
              stroke=2,
              col="orange3",
              fill="white",
              shape=21,
              alpha=1) + 
  
  # taxonomic diversity
  geom_ribbon( aes(ymax=pred_remain_SR_all),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.3,
               fill="gray50") + 
  geom_line( aes(y=pred_remain_SR_all),
             linetype=1,size=1.2,col="gray50") + 
  geom_point( aes(y=remain_all),size=3,
              stroke=2,
              col="gray50",
              fill="white",
              shape=23,
              alpha=1) + 
  
  #geom_line( aes(y=remain_fish_associated*100),linetype=1,size=1.2,col="orange") +
  geom_line( aes(y=pred_remain_SR_secondary),linetype=1,size=1.2,col="black") +
  geom_ribbon( aes(ymax=pred_remain_SR_secondary),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.8,
               fill="black") + 
  geom_point( aes(y=remain_SR_secondary),
              size=3,
              stroke=2,
              col="black",
              fill="white",
              shape=21,
              alpha=1) +
  
  
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Remaining Fish Taxonomic Diversity",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~.*1, name="Remaining Fish Functional Diversity")
  ) + 
  xlab ("") + 
  theme_bw() + 
  ggtitle ("(c) Bleaching vulnerability")+
  theme (axis.title=element_text(size=14),
         axis.text=element_text(size=12),
         plot.title = element_text(size=18))


# plot ----------
#  SR direct -----------------------------
SR1 <- ggplot(analysis_dataset_random, 
       aes (x= (1-remain_coral), y = remain_fish_associated,group=run)) +
  
  # taxonomic diversity
  geom_ribbon( aes(ymax=pred_remain_SR_all),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.01,
               fill="orange") + 
  geom_line( aes(y=pred_remain_SR_all),linetype=1,size=1.2, col= "orange") + 
  
  geom_point( aes(y=remain_all),  
              size=3,
              stroke=2,
              col="orange",
              fill="white",
              shape=23,
              alpha=1
  ) + 
  

  # set third axis
  scale_y_continuous(
    
    # Features of the first axis
    name = "",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~.*1, name="")
  ) + 
  xlab ("") + 
  theme_bw() + 
  ggtitle ("(b) Random removal")+
  labs(subtitle = "Direct loss of Taxonomic Diversity")+
  theme (axis.title=element_text(size=14),
         axis.text = element_text(size=12),
         plot.title = element_text(size=18),
         plot.subtitle = element_text(size=14))



# SR indirect ---------------------------
SR2 <- ggplot(analysis_dataset_random, 
              aes (x= (1-remain_coral), y = remain_fish_associated,group=run)) +
  
  # taxonomic diversity
  geom_line( aes(y=pred_remain_SR_secondary),linetype=1,size=1.2,col="orange3",alpha=0.3) +
  
  geom_ribbon( aes(ymax=pred_remain_SR_secondary),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.01,
               fill="orange3") + 
  geom_point( aes(y=remain_SR_secondary),
              size=3,
              stroke=2,
              col="orange3",
              fill="white",
              shape=21,
              alpha=1) +
  # set third axis
  scale_y_continuous(
    
    # Features of the first axis
    name = "",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~.*1, name="")
  ) + 
  xlab ("") + 
  theme_bw() + 
  labs(subtitle = "Indirect loss of Taxonomic Diversity")+
  ggtitle ("")+
  theme (axis.title=element_text(size=14),
         axis.text = element_text(size=12),
         plot.title = element_text(size=18),
         plot.subtitle = element_text(size=14))


# FD direct -----------------------------

FD1 <- ggplot(analysis_dataset_random, 
              aes (x= (1-remain_coral), y = remain_fish_associated,group=run)) +
  
  # functional diversity
  # direct loss
  geom_line( aes(y=pred_remain_RFS_associated),linetype=1,size=1.2,col="gray50",alpha=0.2) + 
  
  

    geom_ribbon( aes(ymax=pred_remain_RFS_associated),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.01,
               fill="gray50") + 
  
  geom_point( aes(y=remain_FD),
              size=3,
              stroke=2,
              col="gray50",
              fill="white",
              shape=23,
              alpha=1) +

  # third axis
  scale_y_continuous(
    
    # Features of the first axis
    name = "",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~.*1, name="")
  ) + 
  xlab ("") + 
  theme_bw() + 
  labs(subtitle = "Direct loss of Functional Diversity")+
  ggtitle ("")+
  theme (axis.title=element_text(size=14),
         axis.text = element_text(size=12),
         plot.title = element_text(size=10),
         plot.subtitle = element_text(size=14))


# FD indirect -----------------------------

FD2 <- ggplot(analysis_dataset_random, 
              aes (x= (1-remain_coral), y = remain_fish_associated,group=run)) +
  
  # functional diversity
  geom_line( aes(y=pred_remain_RFS_secondary),linetype=1,size=1.2, col = "black") + 
  geom_ribbon( aes(ymax=pred_remain_RFS_secondary),
               ymin=0,
               linetype=1,
               size=1.2, 
               col = NA,
               alpha=0.01,
               fill="black") + 
  geom_point( aes(y=remain_RFS_secondary),
              size=3,
              stroke=2,
              col="black",
              fill="white",
              shape=21,
              alpha=1) + 
  
  # add third axis
  scale_y_continuous(
    
    # Features of the first axis
    name = "",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~.*1, name="")
  ) + 
  xlab ("") + 
  theme_bw() + 
  labs(subtitle = "Indirect loss of fish Functional Diversity")+
  ggtitle ("")+
  theme (axis.title=element_text(size=14),
         axis.text = element_text(size=12),
         plot.title = element_text(size=10),
         plot.subtitle = element_text(size=14))


pdf(here ("output", "robustness_curvesPDF.pdf"),width = 20,height = 7)
#png(here ("output", "robustness_curves.png"),width = 50,height = 15,units = "cm",res=300)

# arrange plots
grid.arrange(main_plot,
             SR1,SR2,
             FD1,FD2,
             vuln_plot,
             nrow=2,ncol=4,
             layout_matrix = rbind (c(1,2,3,6),
                                    c(1,4,5,6))
             )

dev.off()


# comparison --------------------------------------
# violin plot

# plot
do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))  %>%
  
  melt() %>%
  
  mutate(scenario = "random") %>% 
  ggplot (aes(x=Var2,y=value)) + 
  geom_violin(col="gray40",) +
  
  # random
  geom_point (data= do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))  %>%
                melt() %>%
                group_by(Var2) %>%
                summarize(value=mean(value)) %>%
                mutate (Loss = c("Direct","Indirect","Direct","Indirect")),
              aes(x=Var2,y=(value),shape=Loss),
              stroke=2,
              col="gray40",
              size=4) +

  geom_errorbar (data= do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))  %>%
                melt() %>%
                group_by(Var2) %>%
                summarize(LCI=quantile(value,0.025),
                          UCI=quantile(value,0.975)) %>%
                mutate (Loss = c("Direct","Indirect","Direct","Indirect")),
              aes(x=Var2,ymin=LCI,ymax=UCI),
              stroke=2,
              inherit.aes = F,
              size=1,
              col="gray40",
              width = 0.1)+


  scale_shape_manual(values = c("Direct" = 23, "Indirect"=21)) +
  geom_line (data= do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))  %>%
               melt() %>%
               group_by(Var2) %>%
               summarize(value=mean(value)) %>%
               mutate (gr = substr(Var2,1,2)),
             aes(x=Var2,y=(value),group=gr),
             col="gray40",
             size=1) + 
  
  # degree 
  geom_point (data= rob_analysis_degree %>%
                mutate (Loss = c("Direct","Indirect","Direct","Indirect")),
              aes(x=variable,y=(value),shape=Loss),
              size=4,
              stroke=2,
              col="orange")+
  scale_shape_manual(values =c("Direct" = 23, "Indirect"=21)) +
  geom_line(data= rob_analysis_degree %>%
              mutate (gr = substr(variable,1,2)),
            aes(x=variable,y=(value),group=gr),
            position = position_dodge(),
            size=1,
            col="orange")+
  
  # vulnerability
  # degree 
  geom_point (data= rob_analysis_vulnerability%>%
                mutate (Loss = c("Direct","Indirect","Direct","Indirect")),
              aes(x=variable,y=(value),shape=Loss),
              position = position_dodge(),
              size=4,
              stroke=2,
              col="blue")+
  scale_shape_manual(values = c("Direct" = 23, "Indirect"=21)) +
    geom_line(data= rob_analysis_vulnerability %>%
              mutate (gr = substr(variable,1,2)),
            aes(x=variable,y=(value),group=gr),
            position = position_dodge(),
            size=1,
            col="blue")+
  
  # general configs
  labs(x="",y="Network robustness (R)")+
  scale_x_discrete("Dimension",labels=c("Taxonomic Diversity","","Functional Diversity",""),
                   position = "bottom")+
  theme(axis.title = element_text(size=15),
      legend.position = c(0.085,0.9),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(hjust=0.15,size=12),
      axis.text.y = element_text(size=12))
  

ggsave (file = here ("output", "comparison.png"),width=5,height = 5)

# percentages
# random vs degree
table(do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) [,1] > rob_analysis_degree [1,2])/length(do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))[,1])
table(do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) [,2] > rob_analysis_degree [2,2])/length(do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))[,1])
table(do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) [,3] > rob_analysis_degree [3,2])/length(do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))[,1])
table(do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) [,4] > rob_analysis_degree [4,2])/length(do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))[,1])


# random vs vulnerability
table(do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) [,1] > rob_analysis_vulnerability [1,2])/length(do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))[,1])
table(do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) [,2] > rob_analysis_vulnerability [2,2])/length(do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))[,1])
table(do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) [,3] > rob_analysis_vulnerability [3,2])/length(do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))[,1])
table(do.call(rbind, sapply (rob_analysis_random, "[", "robustness")) [,4] > rob_analysis_vulnerability [4,2])/length(do.call(rbind, sapply (rob_analysis_random, "[", "robustness"))[,1])



