library(multifunc)
library(tidyverse)
library(patchwork)
library(purrr)
library(Hmisc)
library(pheatmap)
library(forcats)
library(ggpubr)
library(car)
library(scales)

enzyme_data_22 <- enzyme_22 %>%
  mutate(fertilization2 = factor(fertilization,levels = c('CK','HO','HF','MF')),
         stage = factor(stage,levels = c('V12','R2','mature'))) %>%
  select(u_BG,u_BX,u_CBH,stage,residue,fertilization2) %>%
  arrange(stage,residue,fertilization2)

fungal_traits_relative <- fungi_relative_table %>%
  cbind(taxa_fungi) %>%
  mutate(GENUS = sapply(strsplit(Genus,split = '_'),'[',3)) %>%
  left_join(fungi_traits , by = 'GENUS') %>%
  select(1:144,primary_lifestyle) %>%
  filter(primary_lifestyle %in% c('plant_pathogen','arbuscular_mycorrhizal','ectomycorrhizal')) %>%
  pivot_longer(!primary_lifestyle,names_to = 'sample_name',values_to = 'relative_abundance') %>%
  group_by(sample_name,primary_lifestyle) %>%
  mutate(relative_abundance = sum(relative_abundance)) %>%
  unique() %>%
  left_join(group,by = 'sample_name') %>%
  filter(str_detect(sample_name,'1|2|3|4')) %>%
  arrange(stage,residue,fertilization) %>%
  pivot_wider(names_from = 'primary_lifestyle',values_from = 'relative_abundance')

alpha_total_rep4 <- cbind(alpha_bac_plot[,1:3],
                          fungi_Shannon = alpha_fungi_plot$Shannon,fungi_richness = alpha_fungi_plot$Richness) %>%
  filter(str_detect(sample_name,'1|2|3|4')) %>%
  left_join(alpha_proto_plot[,1:3],by = 'sample_name')

emf_data <- cbind(physical22_sem[,-c(1:3)],CO2 = respiration_co2$CO2,
                  BG = enzyme_data_22$u_BG,BX = enzyme_data_22$u_BX,CBH = enzyme_data_22$u_CBH,
                  pathogen = c(1-fungal_traits_relative$plant_pathogen),
                  mycorrhizoa = c(fungal_traits_relative$ectomycorrhizal+fungal_traits_relative$arbuscular_mycorrhizal)) %>%
  as.data.frame() %>%
  select(MBN,MBC,CO2,BG,BX,CBH,pathogen,mycorrhizoa) %>%
  cbind(residue = physical22_sem$residue,.,bac_shannon =alpha_total_rep4$Shannon.x,
        fungi_shannon =alpha_total_rep4$fungi_Shannon,
        protists_shannon  = alpha_total_rep4$Shannon.y)

emf_heatmap_function <- function(residueeeeeee){
  emf_heatmap_data <- emf_data %>%
    filter(residue ==residueeeeeee) %>%
    select(!residue) %>%
    as.matrix()
  spearman_remove_emf<-rcorr(emf_heatmap_data,type = "spearman")
  spearman_remove_emf_padj <- p.adjust(spearman_remove_emf$P,method = "BH") %>% matrix(ncol = 11)
  r_matrixxx <- spearman_remove_emf$r
  r_matrix <- r_matrixxx[-c(9:11),9:11]
  print(r_matrix)
  print(spearman_remove_emf_padj[-c(9:11),9:11])
  star_matrix <- ifelse(spearman_remove_emf_padj < 0.001, "***",
                        ifelse(spearman_remove_emf_padj < 0.01, "**",
                               ifelse(spearman_remove_emf_padj < 0.05, "*", "")))
  pheatmap(t(r_matrix),
           display_numbers = t(star_matrix[-c(9:11),9:11]), 
           fontsize_number  = 16,
           number_format = "",            
           cluster_rows = FALSE,           
           cluster_cols = FALSE,          
           na_col = "white",              
           color = colorRampPalette(c("#7D0B5D", "white", "#428949"))(100),
           filename = paste0('caoxinzhuang_results/6.0emf/',residueeeeeee,'_emf_total.pdf'),
           width = 7, height = 2.18)
}

emf_heatmap_function('Turnover')
emf_heatmap_function('Remove')
while (!is.null(dev.list())) dev.off()
allVars<-c("MBN",'BX','CBH','BG','CO2','pathogen','mycorrhizoa', "MBC")

german_vars_total <- whichVars(emf_data, allVars)
species <- relevantSp(emf_data,5:ncol(emf_data))

ret <- dplyr::mutate(data, dplyr::across(vars, standardizeUnitScale, .names = "{.col}.std")) %>%
  dplyr::select(paste0(vars, ".std"))

# emf_test <- emf_data %>%
#   select(allVars) %>%
#   apply(.,2,scale) %>%
#   as.data.frame() %>%
#   rowwise() %>%
#   mutate(meanFunction = mean(c_across(everything()), na.rm = TRUE))

emf_total<-cbind(emf_data, getStdAndMeanFunctions(emf_data, german_vars))

unique(physical22_sem$stage)
emf_anova <- emf_total[,-1] %>%
  cbind(physical22_sem[,1:3])  %>%
  mutate(randf = paste(residue,fertilization,sep = '_')) %>%
  arrange(stage,residue,fertilization)

emf_anova %>%
  group_by(stage,residue) %>%
  summarise(mean_emf = mean(meanFunction)) %>%
  pivot_wider(names_from = residue,values_from = mean_emf) %>%
  mutate(proportion = (Turnover-Remove)/Remove)

emf_anova %>%
  cbind(weight = sem_data$weight,height = sem_data$height) %>%
  filter(stage == 'Maturing'&residue=='Remove') %>%
  lmer(height~meanFunction+(1|fertilization),data = .) %>%
  summary()

for (i in unique(emf_anova$stage)) {
  dataaa <- filter(emf_anova,stage ==i )
  alpha_proto_aov <- aov(meanFunction~residue*fertilization,data = dataaa)
  alpha_tukey_aov <- aov(meanFunction~randf,data = dataaa)
  anova_summary <- summary(alpha_proto_aov)
  leveneTestaaa <- leveneTest(meanFunction~residue*fertilization,data = dataaa, center = mean)
  residues <- residuals(alpha_proto_aov) %>% shapiro.test()
  kruskal_fertilization <- kruskal.test(meanFunction ~ fertilization, data = dataaa)
  wilcox_residue <- wilcox.test(meanFunction ~ residue, data = dataaa)
  print('------------------------------------')
  print(i)
  ifelse(residues$p.value<0.05|leveneTestaaa[[3]]<0.05,print(kruskal_fertilization),print(anova_summary))
  ifelse(residues$p.value<0.05|leveneTestaaa[[3]]<0.05,print(wilcox_residue),print(i))
}
emf_anova
emf_total_p <- emf_anova %>%
  ggplot(aes(x = fertilization, y = meanFunction)) + 
  geom_boxplot(aes(  fill = residue),
               position = position_dodge(0.9),alpha = 0.8,outlier.colour = 'grey70')+
  scale_fill_manual(values = CAOXINZHUANG2)+
  scale_color_manual(values = CAOXINZHUANG2)+
  scale_shape_manual(values = c(17,15,16,21))+
  labs(y = 'SMF' , x = NULL)+ 
  facet_grid(.~stage)+
  theme_NMDS2+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.protokground = element_blank(),
        strip.text = element_blank())
emf_total_p

ggsave(emf_total_p,file = 'caoxinzhuang_results/6.0emf/emf_total_p_review.pdf',width = 5,height = 3.09)


scale_01_function <- function(x){
  (x-min(x))/(max(x)-min(x))
}

remove_input_data <-  emf_data %>%
  filter(residue== 'Remove') %>%
  mutate(bac_shannon_scale = scale_01_function(bac_shannon),
         fungi_shannon_scale = scale_01_function(fungi_shannon),
         protists_shannon_scale = scale_01_function(protists_shannon),
         total = scale_01_function(c(bac_shannon_scale+fungi_shannon_scale+protists_shannon_scale))) %>%
  select(total,bac_shannon_scale,fungi_shannon_scale,protists_shannon_scale,
         "MBN",'BX','CBH','BG','CO2','pathogen','mycorrhizoa', "MBC")
turnover_input_data <-  emf_data %>%
  filter(residue== 'Turnover') %>%
  mutate(bac_shannon_scale = scale_01_function(bac_shannon),
         fungi_shannon_scale = scale_01_function(fungi_shannon),
         protists_shannon_scale = scale_01_function(protists_shannon),
         total = scale_01_function(c(bac_shannon_scale+fungi_shannon_scale+protists_shannon_scale))) %>%
  select(total,bac_shannon_scale,fungi_shannon_scale,protists_shannon_scale,
         "MBN",'BX','CBH','BG','CO2','pathogen','mycorrhizoa', "MBC")

emf_threshold_function_mean_value <- function(input_data,residueeeee,ressssidue){
  german_vars <- whichVars(input_data, allVars)
  species <- relevantSp(input_data,5:ncol(input_data))
  input_data<-cbind(input_data, getStdAndMeanFunctions(input_data, german_vars))
  write.csv(input_data,file = paste0('caoxinzhuang_results/5.0/intermediate/',residueeeee,'EMF.csv'))
  lm_data <-input_data %>%
    select(meanFunction,total,bac_shannon_scale,fungi_shannon_scale,protists_shannon_scale) %>%
    pivot_longer(!c(meanFunction),names_to = 'Diversity')
  for (diversssity in c('total','bac_shannon_scale','fungi_shannon_scale','protists_shannon_scale')) {
    ave_fit <- lm(as.formula(paste('meanFunction~',diversssity)),
                  data =  input_data)
    print(summary(ave_fit))
  }
  print(lm_data)
  ggplot(aes(x = value, y = meanFunction,col = Diversity), data = lm_data) +
    theme_bw(base_size = 15) + 
    stat_smooth(method = "lm", size = 1,se=FALSE) +
    scale_color_manual(values = sunsx)+
    labs(x = 'Diversity',y = 'SMF',subtitle = ressssidue)+
    theme_NMDS2
}

harvest_mean_plot <- emf_threshold_function_mean_value(remove_input_data ,'Remove','Harvest')
retention_mean_plot <- emf_threshold_function_mean_value(turnover_input_data,'Turnover','Retention')
mean_emf_plot_total <- ggarrange(harvest_mean_plot,retention_mean_plot,ncol = 2)
mean_emf_plot_total
ggsave(mean_emf_plot_total,file = 'caoxinzhuang_results/6.0emf/diversity_emf_total.pdf',width = 10,height = 5)

germanyThresh_data <- data.frame()
single_threshold_emf <- function(input_dataaaaa,alpha_diversity){
  input_data <- input_dataaaaa %>%
    rstatix::mutate(Diversity = input_dataaaaa[,alpha_diversity]) %>%
    dplyr::select(Diversity,"MBN",'BX','CBH','BG','CO2','pathogen','mycorrhizoa', "MBC")
  german_vars <- whichVars(input_data, allVars)
  germanyThresh<-getFuncsMaxed(input_data, german_vars, threshmin=0.05, threshmax=0.99,
                               prepend=c("plot","Diversity"), maxN=4)
  germanyThresh_data <<- germanyThresh
    for (iiii in c(0.1,0.25,0.5,0.65)) {
      print(iiii)
      dataaaadata <- filter(germanyThresh,thresholds == iiii)
      ifelse(nrow(dataaaadata)==0,print('Bucheng________________'),
             print(summary(mfuncGermanyLinear<-glm(funcMaxed ~ Diversity, 
                                                   data=dataaaadata, 
                                                   family=quasipoisson(link="identity")))))
      print('__________________________________')
    }
  gcPlot<-subset(germanyThresh, germanyThresh$thresholds %in% qw( 0.1,0.25,0.5,0.65))
  gcPlot$percent<-paste(100*gcPlot$thresholds, "%", sep="")
  gcPlot %>%
    ggplot(aes(x = Diversity,y = funcMaxed,col = percent))+
    scale_y_continuous(limits = c(1,8))+
    stat_smooth(method="glm",method.args = list(family=quasipoisson(link="identity")),se=FALSE)
}

germanyThresh_data <- data.frame()
col_icamp
remove_bacteria_single_threshold <- single_threshold_emf(remove_input_data,2)+
  theme_NMDS2+
  labs(y = 'Number of functions above the threshold',subtitle = 'Bacteria')+
  scale_color_manual(values =col_icamp)
remove_fungi_single_threshold <- single_threshold_emf(remove_input_data,3)+
  theme_NMDS2+
  labs(y = 'Number of functions above the threshold',subtitle = 'Fungi')+
  scale_color_manual(values =col_icamp)
remove_protists_single_threshold <- single_threshold_emf(remove_input_data,4)+
  theme_NMDS2+
  labs(y = 'Number of functions above the threshold',subtitle = 'Protists')+
  scale_color_manual(values =col_icamp)

remove_protists_single_threshold
plot(remove_protists_single_threshold)
remove_fungi_single_threshold
remove_bacteria_single_threshold

turnover_bacteria_single_threshold <- single_threshold_emf(turnover_input_data,2)+
  theme_NMDS2+
  labs(y = 'Number of functions above the threshold',subtitle = 'Bacteria')+
  scale_color_manual(values =col_icamp)
turnover_fungi_single_threshold <- single_threshold_emf(turnover_input_data,3)+
  theme_NMDS2+
  labs(y = 'Number of functions above the threshold',subtitle = 'Fungi')+
  scale_color_manual(values =col_icamp)
turnover_protists_single_threshold <- single_threshold_emf(turnover_input_data,4)+
  theme_NMDS2+
  labs(y = 'Number of functions above the threshold',subtitle = 'Protists')+
  scale_color_manual(values =col_icamp)
turnover_protists_single_threshold
turnover_fungi_single_threshold
turnover_bacteria_single_threshold



single_threshold_total_p <- ggarrange(remove_bacteria_single_threshold,remove_fungi_single_threshold,
                                      remove_protists_single_threshold,turnover_bacteria_single_threshold,
                                      turnover_fungi_single_threshold,turnover_protists_single_threshold,ncol = 3,
                                      nrow = 2)
ggsave(single_threshold_total_p,file = 'caoxinzhuang_results/6.0emf/diversity_emf_multithreshold.pdf',width = 10,height =6.18 )

EMF_supplementary_data <- function(data){
  germanyThresh <- data
  germanyThresh$percent <- 100*germanyThresh$thresholds
  germanyThresh
  ggplot(data = germanyThresh, aes(x = Diversity, y = funcMaxed, group = percent)) +
    ylab(expression("Number of Functions >= Threshold")) +
    xlab("Diversity") +  
    stat_smooth(method = "glm",
                method.args = list(family = quasipoisson(link = "identity")), lwd = 0.8, fill = NA,
                aes(color = percent)) +  
    theme_bw(base_size = 14) +
    scale_color_gradientn(name = "Percent of \nMaximum",
                          colours = rev(c("red", "orange", "yellow", "green", "cyan", "blue", "purple")), 
                          breaks = seq(0, 1, by = 0.1), labels = scales::percent(seq(0, 1, by = 0.1)),
                          guide = guide_colorbar(reverse = FALSE))
}

EMF_supplementary_data(germanyThresh_data)

germanyLinearSlopes<-getCoefTab(funcMaxed ~ Diversity, data=germanyThresh_data,
                                coefVar="Diversity", family=quasipoisson(link="identity"))
germanyLinearSlopes
####### Plot the values of the diversity slope at# different levels of the threshold######
germanSlopes <- ggplot(germanyLinearSlopes, aes(x = thresholds * 100,
                                                y = estimate,
                                                ymax = estimate + 1.96 * std.error, 
                                                ymin = estimate - 1.96 * std.error)) + 
  geom_ribbon(fill = "#FFCC99") +  
  geom_point(color = "#993399") + 
  ylab("Change in Number of Functions\n") + 
  xlab("\nThreshold (%)") +  
  geom_abline(intercept = 0, slope = 0, lwd = 1, linetype = 2) + 
  theme_bw(base_size = 14)

germanSlopes

germanIDX <- getIndices(germanyLinearSlopes, germanyThresh_data, funcMaxed ~ Diversity)
germanIDX

