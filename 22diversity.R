#load packages
library(tidyverse)
library(ggpubr)
# library(ggbeeswarm)
library(gghalves)
library(car)
library(vegan)
library(rstatix)
library(reshape2)

# prepare data
asv_bac <- read.csv('sequence/bac_ASV.csv' , row.names = 1 , header = T)
asv_fungi <- read.csv('sequence/fungi_ASV.csv',row.names = 1,header = T)

asv_proto <- read.csv('sequence/proto_ASVs.csv',row.names = 1,header = T) %>%
  filter(Kingdom == "Eukaryota" &(Phylum %in% c("Amoebozoa",'Excavata','Obazoa', "TSAR"))) %>%
  select(!c(v.R.C.R.5,R.R.D.R.5))
dim(asv_fungi)
group <- read.delim('sequence/caoxinzhuang_group.txt' ,row.names = 1 , header = T)
taxa_bac <- asv_bac[,145:151] 
taxa_fungi <- asv_fungi[,145:151]

taxa_protists <- asv_proto[,143:149]
bac_asv <- asv_bac[,1:144] 
fungi_asv <- asv_fungi[,1:144]
protists_asv <- asv_proto[,1:142]

#rare curve
rare_curve_plot <- function(asv,asadsdadsa){
  rare_curve_bac <- rarecurve(t(asv), step = 500,  xlab = "Sample Size", ylab = "Species",label = F)
  rare_curve_bac
  
  rare_curve_df_bac <- data.frame(
    Sample = character(),
    Reads = numeric(),
    Species = numeric(),
    stringsAsFactors = FALSE
  )
  for(i in 1:length(rare_curve_bac)) {

    sample_data <- rare_curve_bac[[i]]
    sample_name <- colnames(bac_asv)[i]
    temp_df <- data.frame(
      Sample = sample_name,  
      Reads = attr(sample_data, "Subsample"),
      Species = as.vector(sample_data),
      stringsAsFactors = FALSE
    )
    

    rare_curve_df_bac <- rbind(rare_curve_df_bac, temp_df)
  }
  head(rare_curve_df_bac)
  rare_curve_df_bac <- rare_curve_df_bac %>%
    left_join(group,by = c('Sample' = "sample_name")
    )
  
  ggplot(rare_curve_df_bac, aes(x = Reads, y = Species, group = Sample, color = residue)) +
    geom_line() +
    theme_bw() +
    labs(title = asadsdadsa) + 
    theme(legend.position = "right") +
    scale_color_manual(values = CAOXINZHUANG2)+
    theme_NMDS2
}
bac_rare_curve <- rare_curve_plot(bac_asv,'Bacteria')
fungi_rare_curve <- rare_curve_plot(fungi_asv,'Fungi')
protists_rare_curve <- rare_curve_plot(protists_asv,'Protists')

rare_curve_total <- ggarrange(bac_rare_curve,fungi_rare_curve,protists_rare_curve,ncol = 2,nrow = 2)
rare_curve_total
ggsave(filename = 'caoxinzhuang_results/6.0emf/supplementary/rare_curve.pdf',rare_curve_total,height = 6.18,width = 10)

bac_asv_rari <- as.data.frame(t(rrarefy(t(bac_asv),min(colSums(bac_asv))))) %>%
  filter(.,rowSums(.)>0)
colSums(bac_asv_rari)
fungi_asv_rari <- as.data.frame(t(rrarefy(t(fungi_asv),min(colSums(fungi_asv))))) %>%
  filter(.,rowSums(.)>0)
colSums(fungi_asv_rari)
proto_asv_rari <- as.data.frame(t(rrarefy(t(proto_asv),min(colSums(proto_asv))))) %>%
  filter(.,rowSums(.)>0)
colSums(proto_asv_rari)
# arrange group
group <- rownames_to_column(group , var = 'sample_name')
group$stage <- factor(group$stage , levels = c('Jointing','Filling','Maturity'))
group$residue <- factor(group$residue,levels = c('Remove','Turnover'))
group$fertilization <- factor(group$fertilization , 
                              levels = c('ck','high_organic','high_fertilizer','mineral_fertilizer'))

#### useful things####
alpha <- function(x,tree=NULL){
  est <- estimateR(x)
  Richness <- est[1,]
  Shannon <- diversity(x,index='shannon',base=exp(1))
  result <- data.frame(Richness,Shannon)
  result
}
CAOXINZHUANG2 <- c('#B44847','#00A1D6')
theme_NMDS <- theme_bw(base_size = 15)+
  theme(legend.position = 'none',panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 16,color="black"),
        axis.text.y = element_text(size = 16,color="black"),
        axis.title.y = element_text(size = 18,color="black"),
        axis.title.x = element_text(size = 18,color="black"))

####alpha####
#####bacteria#####

alpha_bac <- bac_asv_rari%>%
  t() %>%
  alpha() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name')

alpha_bac_plot <- left_join(alpha_bac,group,by = "sample_name") %>%
  mutate(randf = paste(residue,fertilization,sep = '_')) %>%
  arrange(stage,residue,fertilization)
head()
alpha_proto_plot %>%
  group_by(stage,residue)%>%
  summarise(mean_shannon = mean(Shannon))


for (i in unique(alpha_bac_plot$stage)) {
  dataaa <- filter(alpha_bac_plot,stage ==i )
  alpha_proto_aov <- aov(Shannon~residue*fertilization,data = dataaa)
  alpha_tukey_aov <- aov(Shannon~randf,data = dataaa)
  tukey <- tukey_hsd(alpha_tukey_aov)
  tukey_sig <- tukey %>% filter(p.adj < 0.05)
  anova_summary <- summary(alpha_proto_aov)
  leveneTestaaa <- leveneTest(Shannon~residue*fertilization,data = dataaa, center = mean)
  residues <- residuals(alpha_proto_aov) %>% shapiro.test()
  kruskal_fertilization <- kruskal.test(Shannon ~ fertilization, data = dataaa)
  wilcox_residue <- wilcox.test(Shannon ~ residue, data = dataaa)
  print('------------------------------------')
  print(i)
  ifelse(residues$p.value<0.05|leveneTestaaa[[3]]<0.05,print(kruskal_fertilization),print(anova_summary))
  ifelse(residues$p.value<0.05|leveneTestaaa[[3]]<0.05,print(wilcox_residue),print(tukey_sig))
}

alpha_stage_bac_fertilizer <- ggplot(alpha_bac_plot,aes(x = fertilization, y = Shannon)) + 
  geom_boxplot(aes( fill = residue),
               position = position_dodge(0.9),alpha = 0.8,outlier.colour = 'grey70')+
  scale_fill_manual(values = CAOXINZHUANG2)+
  scale_color_manual(values = CAOXINZHUANG2)+
  scale_y_continuous(limits = c(3.4,6.5),breaks = seq(3.5, 6.5, by = 0.5))+
  scale_shape_manual(values = c(17,15,16,21))+
  labs(y = 'Bacteria Shannon Diversity' , 
       x = NULL)+ 
  facet_grid(.~stage)+
  theme_NMDS2+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_blank(),
        strip.text = element_blank())

alpha_stage_bac_fertilizer
#####fungi#####
alpha_fungi <- fungi_asv_rari %>%
  t() %>%
  alpha() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name')
group 
alpha_fungi_plot <- left_join(alpha_fungi,group,by = "sample_name") %>%
  mutate(randf = paste(residue,fertilization,sep = '_')) %>%
  arrange(stage,residue,fertilization)

for (i in unique(alpha_fungi_plot$stage)) {
  dataaa <- filter(alpha_fungi_plot,stage ==i )
  alpha_proto_aov <- aov(Shannon~residue*fertilization,data = dataaa)
  alpha_tukey_aov <- aov(Shannon~randf,data = dataaa)
  anova_summary <- summary(alpha_proto_aov)
  leveneTestaaa <- leveneTest(Shannon~residue*fertilization,data = dataaa, center = mean)
  residues <- residuals(alpha_proto_aov) %>% shapiro.test()
  kruskal_fertilization <- kruskal.test(Shannon ~ fertilization, data = dataaa)
  wilcox_residue <- wilcox.test(Shannon ~ residue, data = dataaa)
  print('------------------------------------')
  print(i)
  ifelse(residues$p.value<0.05|leveneTestaaa[[3]]<0.05,print(kruskal_fertilization),print(anova_summary))
  ifelse(residues$p.value<0.05|leveneTestaaa[[3]]<0.05,print(wilcox_residue),print(i))
}

alpha_stage_fungi_fertilizer <- ggplot(alpha_fungi_plot,aes(x = fertilization, y = Shannon)) + 
  geom_boxplot(aes(  fill = residue),
               position = position_dodge(0.9),alpha = 0.8,outlier.colour = 'grey70')+
  scale_fill_manual(values = CAOXINZHUANG2)+
  scale_color_manual(values = CAOXINZHUANG2)+
  scale_shape_manual(values = c(17,15,16,21))+
  labs(y = 'Fungi Shannon Diversity' ,
       x = NULL)+ 
  facet_grid(.~stage)+
  theme_NMDS2+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_blank(),
        strip.text = element_blank())

alpha_stage_fungi_fertilizer
#####protozoa#####
alpha_proto <- proto_asv_rari%>%
  t() %>%
  alpha() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name')

alpha_proto_plot <- left_join(alpha_proto,group,by = "sample_name") %>%
  mutate(randf = paste(residue,fertilization,sep = '_')) %>%
  arrange(stage,residue,fertilization)

lm(Richness~Shannon,data = alpha_proto_plot) %>% summary()

for (i in unique(alpha_fungi_plot$stage)) {
  dataaa <- filter(alpha_proto_plot,stage ==i )
  alpha_proto_aov <- aov(Shannon~residue*fertilization,data = dataaa)
  alpha_tukey_aov <- aov(Shannon~randf,data = dataaa)
  anova_summary <- summary(alpha_proto_aov)
  leveneTestaaa <- leveneTest(Shannon~residue*fertilization,data = dataaa, center = mean)
  residues <- residuals(alpha_proto_aov) %>% shapiro.test()
  kruskal_fertilization <- kruskal.test(Shannon ~ fertilization, data = dataaa)
  wilcox_residue <- wilcox.test(Shannon ~ residue, data = dataaa)
  print('------------------------------------')
  print(i)
  ifelse(residues$p.value<0.05|leveneTestaaa[[3]]<0.05,print(kruskal_fertilization),print(anova_summary))
  ifelse(residues$p.value<0.05|leveneTestaaa[[3]]<0.05,print(wilcox_residue),print(i))
}

alpha_stage_proto_fertilizer <-
  ggplot(alpha_proto_plot,aes(x = fertilization, y = Shannon)) + 
    geom_boxplot(aes(  fill = residue),
                 position = position_dodge(0.9),alpha = 0.8,outlier.colour = 'grey70')+
  scale_fill_manual(values = CAOXINZHUANG2)+
  scale_color_manual(values = CAOXINZHUANG2)+
  scale_shape_manual(values = c(17,15,16,21))+
  labs(y = 'Protozoa Shannon Diversity' , 
       x = NULL)+ 
  facet_grid(.~stage)+
  theme_NMDS2+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.protokground = element_blank(),
        strip.text = element_blank())
alpha_stage_proto_fertilizer
#####output#####


alpha_proto_plot %>%
  select(Richness:Shannon) %>%
  as.matrix() %>%
  rcorr(. , type = 'spearman')


alpha_p <- ggarrange(emf_total_p,alpha_stage_bac_fertilizer,alpha_stage_fungi_fertilizer,
                     alpha_stage_proto_fertilizer,nrow = 2,ncol = 2,align = 'h')
alpha_p
ggsave(alpha_p,filename='caoxinzhuang_results/6.0emf/alpha2.0.pdf',width = 10,height = 6.18)

# ggsave(alpha_stage_fungi_fertilizer,
#        filename='caoxinzhuang_results/2.0/test/alpha_stage_fungi_fertilizer.pdf',
#        width = 8,height = 4)



#### beta####
#####bacteria#####
# Jointing Stage

bac_distance_jointing <- bac_asv_rari %>%
  select(starts_with('V')) %>%
  t() %>%
  vegdist(method = 'bray')

bac_data_jointing <- bac_asv_rari %>%
  select(starts_with('V')) %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by = 'sample_name')
PERMANOVA_bac_jointing <- adonis2(bac_distance_jointing ~ residue * fertilization , 
                                  bac_data_jointing , permutations = 9999)
PERMANOVA_bac_jointing
bac_NMDS_jointing <- metaMDS(bac_distance_jointing, k = 2)
bac_NMDS_jointing$stress
bac_NMDS_jointing_points <- as.data.frame(bac_NMDS_jointing$points) %>% 
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by='sample_name')

bac_NMDS_jointing_p <-ggplot(bac_NMDS_jointing_points ,aes(x=MDS1,y=MDS2))+
  geom_point(aes(shape=fertilization,color=residue),size=3.5)+
  scale_color_manual(values = CAOXINZHUANG2)+
  scale_shape_manual(values = c(17,15,16,21))+
  labs(x=NULL,y=NULL)+
  geom_text(x = min(bac_NMDS_filling_points$MDS1) + 0.05,
            y = min(bac_NMDS_filling_points$MDS2) +0.2 ,size = 6,
            label = paste0('Stress = ',round(bac_NMDS_jointing$stress , digits = 3)))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme_NMDS
bac_NMDS_jointing_p
#Filling Stage
bac_distance_filling <- bac_asv_rari %>%
  select(starts_with('R')) %>%
  t() %>%
  vegdist(method = 'bray')

bac_data_filling <- bac_asv_rari %>%
  select(starts_with('R')) %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by = 'sample_name')

PERMANOVA_bac_filling <- adonis2(bac_distance_filling ~ residue * fertilization , 
                                  bac_data_filling , permutations = 9999)
PERMANOVA_bac_filling
bac_NMDS_filling <- metaMDS(bac_distance_filling, k = 2)
bac_NMDS_filling$stress
bac_NMDS_filling_points <- as.data.frame(bac_NMDS_filling$points) %>% 
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by='sample_name')
bac_NMDS_filling_p <- ggplot(bac_NMDS_filling_points ,aes(x=MDS1,y=MDS2))+
  geom_point(aes(shape=fertilization,color=residue),size=3.5)+
  scale_color_manual(values = CAOXINZHUANG2)+
  scale_shape_manual(values = c(17,15,16,21))+
  labs(x=NULL,y=NULL)+
  geom_text(x = min(bac_NMDS_filling_points$MDS1) + 0.2,
            y = min(bac_NMDS_filling_points$MDS2) +0.05 ,size = 6,
            label = paste0('Stress = ',round(bac_NMDS_filling$stress , digits = 3)))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme_NMDS
bac_NMDS_filling_p
#Maturing Stage
bac_distance_maturing <- bac_asv_rari %>%
  select(starts_with('m')) %>%
  t() %>%
  vegdist(method = 'bray')

bac_data_maturing <- bac_asv_rari %>%
  select(starts_with('m')) %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by = 'sample_name')
PERMANOVA_bac_maturing <- adonis2(bac_distance_maturing ~ residue * fertilization , 
                                 bac_data_maturing , permutations = 9999)
PERMANOVA_bac_maturing
bac_NMDS_maturing <- metaMDS(bac_distance_maturing, k = 2)
bac_NMDS_maturing$stress
bac_NMDS_maturing_points <- as.data.frame(bac_NMDS_maturing$points) %>% 
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by='sample_name')
bac_NMDS_maturing_p <- ggplot(bac_NMDS_maturing_points ,aes(x=MDS1,y=MDS2))+
  geom_point(aes(shape=fertilization,color=residue),size=3.5)+
  scale_color_manual(values = CAOXINZHUANG2)+
  scale_shape_manual(values = c(17,15,16,21))+
  labs(x=NULL,y=NULL)+
  geom_text(x = min(bac_NMDS_maturing_points$MDS1) + 0.25,
            y = min(bac_NMDS_maturing_points$MDS2) +0.05 ,size = 6,
            label = paste0('Stress = ',round(bac_NMDS_maturing$stress , digits = 3)))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme_NMDS
bac_NMDS_maturing_p
#####fungi#####
# Jointing Stage
fungi_distance_jointing <- fungi_asv_rari %>%
  select(starts_with('V')) %>%
  t() %>%
  vegdist(method = 'bray')

fungi_data_jointing <- fungi_asv_rari %>%
  select(starts_with('V')) %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by = 'sample_name')
PERMANOVA_fungi_jointing <- adonis2(fungi_distance_jointing ~ residue * fertilization , 
                                  fungi_data_jointing , permutations = 9999)
PERMANOVA_fungi_jointing
fungi_NMDS_jointing <- metaMDS(fungi_distance_jointing, k = 2)
fungi_NMDS_jointing$stress
fungi_NMDS_jointing_points <- as.data.frame(fungi_NMDS_jointing$points) %>% 
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by='sample_name')
fungi_NMDS_jointing_p <- ggplot(fungi_NMDS_jointing_points ,aes(x=MDS1,y=MDS2))+
  geom_point(aes(shape=fertilization,color=residue),size=3.5)+
  scale_color_manual(values = CAOXINZHUANG2)+
  scale_shape_manual(values = c(17,15,16,21))+
  labs(x= 'Jointing Stage',y=NULL)+
  geom_text(x = min(fungi_NMDS_jointing_points$MDS1) + 0.35,
            y = min(fungi_NMDS_jointing_points$MDS2) +0.05 ,size = 6,
            label = paste0('Stress = ',round(fungi_NMDS_jointing$stress , digits = 3)))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme_NMDS
fungi_NMDS_jointing_p
#Filling Stage
fungi_distance_filling <- fungi_asv_rari %>%
  select(starts_with('R')) %>%
  t() %>%
  vegdist(method = 'bray')

fungi_data_filling <- fungi_asv_rari %>%
  select(starts_with('R')) %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by = 'sample_name')
PERMANOVA_fungi_filling <- adonis2(fungi_distance_filling ~ residue * fertilization , 
                                 fungi_data_filling , permutations = 9999)
PERMANOVA_fungi_filling
fungi_NMDS_filling <- metaMDS(fungi_distance_filling, k = 2)
fungi_NMDS_filling$stress
fungi_NMDS_filling_points <- as.data.frame(fungi_NMDS_filling$points) %>% 
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by='sample_name')
fungi_NMDS_filling_p <- ggplot(fungi_NMDS_filling_points ,aes(x=MDS1,y=MDS2))+
  geom_point(aes(shape=fertilization,color=residue),size=3.5)+
  scale_color_manual(values = CAOXINZHUANG2)+
  scale_shape_manual(values = c(17,15,16,21))+
  labs(x='Filling Stage',y=NULL)+
  geom_text(x = min(fungi_NMDS_filling_points$MDS1) + 0.35,
            y = min(fungi_NMDS_filling_points$MDS2) +0.05 ,size = 6,
            label = paste0('Stress = ',round(fungi_NMDS_filling$stress , digits = 3)))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme_NMDS
fungi_NMDS_filling_p
#Maturing Stage
fungi_distance_maturing <- fungi_asv_rari %>%
  select(starts_with('m')) %>%
  t() %>%
  vegdist(method = 'bray')

fungi_data_maturing <- fungi_asv_rari %>%
  select(starts_with('m')) %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by = 'sample_name')
PERMANOVA_fungi_maturing <- adonis2(fungi_distance_maturing ~ residue * fertilization , 
                                  fungi_data_maturing , permutations = 9999)
PERMANOVA_fungi_maturing
fungi_NMDS_maturing <- metaMDS(fungi_distance_maturing, k = 2)
fungi_NMDS_maturing$stress
fungi_NMDS_maturing_points <- as.data.frame(fungi_NMDS_maturing$points) %>% 
  rownames_to_column(var = 'sample_name') %>%
  left_join(.,group,by='sample_name')
fungi_NMDS_maturing_p <- ggplot(fungi_NMDS_maturing_points ,aes(x=MDS1,y=MDS2))+
  geom_point(aes(shape=fertilization,color=residue),size=3.5)+
  scale_color_manual(values = CAOXINZHUANG2)+
  scale_shape_manual(values = c(17,15,16,21))+
  labs(x = 'Maturing Stage',y=NULL)+
  geom_text(x = min(fungi_NMDS_maturing_points$MDS1) + 0.25,
            y = min(fungi_NMDS_maturing_points$MDS2) +0.05 ,size = 6,
            label = paste0('Stress = ',round(fungi_NMDS_maturing$stress , digits = 3)))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme_NMDS
fungi_NMDS_maturing_p

#####output#####
PERMANOVA_rp
PERMANOVA_rp <- rbind(PERMANOVA_bac_jointing,PERMANOVA_bac_filling,PERMANOVA_bac_maturing,
                      PERMANOVA_fungi_jointing,PERMANOVA_fungi_filling,PERMANOVA_fungi_maturing) %>%
  select(R2,`Pr(>F)`) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'factor') %>%
  filter(str_detect(factor, 'e'),!str_detect(factor, 'Residual')) %>%
  cbind(.,data.frame('Species' = rep(c('Bacteria','Fungi'),each = 9),
                     'Stage' = rep(rep(c('Jointing','Filling','Maturing'),each = 3),n = 2),
                     'Factor' = rep(c('R','F','R*F' ), n = 6))) %>%
  .[,-1] %>%
  mutate(label = if_else(`Pr(>F)` < 0.001,"***",
                         if_else(`Pr(>F)` <0.01,"**",
                                 if_else(`Pr(>F)`<0.05,"*",""))))
PERMANOVA_rp$Stage <- factor(PERMANOVA_rp$Stage,levels = c('Jointing','Filling','Maturing'))
PERMANOVA_rp$Factor <- factor(PERMANOVA_rp$Factor,levels = c('R','F','R*F' ))
PERMANOVA_rp_p <- ggplot(PERMANOVA_rp)+
  geom_col(aes(x = Stage, y = R2,fill = Factor),color = 'grey',
           position = 'dodge')+
  scale_fill_manual(values =  c("#6BB88A","#99749F","#F0BB47"))+
  geom_text(aes(x = Stage,y = R2 + 0.01,label = label),size = 5,
            position = position_dodge2(width = 1))+
  labs(y = NULL , title = 'PERMANOVA R^2')+
  facet_wrap(Species~.,nrow = 2,strip.position = 'left',scales = 'free_y')+
  theme_NMDS+
  theme(axis.text.x = element_text(size = 17),axis.title = element_text(size = 25),
        legend.position = 'none')
PERMANOVA_rp_p
nmds_p <- ggarrange(ggarrange(bac_NMDS_jointing_p,bac_NMDS_filling_p,bac_NMDS_maturing_p,
                              fungi_NMDS_jointing_p,fungi_NMDS_filling_p,
                              fungi_NMDS_maturing_p,nrow = 2,ncol = 3),
                                ncol = 2,widths = c(12,3.5),PERMANOVA_rp_p)
nmds_p
ggsave(PERMANOVA_rp_p,filename='caoxinzhuang_results/2.0/PERMANOVA_rp_p.pdf',width = 5,height = 9)
# ggsave(nmds_p,filename='caoxinzhuang_results/2.0/beta.jpg',width = 18,height = 9)
# ggsave(nmds_p,filename='caoxinzhuang_results/2.0/beta.pdf',width = 18,height = 9)

####community composition####
#####bacteria#####
bac_relative_table_t <- t(bac_asv) / colSums(bac_asv) 
bac_relative_table <- t(bac_relative_table_t) %>% as.data.frame()

bac_composition_raw <- bac_relative_table %>%
  cbind(phylum = taxa_bac$Phylum) %>%
  pivot_longer(!phylum , names_to = 'sample_name' , values_to = 'relative_abundance') %>%
  group_by(sample_name , phylum) %>%
  summarise(relative = sum(relative_abundance)) %>%
  left_join(group , 'sample_name') 

bac_composition_stage_residue <- bac_composition_raw %>%
  group_by(stage , residue, phylum) %>%
  summarise(mean_abundance = mean(relative ),
            sd_abundance = sd(relative))
bac_composition_stage_residue
top_bac_phylum <- bac_composition_stage_residue %>%
  group_by(phylum) %>%
  summarise(abundance = sum(mean_abundance),
            mean = mean(mean_abundance)) %>%
  arrange(desc(mean))
tail(top_bac_phylum)
for(i in 1:12) {
  p = 1 - top_bac_phylum$mean[[i]]
  print(p)
  print(i)
}

bac_phylum_ttest <- bac_composition_raw %>%
  filter(phylum %in% top_bac_phylum$phylum[1:10] ) %>%
  group_by(stage, phylum) %>%
  t_test(relative~residue) %>%
  mutate(padj = p.adjust(.$p,method = 'BH')) %>%
  mutate(label = if_else(padj < 0.001,"***",
                         if_else(padj <0.01,"**",
                                 if_else(padj<0.05,"*",""))),
         sign = paste0(phylum,stage)) %>%
  select('padj','label','sign')
bac_phylum_ttest
bac_community_composition_phylum_p <- bac_composition_stage_residue %>%
  filter(phylum %in% top_bac_phylum$phylum[1:10]) %>%
  mutate(phylum2 = factor(phylum,levels = top_bac_phylum$phylum[1:10]),
         sign = paste0(phylum,stage)) %>%
  left_join(bac_phylum_ttest,by = 'sign') %>%
  ggplot(aes(x = phylum2, y = mean_abundance ))+
  geom_col(aes(fill = residue),position = 'dodge',col = 'grey',width = 0.7)+
  labs(x = 'Bacteria Phylum' ,y = 'Relative Abundance')+
  geom_errorbar(aes(ymin = mean_abundance, ymax = mean_abundance + sd_abundance),
                position = position_dodge2(1),width = 0.7,linewidth = 0.8)+
  scale_fill_manual(values = CAOXINZHUANG2)+
  geom_text(aes(x = phylum , y = 0.65,label = label),size = 5)+
  facet_grid(stage~.)+
  theme_NMDS+
  theme(axis.text.x = element_text(size = 12,angle = 340,vjust = 0.1,hjust = 0.2))

bac_community_composition_phylum_p2 <- bac_composition_stage_residue %>%
  filter(phylum %in% top_bac_phylum$phylum[5:10]) %>%
  mutate(phylum2 = factor(phylum,levels = top_bac_phylum$phylum[5:10]),
         sign = paste0(phylum,stage)) %>%
  left_join(bac_phylum_ttest,by = 'sign') %>%
  ggplot(aes(x = phylum2, y = mean_abundance ))+
  geom_col(aes(fill = residue),position = 'dodge',col = 'grey',width = 0.7)+
  labs(x = 'Bacteria Phylum' ,y = 'Relative Abundance')+
  geom_errorbar(aes(ymin = mean_abundance, ymax = mean_abundance + sd_abundance),
                position = position_dodge2(1),width = 0.7,linewidth = 0.8)+
  scale_fill_manual(values = CAOXINZHUANG2)+
  geom_text(aes(x = phylum , y = 0.08,label = label),size = 5)+
  facet_grid(stage~.)+
  theme_NMDS+
  theme(axis.text.x = element_text(size = 12,angle = 340,vjust = 0.1,hjust = 0.2))
bac_community_composition_phylum_p2
#####fungi#####
fungi_relative_table_t <- t(fungi_asv) / colSums(fungi_asv)
fungi_relative_table <- t(fungi_relative_table_t) %>% as.data.frame()

fungi_composition_raw <- fungi_relative_table %>%
  cbind(phylum = taxa_fungi$Phylum) %>%
  pivot_longer(!phylum , names_to = 'sample_name' , values_to = 'relative_abundance') %>%
  group_by(sample_name , phylum) %>%
  summarise(relative = sum(relative_abundance)) %>%
  left_join(group , 'sample_name') 

fungi_composition_stage_residue <- fungi_composition_raw %>%
  group_by(stage , residue, phylum) %>%
  summarise(mean_abundance = mean(relative ),
            sd_abundance = sd(relative))

top_fungi_phylum <- fungi_composition_stage_residue %>%
  group_by(phylum) %>%
  summarise(abundance = sum(mean_abundance),
            mean = mean(mean_abundance)) %>%
  arrange(desc(mean))
top_fungi_phylum
for(i in 1:10) {
  p = 1 - top_fungi_phylum$mean[[i]]
  print(p)
  print(i)
}

fungi_phylum_ttest <- fungi_composition_raw %>%
  filter(phylum %in% top_fungi_phylum$phylum[1:3] ) %>%
  group_by(stage, phylum) %>%
  t_test(relative~residue) %>%
  mutate(padj = p.adjust(.$p,method = 'BH')) %>%
  mutate(label = if_else(padj < 0.001,"***",
                         if_else(padj <0.01,"**",
                                 if_else(padj<0.05,"*",""))),
         sign = paste0(phylum,stage)) %>%
  select('padj','label','sign')
top_fungi_phylum$phylum
fungi_community_composition_phylum_p <- fungi_composition_stage_residue %>%
  filter(phylum %in% top_fungi_phylum$phylum[1:3]) %>%
  mutate(phylum2 = factor(phylum,levels = top_fungi_phylum$phylum[1:3]),
         sign = paste0(phylum,stage)) %>%
  left_join(fungi_phylum_ttest,by = 'sign') %>%
  mutate(phylum3 = sapply(strsplit(phylum,split = '_'),'[',3)) %>%
  ggplot(aes(x = phylum3, y = mean_abundance )) +
  geom_col(aes(fill = residue),position = position_dodge2(),col = 'grey',width = 0.7) +
  labs(x = 'Fungi Phylum' ,y = NULL) +
  geom_errorbar(aes(ymin = mean_abundance, ymax = mean_abundance + sd_abundance),
                position = position_dodge2(1),width = 0.7,linewidth =0.8) +
  scale_fill_manual(values = CAOXINZHUANG2) +
  geom_text(aes(x = phylum3 , y = 0.65,label = label),size = 6) +
  facet_grid(stage~.) +
  theme_NMDS +
  theme(axis.text.x = element_text(size = 12,angle = 340,vjust = 0.1,hjust = 0.2))
fungi_community_composition_phylum_p

#####protists#####
head(protists_asv)
protists_relative_table_t <- t(protists_asv) / colSums(protists_asv)
protists_relative_table <- t(protists_relative_table_t) %>% as.data.frame()
unique(taxa_protists$Genus)
protists_composition_raw <- protists_relative_table %>%
  cbind(genus = taxa_protists$Genus) %>%
  mutate(genus = ifelse(is.na(genus) == T,'Na',genus)) %>%
  pivot_longer(!genus , names_to = 'sample_name' , values_to = 'relative_abundance') %>%
  group_by(sample_name , genus) %>%
  summarise(relative = sum(relative_abundance)) %>%
  left_join(group , 'sample_name') 
protists_composition_raw
protists_composition_stage_residue <- protists_composition_raw %>%
  group_by(stage , residue, genus) %>%
  summarise(mean_abundance = mean(relative ),
            sd_abundance = sd(relative))

top_protists_genus <- protists_composition_stage_residue %>%
  group_by(genus) %>%
  summarise(abundance = sum(mean_abundance),
            mean = mean(mean_abundance)) %>%
  arrange(desc(mean))
top_protists_genus
for(i in 1:20) {
  p = 1 - top_protists_genus$mean[[i]]
  print(p)
  print(i)
}

protists_genus_ttest <- protists_composition_raw %>%
  filter(genus %in% top_protists_genus$genus[1:20] ) %>%
  group_by(stage, genus) %>%
  t_test(relative~residue) %>%
  mutate(padj = p.adjust(.$p,method = 'BH')) %>%
  mutate(label = if_else(padj < 0.001,"***",
                         if_else(padj <0.01,"**",
                                 if_else(padj<0.05,"*",""))),
         sign = paste0(genus,stage)) %>%
  select('padj','label','sign')
protists_genus_ttest
top_protists_genus$genus
protists_community_composition_genus_p <- protists_composition_stage_residue %>%
  filter(genus %in% top_protists_genus$genus[1:20]) %>%
  mutate(genus2 = factor(genus,levels = top_protists_genus$genus[1:20]),
         sign = paste0(genus,stage)) %>%
  left_join(protists_genus_ttest,by = 'sign') %>%
  ggplot(aes(x = genus2, y = mean_abundance )) +
  geom_col(aes(fill = residue),position = position_dodge2(),col = 'grey',width = 0.7) +
  labs(x = 'protists genus' ,y = NULL) +
  geom_errorbar(aes(ymin = mean_abundance, ymax = mean_abundance + sd_abundance),
                position = position_dodge2(1),width = 0.7,linewidth =0.8) +
  scale_fill_manual(values = CAOXINZHUANG2) +
  geom_text(aes(x = genus2 , y = 0.28,label = label),size = 6) +
  facet_grid(stage~.) +
  theme_NMDS +
  theme(axis.text.x = element_text(size = 12,angle = 340,vjust = 0.1,hjust = 0.2))
protists_community_composition_genus_p
#####output#####
community_composition_phylum_p <- ggarrange(bac_community_composition_phylum_p,
                                            fungi_community_composition_phylum_p,
                                            heights = 6,widths = c(9,4))
community_composition_phylum_p
ggsave(bac_community_composition_phylum_p2,filename='caoxinzhuang_results/2.0/bac_community_composition_phylum_p2.pdf',width = 6,height = 6)

# ggsave(community_composition_phylum_p,filename='caoxinzhuang_results/2.0/community_composition_phylum_p.pdf',width = 16,height = 6)
# ggsave(community_composition_phylum_p,filename='caoxinzhuang_results/2.0/community_composition_phylum_p.jpg',width = 16,height = 6)
