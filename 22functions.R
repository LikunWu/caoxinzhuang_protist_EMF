library(tidyverse)
library(DESeq2)
library(ggvenn)
library(ggpubr)
library(pheatmap)
library(readxl)
library(rstatix)

####FAPROTAX####
#####prepare data#####
bac_faprotax_raw <- asv_bac %>%
  rownames_to_column(var = 'asv') %>%
  mutate(Taxonomy = paste(Kingdom,Phylum,Class,Order,Family,Genus,Species,sep = ';')) %>%
  .[,-c(146:152)]

write.csv(bac_faprotax_raw,file = 'caoxinzhuang_results/1.0/intermediate/bac_faprotax.csv',quote = F,row.names = F)
# predict with FAPROTAX
bac_faprotax <- read.delim('caoxinzhuang_results/1.0/intermediate/caoxz.faprotax.txt',row.names = 1)

bac_faprotax_prepare <- bac_faprotax[,-1] %>%
  .[rowSums(.)>0,] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name') %>%
  left_join(group,by = 'sample_name') %>%
  pivot_longer(!colnames(group),names_to = 'process',values_to = 'abundance') %>%
  filter(process %in% c('methanotrophy','methanogenesis','nitrification','anammox','nitrate_denitrification',
                        'nitrite_denitrification','nitrous_oxide_denitrification','denitrification',
                        'cellulolysis','xylanolysis','chitinolysis','nitrogen_fixation','nitrate_ammonification',
                        'nitrite_respiration','ligninolysis','nitrate_respiration',
                        'nitrate_reduction','nitrogen_respiration','ureolysis')) 

bac_faprotax_p1 <- bac_faprotax_prepare%>%
  group_by(stage,residue,process) %>%
  summarise(mean_abundance = mean(abundance)) %>%
  mutate(group = paste(stage,residue,sep = '_'))
bac_faprotax_ttest <- bac_faprotax_prepare %>%
  group_by(stage,process) %>%
  t_test(abundance~residue,data = .) %>%
  mutate(padj = p.adjust(p,method = 'BH')) %>%
  mutate(sign = if_else(padj< 0.001,"***",
                        if_else(padj <0.01,"**",
                                if_else(padj<0.05,"*",""))))

bac_faprotax_p1$group <- factor(bac_faprotax_p1$group,levels = c(
  'Jointing_Remove','Jointing_Turnover','Filling_Remove',
  'Filling_Turnover','Maturing_Remove','Maturing_Turnover'
))
# 
# bac_faprotax2
# 
bac_enriched_faprotax_p <-
  ggplot(bac_faprotax2)+
  geom_point(aes(x=process ,y= group , size = mean_abundance,
                 fill = log10(mean_abundance+1)),shape = 21,color = 'grey') +
  labs(x=NULL,y=NULL) +
  scale_size_continuous(range = c(0,20))+
  theme_bw()+
  scale_fill_gradient(low = '#F3F1F4', high = '#095527')+
  theme(axis.text.x=element_text(angle=45,hjust = 1,vjust=1,size = 12),
        axis.text.y = element_text(size = 20),
        legend.position='right')
#### use wider table as input !!!
# !!!start from 'mutate...' %>%
#   .[,-c(1:2)] %>%
#   pivot_wider(names_from = 'group',values_from = 'mean_abundance')
# bac_faprotax3 <- bac_faprotax2[,-1]
# rownames(bac_faprotax3) <- bac_faprotax2[[1]]
# annotation_col <- data.frame(
#   Stage = rep(c('Jointing', 'Filling','Maturing'),each = 2),
#   Residue = rep(c('Remove','Turnover'),n = 3))
# rownames(annotation_col) <- colnames(bac_faprotax3)
# ann_colors <- list(Stage =c(Jointing="#f89820",Filling="#e32a49",Maturing ="#5cbcc8"),
#                    Residue = c(Remove = '#B44847',Turnover = '#00A1D6'))
# bac_faprotax_p <- log10(bac_faprotax3+1) %>%
#   pheatmap(cluster_rows = F,
#            cluster_cols = F,
#            color = colorRampPalette(c('#F3F1F4', '#095527'))(100),
#            show_colnames = F,
#            show_rownames = T,
#            angle_col = 0,
#            border_color = 'gray',
#            fontsize = 12,
#            legend = T,
#            annotation_colors = ann_colors,
#            annotation_col = annotation_col,
#            annotation_legend = F)
# bac_faprotax_p

#####fungal traits#####
# import database
fungi_traits <- read_excel('C:/a_usefulness/nwafu/database/fungaltraits-master/13225_2020_466_MOESM4_ESM.xlsx')
# select enriched fungi ASVs
fungi_enrich_fungltrait_data <- asv_fungi %>%
  rownames_to_column(var = 'asv') %>%
  mutate(GENUS = sapply(strsplit(Genus,split = '_'),'[',3))

left_join(fungi_enrich_fungltrait_data,fungi_traits , by = 'GENUS') %>%
  select(2:145,primary_lifestyle) %>%
  na.omit() %>%
  filter(primary_lifestyle != 'animal_parasite') %>%
  filter(primary_lifestyle != 'lichen_parasite') %>%
  filter(primary_lifestyle != 'algal_parasite') %>%
  filter(primary_lifestyle != 'nectar/tap_saprotroph') %>%
  select(!primary_lifestyle) %>%
  t() %>%
  as.data.frame() %>%
  rowSums(.)/rowSums(t(asv_fungi[,1:144])) 
fungi_enrich_fungltrait_data

fungi_enrich_fungltrait_prepare <- left_join(fungi_enrich_fungltrait_data,fungi_traits , by = 'GENUS') %>%
  select(2:145,primary_lifestyle) %>%
  pivot_longer(!primary_lifestyle,names_to = 'sample_name',values_to = 'abundance') %>%
  na.omit() %>%
  filter(primary_lifestyle != 'animal_parasite') %>%
  filter(primary_lifestyle != 'lichen_parasite') %>%
  filter(primary_lifestyle != 'algal_parasite') %>%
  filter(primary_lifestyle != 'nectar/tap_saprotroph') %>%
  left_join(group,by = 'sample_name')

fungi_enrich_fungltrait_prepare$primary_lifestyle %>% unique()

fungi_enrich_fungltrait_ttest <- fungi_enrich_fungltrait_prepare %>%
  group_by(primary_lifestyle, stage ) %>%
  t_test(abundance~residue,data=.)%>%
  mutate(padj = p.adjust(p,method = 'BH')) %>%
  mutate(sign = if_else(padj< 0.001,"***",
                        if_else(padj <0.01,"**",
                                if_else(padj<0.05,"*",""))))

fungi_enrich_fungltrait_1 <- fungi_enrich_fungltrait_prepare%>%
  group_by(primary_lifestyle,stage,residue) %>%
  summarise(abundance = sum(abundance)/24)%>%
  mutate(group = paste(stage,residue,sep = '_'))

fungi_enrich_fungltrait_1$group <- factor(fungi_enrich_fungltrait_1$group,levels = c(
  'Jointing_Remove','Jointing_Turnover','Filling_Remove',
  'Filling_Turnover','Maturing_Remove','Maturing_Turnover'
))

fungal_traits_p <-  ggplot(fungi_enrich_fungltrait_1)+
  geom_point(aes(x=primary_lifestyle,y=group,size = abundance,
                 fill = log10(abundance+1)),shape = 21,color = 'grey') +
  labs(x=NULL,y=NULL) +
  scale_size_continuous(range = c(0,20))+
  theme_bw()+
  scale_fill_gradient(low = '#F3F1F4', high = '#095527')+
  theme(axis.text.x=element_text(angle=45,hjust = 1,vjust=1,size = 12),
        axis.text.y = element_text(size = 20),
        legend.position='right')

fungal_traits_p
#####output#####
ggsave('caoxinzhuang_results/1.0/bac_enriched_faprotax_p.pdf',bac_enriched_faprotax_p,width = 12,height = 6)
ggsave('caoxinzhuang_results/1.0/fungal_traits_p.pdf',fungal_traits_p,width = 18,height = 6)
