library(tidyverse)
library(ggpubr)
# library(ggbeeswarm)
library(gghalves)
library(car)
library(vegan)
library(rstatix)
library(reshape2)
library(Hmisc)

####data input------------------------------
asv_bac <- read.csv('sequence/bac_ASV.csv' , row.names = 1 , header = T)
asv_fungi <- read.csv('sequence/fungi_ASV.csv',row.names = 1,header = T)
asv_proto <- read.csv('sequence/proto_ASVs.csv',row.names = 1,header = T) %>%
  filter(Kingdom == "Eukaryota" &(Phylum %in% c("Amoebozoa",'Excavata','Obazoa', "TSAR"))) %>%
  select(!c(v.R.C.R.5,R.R.D.R.5))

group <- read.delim('sequence/caoxinzhuang_group.txt' ,row.names = 1 , header = T)
taxa_bac <- asv_bac[,145:151] 
taxa_fungi <- asv_fungi[,145:151]
taxa_protists <- asv_proto[,143:149]
bac_asv <- asv_bac[,1:144] 
fungi_asv <- asv_fungi[,1:144]
protists_asv <- asv_proto[,1:142]

cols_composition <- c('#3AA9AE', '#E67D3A', '#EC9A5F', '#D46A2E', '#F2B075', '#C2581F', '#FFC38D',
                      '#B688CE', '#A274C0', '#C79EDD', '#9465B2', '#D8B5EB', '#8053A4',
                      '#E9CDFF', '#6E4196', '#B28FD6', '#8B5FB9', '#D1A3F0', '#9E77CC','grey80')
####genus level relative abundance--------------------------------------------
#####protist function annotation----------------------------------------------
ref_df <- read.delim("C:/a_usefulness/nwafu/database/protist/protist_functions.txt", sep = "\t", stringsAsFactors = FALSE)
head(ref_df)
ref_table <- ref_df %>%
  filter( Function.groups!= "#N/A" & !is.na(Function.groups)) %>%
  select(taxon, Function.groups) %>%
  distinct()
head(taxon_function)
taxon_function <- setNames(ref_table$Function.groups, ref_table$taxon)
hierarchy <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")

annotate_function <- function(row) {
  for(level in hierarchy) {
    taxon <- row[[level]]
    if(!is.na(taxon) && taxon %in% names(taxon_function)) {
      return(taxon_function[[taxon]])
    }
  }
  return(NA) 
}
protist_table_functions <- asv_df %>%
  mutate(Function = apply(., 1, annotate_function))
protist_table_functions

protist_table_functions_asv <- protist_table_functions %>%
  select(Kingdom:ncol(.)) %>%
  rownames_to_column(var = 'asv')
head(protist_table_functions_asv)
unique(protist_table_functions_asv$Family)
#####protist_function_relative_abundance----------------------------------------
protist_relative_table_t <- t(protist_table_functions2[,-ncol(protist_table_functions2)]) / colSums(protist_table_functions2[,-ncol(protist_table_functions2)]) 

protist_relative_table_functions <- t(protist_relative_table_t) %>%
  as.data.frame() %>%
  cbind(functions = protist_table_functions$Function) %>%
  pivot_longer(!functions,names_to = 'sample_name') %>%
  group_by(sample_name,functions) %>%
  summarise(sum_value = sum (value))

protist_genus_relative_abundance <- t(protist_relative_table_t) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'asv') %>%
  pivot_longer(!asv,names_to = 'sample_name') %>%
  left_join(group,by = 'sample_name') %>%
  group_by(residue,asv) %>%
  summarise(mean_relative = mean(value)) %>%
  pivot_wider(names_from = residue,values_from = mean_relative) %>%
  arrange(desc(Turnover)) %>%
  left_join(protist_table_functions_asv,by = 'asv')

protist_genus_relative_abundance %>%
  filter(asv %in% enrich_protist_asv_id$asv_id)
enrich_protist_asv_id

protist_relative_table_Phagotrophs <- protist_relative_table_functions %>%
  filter(functions=='Phagotrophs') %>%
  left_join(group,by = 'sample_name') %>%
  arrange(stage,residue,fertilization)

 protist_relative_table_Phagotrophs %>%
   group_by(stage) %>%
  t_test(sum_value ~ residue,data = .)

protist_function_composition <- protist_relative_table_functions %>%
  left_join(group,by = 'sample_name') %>%
  arrange(stage,residue,fertilization) %>%
  group_by(functions,stage,residue,fertilization) %>%
  summarise(mean_value = mean(sum_value)) %>%
  mutate(group222 = paste(fertilization,residue,sep = '_'),
         functions222 = if_else(functions=='','NA',functions)) %>%
  ggplot(aes(x = group222, y = mean_value, fill = functions)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c( '#E67D3A','#B688CE','#3AA9AE', '#F2B075', '#B28FD6', '#8053A4'))+
  labs(x = NULL,y = 'Mean Relative Abundance')+
  facet_grid(. ~ stage)+
  theme_NMDS2+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_blank(),
        strip.text = element_blank())
ggsave(protist_function_composition,filename = 'caoxinzhuang_results/6.0emf/protist_function_composition.pdf',
       height = 6.18,width = 10)
#####bacteria Genus relative abundance protist----------------------------------
bac_relative_table_t <- t(bac_asv) / colSums(bac_asv) 
bac_relative_table <- t(bac_relative_table_t) %>% as.data.frame()

bac_composition_raw <- bac_relative_table %>%
  cbind( taxa_bac) %>%
  mutate(Genussss = paste(Phylum,Family,Genus,sep = '_aaa_'),
         Genus2= ifelse(is.na(Genus)==TRUE,paste('AA_aaa_AA_aaa',Phylum,'Unknown',sep ='_'),Genussss)) %>%
  select(!c(ASV_ID:Genussss)) %>%
  pivot_longer(!Genus2, names_to = 'sample_name' , values_to = 'relative_abundance') %>%
  group_by(sample_name,Genus2) %>%
  summarise(relative = sum(relative_abundance)) 

top_bac_Genus <- bac_composition_raw %>%
  group_by(Genus2) %>%
  summarise(abundance = sum(relative),
            mean = mean(relative)) %>%
  arrange(desc(mean))

bac_composition_wider <- bac_composition_raw %>%
  mutate(Genus2 = ifelse(Genus2 %in% top_bac_Genus$Genus2[1:19],Genus2,'Others'),
         Genus3 = ifelse(Genus2 %in% c('Others','No'),'Others',Genus2)) %>%
  group_by(sample_name,Genus3) %>%
  summarise(relative_abundance = sum(relative)) %>%
  as.data.frame() %>%
  pivot_wider(names_from = 'Genus3',values_from = 'relative_abundance') %>%
  filter(!sample_name %in% c('v.R.C.R.5','R.R.D.R.5'))

bac_Genus_relative_abundance_Phagotrophs <- bac_composition_wider %>%
  left_join(protist_relative_table_Phagotrophs,by = 'sample_name')
bac_Genus_relative_abundance_Phagotrophs
bac_composition_wider
 bac_Genus_relative_abundance_Phagotrophs_genus_relative_abundance_long <- 
  bac_composition_wider %>%
  left_join(group,by = 'sample_name')  %>%
  pivot_longer(!colnames(group)) %>%
  group_by(stage,residue,fertilization,name) %>%
  summarise(mean_relative_abundance  = mean(value)) %>%
    mutate(group222 = paste(fertilization,residue,sep = '_'),
           final_col = ifelse(name == "Others", "Others",str_split(name, "_aaa_",simplify = TRUE)[,3]),
           namesss = final_col %>% 
             fct_inorder() %>%   
             fct_relevel("Others", after = Inf) )

  ggplot(bac_Genus_relative_abundance_Phagotrophs_genus_relative_abundance_long, 
       aes(x = group222, y = mean_relative_abundance, fill = namesss)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(cols_composition,'red'))+
  labs(x = NULL,y = 'Mean Relative Abundance')+
  facet_grid(. ~ stage)+
   theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_blank(),
         strip.text = element_blank())
bac_Genus_relative_abundance_Phagotrophs_genus_relative_abundance_long
bac_genus_relative_abundance_p
 #####fungi Genus relative abundance protist----------------------------------
 fungi_relative_table_t <- t(fungi_asv) / colSums(fungi_asv) 
 fungi_relative_table <- t(fungi_relative_table_t) %>% as.data.frame()

 fungi_composition_raw <- fungi_relative_table %>%
   cbind( taxa_fungi) %>%
   mutate(Genussss = paste(Phylum,Family,Genus,sep = '_a_'),
          Genus2= ifelse(is.na(Genus)==TRUE,'No',Genussss)) %>%
   select(!c(Kingdom:Genussss)) %>%
   pivot_longer(!Genus2, names_to = 'sample_name' , values_to = 'relative_abundance') %>%
   group_by(sample_name , Genus2) %>%
   summarise(relative = sum(relative_abundance)) 
 
 top_fungi_Genus <- fungi_composition_raw %>%
   group_by(Genus2) %>%
   summarise(abundance = sum(relative),
             mean = mean(relative)) %>%
   arrange(desc(mean))
 top_fungi_Genus[1:19,]
 fungi_composition_wider <- fungi_composition_raw %>%
   mutate(Genus2 = ifelse(Genus2 %in% top_fungi_Genus$Genus2[1:19],Genus2,'Others'),
          Genus3 = ifelse(Genus2 %in% c('Others','No'),'Others',Genus2)) %>%
   group_by(sample_name,Genus3) %>%
   summarise(relative_abundance = sum(relative)) %>%
   as.data.frame() %>%
   pivot_wider(names_from = 'Genus3',values_from = 'relative_abundance') %>%
   filter(!sample_name %in% c('v.R.C.R.5','R.R.D.R.5'))
 
 fungi_Genus_relative_abundance_Phagotrophs <- fungi_composition_wider %>%
   left_join(protist_relative_table_Phagotrophs,by = 'sample_name')
 
 fungi_Genus_relative_abundance_Phagotrophs_genus_relative_abundance_long <- 
   fungi_composition_wider %>%
   left_join(group,by = 'sample_name')  %>%
   pivot_longer(!colnames(group)) %>%
   group_by(stage,residue,fertilization,name) %>%
   summarise(mean_relative_abundance  = mean(value)) %>%
   mutate(group222 = paste(fertilization,residue,sep = '_'),
          final_col = ifelse(name == "Others", "Others",str_split(name, "_a_",simplify = TRUE)[,3]),
          namesss = final_col %>% 
            fct_inorder() %>%   
            fct_relevel("Others", after = Inf) )

a <-  fungi_Genus_relative_abundance_Phagotrophs_genus_relative_abundance_long %>%
   group_by(residue,namesss) %>%
   summarise(mean = mean(mean_relative_abundance))

bac_Genus_relative_abundance_Phagotrophs_genus_relative_abundance_long %>%
   group_by(residue,namesss) %>%
 summarise(mean_va = mean(mean_relative_abundance)) %>%
  head(n = 20)

ggplot(fungi_Genus_relative_abundance_Phagotrophs_genus_relative_abundance_long, 
        aes(x = group222, y = mean_relative_abundance, fill = namesss)) +
   geom_col(position = "stack") +
   scale_fill_manual(values = c(
     '#B688CE', '#A274C0', '#C79EDD', '#9465B2','#D8B5EB', '#8053A4', '#E9CDFF',
     '#6E4196', '#B28FD6', '#8B5FB9', '#D1A3F0', '#9E77CC', '#E3BAF5', '#734AA8',
     '#C08AE2', '#A168BD','#E67D3A', '#EC9A5F','#3AA9AE', 
     'grey80'))+
   labs(x = NULL,y = 'Mean Relative Abundance')+
   facet_grid(. ~ stage)
 fungi_genus_relative_abundance_p
 bac_genus_relative_abundance_p

#####heatmap---------------------------------------------------------------------
genus_relative_abundance_phagotroph_function <- function(datadata,residueeeeee){
  dataaaa <- datadata %>%
    filter(residue==residueeeeee)
  genus_relative_abundance_Phagotrophs <- dataaaa %>%
    select(!c(sample_name,functions,stage,residue,fertilization)) %>%
    scale() %>%  as.matrix() %>% rcorr(. , type = 'spearman')
  dim(genus_relative_abundance_Phagotrophs$r)
  Phagotrophs__padj <- p.adjust(genus_relative_abundance_Phagotrophs$P,method = "BH") %>% matrix(ncol = 21)
  Phagotrophs_matrixxx <- genus_relative_abundance_Phagotrophs$r
  
  Phagotrophs_matrixxx[-21,21]
  star_Phagotrophs<- ifelse(Phagotrophs__padj < 0.001, "***",
                            ifelse(Phagotrophs__padj < 0.01, "**",
                                   ifelse(Phagotrophs__padj < 0.05, "*", "")))
  star_Phagotrophs[-21,21]
  final_matrix <- cbind(R2 = Phagotrophs_matrixxx[-21,21],star = star_Phagotrophs[-21,21]) %>%
    as.data.frame()%>%
    rownames_to_column(var = 'name') %>%
    mutate(final_col = ifelse(name == "Others", "Others",str_split(name, "_aaa_",simplify = TRUE)[,3]),
           namesss = final_col %>% 
             fct_inorder() %>%   
             fct_relevel("Others", after = Inf),residue=residueeeeee) %>%
    arrange(namesss)
  print(final_matrix)
}

bac_remove_heatmap_data <- genus_relative_abundance_phagotroph_function(bac_Genus_relative_abundance_Phagotrophs,'Remove') %>%
  mutate(speciess = 'Bacteria')
bac_turnover_heatmap_data <- genus_relative_abundance_phagotroph_function(bac_Genus_relative_abundance_Phagotrophs,'Turnover')%>%
  mutate(speciess = 'Bacteria')
fungi_remove_heatmap_data <- genus_relative_abundance_phagotroph_function(fungi_Genus_relative_abundance_Phagotrophs,'Remove')%>%
  mutate(speciess = 'Fungi')
fungi_turnover_heatmap_data <- genus_relative_abundance_phagotroph_function(fungi_Genus_relative_abundance_Phagotrophs,'Turnover')%>%
  mutate(speciess = 'Fungi')

bac_turnover_heatmap_data
bac_heatmap_Phagotrophs <- rbind(bac_remove_heatmap_data,bac_turnover_heatmap_data) %>%
  as.data.frame() %>%
  mutate(R222 = as.numeric(R2),namesss2 = namesss %>% 
           fct_inorder() %>%   
           fct_relevel("Others", after = Inf),
         genus = factor(namesss2,levels = rev(levels(namesss2)))) %>%
  ggplot(aes(x = residue, y = genus, fill = R222)) + 
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(round(R222, 2), star)), 
            size = 4, color = "black") + 
  scale_fill_gradient2(
    low = "#7D0B5D", mid = "white", high = "#428949",
    midpoint = 0)+
  labs(x = 'Residue',y = 'Genus')+
  theme_NMDS2+
  theme(legend.position = 'right')
bac_heatmap_Phagotrophs
fungi_heatmap_Phagotrophs <- rbind(fungi_remove_heatmap_data, fungi_turnover_heatmap_data) %>%
  as.data.frame() %>%
  mutate(R222 = as.numeric(R2),
    namesss2 = namesss %>%
      fct_inorder() %>% 
      fct_relevel("Others", after = Inf),
    genus = factor(namesss2, levels = rev(levels(namesss2))),
    genus222 = ifelse(
      genus == "Others",
      "Others",  
      str_split(genus, 'g__', simplify = TRUE)[,2]),
    genusss = factor(genus222,levels = c("Others",rev(unique(genus222[genus222 != "Others"]))))) %>%
  ggplot(aes(x = residue, y = genusss, fill = R222)) + 
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(round(R222, 2), star)), 
            size = 4, color = "black") + 
  scale_fill_gradient2(
    low = "#7D0B5D", mid = "white", high = "#428949",
    midpoint = 0)+
  labs(x = 'Residue',y = 'Genus')+
  theme_NMDS2+
  theme(legend.position = 'right')



relative_abundance_heatmap <- ggarrange(bac_genus_relative_abundance_p,bac_heatmap_Phagotrophs,
                                        fungi_genus_relative_abundance_p,fungi_heatmap_Phagotrophs,
                                        ncol = 2,nrow = 2,widths = c(10,7))
ggsave(relative_abundance_heatmap,filename = 'caoxinzhuang_results/6.0emf/genus_heatmap2.pdf',
       height = 9.27,width = 15)
