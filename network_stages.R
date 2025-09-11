library(Hmisc)
library(igraph)
library(tidyverse)
library(vegan)
library(ggpubr)
library(rstatix)
library(randomForest)
library(rfPermute)
# rm(list = ls())

# prepare data
asv_bac <- read.csv('sequence/bac_ASV.csv' , row.names = 1 , header = T)
asv_fungi <- read.csv('sequence/fungi_ASV.csv',row.names = 1,header = T)
# asv_proto <- read.csv('sequence/proto_ASVs.csv',row.names = 1,header = T) %>%
#   filter(Kingdom == "Eukaryota" &(Phylum %in% c("Amoebozoa", "Excavata", "TSAR"))) %>%
#   select(!c(v.R.C.R.5,R.R.D.R.5))
group <- read.delim('sequence/caoxinzhuang_group.txt' ,row.names = 1 , header = T)
taxa_bac <- asv_bac[,145:151] 
taxa_fungi <- asv_fungi[,145:151]
taxa_proto <- asv_proto[,143:149]
bac_asv <- asv_bac[,1:144] 
fungi_asv <- asv_fungi[,1:144]
proto_asv <- asv_proto[,1:142]
group
group <- rownames_to_column(group , var = 'sample_name')
group$stage <- factor(group$stage , levels = c('Jointing','Filling','Maturity'))
group$fertilization <- factor(group$fertilization , 
                              levels = c('high_organic','high_fertilizer','mineral_fertilizer','ck'))
group
CAOXINZHUANG2 <- c('#B44847','#00A1D6')
species_col <- c('#8A37A6','#F2CF80','#50F2B0')
theme_NMDS <- theme_bw(base_size = 15)+
  theme(legend.position = 'none',panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22))

bacteria_network_data <- bac_asv %>%
  select(!c(v.R.C.R.5,R.R.D.R.5))%>%
  rownames_to_column(var = 'asv') %>%
  mutate(asv2 = paste0('B11b111',asv)) %>%
  select(!asv) %>%
  column_to_rownames(var = 'asv2')

bacteria_network_data_01_detect <- bacteria_network_data %>%
  apply(.,2,function(x)ifelse(x>0,1,0)) %>%
  as.data.frame()

fungi_network_data <- fungi_asv %>%
  select(!c(v.R.C.R.5,R.R.D.R.5))%>%
  rownames_to_column(var = 'asv') %>%
  mutate(asv2 = paste0('F11b111',asv)) %>%
  select(!asv) %>%
  column_to_rownames(var = 'asv2')

fungi_network_data_01_detect <- fungi_network_data %>%
  apply(.,2,function(x)ifelse(x>0,1,0)) %>%
  as.data.frame()

protozoa_network_data <- proto_asv %>%
  rownames_to_column(var = 'asv') %>%
  mutate(asv2 = paste0('P11b111',asv)) %>%
  select(!asv) %>%
  column_to_rownames(var = 'asv2')

protozoa_network_data_01_detect <- protozoa_network_data %>%
  apply(.,2,function(x)ifelse(x>0,1,0)) %>%
  as.data.frame()

# functiona codes-----------------------------------------------------------------------------------
#filter function for network
fil_low_ric_inner <- function(residueeeee){
  bacteria_network_data_01_detect_depth <- bacteria_network_data_01_detect %>%
    select(starts_with(residueeeee))  %>%
    mutate(detection = rowSums(.)/ncol(.)) %>%
    filter(detection>=0.3)
  bac_asv_network <- bacteria_network_data %>%
    select(starts_with(residueeeee))  %>%
    filter(rownames(.) %in% rownames(bacteria_network_data_01_detect_depth))
  fungi_network_data_01_detect_depth <- fungi_network_data_01_detect %>%
    select(starts_with(residueeeee))  %>%
    mutate(detection = rowSums(.)/ncol(.)) %>%
    filter(detection>=0.3)
  fungi_asv_network <- fungi_network_data %>%
    select(starts_with(residueeeee))  %>%
    filter(rownames(.) %in% rownames(fungi_network_data_01_detect_depth))
  protozoa_network_data_01_detect_depth <- protozoa_network_data_01_detect %>%
    select(starts_with(residueeeee))  %>%
    mutate(detection = rowSums(.)/ncol(.)) %>%
    filter(detection>=0.3)
  protozoa_asv_network <- protozoa_network_data %>%
    select(starts_with(residueeeee))  %>%
    filter(rownames(.) %in% rownames(protozoa_network_data_01_detect_depth))
  network_data_total <- rbind(bac_asv_network,fungi_asv_network,protozoa_asv_network)
  return(network_data_total)
}

r.cut <- 0.7
p.cut <- 0.001

network_calculateeee <- function(asv_tabbble){
  network_a <- asv_tabbble %>% t()  %>% as.matrix() %>% rcorr(. , type = 'spearman')
  network_ar<- network_a$r
  network_a_p.adj<- p.adjust(network_a$P,method = "BH") 
  network_ar[network_a_p.adj>p.cut|abs(network_ar)<r.cut] <- 0
  diag(network_ar) <- 0
  g_network<- graph.adjacency(network_ar,weighted=TRUE,mode="undirected") %>% simplify()
  g_network_ar.bulk<-delete.vertices(g_network,names(degree(g_network)[degree(g_network)==0]))
}

sub_network_properties <- list()
network_properties <- data.frame()
robust_total <- data.frame(timessss = c(1:100))

caoxinzhuang_network_calculate_stages <- function(residuessssss,residuess_ssss,stageeeee){
  asvvvvv <- fil_low_ric_inner(residuessssss)
  networkkk_total <- network_calculateeee(asv_tabbble = asvvvvv)
  left_nodes_number_total <- length(V(networkkk_total)$name)
  robust <- vector('numeric',100)
  START_TIME1 <- Sys.time()
  print(START_TIME1)
  for( i in 1:100) {
    asv_robust <- asvvvvv %>%
      filter(rownames(.) %in% sample_frac(data.frame(V(networkkk_total)$name),0.5)[[1]])
    network_rou<- network_calculateeee(asv_robust)
    left_nodes_number <- length(V(network_rou)$name)
    robust[[i]] <- left_nodes_number/left_nodes_number_total
  }
  robust2 <- as.data.frame(robust)
  colnames(robust2) <- paste(residuess_ssss,stageeeee,sep = '_a_a_a_')
  robust_total <<- cbind(robust_total,robust2)
  START_TIME2 <- Sys.time()
  print(START_TIME2)
  
  V(networkkk_total)$degree <- degree(networkkk_total)
  E(networkkk_total)$correlation <- E(networkkk_total)$weight
  E(networkkk_total)$cor <- ifelse(E(networkkk_total)$correlation >0,1,-1)
  E(networkkk_total)$weight <- abs(E(networkkk_total)$weight)
  #edge_list
  
  edge_g_network_ar.bulk <- data.frame(as_edgelist(networkkk_total))
  edge_g_network_ar.bulk_list <- data.frame(
    source = edge_g_network_ar.bulk[[1]],
    target = edge_g_network_ar.bulk[[2]],
    weight = E(networkkk_total)$weight,
    correlation = E(networkkk_total)$correlation,
    cor = E(networkkk_total)$cor) 
  source_species <- separate(edge_g_network_ar.bulk_list,
                             col = source,into = c('Species1','asv1'),sep = '11b111')
  target_species <- separate(edge_g_network_ar.bulk_list,
                             col = target,into = c('Species2','asv2'),sep = '11b111')
  linkages_composition <- cbind(source_species,target_species) %>%
    as.data.frame() %>%
    select(Species1,Species2)%>%
    mutate(linkage_type = ifelse(Species1==Species2,'WT','CT'),
           linkage_protist = ifelse(c(Species1!='P'&Species2=='P'),'Protist','other'))
  linkages_composition
  properties_value <- data.frame(nodes = vcount(networkkk_total),
                                 edges = ecount(networkkk_total),
                                 Cross_trophic = table(linkages_composition$linkage_type)[[1]],
                                 Within_trophic = table(linkages_composition$linkage_type)[[2]],
                                 protist_cross = table(linkages_composition$linkage_protist)[[2]],
                                 Modularity = modularity(cluster_fast_greedy(networkkk_total,modularity = T)),
                                 Clustering_coefficient = transitivity(networkkk_total),
                                 Average_path_length = mean_distance(networkkk_total,weight = NULL),
                                 Network_diameter = diameter(networkkk_total,weights=NA),
                                 Average_degree = mean(degree(networkkk_total)),
                                 Graph_density = edge_density(networkkk_total),
                                 residue = residuess_ssss,stage =stageeeee)
  
  modules <- cluster_fast_greedy(networkkk_total) 
  
  nodes_module <- membership(modules) %>% as.data.frame()
  
  node_g_network_ar.bulk_list <- data.frame(
    Id = V(networkkk_total)$name, degree = V(networkkk_total)$degree)
  
  node_g_network_ar.bulk_list2 <-  separate(node_g_network_ar.bulk_list,
                                            col = Id,into = c('Species','asv'),sep = '11b111')
  node_g_network_ar.bulk <- cbind(node_g_network_ar.bulk_list,module = nodes_module[[1]],
                                  species = node_g_network_ar.bulk_list2$Species) %>%
    as.data.frame() %>%
    mutate(Polygon = if_else(species=='B',0,if_else(species=='F',3,5)))
  degree1 <- data.frame(degree = V(networkkk_total)$degree,residue = residuess_ssss,stage=stageeeee)
  network_properties <<- rbind(network_properties,properties_value)
  write.table(edge_g_network_ar.bulk_list,
              file = paste0('caoxinzhuang_results/6.0emf/network/networkdata/',residuess_ssss,stageeeee,'_edge.csv'),
              sep = '\t', row.names = F, quote = F)
  write.table(node_g_network_ar.bulk,
              file = paste0('caoxinzhuang_results/6.0emf/network/networkdata/',residuess_ssss,stageeeee,'_node.csv'),
              sep = '\t', row.names = F, quote = F)
  # # sub networks calculate
  module_arrange <- node_g_network_ar.bulk$module %>% table()%>% as.data.frame() %>% arrange(desc(Freq))
  module_1 <- node_g_network_ar.bulk %>% filter(module==module_arrange$.[[1]])
  module_2 <- node_g_network_ar.bulk %>% filter(module==module_arrange$.[[2]])
  module_3 <- node_g_network_ar.bulk %>% filter(module==module_arrange$.[[3]])
  sample_names <- colnames(asvvvvv)
  print(sample_names)
  for (samppppleeee in sample_names) {
    sub_asvsss <- asvvvvv %>% filter(.[[samppppleeee]] > 0) %>% rownames()
    sub_nodes <- V(networkkk_total)[name %in% sub_asvsss]
    if (length(sub_nodes) < 2) {
      warning(paste("Skipping sample", samppppleeee, ": insufficient nodes (<2)"))
      next
    }
    sub_network <- induced_subgraph(networkkk_total, sub_nodes)
    E(sub_network)$weight <- abs(E(sub_network)$weight)
    sub_network <- induced_subgraph(networkkk_total,sub_nodes)
    E(sub_network)$weight <- abs(E(sub_network)$weight)
    edge_g_network_sub <- data.frame(as_edgelist(sub_network))
    edge_g_network_sub_list <- data.frame(
      source = edge_g_network_sub[[1]],
      target = edge_g_network_sub[[2]],
      weight = E(sub_network)$weight,
      correlation = E(sub_network)$correlation,
      cor = E(sub_network)$cor) 
    source_species_sub <- separate(edge_g_network_sub_list,
                                   col = source,into = c('Species1','asv1'),sep = '11b111')
    target_species_sub <- separate(edge_g_network_sub_list,
                                   col = target,into = c('Species2','asv2'),sep = '11b111')
    linkages_composition_sub <- cbind(source_species_sub,target_species_sub) %>%
      as.data.frame() %>%
      select(Species1,Species2)%>%
      mutate(linkage_type = ifelse(Species1==Species2,'WT','CT'),
             linkage_protist = ifelse(c(Species1!='P'&Species2=='P'),'Protist','other'))
    table_link_type <- table(linkages_composition_sub$linkage_type)
    WT_count <- ifelse("WT" %in% names(table_link_type), table_link_type[["WT"]], 0)
    CT_count <- ifelse("CT" %in% names(table_link_type), table_link_type[["CT"]], 0)
    table_link_protist <- table(linkages_composition_sub$linkage_protist)
    protist_count <- ifelse("Protist" %in% names(table_link_protist), table_link_protist[["Protist"]], 0)
    
    node_s_subnetwork <- data.frame(Id_sub = V(sub_network)$name)
    module_1_diversity <- node_s_subnetwork %>% filter(Id_sub %in% module_1$Id) %>% nrow()
    module_2_diversity <- node_s_subnetwork %>% filter(Id_sub %in% module_2$Id) %>% nrow()
    module_3_diversity <- node_s_subnetwork %>% filter(Id_sub %in% module_3$Id) %>% nrow()
    properties_value <- data.frame(
      nodes = vcount(sub_network),
      edges = ecount(sub_network),
      module1_diver = module_1_diversity,
      module2_diver = module_2_diversity,
      module3_diver = module_3_diversity,
      Cross_trophic = CT_count,  
      Within_trophic = WT_count, 
      protist_cross = protist_count,  
      Modularity = modularity(cluster_fast_greedy(sub_network, modularity = T)),
      Clustering_coefficient = transitivity(sub_network),
      Average_path_length = mean_distance(sub_network, weight = NULL),
      Network_diameter = diameter(sub_network, weights = NA),
      Average_degree = mean(degree(sub_network)),
      Graph_density = edge_density(sub_network),
      sample = samppppleeee
    )
    sub_network_properties[[samppppleeee]] <<- properties_value
  }
}
####network calculate-----------------------------------------------------------
sub_network_properties
caoxinzhuang_network_calculate_stages('v.R','harvest','Jointing')
caoxinzhuang_network_calculate_stages('v.T','return','Jointing')
caoxinzhuang_network_calculate_stages('R.R','harvest','Filling')
caoxinzhuang_network_calculate_stages('R.T','return','Filling')
caoxinzhuang_network_calculate_stages('m.R','harvest','Maturity')
caoxinzhuang_network_calculate_stages('m.T','return','Maturity')
sub_network_properties

sub_network_properties_total <- bind_rows(sub_network_properties)

sub_network_properties_total
# robust_total_emf <- robust_total
#####network robust-----------------------------------------------------------
# robust_p <- 
  robust_total_emf %>%
  pivot_longer(!timessss) %>%
  separate(col = name,into = c('straw','stage'),sep = '_a_a_a_') %>%
  group_by(straw,stage)%>%
  summarise(mean_value = mean(value),sd_value = sd(value)) %>%
  ggplot(aes(x = stage, y = mean_value,fill = straw))+
  geom_bar(stat = 'identity',width = 0.5,position = position_dodge2())+
  geom_errorbar(aes(x = stage, ymin = mean_value,ymax = mean_value+sd_value),position = position_dodge2()
                ,width = 0.3,size = 0.8)+
  scale_fill_manual(values = CAOXINZHUANG2)+
  labs(x = NULL,y = 'Robust')+
  theme_NMDS2+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
robust_p
robust_total_emf
#####network properties and EMF-----------------------------------------------------------
colnames(sub_network_properties_total)
network_emf <- sub_network_properties_total %>%
  as.data.frame() %>%
  left_join(group,by = c('sample'='sample_name')) %>%
  arrange(stage,residue,fertilization) %>%
  filter(str_detect(sample, "\\.(1|2|3|4)$"))  %>%
  mutate(emf = emf_total$meanFunction)
colnames(network_emf)

network_emf %>%
  group_by(stage) %>%
  t_test(nodes~residue,data =.)

######Random Forest-----------------------------------------------------------
network_emf_rfm_harvest <- network_emf %>%
  select(!c(nodes,edges,sample:fertilization,Clustering_coefficient)) %>%
  apply(., 2, scale_01_function)
set.seed(20)

rf_model1111111111111 <- rfPermute(
  emf~.,data = network_emf_rfm_harvest,ntree = 500,
  num.rep = 500 
)
importance(rf_model1111111111111, scale = TRUE)
summary(rf_model1111111111111)

importance(rf_model, scale = TRUE)  # 重要性得分
summary(rf_model)   

network_rfm_harvest_model <- randomForest(emf~.,data = network_emf_rfm_harvest,ntree = 5000,importance = TRUE)
print(network_rfm_harvest_model)
randomForest::importance(network_rfm_harvest_model) %>%
  as.data.frame() %>%rownames_to_column(var = 'indicator') %>%
  arrange(desc(`%IncMSE`))
rmf_network_total_importance <- randomForest::importance(network_rfm_harvest_model) %>%
  as.data.frame() %>%rownames_to_column(var = 'indicator') %>%
  arrange(desc(`%IncMSE`))%>%
  mutate(indicators = factor(indicator,levels = indicator),
         types = ifelse(indicators %in% c('protist_cross','Within_trophic','Cross_trophic'),'Interactions',
                        ifelse(indicators%in% c('module1_diver','module2_diver','module3_diver'),
                               'Module_Diversity','Complexity'))) 
rmf_network_total_importance
network_rfm_p <- ggplot(rmf_network_total_importance)+
  geom_bar(aes(x = indicators,y = `%IncMSE`,fill = types),stat = 'identity')+
  scale_fill_manual(values = c('#A8DDE3','#8D3C90','#FAEBBD'))+
  theme_NMDS2

Protist_corre_interactions_p <- 
  ggplot(network_emf,aes(x = stage,y = protist_cross,fill = residue))+
  geom_boxplot(position = position_dodge(0.9),alpha = 0.8,outliers = F,size = 0.6)+
  labs(x = 'Stage',y = 'Protist determined \nmulti-trophic interactions')+
  scale_fill_manual(values = CAOXINZHUANG2)+
  theme_NMDS2

Protist_corre_regression_p <- network_emf %>%
  ggplot(aes(x = protist_cross ,y = emf,col = residue))+
  geom_point()+
  stat_smooth(method = "lm", size = 1,se=FALSE)+
  scale_color_manual(values = CAOXINZHUANG2)+
  labs(x = 'Protist determined multi-trophic interactions',y = 'EMF')+
  theme_NMDS2
 
network_sub_plot <- ggarrange(network_rfm_p,Protist_corre_regression_p,Protist_corre_interactions_p,
                              ncol = 3,nrow = 1,widths = c(2,1,1),heights = 1)
ggsave(network_sub_plot,file = 'caoxinzhuang_results/6.0emf/network/network_emf.pdf',width = 12,height = 3)

Module1_emf <- network_emf %>%
  ggplot(aes(x = module1_diver,y = emf,col = residue))+
  geom_point()+
  stat_smooth(method = "lm", size = 1,se=FALSE)+
  scale_color_manual(values = CAOXINZHUANG2)+
  labs(x = 'Biodiversity',y = 'EMF',subtitle = 'Module 1')
Module1_emf
Module2_emf <- ggplot(network_emf,aes(x = module2_diver,y = emf,col = residue))+
  geom_point()+
  stat_smooth(method = "lm", size = 1,se=FALSE)+
  scale_color_manual(values = CAOXINZHUANG2)+
  labs(x = 'Biodiversity',y = 'EMF',subtitle = 'Module 2')+
  theme_NMDS2
Module2_emf
Module3_emf <- ggplot(network_emf,aes(x = module3_diver,y = emf,col = residue))+
  geom_point()+
  stat_smooth(method = "lm", size = 1,se=FALSE)+
  scale_color_manual(values = CAOXINZHUANG2)+
  labs(x = 'Biodiversity',y = 'EMF',subtitle = 'Module 3')+
  theme_NMDS2
Module3_emf

module_emf_p_total <- ggarrange(Module1_emf,Module2_emf,Module3_emf,robust_p,ncol = 4)
module_emf_p_total
ggsave(module_emf_p_total,file = 'caoxinzhuang_results/5.0/network/network_diversity_emf.pdf',width = 16,height = 4)
sub_network_properties_total[,-c(3:5)] %>%
  as.data.frame() %>%
  left_join(group,by = c('sample'='sample_name')) %>%
  arrange(stage,residue,fertilization) %>%
  pivot_longer(!c(sample,stage,residue,fertilization),names_to = 'indicator') %>%
  group_by(indicator) %>%
  t_test(value~residue,data = .)

network_properties_total <- 
  sub_network_properties_total[,-c(3:5)] %>%
  as.data.frame() %>%
  left_join(group,by = c('sample'='sample_name')) %>%
  arrange(stage,residue,fertilization) %>%
  pivot_longer(!c(sample,stage,residue,fertilization),names_to = 'indicator') %>%
  group_by(residue,indicator)%>%
  summarise(mean_value = mean(value),sd_value = sd(value)) %>%
  filter(indicator!='protist_cross') %>%
  ggplot(aes(x = residue, y = mean_value,fill = residue))+
  geom_bar(stat = 'identity',width = 0.5)+
  geom_errorbar(aes(x = residue, ymin = mean_value,ymax = mean_value+sd_value),width = 0.3,size = 0.8)+
  scale_fill_manual(values = CAOXINZHUANG2)+
  theme_NMDS2+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  facet_wrap(.~indicator,scales = 'free_y')
ggsave(network_properties_total,filename = 'caoxinzhuang_results/6.0emf/network_properties_total.pdf',
       height = 6.18,width = 10)
####nodes composition------------------------------------------------------------------------------
harvest_nodes <- read.csv('caoxinzhuang_results/5.0/intermediate/harvest_node.csv',sep = '\t')
retention_nodes <- read.csv('caoxinzhuang_results/5.0/intermediate/retention_node.csv',sep = '\t')

module_arrange_retention <- retention_nodes$module %>% table()%>% as.data.frame() %>% arrange(desc(Freq))
module_arrange_harvest<- harvest_nodes$module %>% table()%>% as.data.frame() %>% arrange(desc(Freq))
retention_module1_composition <- retention_nodes %>%
  filter(module==module_arrange_retention$.[[1]]) %>% 
  count(species) %>%
  as.data.frame() %>% mutate(n_pro  = n/sum(.[,2]),module = 1)
retention_module1_composition
retention_module2_composition <- retention_nodes %>%
  filter(module==module_arrange_retention$.[[2]]) %>% 
  count(species) %>%
  as.data.frame() %>% mutate(n_pro  = n/sum(.[,2]),module =2)
retention_module3_composition <- retention_nodes %>%
  filter(module==module_arrange_retention$.[[3]]) %>% 
  count(species) %>%
  as.data.frame() %>% mutate(n_pro  = n/sum(.[,2]),module = 3)
retention_module3_composition
harvest_module1_composition <- harvest_nodes %>%
  filter(module==module_arrange_harvest$.[[1]]) %>% 
  count(species) %>%
  as.data.frame() %>% mutate(n_pro  = n/sum(.[,2]),module = 1)
harvest_module1_composition
harvest_module2_composition <- harvest_nodes %>%
  filter(module==module_arrange_harvest$.[[2]]) %>% 
  count(species) %>%
  as.data.frame() %>% mutate(n_pro  = n/sum(.[,2]),module =2)
harvest_module2_composition
harvest_module3_composition <- harvest_nodes %>%
  filter(module==module_arrange_harvest$.[[3]]) %>% 
  count(species) %>%
  as.data.frame() %>% mutate(n_pro  = n/sum(.[,2]),module = 3)

retention_module_composition <- rbind(retention_module1_composition,retention_module2_composition,retention_module3_composition)
harvest_module_composition <- rbind(harvest_module1_composition,harvest_module2_composition,harvest_module3_composition)

ggplot(retention_module1_composition,aes(x = module,y = n_pro,fill = species))+
  geom_bar(stat="identity",position = "stack")