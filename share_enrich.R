#load package
library(tidyverse)
library(DESeq2)
library(ggvenn)
library(ggpubr)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(rstatix)
library(ggsankey)
cols_venn <- c(Remove = '#B44847' , Turnover = '#00A1D6')
CAOXINZHUANG2 <- c('#B44847','#00A1D6')
sunsx <- c('#DCB46C','#EEBCBF','#A8DDE3','#F7E1AF','#466E88','#D6585B')
cols_sankey <- c('#D6585B','#DCB46C','#EEBCBF','#F7E1AF','#466E88','#2A9E89','#E56F51'
                 ,'#B44847','#00A1D6','#A8DDE3')

# prepare NULL# prepare data
asv_bac <- read.csv('sequence/bac_ASV.csv' , row.names = 1 , header = T)
asv_fungi <- read.csv('sequence/fungi_ASV.csv',row.names = 1,header = T)
group <- read.delim('sequence/caoxinzhuang_group.txt' ,row.names = 1 , header = T)

taxa_bac <- asv_bac[,145:151] 
taxa_fungi <- asv_fungi[,145:151]
bac_asv <- asv_bac[,1:144]
fungi_asv <- asv_fungi[,1:144]
# arrange group
group <- rownames_to_column(group , var = 'sample_name')
group$stage <- factor(group$stage , levels = c('Jointing','Filling','Maturity'))
group$fertilization <- factor(group$fertilization , 
                              levels = c('high_organic','high_fertilizer','mineral_fertilizer','ck'))
group$residue <- factor(group$residue,levels = c('Remove','Turnover'))
group
####select####
#bacteria
bac_relative_table_t <- t(bac_asv) / colSums(bac_asv) 
bac_relative_table <- t(bac_relative_table_t) %>% as.data.frame()

bac_asv_0.0001_list <- bac_relative_table %>%
  rownames_to_column(var = 'asv') %>%
  pivot_longer(!asv , names_to = 'sample_name' , values_to = 'abundance') %>%
  group_by(asv) %>%
  summarise(mean_asv = mean(abundance)) %>%
  filter(mean_asv >= 0.0001) # select ASVs with mean relative abundance >= 0.0001
bac_asv_0.0001_list
#fungi
fungi_relative_table_t <- t(fungi_asv) / colSums(fungi_asv)
fungi_relative_table <- t(fungi_relative_table_t) %>% as.data.frame()  

fungi_asv_0.0001_list <- fungi_relative_table %>%
  rownames_to_column(var = 'asv') %>%
  pivot_longer(!asv , names_to = 'sample_name' , values_to = 'abundance') %>%
  group_by(asv) %>%
  summarise(mean_asv = mean(abundance)) %>%
  filter(mean_asv >= 0.0001) # select ASVs with mean relative abundance >= 0.0001
fungi_asv_0.0001_list


####venn####
bac_fungi_venn <- function(asv1,stage){
  venn_data <- asv1 %>%
    rownames_to_column(var = 'asv') %>%
    filter(asv %in% bac_asv_0.0001_list[[1]]) %>%
    column_to_rownames(var = 'asv') %>%
    select( starts_with(stage)) %>%
    filter(.,rowSums(.)>0) %>%
    rownames_to_column(var = 'asv') %>%
    pivot_longer(!asv , names_to = 'sample_name',values_to = 'abundance') %>%
    left_join(group, by = 'sample_name') %>%
    group_by(asv,residue) %>%
    summarise(exi = sum(abundance)) %>%
    pivot_wider(names_from = residue,values_from = exi) 
  list(Remove = venn_data %>% filter(Remove!=0) %>% select(asv) %>% unlist(),
       Turnover = venn_data %>% filter(Turnover!=0) %>% select(asv) %>% unlist()) %>%
    ggvenn(fill_color = CAOXINZHUANG2,
           fill_alpha = 0.85,
           stroke_alpha = 0.1,
           stroke_color = 'grey',
           stroke_size = 1.5,
           text_size = 6.5,
           set_name_size = 0
    )
}
venn_bac1 <- bac_fungi_venn(bac_asv,'V')
venn_bac2 <- bac_fungi_venn(bac_asv,'R')
venn_bac3 <- bac_fungi_venn(bac_asv,'m')
venn_fungi1 <- bac_fungi_venn(fungi_asv,'V')
venn_fungi2 <- bac_fungi_venn(fungi_asv,'R')
venn_fungi3 <- bac_fungi_venn(fungi_asv,'m')
venn_total_p <- ggarrange(venn_bac1,venn_bac2,venn_bac3,venn_fungi1,venn_fungi2,venn_fungi3,ncol = 3,nrow = 2)
venn_total_p
# enrich
enrich_asv <- function(asv11,ss,sstage){
  asv11 %>%
    rownames_to_column(var = 'asv') %>%
    filter(asv %in% bac_asv_0.0001_list[[1]]) %>%
    column_to_rownames(var = 'asv') %>%
    select( starts_with(ss)) %>%
    DESeqDataSetFromMatrix(colData = filter(group , stage == sstage), design = ~ residue) %>%
    DESeq(fitType = 'mean',minReplicatesForReplace = 6, parallel = T) %>%
    results(contrast = c('residue', 'Turnover', 'Remove')) %>%
    as.data.frame()  %>%
    filter(abs(log2FoldChange)>=2 & padj <=0.05) %>% 
    mutate(stage = sstage)
}


#####bacteria#####
enrich_bac_jointing <- enrich_asv(bac_asv,'V','Jointing') 
enrich_bac_filling <- enrich_asv(bac_asv,'R','Filling')
enrich_bac_maturing <- enrich_asv(bac_asv,'m','Maturing')
enrich_bac <- rbind(rownames_to_column(enrich_bac_jointing,var = 'asv') , 
                    rownames_to_column(enrich_bac_filling,var = 'asv'), 
                    rownames_to_column(enrich_bac_maturing,var = 'asv')) 
enrich_bac$stage <- factor(enrich_bac$stage , 
                                    levels = c('Jointing','Filling','Maturing'))
enrichd_bac_ASVs_p <- enrich_bac %>%
  filter(log2FoldChange<10) %>%
  ggplot()+
  geom_jitter(aes(x = stage , y = log2FoldChange,
                  color = log2FoldChange),size = 3 ,
              width = 0.3,shape = 16 ,fill = 'orange') +
  labs(x = NULL , y = 'log2 FC' , title = 'Bacteria')+
  scale_color_gradient2(low = '#B44847',mid = '#F3F1F4',high = '#00A1D6')+
  theme_NMDS 

enrichd_bac_ASVs_p
enriched_bac_asv_name <- enrich_bac$asv %>%
  unique()
enriched_bac_asv_name

write.table(enriched_bac_asv_name,file = 'caoxinzhuang_results/2.0/intermediate/enrich_bac_count.txt',
            quote = F,row.names = F)
enrich_bac_count <- enrich_bac %>%
  mutate(enrich = if_else(log2FoldChange >=0 , 'Turnover','Remove'))
enriched_bac_asv_name

enrich_bac_count
bac_enrich_type <- enrich_bac_count %>%
  select('asv','stage','enrich') %>%
  unique()

bac_enrich_taxa <- taxa_bac %>%
  rownames_to_column(var = 'asv') %>%
  filter(asv %in% enriched_bac_asv_name) %>%
  left_join(bac_enrich_type,by = 'asv') %>%
  left_join(bac_asv_0.0001_list,by = 'asv') %>%
  select(!c(ASV_ID,Kingdom)) %>%
  arrange(stage,enrich)
bac_enrich_taxa %>%
  select(asv,mean_asv,enrich) %>%
  unique() %>%
  summarise(abundance = sum(mean_asv))

# bac_sankey_data <- bac_enrich_taxa %>%
#   select(enrich,Phylum,Genus) 
# write.csv(bac_enrich_taxa,file = 'caoxinzhuang_results/4.0/intermediate/bac_enrich_taxa.csv',quote = F)

# sankey_genus_count_bac <- table(bac_sankey_data$Genus) %>%
#   as.data.frame() %>%
#   left_join(bac_sankey_data,by = c('Var1' = 'Genus')) %>%
#   unique() %>%
#   relocate(Freq,.after = ncol(.))

sankey_genus_count_bac<- bac_enrich_taxa %>%
  na.omit() %>%
  make_long(enrich,Genus,Phylum)

sankey_genus_bac_p <-
  ggplot(sankey_genus_count_bac, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +
  geom_sankey(flow.alpha = 0.5, #条带不透明度
              smooth = 8, #条带弯曲度
              width = 0.12) + #节点宽度
  geom_sankey_text(size = 3.2, #标签字号
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none') #隐藏图例

sankey_genus_bac_p
#####fungi#####

enrich_asv_fungi <- function(asv11,ss,sstage){
    asv11 %>%
    rownames_to_column(var = 'asv') %>%
    filter(asv %in% fungi_asv_0.0001_list[[1]]) %>%
    column_to_rownames(var = 'asv') %>%
    select(starts_with(ss)) %>%
    DESeqDataSetFromMatrix(colData = filter(group ,stage == sstage), design = ~ residue) %>%
    DESeq(fitType = 'mean',minReplicatesForReplace = 6, parallel = F) %>%
    results(contrast = c('residue', 'Turnover', 'Remove')) %>%
    as.data.frame()  %>%
    filter(abs(log2FoldChange)>=2 & padj <=0.05) %>% 
    mutate(stage = sstage)
}

enrich_fungi_jointing <- enrich_asv_fungi(fungi_asv,'V','Jointing')
enrich_fungi_filling <- enrich_asv_fungi(fungi_asv,'R','Filling')
enrich_fungi_maturing <- enrich_asv_fungi(fungi_asv,'m','Maturing')

enrich_fungi <- rbind(rownames_to_column(enrich_fungi_jointing,var = 'asv') , 
                      rownames_to_column(enrich_fungi_filling,var = 'asv'), 
                      rownames_to_column(enrich_fungi_maturing,var = 'asv')) 
enrich_fungi$stage <-  factor(enrich_fungi$stage , 
                              levels = c('Jointing','Filling','Maturing'))
enriched_fungi_asv_name <- enrich_fungi$asv %>%
  unique()
str(enriched_fungi_asv_name)
enrich_fungi_count <- enrich_fungi %>%
  mutate(enrich = if_else(log2FoldChange >=0 , 'Turnover','Remove'))

enrichd_fungi_ASVs_p <- ggplot(enrich_fungi)+
  geom_jitter(aes(x = stage , y = log2FoldChange,
                  color = log2FoldChange),size = 3,
              width = 0.3,shape = 16 ,fill = 'orange') +
  labs(x = NULL , y = NULL , title = 'Fungi')+
  scale_color_gradient2(low = '#B44847',mid = '#F3F1F4',high = '#00A1D6')+
  theme_NMDS
enrichd_fungi_ASVs_p

enrichd_ASVs_p <- ggarrange(enrichd_bac_ASVs_p,enrichd_fungi_ASVs_p,widths = c(5.1,5),heights = 5,
                            nrow = 1,ncol = 2,align = 'v')
enrichd_ASVs_p
write.table(enriched_fungi_asv_name,file = 'caoxinzhuang_results/2.0/intermediate/enrich_fungi_count.txt',
            quote = F,row.names = F)


enrich_fungi_count
fungi_enrich_type <- enrich_fungi_count %>%
  select('asv','stage','enrich') %>%
  unique()

fungi_enrich_taxa <- taxa_fungi%>%
  rownames_to_column(var = 'asv') %>%
  filter(asv %in% enriched_fungi_asv_name) %>%
  left_join(fungi_enrich_type,by = 'asv') %>%
  left_join(fungi_asv_0.0001_list,by = 'asv') %>%
  arrange(stage,enrich) %>%
  select(!c(Kingdom,Species))

fungi_enrich_taxa %>%
  select(asv,mean_asv) %>%
  unique() %>%
  summarise(abun = sum(mean_asv))

write.csv(fungi_enrich_taxa,file = 'caoxinzhuang_results/4.0/intermediate/fungi_enrich_taxa.csv',quote = F)

fungi_enrich_taxa
# sankey_genus_count_fungi <- table(fungi_sankey_data$Genus) %>%
#   as.data.frame() %>%
#   left_join(fungi_sankey_data,by = c('Var1' = 'Genus')) %>%
#   unique() %>%
#   relocate(Freq,.after = ncol(.))

sankey_genus_count_fungi<- fungi_enrich_taxa %>%
  na.omit() %>%
  make_long(enrich,Genus,Phylum)

sankey_genus_fungi_p <-
  ggplot(sankey_genus_count_fungi, aes(x = x, next_x = next_x,
                                                         node = node, next_node = next_node,
                                                         fill = node,
                                                         label = node)) +
  geom_sankey(flow.alpha = 0.5, #条带不透明度
              smooth = 8, #条带弯曲度
              width = 0.12) + #节点宽度
  geom_sankey_text(size = 3.2, #标签字号
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none') #隐藏图例

sankey_genus_fungi_p
sankey_genus_bac_p
#####save results#####

ggsave(venn_total_p,filename = 'caoxinzhuang_results/2.0/venn_total_p.pdf',width = 16,height = 6)
ggsave(enrichd_ASVs_p,filename = 'caoxinzhuang_results/2.0/enrichd_ASVs_p.pdf',width = 9,height = 6)
ggsave(sankey_genus_bac_p,filename = 'caoxinzhuang_results/4.0/intermediate/sankey_genus_bac_p.2.pdf',width = 6,height = 6)
ggsave(sankey_genus_fungi_p,filename = 'caoxinzhuang_results/4.0/intermediate/sankey_genus_fungi_p.2.pdf',width = 6,height = 6)

####phylogenetic tree####
#####bacteria#####
# enriched_bac_asv_name as list for selecting ASVs sequences 
# qiime2 select represent sequences and build phylogenetic tree
enrich_bac_tree <- read.tree('caoxinzhuang_results/2.0/intermediate/enrich_bacteria_tree.nwk')
numtip <- length(enrich_bac_tree$tip.label)
numtip
enrich_bac_tree$tip.label <- str_replace_all(enrich_bac_tree$tip.label,'[^ASV_\\d]','')
enrich_bac_tree$tip.label
enrich_bac_tree
#tree trunk part
enrich_bac_asv_phylum <- asv_bac %>%
  rownames_to_column(var = 'asv') %>%
  filter(asv %in% enriched_bac_asv_name) %>%
  select(asv,Phylum)
list(enrich_bac_asv_phylum$Phylum) %>% table() # list phylum composition

enrich_bac_asv_phylum_group <- list(Actinobacteriota = filter(enrich_bac_asv_phylum,
                                                            Phylum == 'Actinobacteriota') %>% select(asv) %>% as.list()%>% .[[1]],
                                    Chloroflexi = filter(enrich_bac_asv_phylum,
                                                            Phylum == 'Chloroflexi') %>% select(asv)%>% as.list() %>% .[[1]],
                                    Proteobacteria  = filter(enrich_bac_asv_phylum,
                                                              Phylum == 'Proteobacteria') %>% select(asv)%>% as.list()%>% .[[1]],
                                    Others = filter(enrich_bac_asv_phylum,
                                                    Phylum %in% c('Abditibacteriota','Acidobacteriota',
                                                                'Armatimonadota','Cyanobacteria',
                                                                'Firmicutes','Planctomycetota',
                                                                'Verrucomicrobiota','Patescibacteria')) %>% select(asv)%>% as.list()%>% .[[1]])

enrich_bac_tree_group <- groupOTU(enrich_bac_tree,enrich_bac_asv_phylum_group) # add phylum data to tree


dat1_enrich_bac <- data.frame(ID = enrich_bac_tree$tip.label) %>%
  left_join(rownames_to_column(asv_bac,var = 'asv'), by =c ('ID'='asv')) %>%
  select(ID,Phylum) %>%
  mutate(phylum2 = if_else(Phylum %in% c('Abditibacteriota','Acidobacteriota',
                                         'Armatimonadota','Cyanobacteria',
                                         'Firmicutes','Planctomycetota',
                                         'Verrucomicrobiota','Patescibacteria'),'Others',Phylum)) %>%
  left_join(bac_asv_0.0001_list,by = c('ID' = 'asv')) %>%
  left_join(enrich_bac_count,by =  c('ID' = 'asv')) 

dat2_enrich_bac <- data.frame(ID = enrich_bac_tree$tip.label)%>%
  mutate(Jointing = ifelse(ID %in% rownames(enrich_bac_jointing),'1','0')) %>%
  mutate(Filling = ifelse(ID %in% rownames(enrich_bac_filling),'1','0')) %>%
  mutate(Maturing = ifelse(ID %in% rownames(enrich_bac_maturing),'1','0')) %>%
   pivot_longer(!ID,names_to = 'stage',values_to = 'Existence') %>%
  mutate(stage2 = str_replace_all(.$stage,c('Jointing'='1','Filling'='2','Maturing'='3')))
dat2_enrich_bac$stage2 <- as.numeric(dat2_enrich_bac$stage2)

dat3_enrich_bac <-
  data.frame(ID = enrich_bac_tree$tip.label) %>%
  left_join(.,rownames_to_column(bac_asv,var = 'ID') , by = 'ID') %>%
  pivot_longer(!ID,names_to = 'sample_name',values_to = 'abundance') %>%
  left_join(group,by = 'sample_name') %>%
    group_by(ID,residue) %>%
    summarise(mean_abundance = mean(abundance)) %>%
    pivot_wider(names_from = 'residue' , values_from = 'mean_abundance') %>%
    mutate(abundance_change = (Turnover - Remove) / Remove,
           labels = if_else(abundance_change>0,'1','-1'),
           abs_abun = abs(abundance_change),
           log_abun = log10(100*abs_abun),
           minus = as.numeric(labels),
           real_abun = log_abun*minus) %>%
left_join(rownames_to_column(asv_bac,var = 'asv'), by =c ('ID'='asv')) %>%
  select(ID,real_abun,Phylum) %>%
  mutate(phylum2 = if_else(Phylum %in% c('Abditibacteriota','Acidobacteriota',
                                         'Armatimonadota','Cyanobacteria',
                                         'Firmicutes','Planctomycetota',
                                         'Verrucomicrobiota','Patescibacteria'),'Others',Phylum)) 

bac_enrich_tree <-
  ggtree(enrich_bac_tree_group,aes(col = group),
       layout="cir", branch.length='none', size=0.5) + 
  geom_treescale(x=-8, y=0, fontsize=0.1, linesize=0.1)+
  # geom_tiplab(size=3, color="black", aes(angle=0))+
  scale_color_manual(values=sunsx)+
  new_scale_fill()+
  geom_fruit(data=dat1_enrich_bac,
             geom=geom_point,
             mapping = aes(y=ID, fill= enrich),
             shape = 21,color = 'grey',size = 4,alpha = 0.75,
             position=position_identityx()) +
  scale_fill_manual(values = CAOXINZHUANG2)+
  new_scale_fill()+
  geom_fruit(data=dat2_enrich_bac,
             geom=geom_tile,
             mapping=aes(y=ID, x=stage2, fill=Existence),
             color = 'grey',alpha = 1,offset = 0.11,
             pwidth=0.2,
             position=position_identityx()) + 
  scale_fill_manual(values = c('white','#E4DBDE'))+
  new_scale_fill() +
  geom_fruit(data=dat3_enrich_bac,
             geom=geom_bar(),
             mapping = aes(x = real_abun,y=ID, fill= phylum2),
             stat="identity",offset = 0.35,
             orientation="y",
             pwidth = 0.8) +
  scale_fill_manual(values = sunsx)+
  theme(legend.position=c(0.95, 0.5),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        legend.spacing.y = unit(0.01, "cm"))
bac_enrich_tree 
ggsave(bac_enrich_tree,filename = 'caoxinzhuang_results/2.0/bacenrich_tree_2.0.pdf',width = 8,height = 8)
#####fungi#####
# enriched_fungi_asv_name as list for selecting ASVs sequences 
# qiime2 select represent sequences and build phylogenetic tree
enrich_fungi_tree <- read.tree('caoxinzhuang_results/2.0/intermediate/enrich_fungi_tree.nwk')
numtip <- length(enrich_fungi_tree$tip.label)
numtip
enrich_fungi_tree$tip.label <- str_replace_all(enrich_fungi_tree$tip.label,'[^ASV_\\d]','')
enrich_fungi_tree$tip.label

#tree trunk part
enriched_fungi_asv_name
enrich_fungi_asv_phylum <- asv_fungi %>%
  rownames_to_column(var = 'asv') %>%
  filter(asv %in% enriched_fungi_asv_name) %>%
  select(asv,Phylum)
list(enrich_fungi_asv_phylum$Phylum) %>% table() # list phylum composition
enrich_fungi_asv_phylum_group <- list(Ascomycota = filter(enrich_fungi_asv_phylum,
                                                              Phylum == 'p__Ascomycota') %>% select(asv) %>% as.list()%>% .[[1]],
                                      Basidiomycota = filter(enrich_fungi_asv_phylum,
                                                         Phylum == 'p__Basidiomycota') %>% select(asv)%>% as.list() %>% .[[1]],
                                      Mortierellomycota = filter(enrich_fungi_asv_phylum,
                                                              Phylum == 'p__Mortierellomycota') %>% select(asv)%>% as.list()%>% .[[1]],
                                      Others = filter(enrich_fungi_asv_phylum,
                                                    Phylum == 'p__Chytridiomycota') %>% select(asv)%>% as.list()%>% .[[1]])
enrich_fungi_asv_phylum_group
enrich_fungi_tree_group <- groupOTU(enrich_fungi_tree,enrich_fungi_asv_phylum_group) # add phylum data to tree


dat1_enrich_fungi <- data.frame(ID = enrich_fungi_tree$tip.label) %>%
  left_join(rownames_to_column(asv_fungi,var = 'asv'), by =c ('ID'='asv')) %>%
  select(ID,Phylum) %>%
  mutate(phylum2 = if_else(Phylum %in%  c('p__Ascomycota','p__Basidiomycota','p__Mortierellomycota'),Phylum,'Others')) %>%
  left_join(fungi_asv_0.0001_list,by = c('ID' = 'asv')) %>%
  left_join(enrich_fungi_count,by =c( 'ID' = 'asv')) %>%
  mutate(phylum3 = str_replace_all(phylum2,'p__','')) 
asv_fungi

dat2_enrich_fungi <- data.frame(ID = enrich_fungi_tree$tip.label)%>%
  mutate(Jointing = ifelse(ID %in% rownames(enrich_fungi_jointing),'1','0')) %>%
  mutate(Filling = ifelse(ID %in% rownames(enrich_fungi_filling),'1','0')) %>%
  mutate(Maturing = ifelse(ID %in% rownames(enrich_fungi_maturing),'1','0')) %>%
  pivot_longer(!ID,names_to = 'stage',values_to = 'Existence') %>%
  mutate(stage2 = str_replace_all(.$stage,c('Jointing'='1','Filling'='2','Maturing'='3')))

dat2_enrich_fungi$stage2 <- as.numeric(dat2_enrich_fungi$stage2)

dat3_enrich_fungi <-
  data.frame(ID = enrich_fungi_tree$tip.label) %>%
  left_join(.,rownames_to_column(fungi_asv,var = 'ID') , by = 'ID') %>%
  pivot_longer(!ID,names_to = 'sample_name',values_to = 'abundance') %>%
  left_join(group,by = 'sample_name') %>%
  group_by(ID,residue) %>%
  summarise(mean_abundance = mean(abundance+0.1)) %>%
  pivot_wider(names_from = 'residue' , values_from = 'mean_abundance') %>%
  mutate(abundance_change = (Turnover - Remove) / Remove,
         labels = if_else(abundance_change>0,'1','-1'),
         abs_abun = abs(abundance_change),
         log_abun = log10(100*abs_abun),
         minus = as.numeric(labels),
         real_abun = log_abun*minus) %>%
  left_join(rownames_to_column(asv_fungi,var = 'asv'), by =c ('ID'='asv')) %>%
  select(ID,real_abun,Phylum) %>%
  mutate(phylum2 = if_else(Phylum %in% 
                             c('p__Ascomycota','p__Basidiomycota','p__Mortierellomycota'),Phylum,'Others')) %>%
  mutate(phylum3 = str_replace_all(phylum2,'p__',''))

fungi_enrich_tree <-
  ggtree(enrich_fungi_tree_group,aes(col = group),
         layout="cir", branch.length='none', size=0.5) + 
  geom_treescale(x=-8, y=0, fontsize=0.1, linesize=0.1)+
  # geom_tiplab(size=3, color="black", aes(angle=0))+
  scale_color_manual(values=sunsx)+
  new_scale_fill()+
  geom_fruit(data=dat1_enrich_fungi,
             geom=geom_point,
             mapping = aes(y=ID, fill= enrich),
             shape = 21,color = 'grey',size = 4,alpha = 0.75,
             position=position_identityx()) +
  scale_fill_manual(values = CAOXINZHUANG2)+
  new_scale_fill()+
  geom_fruit(data=dat2_enrich_fungi,
             geom=geom_tile,
             mapping=aes(y=ID, x=stage2, fill=Existence),
             color = 'grey',alpha = 1,offset = 0.11,
             pwidth=0.2,
             position=position_identityx()) + 
  scale_fill_manual(values = c('white','#E4DBDE'))+
  new_scale_fill() +
  geom_fruit(data=dat3_enrich_fungi,
             geom=geom_bar(),
             mapping = aes(x = real_abun,y=ID, fill= phylum3),
             stat="identity",offset = 0.35,
             orientation="y",
             pwidth = 0.8) +
  scale_fill_manual(values = sunsx)+
  theme(legend.position=c(0.95, 0.5),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        legend.spacing.y = unit(0.01, "cm"))

ggsave(fungi_enrich_tree,filename = 'caoxinzhuang_results/2.0/fungienrich_tree2.0.pdf',width = 8,height = 8)


