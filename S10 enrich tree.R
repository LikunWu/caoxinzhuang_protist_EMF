library(tidyverse)
library(DESeq2)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(reshape2)
library(Hmisc)

asv_proto <- read.csv('sequence/proto_ASVs.csv',row.names = 1,header = T) %>%
  filter(Kingdom == "Eukaryota" &(Phylum %in% c("Amoebozoa",'Excavata','Obazoa', "TSAR"))) %>%
  select(!c(v.R.C.R.5,R.R.D.R.5))
taxa_protists <- asv_proto[,143:149]
protists_asv <- asv_proto[,1:142]
view(protists_asv)
group_rownames <- group %>%
  filter(!sample_name %in% c('v.R.C.R.5','R.R.D.R.5'))


protist_dds <- DESeqDataSetFromMatrix(data.frame(protists_asv),colData = group_rownames,design = ~ residue)
protist_deseq2 <- DESeq(protist_dds)

protist_desq2_results <- results(protist_deseq2) %>%
  as.data.frame() %>%
  filter(abs(log2FoldChange)>=2&padj<0.05)

protist_desq2_results

write.csv(protist_desq2_results,file = 'caoxinzhuang_results/6.0emf/supplementary/enrich_protist.csv')

####phylogenetric tree----------------------------------------------------------
enrich_tree <- read.tree('caoxinzhuang_results/6.0emf/supplementary/intermediate/tree.nwk')
enrich_tree$tip.label <- gsub("['\"]", "", enrich_tree$tip.label) 

enrich_protist_asv_id <- protist_desq2_results %>%
  rownames_to_column(var = 'asv_id') %>%
  mutate(enrich_site = if_else(log2FoldChange>1,1,0))

enrich_protist_asv_taxa <- taxa_protists %>%
  rownames_to_column(var = 'asv_id') %>%
  filter(asv_id %in% enrich_protist_asv_id$asv_id)

enrich_protist_asv_group <- list(Cercozoa = filter(enrich_protist_asv_taxa,
                                                             Order == 'Cercozoa') %>% select(asv_id) %>% as.list()%>% .[[1]],
                                    Ciliophora = filter(enrich_protist_asv_taxa,
                                                        Order == 'Ciliophora') %>% select(asv_id) %>% as.list()%>% .[[1]],
                                    Amoebozoa = filter(enrich_protist_asv_taxa,
                                                          Phylum == 'Amoebozoa') %>% select(asv_id)%>% as.list()%>% .[[1]],
                                    Others = filter(enrich_protist_asv_taxa,Phylum=='TSAR'&
                                                    !Order %in% c('Cercozoa','Ciliophora')) %>% select(asv_id)%>% as.list()%>% .[[1]])
enrich_protist_asv_group
enrich_protist_asv_group_tree <- groupOTU(enrich_tree,enrich_protist_asv_group) # add phylum data to tree

enrichd_asv_rela_abun <- t(protist_relative_table_t) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'asv') %>%
  pivot_longer(!asv,names_to = 'sample_name') %>%
  left_join(group,by = 'sample_name') %>%
  filter(str_detect(sample_name,'1|2|3|4')&asv %in% enrich_protist_asv_id$asv_id) %>%
  pivot_wider(names_from = asv,values_from = value) %>%
  arrange(stage,residue,fertilization) %>%
  cbind(EMF = emf_anova$meanFunction)

enrich_protist_asv_group <- list(Cercozoa = filter(enrich_protist_asv_taxa,
                                                   Order == 'Cercozoa') %>% select(asv_id) %>% as.list()%>% .[[1]],
                                 Ciliophora = filter(enrich_protist_asv_taxa,
                                                     Order == 'Ciliophora') %>% select(asv_id) %>% as.list()%>% .[[1]],
                                 Amoebozoa = filter(enrich_protist_asv_taxa,
                                                    Phylum == 'Amoebozoa') %>% select(asv_id)%>% as.list()%>% .[[1]],
                                 Others = filter(enrich_protist_asv_taxa,Phylum=='TSAR'&
                                                   !Order %in% c('Cercozoa','Ciliophora')) %>% select(asv_id)%>% as.list()%>% .[[1]])
head(enrich_protist_asv_taxa)
protist_enrich_asv_taxa <-enrich_protist_asv_taxa %>%
  mutate(taxa = ifelse(Order == 'Cercozoa','Cercozoa',ifelse(
    Order == 'Ciliophora','Ciliophora',ifelse(
      Phylum == 'Amoebozoa','Amoebozoa','Others'))))%>%
  select(asv_id,taxa)

enrichd_asv_rela_abun_ring <- enrichd_asv_rela_abun %>%
  select(starts_with('ASV_'),residue) %>%
  pivot_longer(!residue,names_to = 'asv_id') %>%
  group_by(asv_id,residue) %>%
  summarise(mean_real_abun = mean(value)) %>%
  pivot_wider(names_from = residue,values_from = mean_real_abun) %>%
  left_join(protist_enrich_asv_taxa, by = 'asv_id')

protist_emf <- enrichd_asv_rela_abun %>%
  select(!c(sample_name,stage,residue,fertilization)) %>%
  scale() %>%  as.matrix() %>% rcorr(. , type = 'spearman')

asv_count <- c(nrow(enrich_protist_asv_taxa)+1)
protist_emf__padj <- p.adjust(protist_emf$P,method = "BH") %>% matrix(ncol = asv_count)

protist_emf_matrixxx <- protist_emf$r[-asv_count,asv_count] 
protist_emf_matrixxx
star_protist_emf<- ifelse(protist_emf__padj < 0.001, "***",
                          ifelse(protist_emf__padj < 0.01, "**",
                                 ifelse(protist_emf__padj < 0.05, "*", "")))
star_protist_emf[-asv_count,asv_count] 
protist_emf_correlation$asv_id

protist_emf_correlation <- protist_emf_matrixxx%>%
  as.data.frame()%>%
  rownames_to_column(var = 'asv_id')

colnames(protist_emf_correlation) <- c('asv_id','R')

# enrich_tree <-
  ggtree(enrich_protist_asv_group_tree,aes(col = group),
       layout="cir", size=0.5, right = TRUE, open.angle = 80,branch.length = 'none') + 
  geom_treescale(x=-8, y=0, fontsize=0.1, linesize=0.1)+
  scale_color_manual(values = c('white','#3AA9AE','#E67D3A','#6E4196','grey80'))+
  geom_tiplab(size=3, color="black", aes(angle=0))
  new_scale_fill() +
 geom_fruit(data=enrichd_asv_rela_abun_ring,
                              geom=geom_bar(),
                              mapping = aes(x = -Remove,y=asv_id, fill= taxa),
                             stat="identity",offset = 0.6,
                             orientation="y",
                               pwidth = 0.5) +
  scale_fill_manual(values = c('#3AA9AE','#E67D3A','#6E4196','grey80'))+
  new_scale_fill() +
  geom_fruit(data=enrichd_asv_rela_abun_ring,
             geom=geom_bar(),
             mapping = aes(x = Turnover,y=asv_id, fill= taxa),
             stat="identity",offset = -0.45,
             orientation="y",
             pwidth = 0.5) +
  scale_fill_manual(values = c('#3AA9AE','#E67D3A','#6E4196','grey80'))+
  new_scale_fill() +
  geom_fruit(data=protist_emf_correlation,
                            geom=geom_bar(),
                    mapping = aes(x = 0.1,y=asv_id, fill= R),
                     stat="identity",offset = 0,
                     orientation="y",
                        pwidth = 0.15) +
  scale_fill_gradient2(
      low = "#7D0B5D", mid = "white", high = "#428949",midpoint = 0)
enrich_tree
ggsave(filename = 'caoxinzhuang_results/6.0emf/supplementary/enrich_tree.pdf',enrich_tree,height = 10,width = 10)
