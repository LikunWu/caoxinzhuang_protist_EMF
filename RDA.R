library(vegan)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(openxlsx)
library(rdacca.hp)

theme_rda <- theme_bw()+
  theme(axis.title = element_text(size = 18,colour = "black"))+
  theme(axis.text = element_text(size = 18,color="black"))+
  theme(panel.grid=element_blank())+
  theme(legend.position="right") # 图例设置

physical22 <- read.delim('22_physicochemical/22physicochemical.txt',header = T)
physical22$residue <- factor(physical22$residue,levels = c('Remove','Turnover'))
physical22$stage <- factor(physical22$stage,levels = c('Jointing','Filling','Maturing'))
physical22$fertilization <- factor(physical22$fertilization , 
                                   levels = c('CK','HO','HF','MF')) 
alpha_total <- cbind(alpha_bac_plot,shannon_fungi = alpha_fungi_plot$Shannon)
alpha_total$residue <- factor(alpha_total$residue,levels = c('Remove','Turnover'))
alpha_total$stage <- factor(alpha_total$stage , levels = c('Jointing','Filling','Maturing'))
alpha_total$fertilization <- factor(alpha_total$fertilization ,
                                    levels = c('ck','high_organic','high_fertilizer','mineral_fertilizer'))

physical22_2 <- physical22 %>% arrange(stage,residue,fertilization) %>%
  as.data.frame() %>%
  dplyr::select(!c(stage,fertilization,MBC,MBN)) 
physical22_2
sample_rep <- rep(1:6,times=nrow(alpha_bac_plot)/6)

####BACTERIA####
# bac_total <-   
#   fungi_asv %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = 'sample_name') %>%
#   cbind(rep = sample_rep) %>%
#   left_join(group, by = 'sample_name') %>%
#   arrange(stage,residue,fertilization) %>%
#   filter(rep<=4) %>%
#   column_to_rownames(var = 'sample_name') %>%
#   select(!c(colnames(group)[2:4],rep)) %>%
#   decostand(method = 'hellinger')
# 
# physical22_2_total <- physical22_2 %>%
#   select(!residue) %>%
#   scale()
# rownames(physical22_2_total) <- rownames(physical22_2_total)
# rda_bac_total <- rda(bac_total , physical22_2_total)
# envfit(rda_bac_total , physical22_2_total ,permutations = 9999)

####RENMOVE####
 
remove_group <- group %>%
  cbind(rep = sample_rep) %>%
  arrange(stage,residue,fertilization) %>%
  filter(residue == 'Remove'& rep<=4) %>%
  column_to_rownames(var = 'sample_name') 
RsquareAdj(remove_bac_rda)
anova.cca(remove_bac_rda, permutations = how(nperm = 9999))
remove_env_fit_bac <- envfit(remove_bac_rda , physical22_2_remove,permutations = 9999)
remove_env_fit_bac
vif.cca(remove_bac_rda)

rda_score_bac_remove <- data.frame(remove_bac_rda$CCA$u[,1:2]) %>%
  cbind(remove_group)
rda_env_bac_remove <- data.frame(remove_bac_rda$CCA$biplot[,1:2])

remove_env_fit_bac_data <- cbind(R = remove_env_fit_bac$vectors$r,
                                 p = remove_env_fit_bac$vectors$pvals) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'indicator') %>%
  mutate(label = if_else(p < 0.001,"***",
                         if_else(p <0.01,"**",
                                 if_else(p<0.05,"*",""))))
remove_env_fit_bac_p <- remove_env_fit_bac_data %>%
  arrange(desc(R)) %>%
  ggplot(aes(x = factor(indicator,levels = indicator) , y = R))+
  geom_col(aes(fill =  label))+
  labs(x = 'Indicator' , y = 'R2')+
  geom_text(aes(y = R+0.01,label = label))+
  scale_fill_manual(values = col_icamp)+
  theme_NMDS

bac_remove_rda_p <- ggplot(rda_score_bac_remove, aes(RDA1, RDA2)) +
  geom_point(aes(color = stage),size = 3) + 
  xlab(paste("RDA1 ( ",round(remove_bac_rda$CCA$eig[1]/sum(remove_bac_rda$CCA$eig)*100,2),"%"," )", sep = "")) + 
  ylab(paste("RDA2 ( ",round(remove_bac_rda$CCA$eig[2]/sum(remove_bac_rda$CCA$eig)*100,2),"%"," )", sep = "")) +
  geom_segment(data = rda_env_bac_remove, aes(x = 0, y = 0, xend = rda_env_bac_remove[,1],
                                              yend = rda_env_bac_remove[,2]),
               arrow = arrow(length =unit(0.35,"cm"),type = "closed",angle=22.5),
               linetype=1,colour ="black",size=0.6)+ #添加箭头
  geom_text_repel(data = rda_env_bac_remove, segment.colour = "black",
                  aes(x = rda_env_bac_remove[,1], y = rda_env_bac_remove[,2],
                      label = rownames(rda_env_bac_remove)),size=6) +
  scale_color_manual(values = sunsx)+
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_rda
remove_env_fit_bac_p
bac_remove_rda_p

####turnover####
bac_turnover <-   
  bac_asv %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name') %>%
  cbind(rep = sample_rep) %>%
  left_join(group, by = 'sample_name') %>%
  arrange(stage,residue,fertilization) %>%
  filter(residue =='Turnover'&rep<=4) %>%
  column_to_rownames(var = 'sample_name') %>%
  select(!c(colnames(group)[2:4],rep)) %>%
  decostand(method = 'hellinger')

physical22_2_turnover <- physical22_2 %>%
  filter(residue == 'Turnover') %>% 
  select(!c(residue,TC)) %>%
  scale()
rownames(physical22_2_turnover) <- rownames(bac_turnover)

turnover_bac_rda <- rda(bac_turnover , physical22_2_turnover)
vif.cca(turnover_bac_rda)
mite.cap.hp_bac_turnover<-rdacca.hp(vegdist(bac_turnover,method = "bray"),
                                  physical22_2_turnover,method='RDA',type='R2',scale=FALSE)
mite.cap.hp_bac_turnover
turnover_group <- group %>%
  cbind(rep = sample_rep) %>%
  arrange(stage,residue,fertilization) %>%
  filter(residue == 'Turnover'& rep<=4) %>%
  column_to_rownames(var = 'sample_name') 
RsquareAdj(remove_fungi_rda)
anova.cca(turnover_bac_rda, permutations = how(nperm = 9999))
turnover_env_fit_bac <- envfit(turnover_bac_rda , physical22_2_turnover,permutations = 9999)
turnover_env_fit_bac
vif.cca(turnover_bac_rda)

rda_score_bac_turnover <- data.frame(turnover_bac_rda$CCA$u[,1:2]) %>%
  cbind(turnover_group)
rda_env_bac_turnover <- data.frame(turnover_bac_rda$CCA$biplot[,1:2])
turnover_env_fit_bac_data <- cbind(R = turnover_env_fit_bac$vectors$r,
                                   p = turnover_env_fit_bac$vectors$pvals) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'indicator') %>%
  mutate(label = if_else(p < 0.001,"***",
                         if_else(p <0.01,"**",
                                 if_else(p<0.05,"*",""))))
turnover_env_fit_bac_p <- turnover_env_fit_bac_data %>%
  arrange(desc(R)) %>%
  ggplot(aes(x = factor(indicator,levels = indicator) , y = R))+
  geom_col(aes(fill =  label))+
  labs(x = 'Indicator' , y = 'R2')+
  geom_text(aes(y = R+0.01,label = label))+
  scale_fill_manual(values = c('#2A9E89' , '#F7B7AE'))+
  theme_NMDS
bac_turnover_rda_p <- ggplot(rda_score_bac_turnover, aes(RDA1, RDA2)) +
  geom_point(aes(color = stage),size = 3) + 
  xlab(paste("RDA1 ( ",round(turnover_bac_rda$CCA$eig[1]/sum(turnover_bac_rda$CCA$eig)*100,2),"%"," )", sep = "")) + 
  ylab(paste("RDA2 ( ",round(turnover_bac_rda$CCA$eig[2]/sum(turnover_bac_rda$CCA$eig)*100,2),"%"," )", sep = "")) +
  geom_segment(data = rda_env_bac_turnover, aes(x = 0, y = 0, xend = rda_env_bac_turnover[,1],
                                              yend = rda_env_bac_turnover[,2]),
               arrow = arrow(length =unit(0.35,"cm"),type = "closed",angle=22.5),
               linetype=1,colour ="black",size=0.6)+ #添加箭头
  geom_text_repel(data = rda_env_bac_turnover, segment.colour = "black",
                  aes(x = rda_env_bac_turnover[,1], y = rda_env_bac_turnover[,2],
                      label = rownames(rda_env_bac_turnover)),size=6) +
  scale_color_manual(values = sunsx)+
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_rda
turnover_env_fit_bac_p
bac_turnover_rda_p

####fungi####
####RENMOVE####
fungi_remove <-   
  fungi_asv %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name') %>%
  cbind(rep = sample_rep) %>%
  left_join(group, by = 'sample_name') %>%
  arrange(stage,residue,fertilization) %>%
  filter(residue =='Remove'&rep<=4) %>%
  column_to_rownames(var = 'sample_name') %>%
  select(!c(colnames(group)[2:4],rep)) %>%
  decostand(method = 'hellinger')

physical22_2_remove <- physical22_2 %>%
  filter(residue == 'Remove') %>% 
  select(!c(residue,TC)) %>%
  scale()
rownames(physical22_2_remove) <- rownames(fungi_remove)

remove_fungi_rda <- rda(fungi_remove , physical22_2_remove)
vif.cca(remove_fungi_rda)
mite.cap.hp_fungi_remove<-rdacca.hp(vegdist(fungi_remove,method = "bray"),
                                    physical22_2_remove,method='RDA',type='R2',scale=FALSE)
mite.cap.hp_fungi_remove
remove_group <- group %>%
  cbind(rep = sample_rep) %>%
  arrange(stage,residue,fertilization) %>%
  filter(residue == 'Remove'& rep<=4) %>%
  column_to_rownames(var = 'sample_name') 
RsquareAdj(remove_fungi_rda)
anova.cca(remove_fungi_rda, permutations = how(nperm = 9999))
remove_env_fit_fungi <- envfit(remove_fungi_rda , physical22_2_remove,permutations = 9999)
vif.cca(remove_fungi_rda)
remove_env_fit_fungi
rda_score_fungi_remove <- data.frame(remove_fungi_rda$CCA$u[,1:2]) %>%
  cbind(remove_group)
rda_env_fungi_remove <- data.frame(remove_fungi_rda$CCA$biplot[,1:2])
remove_env_fit_fungi_data <- cbind(R = remove_env_fit_fungi$vectors$r,
                                   p = remove_env_fit_fungi$vectors$pvals) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'indicator') %>%
  mutate(label = if_else(p < 0.001,"***",
                         if_else(p <0.01,"**",
                                 if_else(p<0.05,"*",""))))
remove_env_fit_fungi_p <- remove_env_fit_fungi_data %>%
  arrange(desc(R)) %>%
  ggplot(aes(x = factor(indicator,levels = indicator) , y = R))+
  geom_col(aes(fill =  label))+
  labs(x = 'Indicator' , y = 'R2')+
  geom_text(aes(y = R+0.01,label = label))+
  scale_fill_manual(values = col_icamp)+
  theme_NMDS
fungi_remove_rda_p <- ggplot(rda_score_fungi_remove, aes(RDA1, RDA2)) +
  geom_point(aes(color = stage),size = 3) + 
  xlab(paste("RDA1 ( ",round(remove_fungi_rda$CCA$eig[1]/sum(remove_fungi_rda$CCA$eig)*100,2),"%"," )", sep = "")) + 
  ylab(paste("RDA2 ( ",round(remove_fungi_rda$CCA$eig[2]/sum(remove_fungi_rda$CCA$eig)*100,2),"%"," )", sep = "")) +
  geom_segment(data = rda_env_fungi_remove, aes(x = 0, y = 0, xend = rda_env_fungi_remove[,1],
                                              yend = rda_env_fungi_remove[,2]),
               arrow = arrow(length =unit(0.35,"cm"),type = "closed",angle=22.5),
               linetype=1,colour ="black",size=0.6)+ #添加箭头
  geom_text_repel(data = rda_env_fungi_remove, segment.colour = "black",
                  aes(x = rda_env_fungi_remove[,1], y = rda_env_fungi_remove[,2],
                      label = rownames(rda_env_fungi_remove)),size=6) +
  scale_color_manual(values = sunsx)+
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_rda
remove_env_fit_fungi_p
fungi_remove_rda_p
####turnover####
fungi_turnover <-   
  fungi_asv %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_name') %>%
  cbind(rep = sample_rep) %>%
  left_join(group, by = 'sample_name') %>%
  arrange(stage,residue,fertilization) %>%
  filter(residue =='Turnover'&rep<=4) %>%
  column_to_rownames(var = 'sample_name') %>%
  select(!c(colnames(group)[2:4],rep)) %>%
  decostand(method = 'hellinger')

physical22_2_turnover <- physical22_2 %>%
  filter(residue == 'Turnover') %>% 
  select(!c(residue,TC)) %>%
  scale()
rownames(physical22_2_turnover) <- rownames(fungi_turnover)

turnover_fungi_rda <- rda(fungi_turnover , physical22_2_turnover)
vif.cca(turnover_fungi_rda)
mite.cap.hp_fungi_turnover<-rdacca.hp(vegdist(fungi_turnover,method = "bray"),
                                      physical22_2_turnover,method='RDA',type='R2',scale=FALSE)
mite.cap.hp_fungi_turnover
turnover_group <- group %>%
  cbind(rep = sample_rep) %>%
  arrange(stage,residue,fertilization) %>%
  filter(residue == 'Turnover'& rep<=4) %>%
  column_to_rownames(var = 'sample_name') 
RsquareAdj(turnover_fungi_rda)
anova.cca(turnover_fungi_rda, permutations = how(nperm = 9999))
turnover_env_fit_fungi <- envfit(turnover_fungi_rda , physical22_2_turnover,permutations = 9999)
turnover_env_fit_fungi
vif.cca(turnover_fungi_rda)

rda_score_fungi_turnover <- data.frame(turnover_fungi_rda$CCA$u[,1:2]) %>%
  cbind(turnover_group)
rda_env_fungi_turnover <- data.frame(turnover_fungi_rda$CCA$biplot[,1:2])

fungi_turnover_rda_p <- ggplot(rda_score_fungi_turnover, aes(RDA1, RDA2)) +
  geom_point(aes( color = stage),size = 3) + 
  xlab(paste("RDA1 ( ",round(turnover_fungi_rda$CCA$eig[1]/sum(turnover_fungi_rda$CCA$eig)*100,2),"%"," )", sep = "")) + 
  ylab(paste("RDA2 ( ",round(turnover_fungi_rda$CCA$eig[2]/sum(turnover_fungi_rda$CCA$eig)*100,2),"%"," )", sep = "")) +
  geom_segment(data = rda_env_fungi_turnover, aes(x = 0, y = 0, xend = rda_env_fungi_turnover[,1],
                                                yend = rda_env_fungi_turnover[,2]),
               arrow = arrow(length =unit(0.35,"cm"),type = "closed",angle=22.5),
               linetype=1,colour ="black",size=0.6)+ #添加箭头
  geom_text_repel(data = rda_env_fungi_turnover, segment.colour = "black",
                  aes(x = rda_env_fungi_turnover[,1], y = rda_env_fungi_turnover[,2],
                      label = rownames(rda_env_fungi_turnover)),size=6) +
  scale_color_manual(values = sunsx)+
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_rda

turnover_env_fit_fungi_data <- cbind(R = turnover_env_fit_fungi$vectors$r,
                                     p = turnover_env_fit_fungi$vectors$pvals) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'indicator') %>%
  mutate(label = if_else(p < 0.001,"***",
                         if_else(p <0.01,"**",
                                 if_else(p<0.05,"*",""))))
turnover_env_fit_fungi_p <- turnover_env_fit_fungi_data %>%
  arrange(desc(R)) %>%
  ggplot(aes(x = factor(indicator,levels = indicator) , y = R))+
  geom_col(aes(fill =  label))+
  labs(x = 'Indicator' , y = 'R2')+
  geom_text(aes(y = R+0.01,label = label))+
  scale_fill_manual(values = col_icamp)+
  theme_NMDS
turnover_env_fit_fungi_p

fungi_remove_rda_p
fungi_turnover_rda_p
####output####
remove_env_fit_bac$vectors$r
turnover_env_fit_bac
remove_env_fit_fungi
turnover_env_fit_fungi
rda_total_p <- ggarrange(bac_remove_rda_p,bac_turnover_rda_p,fungi_remove_rda_p,
                         fungi_turnover_rda_p,widths = 8,heights = 8,
                         common.legend = TRUE, legend = "bottom")
rda_total_p
env_fit_total_p <- ggarrange(remove_env_fit_bac_p,turnover_env_fit_bac_p,
                             remove_env_fit_fungi_p,turnover_env_fit_fungi_p,
                             widths = 8,heights = c(2.472,2.472,2.472,2.472),nrow = 2,ncol = 2)
env_fit_total_p

rdacca_hp_plot <- function(aaa,coloooo){
  dataaaa = aaa$Hier.part
  dataaaa %>%
    as.data.frame() %>%
    rownames_to_column(var = 'indicators') %>%
    mutate(indicator = factor(indicators))%>%
    arrange(desc(Individual)) %>%
    ggplot(.)+
    geom_col(aes(x = factor(indicator,levels = indicator),y = `I.perc(%)`),fill = coloooo,alpha = 0.9)+
    labs(x = NULL,y = NULL)+
    theme_NMDS2
}
'#B44847'
bac_harvest_plot <- rdacca_hp_plot(mite.cap.hp_bac_remove,'#B44847')
bac_retention_plot <- rdacca_hp_plot(mite.cap.hp_bac_turnover,'#00A1D6')
fungi_harvest_plot <- rdacca_hp_plot(mite.cap.hp_fungi_remove,'#B44847')
fungi_retention_plot <- rdacca_hp_plot(mite.cap.hp_fungi_turnover,'#00A1D6')

fig_6 <- ggarrange(bac_harvest_plot,bac_retention_plot,fungi_harvest_plot,fungi_retention_plot,ncol = 2,nrow = 2)

# ggsave(rda_total_p, filename = 'caoxinzhuang_results/4.0/intermediate/rda_total_p.pdf',
#        height = 8,width  = 8)
# ggsave(env_fit_total_p, filename = 'caoxinzhuang_results/4.0/intermediate/env_fit_total_p.pdf',
#        height = 4.944,width  = 8)

ggsave(fig_6, filename = 'caoxinzhuang_results/5.0/fig_6.pdf',
       height = 6.18,width  = 10)
