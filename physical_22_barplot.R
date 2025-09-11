library(tidyverse)
library(rstatix)
library(car)

physicalcarphysical_2.0 <- read.delim('22_physicochemical/22_soil_physical.txt',header = T)
colnames(physicalcarphysical_2.0)
physicalcarphysical_2.0 %>% select(stage:fertilization,moisture) %>%
  filter(stage=='Jointing') %>%
  as.data.frame() %>%
  t_test(moisture~residue,data = .)

physical_2.0$stage <- factor(physical_2.0$stage,levels = c('Jointing','Filling','Maturing'))
physical_2.0$fertilization <- factor(physical_2.0$fertilization , 
                                   levels = c('CK','HO','HF','MF')) 
physical_2.0$residue <- factor(physical_2.0$residue,levels = c('Remove','Turnover'))
physical_2.0 <- physical_2.0 %>% arrange(stage,residue,fertilization)

physical_2.0$CO2 <- respiration_co2$CO2
colnames(physical_names)

physical_names <- c('NH4','NO3','SOC','DOC','MBC','AP','AK','CO2','DON','MBN','pH','TC','TN')
theme_NMDS2 <- theme_bw(base_size = 15)+
  theme(legend.position = 'none',panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        axis.title.y = element_text(size = 13,color="black"),
        axis.title.x = element_text(size = 13,color="black"))

physical_long <- physical_2.0 %>% 
  pivot_longer(!c(stage,residue,fertilization), names_to = 'indicator',values_to = 'value') %>%
  mutate(group = paste(fertilization,residue,sep = '_'),
         group2 = factor(group,levels = c('CK_Remove','CK_Turnover','HO_Remove','HO_Turnover',
                                          'HF_Remove','HF_Turnover','MF_Remove','MF_Turnover')))
physical_long_mean_sd <- physical_long %>%
  group_by(stage,residue,fertilization,indicator) %>%
  summarise(mean_value = mean(value),sd_value = sd(value))%>%
  mutate(group = paste(fertilization,residue,sep = '_'),
         group2 = factor(group,levels = c('CK_Remove','CK_Turnover','HO_Remove','HO_Turnover',
                                          'HF_Remove','HF_Turnover','MF_Remove','MF_Turnover')))
physical_long_mean_sd %>%
  filter(indicator=='TC'&stage=='Filling')


physical_plot <- function(indicatorsssss){
  physical_long_mean_sd %>%
    filter(indicator ==indicatorsssss) %>%
    ggplot()+
    geom_bar(aes(x = group2,y = mean_value,fill = residue),stat = 'identity',width = 0.8)+
    geom_errorbar(aes(x = group2,ymax = mean_value+sd_value,ymin = mean_value),size = 0.5,width = 0.8)+
    scale_fill_manual(values = CAOXINZHUANG2)+
    facet_grid(.~stage)+
    labs(x = NULL,y = indicatorsssss)+
    theme_NMDS2+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_blank(),
          strip.text = element_blank())
}
unique(physical_long_mean_sd$indicator)

DOC_plot <- physical_plot('DOC')
SOC_plot <- physical_plot('SOC')
MBC_plot <- physical_plot('MBC')
NO3_plot <- physical_plot('NO3')
pH_plot <- physical_plot('pH')
AK_plot <- physical_plot('AK')
AP_plot <- physical_plot('AP')
DON_plot <- physical_plot('DON')
MBN_plot <- physical_plot('MBN')
TC_plot <- physical_plot('TC')
TN_plot <- physical_plot('TN')

NO3_plot_supple <- physical_long_mean_sd %>%
  filter(indicator =='NO3') %>%
  ggplot()+
  geom_bar(aes(x = group2,y = mean_value,fill = residue),stat = 'identity',width = 0.8)+
  geom_errorbar(aes(x = group2,ymax = mean_value+sd_value,ymin = mean_value),size = 0.5,width = 0.8)+
  scale_fill_manual(values = CAOXINZHUANG2)+
  facet_grid(.~stage)+
  scale_y_continuous(limits = c(0,10))+
  labs(x = NULL,y = 'NO3')+
  theme_NMDS2+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_blank(),
        strip.text = element_blank())

physical_long_mean_sd %>%
  filter(indicator =='pH') %>%
  ggplot()+
  geom_bar(aes(x = group2,y = mean_value,fill = residue),stat = 'identity',width = 0.8)+
  geom_errorbar(aes(x = group2,ymax = mean_value+sd_value,ymin = mean_value),size = 0.5,width = 0.8)+
  scale_fill_manual(values = CAOXINZHUANG2)+
  facet_grid(.~stage)+
  scale_y_continuous(limits = c(0,8.5))+
  labs(x = NULL,y = 'pH')+
  theme_NMDS2+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_blank(),
        strip.text = element_blank())

NO3_plot_supple
mbc_plot_supple <- physical_long_mean_sd %>%
  filter(indicator =='MBC') %>%
  ggplot()+
  geom_bar(aes(x = group2,y = mean_value,fill = residue),stat = 'identity',width = 0.8)+
  geom_errorbar(aes(x = group2,ymax = mean_value+sd_value,ymin = mean_value),size = 0.5,width = 0.8)+
  scale_fill_manual(values = CAOXINZHUANG2)+
  facet_grid(.~stage)+
  scale_y_continuous(limits = c(0,50))+
  labs(x = NULL,y = 'MBC')+
  theme_NMDS2+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_blank(),
        strip.text = element_blank())

physical_supple <- ggarrange(NO3_plot_supple,SOC_plot,mbc_plot_supple,DOC_plot,nrow = 2,ncol = 2)
physical_total <- ggarrange(NO3_plot,SOC_plot,mbc_plot,DOC_plot,nrow = 2,ncol = 2)

ggsave('caoxinzhuang_results/5.0/physical_supple.pdf',physical_supple,width = 10,height = 6)
ggsave('caoxinzhuang_results/5.0/physical_total2.pdf',physical_total,width = 10,height = 6)
physical_total
colnames(physical_2.0)

for (aaa in unique(physical_long$stage)) {
  dataaaa = filter(physical_long,stage == aaa)
  print(aaa)
  for (i in colnames(physical_2.0)[4:16]) {
    data2 = filter(dataaaa,indicator == i)
    anova = aov(value ~ residue*fertilization,data = data2)
    summary__ = summary(anova)
    anova_tukey = aov(value ~ group2,data = data2)
    leveneTest_results <- leveneTest(value ~ residue*fertilization, data = data2)
    tuker = tukey_hsd(anova_tukey,conf.level = 0.95)
    residuals <- residuals(anova)
    shapiro_test <- shapiro.test(residuals)
    kruskal_fertilization <- kruskal.test(value ~ fertilization, data = data2)
    kruskal_residue <- kruskal.test(value ~ residue, data = data2)
    print(i)
    ifelse(shapiro_test$p.value<0.05|leveneTest_results$`Pr(>F)`[[1]]<0.05,
           print(leveneTest_results$`Pr(>F)`),print('Nengcheng________'))
    ifelse(shapiro_test$p.value<0.05|leveneTest_results$`Pr(>F)`[[1]]<0.05,
           print(kruskal_residue),print(summary__))
    ifelse(shapiro_test$p.value<0.05|leveneTest_results$`Pr(>F)`[[1]]<0.05,
           print(kruskal_fertilization),print('Nengcheng________'))
  }
}

final_results_f <- data.frame()
final_results_tukey <- data.frame()

for (aaa in unique(physical_long$stage)) {
  dataaaa <- filter(physical_long, stage == aaa)
  
  for (i in physical_names) {
    data2 <- filter(dataaaa, indicator == i)
    
    anova <- aov(value ~ residue * fertilization, data = data2)
    
    anova_summary <- summary(anova)[[1]]
    
    effect_results <- data.frame(
      stage = aaa,
      indicator = i,
      effect = c("residue", "fertilization", "residue:fertilization"),
      F_value = anova_summary$`F value`[1:3],
      p_value = anova_summary$`Pr(>F)`[1:3]
    )
    print(effect_results)
    data2$group2 <- interaction(data2$residue, data2$fertilization)
    
    anova_tukey <- aov(value ~ group2, data = data2)
    tuker <- tukey_hsd(anova_tukey, conf.level = 0.95)
    
    tukey_results <- tuker %>%
      mutate(
        stage = aaa,
        indicator = i,
        comparison = paste(group1, "-", group2)
      ) %>%
      select(stage, indicator, 
             comparison, estimate, conf.low, conf.high, p.adj,p.adj.signif)
    final_results_f = bind_rows(final_results_f,effect_results)
    final_results_tukey = bind_rows(final_results_tukey,tukey_results)
  }
}

final_results_f %>%
  filter(p_value<0.05)
final_results_f
# &indicator%in% c('NH4','SOC','DOC','MBC','CO2')

supplement_figures <- ggarrange(AK_plot,AP_plot,DON_plot,MBN_plot,TC_plot,TN_plot,NH4_plot,pH_plot,ncol = 2,nrow = 4)

ggsave(supplement_figures,file = 'caoxinzhuang_results/5.0/supplement_figure2.pdf',width = 10,height = 12)
