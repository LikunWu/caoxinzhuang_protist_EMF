library(tidyverse)
library(vegan)
library(rdacca.hp)
library(ggpubr)
library(ggrepel)
library(randomForest)
library(lme4)
library(lmerTest)
library(MuMIn)
library(rstatix)

biomass_22 <-  weight_height_22_long %>%
  group_by(rep , residue , fertilization ,indicator) %>%
  summarise(values = mean(value)) %>%
  pivot_wider(names_from = indicator,values_from = values)  %>%
  mutate(residue = factor(residue,levels = c('R','T')),
         fertilization = factor(fertilization , levels = c('CK','HO','HF','MF'))) %>%
  arrange(residue,fertilization)

respiration_co2 <- respiration22 %>%
  select(!NO2) %>%
  mutate(residue = factor(residue,levels = c('remove','turnover')),
         stage = factor(stage,levels = c('Jointing','Filling','Maturing')),
         fertilization = factor(fertilization,levels = c('ck','high_organic','high_fertilizer','mineral_fertilizer'))) %>%
  arrange(stage,residue,fertilization)

physical22 <- read.delim('22_physicochemical/22physicochemical.txt',header = T)
physical22$residue <- factor(physical22$residue,levels = c('Remove','Turnover'))
physical22$stage <- factor(physical22$stage,levels = c('Jointing','Filling','Maturing'))
physical22$fertilization <- factor(physical22$fertilization , 
                                   levels = c('CK','HO','HF','MF')) 

alpha_total_rep4 <- cbind(alpha_bac_plot[,1:2],fungi_richness = alpha_fungi_plot$Richness) %>%
  filter(str_detect(sample_name,'1|2|3|4')) %>%
  left_join(alpha_proto_plot[,1:2],by = 'sample_name')

physical22_sem <- physical22 %>% 
  arrange(stage,residue,fertilization) %>%
  as.data.frame()
dim(physical22_sem)
sem_data <- cbind(physical22_sem[,-c(1:3)],bac_shannon = alpha_total_rep4$Shannon.x,CO2 = respiration_co2$CO2,
                  fungi_shannon = alpha_total_rep4$fungi_shannon,
                  weight = biomass_22$weight,height = biomass_22$height) %>%
  as.data.frame() %>%
  scale() %>%
  cbind(physical22_sem[,1:3],.)
colnames(sem_data)

remove_data <- sem_data %>% 
  filter(residue =='Remove') 
remove_model <- lmer(weight ~ NH4+NO3+AP+AK+pH+MBC+MBN+DOC+DON+TN+TC+SOC+bac_shannon+fungi_shannon+
                      (1|stage)+(1|fertilization),data = remove_data) 
step_model_remove <- step(remove_model)
final_model_remove <- get_model(step_model_remove)
vif(final_model_remove)
summary(final_model_remove)
remove_model_empty <- lmer(weight ~ (1|fertilization),data = remove_data)
remove_model_AK <- lmer(weight ~ AK+(1|fertilization),data = remove_data) 
remove_model_SOC <- lmer(weight ~ SOC+(1|fertilization),data = remove_data) 

r.squaredGLMM(remove_model_empty)[1]
r.squaredGLMM(remove_model_AK)
r.squaredGLMM(remove_model_SOC)

shapiro.test(residuals(final_model_remove))

ggplot(remove_data,aes(x = SOC,y =weight))+
  geom_point()+
  geom_smooth(method = "lm")

turnover_data <- sem_data %>% 
  filter(residue =='Turnover') 
turnover_model <- lmer(weight ~ NH4+NO3+AP+AK+pH+MBC+MBN+DOC+DON+TN+TC+SOC+bac_shannon+fungi_shannon+
                       (1|stage)+(1|fertilization),data = turnover_data) 
step_model_turnover <- step(turnover_model)
final_model_turnover <- get_model(step_model_turnover)
vif(final_model_turnover)
summary(final_model_turnover)
shapiro.test(residuals(final_model_turnover))

turnover_model_empty <- lmer(weight ~ (1|fertilization),data = turnover_data)
turnover_model_NO3 <- lmer(weight ~ NO3+(1|fertilization),data = turnover_data) 
turnover_model_AK <- lmer(weight ~ AK+(1|fertilization),data = turnover_data) 
turnover_model_DON <- lmer(weight ~ DON+(1|fertilization),data = turnover_data) 
turnover_model_TN <- lmer(weight ~ TN+(1|fertilization),data = turnover_data) 
r.squaredGLMM(turnover_model_empty)
r2_nakagawa(turnover_model_AK)
r2_nakagawa(turnover_model_DON)
r2_nakagawa(turnover_model_TN)

ggplot(turnover_data,aes(x = AK,y =weight))+
  geom_point()+
  geom_smooth(method = "lm")

ggplot(sem_data)+
  geom_boxplot(aes(x = fertilization,y = weight,fill = residue))

# %>%
#   mutate(residuuuue = ifelse(residue=='Remove',0,1),
#          fertilizatiooon = ifelse(fertilization == 'CK',0,
#                                   ifelse(fertilization == 'HO',1,
#                                          ifelse(fertilization == 'HF',2,3)))) %>%
#   filter(stage =='Maturing')
head(sem_data)
rfm_harvest <- sem_data %>%
  filter(residue =='Turnover'&stage=='Maturing') %>%
  select(!c(stage:fertilization,TC,MBC,SOC,AK,AP,DOC))
rfm_harvest
harvest_weight_lm <- lm(height~ fungi_shannon,data=rfm_harvest)
summary(harvest_weight_lm)

randomForest(weight~.,rfm_harvest, importance = TRUE,ntree = 500)

summary(RFMdata_bac)
importance(RFMdata_bac)[,1:2]

bac_shannon_rfs_putdata <- importance(RFMdata_bac, sort.by = NULL, decreasing = TRUE)


bac_shannon_RFM_p <- bac_shannon_rfs_putdata %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*",""))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(desc(X.IncMSE)) %>%
  mutate(group = if_else(label=="","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>%
  ggplot(aes(x = names, y = X.IncMSE))+
  geom_bar(aes(fill = label),stat = "identity")+
  scale_fill_manual(values = cols_rfm)+
  theme_changwu_depth+
  theme(axis.text.x = element_text(angle = 315,hjust = -0.1))+
  geom_text(aes(y = X.IncMSE + 1,label = label),size=10)+
  labs(x = NULL, y = NULL)+
  annotate("text", x = 9, y = 20, label = 'Bacteria\nShannon' , colour = "black",size=6) +
  geom_label( x = 10.7, y = 19, label = expression('R'^2==.62),
              colour = "black",label.size=0,size=5)










ggplot(sem_data)+
  geom_boxplot(aes(x = fertilization,y = weight,fill = residue))

colnames(sem_data)

sem_weight <- lm(weight~NH4+NO3+AP+AK+TN+SOC,data = sem_data)
sem_height <- lm(height~NH4+NO3+AP+AK+TN+SOC,data = sem_data)
sem_weight_microbe <- lm(weight~DOC+DON+MBC+MBN+CO2,data = sem_data)
sem_height_microbe <- lm(height~DOC+DON+MBC+MBN+CO2,data = sem_data)
vif(sem_weight_microbe)
vif(sem_height_microbe)

sem_data2 <- sem_data %>%
  select(!TC)

biomass_retention_path <- list(
  residue = 'residuuuue',
  fertilization = 'fertilizatiooon',
  nutrients = c('AK','TN','SOC'),
  microbial_nutrients = c('DOC','CO2'),
  pH = 'pH',
  bac_diversity = 'bac_shannon',
  fungi_diversity = 'fungi_shannon',
  weight = 'weight',
  height = 'height'
)

modessss_retention <- rep('A',9)
residue_retention <- c(0,0,1,1,1,1,1,1,1)
fertilization_retention <- c(0,0,0,1,1,1,1,1,1)
nutrients_retention <- c(0,0,0,1,1,1,1,1,1)
microbial_nutrients <- c(0,0,0,0,1,1,1,1,1)
pH <- c(0,0,0,0,0,1,1,1,1)
bac_diversity_retention <-  c(0,0,0,0,0,0,0,1,1)
fungi_diversity_retention <-  c(0,0,0,0,0,0,0,1,1)
weight <- c(0,0,0,0,0,0,0,0,0)
height <- c(0,0,0,0,0,0,0,0,0)

path_retention <- cbind(residue_retention,fertilization_retention,nutrients_retention,microbial_nutrients,
                        pH,bac_diversity_retention,fungi_diversity_retention,weight,height)

rownames(path_retention)<-colnames(path_retention)

retention_plspm <- plspm(sem_data2,path_retention,biomass_retention_path,modes = modessss_retention)

summary(retention_plspm)
plot(retention_plspm)


biomass_biomass_path <- list(
  residue = 'residuuuue',
  fertilization = 'fertilizatiooon',
  nutrients = c('AK','TN','SOC'),
  microbial_nutrients = c('DOC','CO2'),
  pH = 'pH',
  bac_diversity = 'bac_shannon',
  fungi_diversity = 'fungi_shannon',
  maize_mass = c('weight','height')
)

modessss_biomass <- rep('A',8)
residue_biomass <- c(0,0,1,1,1,1,1,1)
fertilization_biomass <- c(0,0,0,1,1,1,1,1)
nutrients_biomass <- c(0,0,0,1,1,1,1,1)
microbial_nutrients_biomass <- c(0,0,0,0,1,1,1,1)
pH <- c(0,0,0,0,0,1,1,1)
bac_diversity_biomass <-  c(0,0,0,0,0,0,0,1)
fungi_diversity_biomass <-  c(0,0,0,0,0,0,0,1)
maize_mass <- c(0,0,0,0,0,0,0,0)


path_biomass <- cbind(residue_biomass,fertilization_biomass,nutrients_biomass,microbial_nutrients_biomass,
                        pH,bac_diversity_biomass,fungi_diversity_biomass,maize_mass)

rownames(path_biomass)<-colnames(path_biomass)

maize_mass_plspm <- plspm(sem_data2,path_biomass,biomass_biomass_path,modes = modessss_biomass)

summary(maize_mass_plspm)
plot(retention_plspm)