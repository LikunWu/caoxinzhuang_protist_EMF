library(tidyverse)
library(randomForest)
library(piecewiseSEM)
library(nlme)
library(lme4)

harvest_emf <- read.csv('caoxinzhuang_results/5.0/intermediate/RemoveEMF.csv',header = T)
retention_emf <- read.csv('caoxinzhuang_results/5.0/intermediate/TurnoverEMF.csv',header = T)

emf_bind <- rbind(harvest_emf,retention_emf)
physical22_sem
piecewise_sem_data <- cbind(physical = physical_latent,total = emf_bind$total,bacteria = emf_bind$bac_shannon_scale,
                            fungi = emf_bind$fungi_shannon_scale,protists = emf_bind$protozoa_shannon_scale,
                            emf = emf_bind$meanFunction,cross_kingdom_relation = scale_01_function(network_emf$protist_cross),
                            residue = physical22_2$residue,stage = physical22_sem$stage,
                            fertilization = physical22_sem$fertilization) %>% as.data.frame()

piecewise_sem_data
model <- psem(lme(emf ~ cross_kingdom_relation+protists+bacteria+fungi, random = list(stage = ~1|stage,fertilization = ~1|fertilization),data = piecewise_sem_data),
              lme(cross_kingdom_relation ~ bacteria+protists+residue+fungi, random = list(stage = ~1|stage,fertilization = ~1|fertilization),data = piecewise_sem_data),
              lme(protists ~ residue, random = list(stage = ~1|stage,fertilization = ~1|fertilization),data = piecewise_sem_data),
              lme(bacteria ~ residue, random = list(stage = ~1|stage,fertilization = ~1|fertilization),data = piecewise_sem_data),
              lme(fungi ~ residue, random = list(stage = ~1|stage,fertilization = ~1|fertilization),data = piecewise_sem_data),
              fungi %~~% protists)
summary(model)
sem_summary <- summary(model)

std_coef <- sem_summary$coefficients
std_coef
std_coef <- std_coef[std_coef$Response != "Residual", ]  
std_coef


calc_indirect_effect <- function(path_coefficients, path_chain) {
  effect <- 1
  for (i in 1:(length(path_chain)-1)) {
    from <- path_chain[i]
    to <- path_chain[i+1]
    coef <- path_coefficients$Estimate[path_coefficients$Predictor == from & path_coefficients$Response == to]
    if (length(coef) == 0) return(0)  
    effect <- effect * coef
  }
  return(effect)
}

# 提取所有变量名（排除emf本身）
variables <- c("cross_kingdom_relation", "protists", "bacteria", "fungi", "residue")

# 初始化总效应存储
total_effects <- data.frame(
  Variable = variables,
  Direct = numeric(length(variables)),
  Indirect = numeric(length(variables)),
  Total = numeric(length(variables))
)

# 填充直接效应
for (i in 1:nrow(total_effects)) {
  var <- total_effects$Variable[i]
  direct_coef <- std_coef$Estimate[std_coef$Response == "emf" & std_coef$Predictor == var]
  total_effects$Direct[i] <- ifelse(length(direct_coef) > 0, direct_coef, 0)
}

# 计算间接效应（需要根据模型结构手动定义路径）
# 示例：residue的间接效应
paths <- list(
  residue = list(
    c("residue", "protists", "emf"),                 # residue → protists → emf
    c("residue", "bacteria", "emf"),                 # residue → bacteria → emf
    c("residue", "fungi", "emf"),                    # residue → fungi → emf
    c("residue", "protists", "cross_kingdom_relation", "emf"),  # residue → protists → cross → emf
    c("residue", "bacteria", "cross_kingdom_relation", "emf"),  # residue → bacteria → cross → emf
    c("residue", "fungi", "cross_kingdom_relation", "emf")       # residue → fungi → cross → emf
  ),
  protists = list(
    c("protists", "cross_kingdom_relation", "emf")    # protists → cross → emf
  ),
  bacteria = list(
    c("bacteria", "cross_kingdom_relation", "emf")    # bacteria → cross → emf
  ),
  fungi = list(
    c("fungi", "cross_kingdom_relation", "emf")       # fungi → cross → emf
  )
)

# 计算每个变量的间接效应
for (var in names(paths)) {
  total_indirect <- 0
  for (path in paths[[var]]) {
    total_indirect <- total_indirect + calc_indirect_effect(std_coef, path)
  }
  idx <- which(total_effects$Variable == var)
  total_effects$Indirect[idx] <- total_indirect
}

total_effects$Total <- total_effects$Direct + total_effects$Indirect
total_effects_final <- total_effects %>%
  arrange(desc(Total)) %>%
  mutate(Variable2 = factor(Variable,levels = c('residue','bacteria','fungi','protists','cross_kingdom_relation'))) %>%
  select(!Variable)%>%
  pivot_longer(!Variable2,names_to = 'direction',values_to = 'Effects')
total_effects_final
sem_total_effects_p <- ggplot(total_effects_final)+
  geom_col(aes(x = Variable2,y = Effects,fill = direction),position = position_dodge(width = 0.8), width = 0.7)+
  scale_fill_manual(values = c('#A8DDE3','#8D3C90','#FAEBBD'))+
  labs(x = NULL,y = 'Standard effects')+
  theme_NMDS2

ggsave('caoxinzhuang_results/6.0emf/SEM/sem_total_effects_p.pdf',sem_total_effects_p,width  = 6,height = 3)


####VPA####
vpa_data <- emf_anova %>%
  select(MBN:protists_shannon) %>%
  cbind(cross_trophic = piecewise_sem_data$cross_kingdom_relation) %>%
  apply(.,2,scale_01_function) %>%
  cbind(residue = piecewise_sem_data$residue) %>%
  as.data.frame()
vpa_harvest <- vpa_data %>% filter(residue==1)


??varpart
vpa_harvest_p <- varpart(
  Y = select(vpa_harvest,MBN:mycorrhizoa), 
  X = vpa_harvest$cross_trophic,
  X2 = vpa_harvest$protists_shannon ,      
  X3 = vpa_harvest$fungi_shannon , 
  X4 = vpa_harvest$bac_shannon 
)
rda_model_harvest <- rda(
  formula = select(vpa_harvest, MBN:mycorrhizoa) ~ 
    cross_trophic + protists_shannon + fungi_shannon + bac_shannon,
  data = vpa_harvest
)
set.seed(123) # 设置随机种子保证结果可重复
anova_result_vpa_harvest <- anova(
  object = rda_model_harvest, 
  by = "margin",    # 进行边际效应检验
  permutations = 9999 # 建议置换次数 >= 999
)
vpa_harvest_p
pdf("caoxinzhuang_results/6.0emf/SEM/VPA_Plot_harvest.pdf",     
    width = 4,  
    height = 4,  
    family = "sans") 

plot(vpa_harvest_p,
     bg = c("#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C"),  # 自定义颜色
     Xnames = c("Protist related interactions", "Protist", "Fungi", "Bacteria"),
     main = "Variance Partitioning (4 Groups)"
) 

dev.off()

vpa_retention <- vpa_data %>% filter(residue==2)

vpa_retention_p <- varpart(
  Y = select(vpa_retention,MBN:mycorrhizoa), 
  X = vpa_retention$cross_trophic,
  X2 = vpa_retention$protists_shannon ,      
  X3 = vpa_retention$fungi_shannon , 
  X4 = vpa_retention$bac_shannon 
)

rda_model_retention <- rda(
  formula = select(vpa_retention, MBN:mycorrhizoa) ~ 
    cross_trophic + protists_shannon + fungi_shannon + bac_shannon,
  data = vpa_retention
)
set.seed(123) # 设置随机种子保证结果可重复
anova_result_vpa_retention <- anova(
  object = rda_model_retention, 
  by = "margin",    # 进行边际效应检验
  permutations = 9999 # 建议置换次数 >= 999
)
anova_result_vpa_retention
pdf("caoxinzhuang_results/6.0emf/SEM/VPA_Plot_retention.pdf",     
    width = 4,  
    height = 4,  
    family = "sans") 
plot(vpa_retention_p,
     bg = c("#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C"),  # 自定义颜色
     Xnames = c("Protist related interactions", "Protist", "Fungi", "Bacteria"),
     main = "Variance Partitioning (4 Groups)"
) 

dev.off()
head(vpa_data)

Y_matrix <- vpa_data[, c("MBN", "MBC", "CO2", "BG", "BX", "CBH", "pathogen", "mycorrhizoa")]
head(vpa_data)
