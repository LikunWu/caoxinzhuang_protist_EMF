#调用R包
library(multifunc)
library(tidyverse)
library(patchwork)
library(purrr)
library(forcats)
library(car)
library(scales)

data_caoxinzhuang_turnover <- read.csv('emf_shannon.csv',row.names = 1) %>%
  filter(residue=='Remove') %>%
  select(!c(stage,fertilization,residue,pH,  "AP","NH4", "NO3", "AK",'TC','TN','SOC', "DON", "DOC"))

allVars<-c("MBN",'BX','CBH','BG','CO2','pathogen','mycorrhizoa', "MBC")

bac_retention <- data_caoxinzhuang_turnover %>%
  mutate(Diversity = c(scale(bac_shannon)+scale(fungi_shannon)+scale(protozoa_shannon))) %>%
  select(!c(bac_shannon,fungi_shannon,protozoa_shannon)) %>%
  select(Diversity,everything())

#筛选我们纳入的数据集中，删除掉三分之二以上没有NA的数据行，这一步我们其实也没用，完全就是为了重新定义一下名字
german_vars <- whichVars(bac_retention, allVars)

# 再次查看功能名称
germany <- bac_retention
species <- relevantSp(bac_retention,2:ncol(bac_retention))
#这个数字是根据我们想纳入的功能在原始数据表里的顺序来定的，如，本次是3

#---------------------------------------均值法---------------------------------------------------
germany<-cbind(germany, getStdAndMeanFunctions(germany, german_vars))
germany
# write.csv(germany, 'caoxinzhuang_results/5.0/intermediate/EMF.csv')  

aveFit<-lm(meanFunction ~ Diversity, data=germany)
Anova(aveFit)
summary(aveFit)#看Adjusted R-squared和p-value  就是多样性和EMF的R2和p

# 绘制回归图
ggplot(aes(x = Diversity, y = meanFunction), data = germany) +
  geom_point(size = 3) + 
  theme_bw(base_size = 15) + 
  stat_smooth(method = "lm", colour = "black", size = 2) + 
  xlab("Diversity") +  ylab("multifunctional") + 
  # 添加R²和p值的文本框 
  annotate("text", x = max(germany$Diversity), y = max(germany$meanFunction), 
           label = paste0("R² = ", round(summary(aveFit)$adj.r.squared, 3), "\n",
                          "P  ", 
                          ifelse(summary(aveFit)$coefficients[2, 4] < 0.001, "<0.001",  
                                 ifelse(summary(aveFit)$coefficients[2, 4] < 0.01, "<0.01", 
                                        ifelse(summary(aveFit)$coefficients[2, 4] < 0.05, "<0.05",
                                               round(summary(aveFit)$coefficients[2, 4], 3))))), 
           hjust = 1, vjust = 1, size = 5, fontface = "bold")

#---------------------------------------单阈值法---------------------------------------------------
#单阈值法的所有代码，都是无脑跑，不需要改任何地方#单阈值法不能像均值法一样，提取出一列数据作为多功能性，这个单阈值法的y轴是很多线，应该是很多行的单功能性数值#进行单阈值分析前，需要创建一个新的数据集，如下
germanyThresh<-getFuncsMaxed(germany, german_vars, threshmin=0.05, threshmax=0.99,
                             prepend=c("plot","Diversity"), maxN=7)
#我们仅对阈值0.8的数据进行分析
mfuncGermanyLinear08<-glm(funcMaxed ~ Diversity, data=subset(germanyThresh, 
                           germanyThresh$thresholds=="0.2"), 
                          family=quasipoisson(link="identity"))
mfuncGermanyLinear08
Anova(mfuncGermanyLinear08, test.statistic="F")
summary(mfuncGermanyLinear08)
gcPlot<-subset(germanyThresh, germanyThresh$thresholds %in% qw(0.1,0.25, 0.5, 0.75, 0.9)) #note, using qw as %in% is a string comparison operatorPlot$percent<-paste(100*gcPlot$thresholds, "%", sep="")
gcPlot$percent<-paste(100*gcPlot$thresholds, "%", sep="")

  gcPlot %>%
  ggplot(aes(x = Diversity,y = funcMaxed,col = percent))+
  stat_smooth(method="glm",method.args = list(family=quasipoisson(link="identity")),se=FALSE)

#然后我们看看整体阈值的变化
germanyThresh$percent <- 100*germanyThresh$thresholds
germanyThresh
# 创建一个空白的数据框，用于绘制图例   
# 值得注意的是，图例并不能和下面的实际图在R里做在一起，所以需要手动AI
legend_data <- data.frame(x = 1:7, y = 1, percent = seq(0, 1, length.out = 7))
ggplot(legend_data, aes(x, y, color = percent)) + 
  geom_point(size = 5) + 
  scale_color_gradientn(colours = rev(c("red", "orange", "yellow", "green", "cyan", "blue", "purple")), 
                        breaks = seq(0, 1, by = 0.2), labels = scales::percent(seq(0, 1, by = 0.2)),
                        guide = guide_legend(override.aes = list(size = 5))) + 
  theme_void()

ggplot(data = germanyThresh, aes(x = Diversity, y = funcMaxed, group = percent)) +
  ylab(expression("Number of Functions" >= Threshold)) +
  xlab("MultiNet") +  
  stat_smooth(method = "glm",
              method.args = list(family = quasipoisson(link = "identity")), lwd = 0.8, fill = NA,
              aes(color = percent)) +  theme_bw(base_size = 14) +
  scale_color_gradientn(name = "Percent of \nMaximum",
                        colours = rev(c("red", "orange", "yellow", "green", "cyan", "blue", "purple")), 
                        breaks = seq(0, 1, by = 0.1), labels = scales::percent(seq(0, 1, by = 0.1)),
                        guide = guide_colorbar(reverse = FALSE))



#---------------------------------------多阈值法--------------------------------------------------
germanyLinearSlopes<-getCoefTab(funcMaxed ~ Diversity, data=germanyThresh,
                                 coefVar="Diversity", family=quasipoisson(link="identity"))
germanyLinearSlopes
####### Plot the values of the diversity slope at# different levels of the threshold######
germanSlopes <- ggplot(germanyLinearSlopes, aes(x = thresholds * 100,
                                                y = estimate,
                                                ymax = estimate + 1.96 * std.error, 
                                                ymin = estimate - 1.96 * std.error)) + 
  geom_ribbon(fill = "#FFCC99") +  
  geom_point(color = "#993399") +  # 将深绿色改为紫色（使用16进制颜色代码）
  ylab("Change in Number of Functions\n") + 
  xlab("\nThreshold (%)") +  
  geom_abline(intercept = 0, slope = 0, lwd = 1, linetype = 2) + 
  theme_bw(base_size = 14)

germanSlopes

germanIDX <- getIndices(germanyLinearSlopes, germanyThresh, funcMaxed ~ Diversity)
germanIDX
germanyLinearSlopes$estimate[which(germanyLinearSlopes$thresholds==germanIDX$Tmde)]

germanyThresh$IDX <- 0
germanyThresh$IDX [which(germanyThresh$thresholds %in%c(germanIDX$Tmin, germanIDX$Tmax, germanIDX$Tmde))] <- 1
germanyThresh$funcMaxed

germanyThresh$meanFunction
ggplot(data = germanyThresh, aes(x = Diversity, y = funcMaxed)) + 
  xlab("MultiNet") +
  geom_smooth(method = "glm",method.args = list(family = quasipoisson(link = "identity")),
                  fill = NA, aes(lwd = IDX)  ) + 
  theme_bw(base_size = 14) +
  scale_color_gradientn(   
    name = "Percent of \nMaximum",   
    colours = rev(c("red", "orange", "yellow", "green", "cyan", "blue", "purple")),  
    values = rescale(c(0, 100), c(0, 1))  ) +  scale_size(range = c(0.3, 5), guide = "none")+ 
  annotate(geom = "text", x = -0.7, y = c(0.2, 2, 4.6), label = c("Tmax", "Tmde", "Tmin")) +  
  annotate(geom = "text", x = 1, y = c(germanIDX$Mmin, germanIDX$Mmax, germanIDX$Mmde), 
           label = c("Mmin", "Mmax", "Mmde"))

germanSlopes +
  annotate(geom = "text", y = c(-0.01, -0.01, -0.01, germanIDX$Rmde.linear + 0.02), 
                          x = c(germanIDX$Tmin * 100, germanIDX$Tmde * 100, germanIDX$Tmax * 100,
                                germanIDX$Tmde * 100),      
           label = c("Tmin", "Tmde", "Tmax", "Rmde"), color = "black", vjust = -0.5)

