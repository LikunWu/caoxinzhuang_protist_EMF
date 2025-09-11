library(tidyverse)
library(ggpubr)
library(rstatix)

enzyme_raw23 <- read.delim('23physicochemical/enzyme/enzyme_23.txt',header = T)
enzyme_raw22 <- read.delim('23physicochemical/enzyme/enzyme_22.txt',header = T)
enzyme_group <- read.delim('23physicochemical/enzyme/group_enzyme.txt',header = T)

fs_fr <- read.delim('23physicochemical/enzyme/FS_FR.txt',header = T,row.names = 1) #导入在excel里面计算好均值的ck表

enzyme_22 <- enzyme_raw22 %>%
  pivot_longer(!c(order,enzyme),names_to = 'col',values_to = 'value') %>%
  group_by(order,enzyme,col) %>%
  na.omit() %>% #从这里开始去除NA是因为na.omit命令是去除有NA值的那一行，如果在前面宽表时就去，会误删很多
  summarise(mean_value = mean(value)) %>% #截至到此，均值已计算完毕
  arrange(order,enzyme,col) %>%
  pivot_wider(names_from = 'enzyme',values_from = 'mean_value') %>% #准备数据开始计算酶活
  cbind(enzyme_group) %>% #此时将宽表和分组信息进行合并
  mutate(q_amc = (AMC - BUF)/ fs_fr[[10]], q_mub = (MUB - BUF)/ fs_fr[[9]],
         e_amc = fs_fr[[10]]/(10*0.00005), e_mub = fs_fr[[9]]/(10*0.00005)) %>% # 1.计算几个关键系数,计算公式见protocol  2.mutate命令为在后面加入一列新的数据，当然也可以加入多列，可以用来对宽表数据加入均值等
  mutate(f_ap = (AP - BUF)/q_mub - fs_fr[[1]] , f_app = (APP - BUF)/q_amc - fs_fr[[2]],
         f_asf = (ASF - BUF)/q_mub - fs_fr[[3]] , f_bg = (BG - BUF)/q_mub - fs_fr[[4]],
         f_lap = (LAP - BUF)/q_amc - fs_fr[[5]] , f_bx = (BX - BUF)/q_mub - fs_fr[[6]],
         f_cbh = (CBH - BUF)/q_mub - fs_fr[[7]] , f_nag = (NAG - BUF)/q_amc - fs_fr[[8]]) %>%
  mutate(u_AP = f_ap*50 / (e_mub*0.15*0.5*3),u_APP = f_app*50 / (e_amc*0.15*2*3),
         u_ASF = f_asf*50 / (e_mub*0.15*2*3),u_BG = f_bg*50 / (e_mub*0.15*2*3),
         u_LAP = f_lap*50 / (e_amc*0.15*2*3),u_BX = f_bx*50 / (e_mub*0.15*4*3),
         u_CBH = f_cbh*50 / (e_mub*0.15*4*3),u_NAG = f_nag*50 / (e_amc*0.15*4*3)) %>%
  select(starts_with('u_'),colnames(enzyme_group))

#作图代码，仅供参考
# enzyme_22 %>%
#   pivot_longer(!c(colnames(enzyme_group),order,rep),names_to = 'enzyme',values_to = 'U') %>%
#   filter(U>0)  %>%
#   filter(enzyme == 'u_AP'&U<0.15|enzyme == 'u_BG'&U<0.05|
#            enzyme == 'u_BX'&U<0.03|enzyme == 'u_CBH'&U<0.05) %>%
#   ggplot()+
#   geom_boxplot(aes(x = residue, y = U, fill = residue))+
#   geom_jitter(aes(x = residue, y = U, fill = residue),shape = 21,size = 1.5,
#               width = 0.5,alpha = 0.5)+
#   labs(x = 'Residue', y = 'Enzyme')+
#   scale_fill_manual(values = CAOXINZHUANG2)+
#   facet_wrap(.~enzyme,scales = 'free_y')+
#   theme_NMDS
enzyme_data_22 <- enzyme_22 %>%
    pivot_longer(!c(colnames(enzyme_group),order,rep),names_to = 'enzyme',values_to = 'U') %>%
  filter(enzyme == 'u_BG'|enzyme == 'u_BX'|enzyme == 'u_CBH'|enzyme == 'u_AP') %>%
  mutate(U2 = ifelse(U<0,0,U))
colnames(enzyme_data_22)

write.csv(enzyme_data_22,file = '22_physicochemical/22enzyme.csv')
unique(enzyme_22$u_AP)
dim(enzyme_data_22)
enzyme_data_22 %>%
  group_by(enzyme,fertilization) %>% 
  t_test(U~residue) 

ggplot(enzyme_data_22,aes(x = residue , y = U))+
  geom_boxplot(aes(fill = residue),alpha = 0.8,width = 0.8,size = 1,
               position = position_dodge2(width = 0.5))+
  geom_jitter(size = 3,alpha = 0.3,col = 'grey',
              position = position_dodge2(width = 0.5))+
  scale_fill_manual(values = CAOXINZHUANG2)+
  labs(x = 'Residue', y = 'U')+
  theme_NMDS+
  theme(axis.title.y = element_text(size = 30))+
  facet_wrap(.~enzyme,scales = 'free_y')

dim(enzyme_data_22)
enzyme_data_23 <- enzyme_23 %>%
  pivot_longer(!c(colnames(enzyme_group),order,rep),names_to = 'enzyme',values_to = 'U') %>%
  filter(U>0)  %>%
  filter(enzyme == 'u_AP'&U<0.15|enzyme == 'u_BG'&U<0.05|
           enzyme == 'u_BX'&U<0.05|enzyme == 'u_CBH'&U<0.05) 

enzyme_total <- rbind(enzyme_data_22,enzyme_data_23)

enzyme_total %>%
  filter(enzyme == 'u_CBH') %>%
  leveneTest(U~fertilization,data = .)


enzyme_total %>%
  filter(enzyme == 'u_AP') %>% 
  group_by(residue)%>%
  summarise(meanvalue = mean(U))


enzyme_total %>%
  filter(enzyme == 'u_CBH') %>% 
  aov(U~fertilization+residue:fertilization,data = .) %>%
  tukey_hsd()


for(indica in unique(enzyme_total$enzyme)) {
  data11 <- filter(enzyme_total,enzyme == indica)
  p1 <- ggplot(data11,aes(x = residue , y = U))+
    geom_boxplot(aes(fill = residue),alpha = 0.8,width = 0.8,size = 1,
                 position = position_dodge2(width = 0.5))+
    geom_jitter(size = 3,alpha = 0.3,col = 'grey',
                position = position_dodge2(width = 0.5))+
    scale_fill_manual(values = CAOXINZHUANG2)+
    labs(x = 'Residue', y = indica)+
    theme_NMDS+
    theme(axis.title.y = element_text(size = 30))
  p2 <- ggplot(data11,aes(x = fertilization , y = U))+
    geom_boxplot(aes(fill =  
                       fertilization),width =0.8,size = 2)+
    labs(x = 'Fertilization', y = NULL)+
    scale_fill_manual(values = sunsx)+
    scale_color_manual(values = sunsx)+
    theme_NMDS+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  pp <- ggarrange(p1,p2,widths = c(4,7),heights = 4)
  plot(pp)
  ggsave(filename = paste0(picture_dir,indica,'2.pdf'),pp,width = 8,height = 3)
}
picture_dir
