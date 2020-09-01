setwd("D:/FTICR")

library(tidyverse)
IN=read.table("Influent.tsv",sep = "\t",header = T, stringsAsFactors = F)
PU=read.table("AGC.tsv",sep = "\t",header = T, stringsAsFactors = F)
A2O_PRE60=read.table("A2O_PRE.tsv",sep = "\t",header = T, stringsAsFactors = F)
A2O_A1=read.table("A2O_A1.tsv",sep = "\t",header = T, stringsAsFactors = F)
A2O_A2=read.table("A2O_A2.tsv",sep = "\t",header = T, stringsAsFactors = F)
A2O_o=read.table("A2O_O.tsv",sep = "\t",header = T, stringsAsFactors = F)
A2O_SED=read.table("A2O_SED.tsv",sep = "\t",header = T, stringsAsFactors = F)
AO_Ab=read.table("AO_A.tsv",sep = "\t",header = T, stringsAsFactors = F)
AO_O=read.table("AO_O.tsv",sep = "\t",header = T, stringsAsFactors = F)
AO_SED=read.table("AO_SED.tsv",sep = "\t",header = T, stringsAsFactors = F)
MIX=read.table("MIX.tsv",sep = "\t",header = T, stringsAsFactors = F)
OUT=read.table("Effluent.tsv",sep = "\t",header = T, stringsAsFactors = F)

cols <- c("CHO" = "#cd4f39",
          "CHNO" = "#1874cd",
          "CHOS" = "#eec900",
          "CHOP" = "#00cd00",
          "CHNOS" = "#ff00ff",
          "CHNOP" = "#ff7f00",
          "CHOSP" = "#43cd80",
          "CHNOSP" = "#ffb5c5",
          "Cl-containing" = "#878787")

# formula classes
for(now in list){
now$group[now$group=="CHNOCl"|now$group=="CHNOPCl"|now$group=="CHNOSCl"|now$group=="CHOPCl"|now$group=="CHOCl"|now$group=="CHOSCl"] = "Cl-containing"
ggplot(now,aes(x=O_C,y=H_C,color=group))+
  geom_rect(xmin = 0, xmax = 0.2, ymin = 1.5, ymax = 2.3,fill="white",alpha=0, color="black")+
geom_rect(xmin = 0.2, xmax = 0.52, ymin = 1.5, ymax = 2.2,fill="white",alpha=0, color="black")+
geom_rect(xmin = 0.52, xmax = 0.7, ymin = 1.5, ymax = 2.2,fill="white",alpha=0, color="black")+
geom_rect(xmin = 0.7, xmax = 1.1, ymin = 1.5, ymax = 2.4,fill="white",alpha=0, color="black")+
geom_rect(xmin = 0.25, xmax = 0.67, ymin = 0.75, ymax = 1.5,fill="white",alpha=0, color="black")+
geom_rect(xmin = 0.67, xmax = 0.97, ymin = 0.53, ymax = 1.5,fill="white",alpha=0, color="black")+
geom_rect(xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1.25,fill="white",alpha=0, color="black")+
  geom_point(size=0.6)+
  theme_bw()+
  scale_color_manual(values = cols)+
  scale_x_continuous(limits = c(0,1.15))+
  scale_y_continuous(limits = c(0,2.5),expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5))+
  theme(rect = element_rect(fill=NULL))

CHO=filter(now,group=="CHO")%>% dim()
CHNO=filter(now,group=="CHNO")%>% dim()
CHOP=filter(now,group=="CHOP")%>% dim()
CHOS=filter(now,group=="CHOS")%>% dim()
CHNOS=filter(now,group=="CHNOS")%>% dim()
CHNOP=filter(now,group=="CHNOP")%>% dim()
CHNOSP=filter(now,group=="CHNOSP")%>% dim()
}
# compound classes
cols <- c("Lipids" = "#cd4f39",
          "Proteins" = "#1874cd",
          "Amino_sugar" = "#eec900",
          "Condensed_aromatics" = "#00cd00",
          "Carbohydrate" = "#ff00ff",
          "Lignin" = "#ff7f00",
          "Tannin" = "#43cd80",
          "Others" = "#878787")

lipid = filter(now,H/C>=1.5&H/C<=2.3&O/C>=0&O/C<=0.2&N/C<=0.04&P/C<=0.03) %>% dim()
protein = filter(now,H/C>=1.5&H/C<=2.2&O/C>=0.2&O/C<=0.52&N/C<=0.44&N/C>=0.178&P/C<=0.06) %>% dim()
amino_sugar = filter(now,H/C>=1.5&H/C<=2.2&O/C>=0.52&O/C<=0.7&N/C<=0.182&N/C>0.07&P/C<0.167) %>% dim()
carbohydrate = filter(now,H/C>=1.5&H/C<=2.4&O/C>=0.7&O/C<=1.1&N==0&P==0) %>% dim()
condensed_aromatics = filter(now,H/C>=0.5&H/C<=1.25&O/C>=0&O/C<=0.25) %>% dim()
lignin = filter(now,H/C>=0.75&H/C<=1.5&O/C>=0.25&O/C<=0.67) %>% dim()
tannin = filter(now,H/C>=0.53&H/C<=1.5&O/C>=0.67&O/C<=0.97) %>% dim()

combine_all = rbind(IN,PU,A2O_A1,A2O_A2,A2O_o,A2O_PRE60,A2O_SED,AO_Ab,AO_O,AO_SED,MIX,OUT) %>%
  select(formula,class,group,C,H,O,N,S,P,Cl,theor_mass,DBE,DBEO) %>%
  unique() %>%
  mutate(compound = case_when(H/C>=1.5&H/C<=2.3&O/C>=0&O/C<=0.2&N/C<=0.04&P/C<=0.03 ~ "Lipids",
                              H/C>=1.5&H/C<=2.2&O/C>=0.2&O/C<=0.52&N/C<=0.44&N/C>=0.178&P/C<=0.06 ~ "Proteins",
                              H/C>=1.5&H/C<=2.2&O/C>=0.52&O/C<=0.7&N/C<=0.182&N/C>0.07&P/C<0.167 ~ "Amino_sugar",
                              H/C>=1.5&H/C<=2.4&O/C>=0.7&O/C<=1.1&N==0&P==0 ~ "Carbohydrate",
                              H/C>=0.5&H/C<=1.25&O/C>=0&O/C<=0.25 ~ "Condensed_aromatics",
                              H/C>=0.75&H/C<=1.5&O/C>=0.25&O/C<=0.67 ~ "Lignin",
                              H/C>=0.53&H/C<=1.5&O/C>=0.67&O/C<=0.97 ~ "Tannin",
                              TRUE ~ "Others")) %>%
  mutate(AI=(1+C-0.5*O-S-0.5*N-0.5*P-0.5*H-0.5*Cl)/(C-0.5*O-N-S-P)) %>%
  mutate(NOSC=4-(4*C+H+Cl-2*O-2*S-3*N+5*P)/C)
write_tsv(combine_all,"formula_info.tsv")

# group compound ploting
formula2compound = select(combine_all,c(formula,compound))
G1 = read_tsv("group/G1.txt") %>% left_join(formula2compound)
G2 = read_tsv("group/G2.txt") %>% left_join(formula2compound)
G3 = read_tsv("group/G3.txt") %>% left_join(formula2compound)
G4 = read_tsv("group/G4.txt") %>% left_join(formula2compound)
G5 = read_tsv("group/G5.txt") %>% left_join(formula2compound)
G6 = read_tsv("group/G6.txt") %>% left_join(formula2compound)
G7 = read_tsv("group/G7.txt") %>% left_join(formula2compound)

now = G7
ggplot(now,aes(x=O_C,y=H_C,color=compound))+
    geom_point(size=0.6)+
    theme_bw()+
    scale_color_manual(values = cols)+
    scale_x_continuous(limits = c(0,1.15))+
    scale_y_continuous(limits = c(0,2.5),expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5))+
    theme(rect = element_rect(fill=NULL),
          legend.key.height = unit(5, "pt"),
          legend.key.width = unit(3, "pt"),
          legend.background = element_blank(),
          legend.title = element_blank(),
          text = element_text(size = 4))
ggsave(filename = "group/G7.compound.pdf",device = "pdf",
         width = 3,height = 2)

# plot compound in each group
g1 = G1 %>% group_by(compound) %>% summarise(n=n()) %>% print() %>% rename(G1=n)
g2 = G2 %>% group_by(compound) %>% summarise(n=n()) %>% print() %>% rename(G2=n)
g3 = G3 %>% group_by(compound) %>% summarise(n=n()) %>% print() %>% rename(G3=n)
g4 = G4 %>% group_by(compound) %>% summarise(n=n()) %>% print() %>% rename(G4=n)
g5 = G5 %>% group_by(compound) %>% summarise(n=n()) %>% print() %>% rename(G5=n)
g6 = G6 %>% group_by(compound) %>% summarise(n=n()) %>% print() %>% rename(G6=n)
g7 = G7 %>% group_by(compound) %>% summarise(n=n()) %>% print() %>% rename(G7=n)
compound.in.groups = full_join(g1,g2) %>%
  full_join(g3) %>% full_join(g4) %>% full_join(g5) %>%
  full_join(g6) %>% full_join(g7)

ggplot(compound.in.groups %>%
         pivot_longer(-compound,
                      names_to = "group",
                      values_to = "number"),
       aes(x=group,y=number,fill=compound)) +
  geom_col(position = "fill")+
  theme_bw()+
  scale_fill_manual(values = cols)+
  theme(rect = element_rect(fill=NULL),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 4))
ggsave(filename = "group/Compounds.in.groups.proportion.pdf",device = "pdf",
       width = 3,height = 2)

# plot group in each group
g1 = G1 %>% group_by(group) %>% summarise(n=n()) %>% print() %>% rename(G1=n)
g2 = G2 %>% group_by(group) %>% summarise(n=n()) %>% print() %>% rename(G2=n)
g3 = G3 %>% group_by(group) %>% summarise(n=n()) %>% print() %>% rename(G3=n)
g4 = G4 %>% group_by(group) %>% summarise(n=n()) %>% print() %>% rename(G4=n)
g5 = G5 %>% group_by(group) %>% summarise(n=n()) %>% print() %>% rename(G5=n)
g6 = G6 %>% group_by(group) %>% summarise(n=n()) %>% print() %>% rename(G6=n)
g7 = G7 %>% group_by(group) %>% summarise(n=n()) %>% print() %>% rename(G7=n)
group.in.groups = full_join(g1,g2) %>%
  full_join(g3) %>% full_join(g4) %>% full_join(g5) %>%
  full_join(g6) %>% full_join(g7) %>% rename(formula_group=group)

cols <- c("CHO" = "#cd4f39",
          "CHNO" = "#1874cd",
          "CHOS" = "#eec900",
          "CHOP" = "#00cd00",
          "CHNOS" = "#ff00ff",
          "CHNOP" = "#ff7f00",
          "CHOSP" = "#43cd80",
          "CHNOSP" = "#ffb5c5",
          "Cl-containing" = "#878787")

ggplot(group.in.groups %>%
         pivot_longer(-formula_group,
                      names_to = "group",
                      values_to = "number"),
       aes(x=group,y=number,fill=formula_group)) +
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = cols)+
  theme(rect = element_rect(fill=NULL),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 4))
ggsave(filename = "group/group.in.groups.absolute.pdf",device = "pdf",
       width = 3,height = 2)

M1.summary = M1 %>% group_by(group) %>% summarise(MW=mean(theor_mass),DBE=mean(DBE),O_C=mean(O_C),H_C=mean(H_C),NOSC=mean(NOSC),n=n(),nclass=length(unique(class)))
for(now_group in unique(M1$group)){
  M1.MLBL = sum(M1$H_C[M1$group==now_group]>1.5)/sum(M1$group==now_group)
  print(paste(now_group,"is",M1.MLBL))
}