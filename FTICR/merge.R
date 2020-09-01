setwd("D:/OneDrive/Study/JXZ/New_Taxonomy/FTICR")

library(tidyverse)
library(vegan)
library(stringr)
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

combine <- select(IN,formula,normalization) %>% rename(In=normalization) %>%
  full_join(select(PU,formula,normalization)%>% rename(PU=normalization)) %>%
  full_join(select(A2O_PRE60,formula,normalization)%>% rename(A2O_PRE=normalization)) %>%
  full_join(select(A2O_A1,formula,normalization)%>% rename(A2O_A1=normalization)) %>%
  full_join(select(A2O_A2,formula,normalization)%>% rename(A2O_A2=normalization)) %>%
  full_join(select(A2O_o,formula,normalization)%>% rename(A2O_O=normalization)) %>%
  full_join(select(A2O_SED,formula,normalization)%>% rename(A2O_SED=normalization)) %>%
  full_join(select(AO_Ab,formula,normalization)%>% rename(AO_A=normalization)) %>%
  full_join(select(AO_O,formula,normalization)%>% rename(AO_O=normalization)) %>%
  full_join(select(AO_SED,formula,normalization)%>% rename(AO_SED=normalization)) %>%
  full_join(select(MIX,formula,normalization)%>% rename(MIX=normalization)) %>%
  full_join(select(OUT,formula,normalization)%>% rename(OUT=normalization))
combine[is.na(combine)] = 0
# write_tsv(combine, "combined.normalization.tsv")
row.names(combine) = combine$formula
combine = combine[,-1]
dist <- vegdist(data.frame(t(combine)))
hc = hclust(dist)
plot(hc)

AOvsA2O <- left_join(select(AO_SED,formula,X)%>%rename(AO=X),select(A2O_SED,formula,X)%>%rename(A2O=X)) %>%
  mutate(status="only") %>%
  filter(!is.na(AO)) %>%
  filter(AO>A2O|is.na(A2O))
AOvsA2O$status[!is.na(AOvsA2O$A2O)] = "more"
AOvsA2O <- AOvsA2O[-c(2,3)]
AOvsA2O <- left_join(AOvsA2O,AO_SED)

pdf("AOvsA2O_CHOS.pdf",wi=6,he=4)
ggplot(AOvsA2O %>% filter(group=="CHOS"),aes(x=O/C,y=H/C,color=status,shape=status))+
  geom_point(size=0.6)+
  theme_bw()+
  #scale_color_manual(values = cols)+
  scale_x_continuous(limits = c(0,1.15))+
  scale_y_continuous(limits = c(0,2.5),expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5))+
  theme(rect = element_rect(fill=NULL))
dev.off()

A2OvsAO <- left_join(select(A2O_SED,formula,X)%>%rename(A2O=X),select(AO_SED,formula,X)%>%rename(AO=X)) %>%
  mutate(status="only") %>%
  filter(!is.na(A2O)) %>%
  filter(A2O>AO|is.na(AO))
A2OvsAO$status[!is.na(A2OvsAO$AO)] = "more"
A2OvsAO <- A2OvsAO[-c(2,3)]
A2OvsAO <- left_join(A2OvsAO,A2O_SED)

pdf("A2OvsAO_CHOS.pdf",wi=6,he=4)
ggplot(A2OvsAO %>% filter(group=="CHOS"),aes(x=O/C,y=H/C,color=status,shape=status))+
  geom_point(size=0.6)+
  theme_bw()+
  #scale_color_manual(values = cols)+
  scale_x_continuous(limits = c(0,1.15))+
  scale_y_continuous(limits = c(0,2.5),expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5))+
  theme(rect = element_rect(fill=NULL))
dev.off()

for(now in list(IN,PU,A2O_PRE60,A2O_A1,A2O_A2,A2O_o,A2O_SED,AO_Ab,AO_O,AO_SED,MIX,OUT)){
  now = mutate(now,AI=(1+C-0.5*O-S-0.5*N-0.5*P-0.5*H-0.5*Cl)/(C-0.5*O-N-S-P))
  now = now[is.finite(now$AI),]
  print(mean(na.omit(now$AI)))
}

for(now in list(IN,PU,A2O_PRE60,A2O_A1,A2O_A2,A2O_o,A2O_SED,AO_Ab,AO_O,AO_SED,MIX,OUT)){
  now = mutate(now,NOSC=4-(4*C+H+Cl-2*O-2*S-3*N+5*P)/C)
  now = now[is.finite(now$NOSC),] 
  print(mean(na.omit(now$NOSC)))
}

# MLBL(%)
for(now in list(IN,PU,A2O_PRE60,A2O_A1,A2O_A2,A2O_o,A2O_SED,AO_Ab,AO_O,AO_SED,MIX,OUT)){
  sum_now = length(t(now))
  now = filter(now,H/C>=1.5) %>% t() %>% length()
  print(now/sum_now)
}



combine1 <- select(PU,formula,X) %>% rename(PU=X) %>%
  full_join(select(A2O_o,formula,X)%>% rename(A2O_O=X)) %>%
  full_join(select(AO_O,formula,X)%>% rename(AO_O=X))
combine1[is.na(combine1)] = 0
info = rbind(PU %>%
               mutate(AI=(1+C-0.5*O-S-0.5*N-0.5*P-0.5*H)/(C-0.5*O-N-S-P)) %>%
               mutate(NOSC=4-(4*C+H+Cl-2*O-2*S-3*N+5*P)/C) %>%
               select(c(formula, class, group, theor_mass,DBE,DBEO, O_C, H_C,AI,NOSC)),
             A2O_o %>%
               mutate(AI=(1+C-0.5*O-S-0.5*N-0.5*P-0.5*H)/(C-0.5*O-N-S-P)) %>%
               mutate(NOSC=4-(4*C+H+Cl-2*O-2*S-3*N+5*P)/C) %>%
               select(c(formula, class, group, theor_mass,DBE,DBEO, O_C, H_C,AI,NOSC)),
             AO_O %>%
               mutate(AI=(1+C-0.5*O-S-0.5*N-0.5*P-0.5*H)/(C-0.5*O-N-S-P)) %>%
               mutate(NOSC=4-(4*C+H+Cl-2*O-2*S-3*N+5*P)/C) %>%
               select(c(formula, class, group, theor_mass,DBE,DBEO, O_C, H_C,AI,NOSC))) %>%
  unique()
combine1 <- left_join(combine1, info)
combine1$group[str_detect(combine1$group,pattern = "Cl")] <- "Cl-containing"
write_tsv(combine1, "PU_O.tsv")
