setwd("D:/Study")

library(tidyverse)
library(stringr)
library(vegan)

brite = read_tsv("0.Module.BRITE.txt") %>% select(-L3)
agg = read_tsv("enrichM/aggregate_output.tsv") %>%
  left_join(brite)
agg$L1[is.na(agg$L1)] = "others"
agg$L2[is.na(agg$L2)] = "others"

# summary the abundance of each category under L2
for(i in c(2:(length(agg)-3))){
  sum = agg[c(i,261)] %>% 
    mutate(now=as_vector(agg[,i])) %>%
    group_by(L2) %>% 
    summarise(now=sum(now))
  names(sum)[2] = names(agg)[i]
  if(i==2){result=sum}
  else{result = full_join(result,sum)}
}
write_tsv(result, "5.KEGG.module.L2.tsv")

# plot
ab <- read_tsv("5.KEGG.module.L2.relativeAbundance.tsv") %>%
  pivot_longer(-L2, names_to = "genome", values_to = "abundance") %>%
  mutate(group = genome)
ab$group[str_detect(ab$group,"A2O")] <- "A2O"
ab$group[str_detect(ab$group,"AO")] <- "AO"

pdf("5.module.meanAbundance.pdf",wi=2,he=2)
ggplot(ab %>%
         group_by(L2,group) %>%
         summarise(mean_abundance=mean(abundance)),
       aes(x=group, y=mean_abundance, fill=L2)) +
  geom_col() +
  theme_bw() +
  scale_fill_brewer(palette = "Set3")+
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.background = element_blank(),
        text = element_text(size = 4),
        panel.border = element_rect(size=.5))
dev.off()

# count unique module in AO or A2O
now = "Energy metabolism"
now = "Carbohydrate metabolism"
now = "Amino acid metabolism"
now = "Biosynthesis of other secondary metabolites"
now = "Biosynthesis of terpenoids and polyketides"
now = "Glycan metabolism"
now = "Lipid metabolism"
now = "Metabolism of cofactors and vitamins"
now = "Nucleotide metabolism"
now = "Xenobiotics biodegradation"
eneg <- filter(agg,L2==now) %>%
  select(-c(name,L1,L2)) %>%
  pivot_longer(-module, names_to = "genome", values_to = "abundance") %>%
  mutate(group = genome)
eneg$group[str_detect(eneg$group,"A2O")] <- "A2O"
eneg$group[str_detect(eneg$group,"AO")] <- "AO"
eneg_unique<- eneg %>%
  group_by(module,group) %>%
  summarise(sum_abundance=sum(abundance)) %>%
  pivot_wider(names_from = group,values_from = sum_abundance) %>%
  filter(A2O!=0|AO!=0)
ao <- eneg_unique$AO[eneg_unique$AO!=0] %>% length()
a2o <- eneg_unique$A2O[eneg_unique$A2O!=0] %>% length()
share <- eneg_unique$module[eneg_unique$A2O!=0&eneg_unique$AO!=0] %>% length()
draw.pairwise.venn(ao,a2o,share,category = c("AO", "A2O"),
                   fill = c("#4f94cd", "#c67171"))
par(mfrow=c(1,3), mar=c(4,4,2,1))