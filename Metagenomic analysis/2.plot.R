setwd("D:/Study")
library(tidyverse)

tax <- read.delim("0.GTDB.tax.tsv", sep = '\t', header = T,
                  stringsAsFactors = FALSE)
tax$phylum[tax$phylum=="Proteobacteria"] <- tax$class[tax$phylum=="Proteobacteria"]
tax$phylum[tax$phylum=="Bdellovibrionota_B"|
               tax$phylum=="Bdellovibrionota_A"] <- "Bdellovibrionota"
tax$phylum[is.na(tax$phylum)] <- "Others"
tax <- select(tax,c(genome,phylum))
ao <- read_tsv("1.AO.genome.rpm.tsv") %>%
  left_join(tax) 
a2o <- read_tsv("1.A2O.genome.rpm.tsv") %>%
  left_join(tax) %>% group_by(phylum)

# summary AO's phylum
i=0
for(now in unique(tax$phylum)){
  kkk = ao %>%  filter(phylum==now) %>%
  summarise(AO_A_1=sum(AO_A_1), AO_A_2=sum(AO_A_2),
            AO_A_3=sum(AO_A_3),AO_O_1=sum(AO_O_1),
            AO_O_2=sum(AO_O_2),AO_O_3=sum(AO_O_3)) %>%
    data.frame()
  row.names(kkk)=now
  if(i==0){aaa = kkk}
  else{
   aaa = rbind(aaa,kkk) 
  }
  i= i+1
}
aaa = mutate(aaa,phylum=row.names(aaa)) %>% as_tibble()
write_tsv(aaa, "2.AO.phylum.tsv")

# summary A2O's phylum
i=0
for(now in unique(tax$phylum)){
  kkk = a2o %>%  filter(phylum==now) %>%
    summarise(A2O_pre_1=sum(A2O_pre_1), A2O_pre_2=sum(A2O_pre_2),
              A2O_pre_3=sum(A2O_pre_3),A2O_A1_1=sum(A2O_A1_1),
              A2O_A1_2=sum(A2O_A1_2),A2O_A1_3=sum(A2O_A1_3),
              A2O_A2_1=sum(A2O_A2_1), A2O_A2_2=sum(A2O_A2_2),
              A2O_A2_3=sum(A2O_A2_3),A2O_O_1=sum(A2O_O_1),
              A2O_O_2=sum(A2O_O_2),A2O_O_3=sum(A2O_O_3)) %>%
    data.frame()
  row.names(kkk)=now
  if(i==0){bbb = kkk}
  else{
    bbb = rbind(bbb,kkk) 
  }
  i= i+1
}
bbb = mutate(bbb,phylum=row.names(bbb)) %>% as_tibble()
write_tsv(bbb, "2.A2O.phylum.tsv")
ab <- left_join(bbb,aaa)
write_tsv(ab, "2.phylum.rpm.tsv")

# abundance
ab <- read_tsv("2.phylum.relativeAbundance.tsv") %>%
  arrange(-A2O_A1_1)
ab$phylum <- factor(ab$phylum,levels = ab$phylum)
ab <- pivot_longer(ab, -phylum, names_to = "sample",
               values_to = "abundance")
ab$sample <- factor(ab$sample, levels = c("A2O_pre_1",
                                          "A2O_pre_2",
                                          "A2O_pre_3",
                                          "A2O_A1_1",
                                          "A2O_A1_2",
                                          "A2O_A1_3",
                                          "A2O_A2_1",
                                          "A2O_A2_2",
                                          "A2O_A2_3",
                                          "A2O_O_1",
                                          "A2O_O_2",
                                          "A2O_O_3",
                                          "AO_A_1",
                                          "AO_A_2",
                                          "AO_A_3",
                                          "AO_O_1",
                                          "AO_O_2",
                                          "AO_O_3"))
# collapse
ab$phylum[ab$abundance<1.7] <- "Others"
ab$phylum[ab$phylum=="Verrucomicrobiota"] <- "Others"

pdf("2.phylum.abundance.pdf",wi=4,he=3)
ggplot(ab, aes(x = sample, y = abundance, fill = phylum)) +
  geom_col(position = "fill") + 
  theme_bw() + 
  scale_y_continuous(expand = c(.01,.01)) +
  theme(text = element_text(size = 4),
        axis.text.x = element_text(angle = 60,
                                   hjust = 1,
                                   vjust = 1),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank())+
  scale_fill_brewer(palette = "Set3")
dev.off()
# foldchange
fc <- read_csv("2.phylumFC.csv") %>%
  select(baseMean,phylum,log2FoldChange,padj) %>%
  mutate(status = log2FoldChange) %>%
  arrange(-log2FoldChange)
fc$phylum <- factor(fc$phylum,levels = fc$phylum)
fc$status[fc$log2FoldChange>0] <- "up"
fc$status[fc$log2FoldChange<0] <- "down"
fc$status[fc$padj>0.05] <- "insignificant"
mycol <- c("up"="red","down"="blue","insignificant"="grey")
pdf("2.phylum.log2fc.pdf", wi=4, he=3)
ggplot(fc, aes(x=log2FoldChange, y=phylum)) +
  geom_vline(xintercept = 0, size=.5)+
  geom_segment(aes(x=0, xend=log2FoldChange, 
                   y=phylum, yend=phylum),
                   size= .5, color = "#ababab", linetype = "dashed")+
  geom_point(aes(color = status,size = log10(baseMean))) +
  theme_bw() +
  scale_color_manual(values = mycol) +
  scale_x_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.background = element_blank(),
        text = element_text(size = 4))+
  scale_size_continuous(limits = c(2,6), range = c(.1,4))
dev.off()
