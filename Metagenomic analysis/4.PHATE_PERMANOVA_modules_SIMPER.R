setwd("D:/Study")
library(vegan)
library(tidyverse)
library(plyr)

#make presence and absence file
pa <- read_tsv("enrichM/module_completeness.tsv") %>%
  select(genome, module, completeness) %>%
  pivot_wider(names_from = genome,
              values_from = completeness)
# write_tsv(pa, "4.genome.module.completeness.tsv")
# use process AO or A2O as group
group <- read.delim("0.genome.group.txt", sep = '\t', header = T,
                    stringsAsFactors = FALSE)
# use tax as group information
group <- read.delim("0.GTDB.tax.tsv", sep = '\t', header = T,
                    stringsAsFactors = FALSE)
group$phylum[group$phylum=="Proteobacteria"] <- group$class[group$phylum=="Proteobacteria"]
group$phylum[group$phylum=="Bdellovibrionota_B"|
               group$phylum=="Bdellovibrionota_A"] <- "Bdellovibrionota"
group <- select(group,c(genome,phylum))
np <- data.frame(table(group$phylum)) %>% arrange(-Freq) %>%
  filter(Freq <= 5) %>% select(Var1)
group$phylum[is.na(group$phylum)] <- "Others"
group$phylum[is.element(group$phylum,np$Var1)] <- "Others"
# PHATE
pa <- read_tsv("4.phate.phylum.tsv") %>% left_join(group)
pdf("4.phate.pdf", wi = 4.2, he = 3)
ggplot(pa, aes(x=PHATE1,y=PHATE2,color=phylum)) +
  geom_point()+
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.background = element_blank(),
        text = element_text(size = 4))+
  scale_color_brewer(palette = "Set3")
dev.off()

# PERMANOVA
pa <- read.delim('4.genome.module.completeness.tsv', 
                       row.names = 1, 
                       sep = '\t', 
                       stringsAsFactors = FALSE, 
                       check.names = FALSE)
pa[pa>0]<-1
pa[is.na(pa)]<-0
pa <- data.frame(t(pa))
pa <- pa[group$genome,]
# AOvsA2O
adonis_result <- adonis(pa~group, 
                        group, 
                        permutations = 1000, 
                        distance = 'bray')
# phylum-wise
target = c("Myxococcota", "Chloroflexota", "Acidobacteriota",
           "Actinobacteriota", "Alphaproteobacteria",
           "Bacteroidota", "Bdellovibrionota", "Gammaproteobacteria",
           "Patescibacteria", "Planctomycetota", "Elusimicrobiota",
           "Zixibacteria", "Campylobacterota", "Nitrospirota")
# or
target = group %>% group_by(phylum) %>% summarise(n=n()) %>%
  filter(n>=2) %>% select(phylum)
target = target$phylum
i = 0
for(phy in target){
  nowgroup = group
  nowgroup$phylum[nowgroup$phylum!=phy] = "Others"
  adonis_result = adonis(pa~phylum, 
                        nowgroup, 
                        permutations = 9999, 
                        distance = 'bray')
  if(i==0){
    output <- data.frame(unlist(data.frame(adonis_result$aov.tab, check.names = FALSE)[1, ]))
    i = i+1
    names(output)[1] = phy
  }
  else{
    output = cbind(output,data.frame(unlist(data.frame(adonis_result$aov.tab, check.names = FALSE)[1, ])))
    names(output)[i+1] = phy
    i = i+1 
  }
}
output <- data.frame(t(output))
output <- cbind(rownames(otuput), otuput)
names(output) <- c('', 'Df', 'Sums of squares', 
                   'Mean squares', 'F.Model', 
                   'Variation (R2)', 'Pr (>F)')
write.table(output, file = paste('4.PERMANOVA.phylum.tsv',sep=""), 
            row.names = FALSE, sep = '\t', 
            quote = FALSE, na = '')
# write results
#otuput <- data.frame(adonis_result$aov.tab, 
#                     check.names = FALSE, 
#                     stringsAsFactors = FALSE)

#write.table(output, file = paste('4.', now_tax,'.PERMANOVA.module.tsv',sep=""), 
#            row.names = FALSE, sep = '\t', 
#            quote = FALSE, na = '')

########################
# pair-wise comparison #
########################
group_name <- unique(group$phylum)
adonis_result_pairwise <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(group, phylum %in% c(group_name[i], group_name[j]))
    otu_ij <- pa[group_ij$genome, ]
    adonis_result_otu_ij <- adonis(otu_ij~phylum, group_ij, permutations = 999, distance = 'bray')
    adonis_result_pairwise <- rbind(adonis_result_pairwise, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', unlist(data.frame(adonis_result_otu_ij$aov.tab, check.names = FALSE)[1, ])))
  }
}
adonis_result_pairwise <- data.frame(adonis_result_pairwise, stringsAsFactors = FALSE)
names(adonis_result_pairwise) <- c('group', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
for (i in 1:nrow(adonis_result_pairwise)) {
  if (adonis_result_pairwise[i, 'Pr (>F)'] <= 0.001) adonis_result_pairwise[i, 'Sig'] <- '***'
  else if (adonis_result_pairwise[i, 'Pr (>F)'] <= 0.01) adonis_result_pairwise[i, 'Sig'] <- '**'
  else if (adonis_result_pairwise[i, 'Pr (>F)'] <= 0.05) adonis_result_pairwise[i, 'Sig'] <- '*'
}
write.table(adonis_result_pairwise, '4.PERMANOVA.pairwise.tsv', row.names = FALSE, sep = '\t', quote = FALSE, na = '')

# SIMPER
module_anno = read_tsv("0.Module.BRITE.txt")
mycol = c("Carbohydrate metabolism" = "#ffaeb9",
          "Lipid metabolism" = "#ee6a50",
          "Nucleotide metabolism" = "#ee3b3b",
          "Amino acid metabolism" = "#7b68ee",
          "Glycan metabolism" = "#ffec8b",
          "Metabolism of cofactors and vitamins" = "#ee9a00",
          "Biosynthesis of terpenoids and polyketides" = "#3cb371",
          "Xenobiotics biodegradation" = "#20b2aa",
          "Energy metabolism" = "#4f94cd",
          "Gene set" = "#b5b5b5",
          "Module set" = "#969696")
# now let's rock
for(now_tax in target){
sir = pa
sir.env = group
sir.env$phylum[sir.env$phylum!=now_tax] = "Others"
sir.env$phylum = factor(sir.env$phylum, levels = c(now_tax, "Others"))
sim <- with(sir.env, simper(sir, phylum, permutations = 999, parallel = 4))
# summary(sim)
sim.result = tibble(ID=c(1:length(sim[[1]]$species)))
sim.result$module = sim[[1]]$species
sim.result$contribution = sim[[1]]$average
sim.result$sd = sim[[1]]$sd
sim.result$Others = round(sim[[1]]$ava,digits = 2)
sim.result$target = round(sim[[1]]$avb,digits = 2) # change
sim.result$order = sim[[1]]$ord
sim.result$cusum = sim[[1]]$cusm
sim.result$p = sim[[1]]$p
print(paste(now_tax, ": ", sim[[1]]$overall,sep=""))
write_tsv(sim.result, paste("SIMPER/4.",now_tax,".SIMPER.tsv",sep=""))
# annotate
sim.result = left_join(sim.result, module_anno)
sim.result = mutate(sim.result, log2fc= log2(target/Others))
sim.result$name = str_replace(sim.result$name,
                              pattern = "PATH.*",
                              replacement = paste(round(sim.result$target*100,0),"%",sep=""))
# sort
sim.result$L2 = factor(sim.result$L2, levels = c("Carbohydrate metabolism",
                                                 "Lipid metabolism",
                                                 "Nucleotide metabolism",
                                                 "Amino acid metabolism",
                                                 "Glycan metabolism",
                                                 "Metabolism of cofactors and vitamins",
                                                 "Biosynthesis of terpenoids and polyketides",
                                                 "Xenobiotics biodegradation",
                                                 "Energy metabolism",
                                                 "Gene set",
                                                 "Module set"))
sim.result = sim.result %>% arrange(L1,L2,L3,log2fc)
sim.result$name = factor(sim.result$name, levels = sim.result$name)
#plot
#pdf(paste("4.",now_tax,".SIMPER.pdf",sep=""), wi=4, he =3)
ggplot(sim.result %>% filter(p<0.05&is.finite(log2fc)&!is.na(name)),
       aes(x=as.factor(name),y=log2fc,fill=L2)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(text = element_text(size=4),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(size=.5,linetype = "dashed")) +
  coord_flip() +
  xlab(paste(now_tax)) +
  ylab("Log2(Foldchange)") +
  scale_fill_manual(values = mycol)
ggsave(paste("SIMPER/4.",now_tax,".SIMPER.pdf",sep=""),
       device = "pdf", width = 4, height = 3)
#dev.off()
}
