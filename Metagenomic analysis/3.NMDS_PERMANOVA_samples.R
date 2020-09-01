setwd("D:/Study")

library(vegan)
library(tidyverse)
library(plyr)
# NMDS
# group information
group <- read.delim("0.sampleGroup.tsv", sep = '\t', header = T,
                    stringsAsFactors = FALSE)
# genome abundance profile
pa <- read.delim('0.Coverm_Relative_Abundance.tsv', 
                 row.names = 1, 
                 sep = '\t', 
                 stringsAsFactors = FALSE, 
                 check.names = FALSE)[-1,]
pa <- data.frame(t(pa))
pa <- pa[group$sample,]
nmds1 <- metaMDS(pa, distance = "bray", k = 4, try = 1000)
nmds1
nmds1.stress <- nmds1$stress
nmds1.point <- data.frame(nmds1$point)
nmds1.species <- data.frame(nmds1$species)
# write out
#write.table(nmds1.point, 'nmds.module.tsv',sep = "\t",
#            col.names = T, row.names = T, quote = F)
nmds_plot <- nmds1
nmds_plot$species <- {nmds_plot$species}[1:10, ]
plot(nmds_plot, type = 't', 
     main = paste('Stress =', round(nmds1$stress, 4)))
sample_site <- nmds1.point[1:2]
sample_site$sample <- rownames(sample_site)
names(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
sample_site <- as_tibble(sample_site) %>%
  left_join(group)
tax <- read_tsv("../MENAP/info_MAG_0.81_combined_tax_2.txt")[,c(1,37)]
names(tax)[1] <- "genome"
sample_site <- left_join(sample_site, tax)
# ggplot
group_border <- ddply(sample_site, 'group', 
                      function(df) df[chull(df[[1]], df[[2]]), ])
pdf("3.NMDS.pdf", wi = 4.2, he = 3)
ggplot(sample_site, aes(x = NMDS1, y = NMDS2, 
                        color = group)) +
  geom_point() +
  geom_polygon(data = group_border, 
               alpha = 0.0, show.legend = F) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2")
dev.off()

# PERMANOVA
adonis_result <- adonis(pa~group, 
                        group, 
                        permutations = 10000, 
                        distance = 'bray')
adonis_result
# write results
otuput <- data.frame(adonis_result$aov.tab, 
                     check.names = FALSE, 
                     stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 
                   'Mean squares', 'F.Model', 
                   'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = '3.PERMANOVA.result.tsv', 
            row.names = FALSE, sep = '\t', 
            quote = FALSE, na = '')

#PCoA
distance <- vegdist(pa, method = 'bray')
pcoa <- cmdscale(distance, k = (nrow(pa) - 1), eig = TRUE)
ordiplot(scores(pcoa)[ ,c(1, 2)], type = "text")
point <- data.frame(pcoa$point)
species_coordinate <- wascores(pcoa$points[,1:2], pa)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_site <- data.frame({pcoa$point})[1:2]
sample_site$sample <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
sample_site <- merge(sample_site, group, by = 'sample', all.x = TRUE)
group_border <- ddply(sample_site, 'group', function(df) df[chull(df[[2]], df[[3]]), ])
#plot
pdf("3.PCoA.pdf", wi = 4.2, he = 3)
ggplot(sample_site, aes(PCoA1, PCoA2)) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        text = element_text(color = "black"), 
        panel.grid.minor = element_line(color = '#EBEBEB', size = 0.1)) + 
  geom_vline(xintercept = 0, color = '#696969', size = 0.4, linetype = "dashed") + 
  geom_hline(yintercept = 0, color = '#696969', size = 0.4, linetype = "dashed") + 
  #geom_polygon(data = group_border, aes(fill = metal)) + 
  geom_point(aes(color = group), size = 2, alpha = 0.8) + 
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%'))
dev.off()