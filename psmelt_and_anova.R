all_data <- readRDS("/Volumes/New Volume/Mexico/maurilia/all_data_phyloseq.RDS")
library(phyloseq)
library(plyr)

all_data
head(sample_data(all_data))
unique(sample_data(all_data)$Treatment)

sample_sums(all_data)
all_data_rela <- transform_sample_counts(all_data, function(x) x / sum(x))
all_data_rela
sample_sums(all_data_rela)

# melt the object to table
all_data_rela.melt <- psmelt(all_data_rela)
head(all_data_rela.melt)

library(plyr)

## group by individual treatments
proteo_rela.melt <- subset(all_data_rela.melt, phylum == "Proteobacteria")
proteo_rela.trt <- ddply(proteo_rela.melt, .(Treatment), summarise, total=sum(Abundance))
proteo_rela.trt

proteo_rela.trt.macro <- ddply(proteo_rela.melt, .(Treatment, Macroalga), summarise, total=sum(Abundance))
proteo_rela.trt.macro

head(proteo_rela.melt)
anova(lm(proteo_rela.melt$Abundance ~ proteo_rela.melt$Treatment))
TukeyHSD(aov(proteo_rela.melt$Abundance ~ proteo_rela.melt$Treatment))

