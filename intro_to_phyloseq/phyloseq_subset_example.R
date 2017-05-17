library(phyloseq)
data <- readRDS("physeq_example.RDS")
data

data_1h <- subset_samples(data, time == "1h")
data_1h

otu_table(data_1h)
tax_table(data_1h)
sample_data(data_1h)

data_bacter <- subset_taxa(data, phylum == "Bacteroidetes")
data_bacter

otu_table(data_bacter)
tax_table(data_bacter)
sample_data(data_bacter)

data_avg50 <- filter_taxa(data, function(x) mean(x) > 50, prune = T)
data_avg50

otu_table(data_avg50)
tax_table(data_avg50)
# filter out taxa that sums to less than 500
data_core <- filter_taxa(data_avg50, function(x) sum(x > 5) > (.8*length(x)), TRUE)
data_core
otu_table(data_core)

data_avg50
otu_table(data_avg50)
dist = phyloseq::distance(data_avg50, "bray")
dist
ord_results = ordinate(data_avg50,"PCoA", "bray")
plot_ordination(data_avg50, ord_results)
#   X month time
plot_ordination(data_avg50, ord_results,color="time",label="X",shape="month")
library(vegan)
adonis(dist~time,perm=9999,as(sample_data(data_avg50),"data.frame"))
