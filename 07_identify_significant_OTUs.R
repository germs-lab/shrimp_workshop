##################################################
##						##
## Find the OTUs that were significantly        ##
## different in abundances using phyloseq-DESeq2##
##						##
##################################################
## Author: Fan Yang

## Before start, navigate to where your phyloseq object was saved at. 
## It should be the same as where your "otu count table", "taxa table", and "sample meta data" are at.
## set it as your working directory. 
## For example, in my case (and my operating system is Mac OSX):
setwd("~/Box sync/Mexico/paola/data_for_R/rk1") 
## this means from now on, all of the files I import are coming directly out of the `rk1` folder

#####################
# load libraries    #
#####################
library(phyloseq)
library(DESeq2)

#####################
# NOTE              #
#####################
# `DESeq2` has data normalizaion implimented.
# Therefore, always use non-standardized data but with potentially errorness OTUs and under sampled samples removed.
# ie. "data.mintax5.excluded.phy"
##
# `DESeq2` was originally designed to identify differentially expressed genes (e.g., microarray data)
# it can only compare 2 things at a time (ie, control vs. treated, or high vs. low; it cannot be, control vs. high vs. low)
# it's under line assumption is "majority of the genes (or OTUs) do not change between treatments"

######################################################################
# load the phyloseq object with potentially errorness OTUs removed   #
######################################################################
# if you recall, in step 2, we save a phyloseq object with OTUs summed across all samples less than 5 removed and low coverage and/or low sequencing depth samples excluded to the working directory using command: saveRDS(data.mintax5.excluded.phy, "otu_sum_min_5_excluded_4samples_phyloseq.RDS")
# now we are going to start where we left at by read in the same phyloseq object

data.mintax5.excluded.phy <- readRDS("otu_sum_min_5_excluded_4samples_phyloseq.RDS")

# do a quick check of the phyloseq object:
data.mintax5.excluded.phy
    #> data.mintax5.excluded.phy        
    #phyloseq-class experiment-level object
    #otu_table()   OTU Table:         [ 4942 taxa and 30 samples ]
    #sample_data() Sample Data:       [ 30 samples by 12 sample variables ]
    #tax_table()   Taxonomy Table:    [ 4942 taxa by 8 taxonomic ranks ]

######################################################################
# identify the experimental levels for comparison 		     #
######################################################################
# as I mentioned in the "NOTE", `DESeq2` can only compare 2 levels. 
# if you recall, from our post-hoc community analysis (step 4), we have identified all significantly different pairs of communities. 
# this can be used to guide this part of the analysis. 
# for this demonstration, I will use the pair " T0::Nutrients"

# 1. subset phyloseq object for "T0" and "Nutrients" only. 
t0.nutrients <- subset_samples(data.mintax5.excluded.phy, grepl("T0|Nutrients", Variable))
# make sure no all 0 OTUs are in the phyloseq object
t0.nutrients <- prune_taxa(taxa_sums(t0.nutrients) > 0, t0.nutrients)
t0.nutrients
	#> t0.nutrients                     
	#phyloseq-class experiment-level object
	#otu_table()   OTU Table:         [ 4029 taxa and 14 samples ]
	#sample_data() Sample Data:       [ 14 samples by 12 sample variables ]
	#tax_table()   Taxonomy Table:    [ 4029 taxa by 8 taxonomic ranks ]

# 2. convert the phylsoeq object to DESeq2 object
t0.nutrients.dds <- phyloseq_to_deseq2(t0.nutrients, ~ Variable)
	# this is what the DESeq2 object looks like:
	#> t0.nutrients.dds
	#class: DESeqDataSet 
	#dim: 4029 14 
	#metadata(0):
	#assays(1): counts
	#rownames(4029): OTU_1 OTU_100 ... OTU_996 OTU_999
	#rowRanges metadata column names(0):
	#colnames(14): lp_15 lp_16 ... lp_8 lp_9
	#colData names(12): CIBNOR_id Sample_id ... Experiment SAMPLES
# we will also save the this DESeq2 object to file for later use:
saveRDS(t0.nutrients.dds, "t0.nutirents.dds.RDS")
# you can load it like this:
t0.nutrients.dds <- readRDS("t0.nutirents.dds.RDS")

