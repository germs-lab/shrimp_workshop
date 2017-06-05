##################################################
##						##
## evaluate how microbial communities differ    ##
## based on different experimental factors      ##
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
library(vegan)

######################################################################
# load the phyloseq object with potentially errorness OTUs removed   #
######################################################################
# if you recall, in step 2, we save a phyloseq object with OTUs summed across all samples less than 5 removed and low coverage and/or low sequencing depth samples excluded to the working directory using command: saveRDS(data.mintax5.excluded.phy, "otu_sum_min_5_excluded_4samples_phyloseq.RDS")
# now we are going to start where we left at by read in the same phyloseq object

data.mintax5.excluded.phy <- readRDS("otu_sum_min_5_excluded_phyloseq.RDS")

# do a quick check of the phyloseq object:
data.mintax5.excluded.phy

######################################################################
# standardizing counts for different sequencing depth                #
######################################################################
# as in step 2, we found that sample sequencing depths varied. 
# therefore, we can use all samples. 
# however, different sampling depths (total number of sequences) do make OTU counts among different samples difficult to compare. 
# therefore, we have to standardize the OTU counts based on the sequencing depth somehow. 

# there are two ways of doing this. Both of them sensitive to the true sequencing coverage.# so excluding samples that were not well covered (or drastically different from the rest of your samples) in step 2 are really important 

######################################################################
# we introduced relative abundance during the workshop
######################################################################
data.mintax5.excluded.rela.phy <- transform_sample_counts(data.mintax5.excluded.phy, function(x) x/sum(x))
    # this converts the total number of sequences in each sample to 1 and all OTU counts as fractions of 1 in each sample. 

# save it for later use:
saveRDS(data.mintax5.excluded.rela.phy, "otu_sum_min_5_excluded_4samples_relative_abundance_phyloseq.RDS")

######################################################################
# another one is rarefaction: which resamples your OTU for the same sequencing depth
# ie., if sample 1 has 10k sequences, sample 2 has 15k sequences,
# sample 2 can be resampled in computer to 10k sequences and gives new abundances for each OTU.
######################################################################
data.mintax5.excluded.rarefy.phy <- rarefy_even_depth(data.mintax5.excluded.phy,, sample.size=min(sample_sums(data.mintax5.excluded.phy)-1), rngseed=15879966)
    # the above command means that we are making a new phyloseq object 
    # resampling all OTUs in each sample for a sampling depth of 13367 sequences
    # and we see an output message of:
    #
    #> data.mintax5.excluded.rarefy.phy <- rarefy_even_depth(data.mintax5.excluded.phy,, sample.size=min(sample_sums(data.mintax5.excluded.phy)-1), rngseed=15879966)
    #`set.seed(15879966)` was used to initialize repeatable random subsampling.
    #Please record this for your records so others can reproduce.
    #Try `set.seed(15879966); .Random.seed` for the full vector
    #...
    #109OTUs were removed because they are no longer 
    #present in any sample after random subsampling
    #
    #...
    #
# the message shows that because of resampling, 1090 OTUs are no longer present in your data. 
# this is expected, because when you resample from 30k sequences to 10k sequences, some OTUs that were low in abundance to begin with are likely no longer detectable. 
# so, if that 1090 OTUs are really important for your study (eg., you know you need to focus on the rare species), I would use relative abundance instead. 

## we can save it for later use:
saveRDS(data.mintax5.excluded.rarefy.phy, "otu_sum_min_5_excluded_4samples_rarefied_phyloseq.RDS")

######################################################################
# NOTES ON STANDARDIZATION #
######################################################################
# if you have removed low coverage and/or low seequencing depth samples (step 2),
# both relative abundance and rarefaction should generate very similar results. 

