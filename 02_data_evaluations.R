##################################################
##						##
## evaluate sequencing depths and coverages     ##
## using existing phyloseq object               ##
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
# load `phyloseq`   #
#####################
library(phyloseq)

##################################################
# load the phyloseq object you saved in step 1   #
##################################################
# if you recall, at the end of step 1, we save a phyloseq object to the working directory using command: saveRDS(data.phy, "raw_data_phyloseq.RDS")
# now we are going to start where we left at by read in the same phyloseq object

data.phy <- readRDS("raw_data_phyloseq.RDS")

# do a quick check of the phyloseq object:
data.phy
    #and it looks like this:
    #> data.phy
    #phyloseq-class experiment-level object
    #otu_table()   OTU Table:         [ 9875 taxa and 34 samples ]
    #sample_data() Sample Data:       [ 34 samples by 12 sample variables ]
    #tax_table()   Taxonomy Table:    [ 9875 taxa by 8 taxonomic ranks ]

##################################################
# removing potentially errorness OTUs            #
##################################################
# although we have already removed many low quality and errorness sequences (e.g., chimeras) during sequence processing. There could still be some left whe we made our OTU count table. 
# these potentially errorness sequences would frequently be clustered into their own OTUs, resulting OTUs appeared to be rare (e.g., OTU_10 has a count of 1 in lp_11 only)
# these OTUs could be really rare organisms. But for organisms this rare, they are not comparable among replicates of samples either. 
# Therefore, I usually remove OTU's sum up to less than 5 across all samples before any downstream analyses

data.mintax5.phy <- filter_taxa(data.phy, function(x) sum(x) >= 5, T)

# if you check the minimal OTU sums across all samples, it should be 5. 
min(taxa_sums(data.mintax5.phy))

# we can also save this phyloseq object for later use as well:
saveRDS(data.mintax5.phy, "otu_sum_min_5_phyloseq.RDS")

##################################################
# evaluating sequencing depth and coverage       #
##################################################
# the number of sequences per sample and how many OTUs were covered by the seuqencing can skew your community analyses. 
# I usually do a quick evaluation of sequencing depth and sequencing coverage to remove any samples that were not adequately sequenced. 
# the below histogram plots the distribution of sample sequencing depths. 
hist(sample_sums(data.mintax5.phy), breaks = nsamples(data.mintax5.phy)/2)
    # and we can see that majority the samples have ~20k sequences, but 2 samples have less than 10K sequences.
    # the samples that were less than 10k may not be adequately sequenced. 

# therefore, we need to double check the sequence coverage as well. 
# We use Good's Estimate of Coverage to evaluate how well each sample's community was covered by sequencing. 
# calculate Good's coverage for each sample
ssum <- data.frame(sample_sums(data.mintax5.phy)) #get sample sums
names(ssum) <- "sample_sums"
totu<-t(data.frame(otu_table(data.mintax5.phy))) #transpose subsetted otu table so that samples are in row

ssum$n1 <- rowSums(totu == 1) #counts of singletons in each sample
ssum$C <- 1-(ssum$n1 / ssum$sample_sums) #calculate Good's coverage

## plot Good's coverage histogram
hist(ssum$C, breaks = nrow(ssum)/2)

# keeping in mind that C is determined based on how many OTUs with only 1 count (singletons) are there in your sample. 
# therefore, when sequencing depth is shallow, it can give you a false good coverage

## therefore, we should also plot Good's coverage VS. sequencing depth
plot(C ~ sample_sums, data = ssum, ylim= c(0, 1))

# as we can see that the 4 samples on the left have lower coverage and/or lower sequencing depth comparing to the other smaples.
# personally, I would exclude these 4 samples in the downstream analyses

##################################################################
# excluding samples with low coverage and/or sequencing depth    #
##################################################################
# first, identify the 4 samples we would like to exclude
# sort by sample_sums then by coverage
ssum <- ssum[order(ssum$sample_sums, ssum$C), ]
    # the first 4 are the samples we would like to exclude
sample_to_exclude <- row.names(ssum)[1:4] 
    #and we can see the samples in "sample_to_exclude" are:
    #> sample_to_exclude
    #[1] "lp_24" "lp_1"  "lp_10" "lp_19"

# create a new phyloseq object that with the above 4 samples removed:
data.mintax5.excluded.phy <- prune_samples(!sample_names(data.mintax5.phy) %in% sample_to_exclude, data.mintax5.phy)

data.mintax5.excluded.phy
    #we can see that the new phyloseq contains:
    #> data.mintax5.excluded.phy
    #phyloseq-class experiment-level object
    #otu_table()   OTU Table:         [ 4944 taxa and 30 samples ]
    #sample_data() Sample Data:       [ 30 samples by 12 sample variables ]
    #tax_table()   Taxonomy Table:    [ 4944 taxa by 8 taxonomic ranks ]

# sometimes, some OTUs may be only found in some of the samples. And when you remove the samples, it would result those OTUs to sum up to 0 across all samples left. 
# phyloseq does not remove these all 0 OTUs automatically. But they could cause issues in later analyses. 
# therefore, we should check it and remove OTUs that are all 0s after excluding some of the samples. 
min(taxa_sums(data.mintax5.excluded.phy))
    #> min(taxa_sums(data.mintax5.excluded.phy))
    #[1] 0
    #
    # as we can see, there are OTUs summed up to 0 across all samples

# to remove all 0 OTUs:
data.mintax5.excluded.phy <- prune_taxa(taxa_sums(data.mintax5.excluded.phy) > 0, data.mintax5.excluded.phy)

data.mintax5.excluded.phy
    #and now the phyloseq object looks like this:
    #> data.mintax5.excluded.phy
    #phyloseq-class experiment-level object
    #otu_table()   OTU Table:         [ 4942 taxa and 30 samples ]
    #sample_data() Sample Data:       [ 30 samples by 12 sample variables ]
    #tax_table()   Taxonomy Table:    [ 4942 taxa by 8 taxonomic ranks ]

# and now we can save this new phyloseq object to file for later use:
saveRDS(data.mintax5.excluded.phy, "otu_sum_min_5_excluded_4samples_phyloseq.RDS")

