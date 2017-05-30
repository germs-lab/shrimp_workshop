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

## plot Good's coverage VS. sequencing depth
plot(C ~ sample_sums, data = ssum, ylim= c(0, 1))

# and we can see that although some of the samples have <10k sequences, they were still well covered. 

