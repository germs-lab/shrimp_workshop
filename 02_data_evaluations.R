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

test <- filter_taxa(data.phy, function(x) sum(x) == 1, T)
head(otu_table(test))
## get rid of any taxa that summing up to be less than 5 across all samples
data.taxmin5.phy<-prune_taxa(taxa_sums(data.phy) >= 5, data.phy)
data.taxmin5.phy<-prune_samples(sample_sums(data.taxmin5.phy) > 0, data.taxmin5.phy) #get rid of samples with all 0's, if it exists
print("saving phyloseq object with taxa sum >= 5 across all samples to current directory ... ")
saveRDS(data.taxmin5.phy, "taxsum_min5_sequence_phyloseq.RDS")

## save the histogram of sample sequencing depth distribution
print("Generating histogram on sample sequencing depth to current directory ... ")
pdf("data.taxmin5.sequencing_depth_hist.pdf")
hist(sample_sums(data.taxmin5.phy), breaks = nsamples(data.taxmin5.phy)/2)
dev.off()

## calculate Good's coverage for each sample
si <- data.frame(sample_sums(data.taxmin5.phy)) #get sample sums
names(si) <- "sample_sums"
totu<-t(data.frame(otu_table(data.phy))) #transpose subsetted otu table so that samples are in row

si$n1 <- rowSums(totu == 1) #counts of singletons in each sample
si$C <- 1-(si$n1 / si$sample_sums) #calculate Good's coverage
si$SAMPLES <- row.names(si)
## save the sample coverage information to file
print("Saving the sample coverage information as a text file in current directory ... ")
write.table(si, "data.taxmin5.sample_goods_coverage.txt", sep = "\t", quote = F, row.names = F)

## plot Good's coverage histogram
print("Generating histogram on Good's estimated coverage to current directory ... ")
pdf("data.taxmin5.goods_coverage_hist.pdf")
hist(si$C, breaks = nrow(si)/2)
dev.off()

