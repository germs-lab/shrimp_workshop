##################################################
##						##
## make a phyloseq object from tables           ##
##						##
##################################################
## Author: Fan Yang

## Before start, navigate to where your "otu count table", "taxa table", and "sample meta data" are and set it as your working directory. 
## For example, in my case (and my operating system is Mac OSX):
setwd("/Users/fanyang/Box Sync/Mexico/paola/data_for_R/rk1") 
## this means from now on, all of the files I import are coming directly out of the `rk1` folder

#####################
# load `phyloseq`   #
#####################
library(phyloseq)

####################
# read in tables   #
####################
# read in the otu count table from the working directory (ie. folder `rk1` in this case)
otu <- read.delim("cdhit_otu_table_wide.txt", row.names=1) #`read.delim` takes tab delimited files and makes the first row as the table header
# check otu table dimension
dim(otu) 
    # output shows (number of rows, number of columns)
    #> dim(otu) 
    #[1] 9875   36

# read in the taxa table from the working directory (ie. folder `rk1` in this case)
tax <- read.delim("cdhit_taxa_table_w_repseq.txt") #`read.delim` takes tab delimited files and makes the first row as the table header. 
# we also need to assign the `tax` table with row names that are the same as the `otu` table. 
# let's take a quick look at the first 6 rows (head) and 5 columns of the `otu` table to see the row names.
head(otu[, 1:5])
    # and we see row names are displayed as "OTU_0", "OTU_1", "OTU_10"... etc
    #> head(otu[, 1:5])
    #         lp_0 lp_1 lp_10 lp_11 lp_12
    #OTU_0       0    0     0     0     0
    #OTU_1     159   11    10    32    69
    #OTU_10      0    0     0     1     0
    #OTU_100     5    2     0     2     4
    #OTU_1000    0    0     0     0     2
    #OTU_1001    0    0     0     0     2
# let's a take a quick look at the first 6 rows (head) of the `tax1 table.
head(tax)


# there is a column in your `tax` table with header "OTUS". Let's use this column for row names. 
row.names(tax)<-tax$OTUS 
# check taxa table dimension
dim(tax) # output shows (number of rows, number of columns)

## the number of rows of `otu` should equal to `tax`. 
## if not, make sure you used the right files

# read in the sample metadata (ie. experimental design, treatments, sample id, etc.)
si <- read.delim("rk1_meta_sample_data.txt")
# check the dimension of the sample metadata
dim(si) 

## the number of rows of `si` should equal to the number of columns in `otu`. 
## if not, make sure you used the right files

# your sample metadata should contain a column that shares the same sample ids as displayed in the `otu` header
# we can take a quick look at the first 6 rows and 5 columns of the `otu` table
head(otu[, 1:5]) 
    # and we see:
    #> head(otu[, 1:5])
    #         lp_0 lp_1 lp_10 lp_11 lp_12
    #OTU_0       0    0     0     0     0
    #OTU_1     159   11    10    32    69
    #OTU_10      0    0     0     1     0
    #OTU_100     5    2     0     2     4
    #OTU_1000    0    0     0     0     2
    #OTU_1001    0    0     0     0     2
# than the first 6 rows of the `si` table
head(si)
    # and we see:
    #> head(si)
    #                       CIBNOR_id Sample_id Culture.media  Variable Replicate
    #1   Sediment_control_72h_E1_rep1     S_293  Sediment 72h   Control         1
    #2   Sediment_control_72h_E1_rep2     S_294  Sediment 72h   Control         2
    #3   Sediment_control_72h_E1_rep3     S_295  Sediment 72h   Control         3
    #4   Sediment_control_72h_E1_rep4     S_296  Sediment 72h   Control         4
    #5   Sediment_control_72h_E1_rep5     S_297  Sediment 72h   Control         5
    #6 Sediment_nutrients_72h_E1_rep1     S_298  Sediment 72h Nutrients         1
    #            plate well DNA.concentration X..of.reads X..of.singletons
    #1 Magallon_Plate2   D9              12.2       32679             1125
    #2 Magallon_Plate2  D10              24.4        7990              781
    #3 Magallon_Plate2  D11              27.8       28417             1434
    #4 Magallon_Plate2  D12              17.8       28447             1044
    #5 Magallon_Plate2   E1              19.8       20787             1188
    #6 Magallon_Plate2   E2              16.7       19106              703
    #  Experiment SAMPLES
    #1        RK1    lp_0
    #2        RK1    lp_1
    #3        RK1    lp_2
    #4        RK1    lp_3
    #5        RK1    lp_4
    #6        RK1    lp_5


data.phy<-phyloseq(otu_table(as.matrix(otu),taxa_are_rows=T), tax_table(as.matrix(tax)))
# check phyloseq object size
print("the phyloseq object contains: ")
data.phy
print("saving raw sequence phyloseq object to current directory ... ")
saveRDS(data.phy, "raw_sequence_phyloseq.RDS")

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
