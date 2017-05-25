##################################################
##						##
## make a phyloseq object from tables           ##
##						##
##################################################
## Author: Fan Yang

## Before start, navigate to where your "otu count table", "taxa table", and "sample meta data" are and set it as your working directory. 
## For example, in my case (and my operating system is Mac OSX):
setwd("~/Box sync/Mexico/paola/data_for_R/rk1") 
## this means from now on, all of the files I import are coming directly out of the `rk1` folder

#####################
# load `phyloseq`   #
#####################
library(phyloseq)

####################
# read in tables   #
####################

## getting the otu table ready
##############################
# read in the otu count table from the working directory (ie. folder `rk1` in this case)
otu <- read.delim("cdhit_otu_table_wide.txt", row.names=1) #`read.delim` takes tab delimited files and makes the first row as the table header
# check otu table dimension
dim(otu) 
    # output shows (number of rows, number of columns)
    #> dim(otu) 
    #[1] 9875   36

## getting the taxonomy table ready
###################################
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
    #and we see column "OTUS" share the same information as the `otu` table row names.
    #> head(tax)
    #      repseq     OTUS   domain                phylum
    #1    lp_0|10 OTU_3178 Bacteria        Proteobacteria
    #2 lp_0|10002 OTU_2682 Bacteria        Proteobacteria
    #3 lp_0|10003 OTU_9535 Bacteria unclassified_Bacteria
    #4 lp_0|10004 OTU_2681 Bacteria unclassified_Bacteria
    #5 lp_0|10029  OTU_505  Archaea        Woesearchaeota
    #6 lp_0|10031 OTU_2233 Bacteria        Proteobacteria
    #                        class                       order
    #1         Gammaproteobacteria             Alteromonadales
    #2         Gammaproteobacteria             Alteromonadales
    #3       unclassified_Bacteria       unclassified_Bacteria
    #4       unclassified_Bacteria       unclassified_Bacteria
    #5 unclassified_Woesearchaeota unclassified_Woesearchaeota
    #6         Alphaproteobacteria            Rhodospirillales
    #                       family                          genus
    #1            Alteromonadaceae  unclassified_Alteromonadaceae
    #2            Alteromonadaceae                     Glaciecola
    #3       unclassified_Bacteria          unclassified_Bacteria
    #4       unclassified_Bacteria          unclassified_Bacteria
    #5 unclassified_Woesearchaeota    unclassified_Woesearchaeota
    #6           Rhodospirillaceae unclassified_Rhodospirillaceae

# let's use column "OTUS" for `tax` table row names. 
row.names(tax)<-tax$OTUS 
head(tax) # to double check the tax table
    #now the tax table should look like this:
    #> head(tax) # to double check the tax table
    #             repseq     OTUS   domain                phylum
    #OTU_3178    lp_0|10 OTU_3178 Bacteria        Proteobacteria
    #OTU_2682 lp_0|10002 OTU_2682 Bacteria        Proteobacteria
    #OTU_9535 lp_0|10003 OTU_9535 Bacteria unclassified_Bacteria
    #OTU_2681 lp_0|10004 OTU_2681 Bacteria unclassified_Bacteria
    #OTU_505  lp_0|10029  OTU_505  Archaea        Woesearchaeota
    #OTU_2233 lp_0|10031 OTU_2233 Bacteria        Proteobacteria
    #                               class                       order
    #OTU_3178         Gammaproteobacteria             Alteromonadales
    #OTU_2682         Gammaproteobacteria             Alteromonadales
    #OTU_9535       unclassified_Bacteria       unclassified_Bacteria
    #OTU_2681       unclassified_Bacteria       unclassified_Bacteria
    #OTU_505  unclassified_Woesearchaeota unclassified_Woesearchaeota
    #OTU_2233         Alphaproteobacteria            Rhodospirillales
    #                              family                          genus
    #OTU_3178            Alteromonadaceae  unclassified_Alteromonadaceae
    #OTU_2682            Alteromonadaceae                     Glaciecola
    #OTU_9535       unclassified_Bacteria          unclassified_Bacteria
    #OTU_2681       unclassified_Bacteria          unclassified_Bacteria
    #OTU_505  unclassified_Woesearchaeota    unclassified_Woesearchaeota
    #OTU_2233           Rhodospirillaceae unclassified_Rhodospirillaceae

# check taxa table dimension
dim(tax) # output shows (number of rows, number of columns)
    #output shows
    #> dim(tax) # output shows (number of rows, number of columns)
    #[1] 9875    8

## the number of rows of `otu` should equal to `tax`. 
## if not, make sure you used the right files

## getting the metadata ready
#############################
# read in the sample metadata (ie. experimental design, treatments, sample id, etc.)
si <- read.delim("rk1_meta_sample_data.txt")
# check the dimension of the sample metadata
dim(si) 
    #we see the dimension of the sample metadata is:
    #> dim(si) 
    #[1] 34 12

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
# then the first 6 rows of the `si` table
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

# we can see that column "SAMPLES" in `si` shares the same information as the headers in `otu`. 
# therefore, to match `si` with `otu`, we need to set the row names of `si` using column "SAMPLES". 
row.names(si) <- si$SAMPLES
# and we can see what `si` look like now:  
head(si)
    #> head(si)
    #                          CIBNOR_id Sample_id Culture.media  Variable Replicate
    #lp_0   Sediment_control_72h_E1_rep1     S_293  Sediment 72h   Control         1
    #lp_1   Sediment_control_72h_E1_rep2     S_294  Sediment 72h   Control         2
    #lp_2   Sediment_control_72h_E1_rep3     S_295  Sediment 72h   Control         3
    #lp_3   Sediment_control_72h_E1_rep4     S_296  Sediment 72h   Control         4
    #lp_4   Sediment_control_72h_E1_rep5     S_297  Sediment 72h   Control         5
    #lp_5 Sediment_nutrients_72h_E1_rep1     S_298  Sediment 72h Nutrients         1
    #               plate well DNA.concentration X..of.reads X..of.singletons
    #lp_0 Magallon_Plate2   D9              12.2       32679             1125
    #lp_1 Magallon_Plate2  D10              24.4        7990              781
    #lp_2 Magallon_Plate2  D11              27.8       28417             1434
    #lp_3 Magallon_Plate2  D12              17.8       28447             1044
    #lp_4 Magallon_Plate2   E1              19.8       20787             1188
    #lp_5 Magallon_Plate2   E2              16.7       19106              703
    #     Experiment SAMPLES
    #lp_0        RK1    lp_0
    #lp_1        RK1    lp_1
    #lp_2        RK1    lp_2
    #lp_3        RK1    lp_3
    #lp_4        RK1    lp_4
    #lp_5        RK1    lp_5
    #

###########i#################
## make a phyloseq object  ##
#############################
# now `otu`, `tax`, and `si` are ready to be put together to make a phyloseq object
data.phy<-phyloseq(otu_table(as.matrix(otu),taxa_are_rows=T), tax_table(as.matrix(tax)))
    # explanation to the above code:  
    # `otu_table` and `tax_table` are phyloseq functions to prepare `otu` and `tax` for phyloseq to understand. 
    # both `otu_table` and `tax_table` can only read "matrix". That's why we used `as.matrix` function to convert data frames `otu` and `tax`.
    # "taxa_are_rows=T" is a special option in function `otu_table`. It tells `otu_table` that the matrix ("as.matrix(otu)") have "OTU_0", "OTU_1", "OTU_2", etc. as row names in table `otu`. 

# now we can add metadata table `si` to the phyloseq object `data.phy` we just created 
sample_data(data.phy) <- si
    # `sample_data` is a phyloseq function 
    # the above command assign table `si` to phyloseq object `data.phy` as metadata using a phyloseq function `sample_data`

# we can see what phyloseq object `data.phy` looks like:
data.phy
    #> data.phy
    #phyloseq-class experiment-level object
    #otu_table()   OTU Table:         [ 9875 taxa and 34 samples ]
    #sample_data() Sample Data:       [ 34 samples by 12 sample variables ]
    #tax_table()   Taxonomy Table:    [ 9875 taxa by 8 taxonomic ranks ]

# now the phyloseq object `data.phy` is ready for downstream analyses!


###########i###########################i############
## save the created phyloseq object for later use ##
#######################################i############
# while it's good to be able to reproduce the phyloseq object `data.phy` from otu, taxonomy, and metadata tables every time, the code gets lengthy to start from the very beginning every single time. 
# we can save phyloseq object `data.phy` to file and we can use it directly later. 
# you can think of this as saving your progress in video games ;)

saveRDS(data.phy, "raw_data_phyloseq.RDS")
    # `saveRDS` is the function
    # it takes your phyloseq object `data.phy` and saves it as a .RDS file in your working directory. 
    # "raw_data_phyloseq.RDS" is the file name I'm using here as an example. You can name it anything you want.
    # you should see a new file named "raw_data_phyloseq.RDS" apears in your working directory after running the command. 





