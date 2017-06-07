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

# 3. identify the significantly different OTUs
t0.nutrients.dds <- DESeq(t0.nutrients.dds, test="Wald", fitType="parametric")
    #it should give outputs like this:
    #> t0.nutrients.dds <- DESeq(t0.nutrients.dds, test="Wald", fitType="parametric")
    #using pre-existing size factors
    #estimating dispersions
    #found already estimated dispersions, replacing these
    #gene-wise dispersion estimates
    #mean-dispersion relationship
    #final dispersion estimates
    #fitting model and testing
    #-- replacing outliers and refitting for 16 genes
    #-- DESeq argument 'minReplicatesForReplace' = 7 
    #-- original counts are preserved in counts(dds)
    #estimating dispersions
    #fitting model and testing
	##
	# the comparison is in the same order how your experimental levels are (in this case, "Variable")
	##
    # we can take a look at the levels in "Variable"
    levels(t0.nutrients.dds$Variable)
        #> levels(t0.nutrients.dds$Variable)
        #[1] "Nutrients" "T0"    
		##
        # therefore,the comparison will be comparing the OTU abundances in "Nutrients" to the OTU abundances in "T0"
	
# get the results into a readable table format:
t0.nutrients.res <- data.frame(results(t0.nutrients.dds, cooksCutoff=F))
# it now looks something like this:
head(t0.nutrients.res)
    #> head(t0.nutrients.res)           
    #            baseMean log2FoldChange     lfcSE       stat       pvalue
    #OTU_1    25.07814166       3.106299 0.8492941  3.6575073 2.546800e-04
    #OTU_100   0.19408364      -2.786412 2.5045233 -1.1125519 2.659009e-01
    #OTU_1003  0.08975577       2.106717 2.5720000  0.8190968 4.127312e-01
    #OTU_1004 10.66694147      -6.394586 1.3484221 -4.7422730 2.113336e-06
    #OTU_1005  0.18932912      -1.608187 2.1831842 -0.7366245 4.613507e-01
    #OTU_1007  0.06469455      -1.633472 2.6683278 -0.6121707 5.404249e-01
    #                 padj
    #OTU_1    2.612289e-03
    #OTU_100            NA
    #OTU_1003           NA
    #OTU_1004 5.578751e-05
    #OTU_1005           NA
    #OTU_1007           NA
	##
	# the significance of the differences for each OTU is displayed in column "padj"
    # the column "log2FoldChange" shows how different the OTU abundances were in "Nutrients" vs. "T0" in log2 format
    # a positive number means it's more abundant in "Nutrients"
    # a negative number means it's more abundant in "T0"

# subset for "padj" less than 0.05:
t0.nutrients.res.sig <- subset(t0.nutrients.res, padj < 0.05) 
# and now it looks like this:
head(t0.nutrients.res.sig)
    #> head(t0.nutrients.res.sig)       
    #           baseMean log2FoldChange     lfcSE      stat       pvalue
    #OTU_1    25.0781417       3.106299 0.8492941  3.657507 2.546800e-04
    #OTU_1004 10.6669415      -6.394586 1.3484221 -4.742273 2.113336e-06
    #OTU_1017  8.1745704       3.504332 1.0595986  3.307226 9.422475e-04
    #OTU_1047  0.9562083       3.728892 1.4264269  2.614148 8.945020e-03
    #OTU_1052 15.3830749       3.017087 0.9068023  3.327171 8.773244e-04
    #OTU_1055 14.4562508       2.845759 1.0439347  2.725994 6.410824e-03
    #                 padj
    #OTU_1    2.612289e-03
    #OTU_1004 5.578751e-05
    #OTU_1017 6.720534e-03
    #OTU_1047 3.304901e-02
    #OTU_1052 6.405955e-03
    #OTU_1055 2.610381e-02

# 4. linking taxonomy information to the signficantly different OTUs. 
# it's great to know which OTU is signficantly more abundant in one condition than the other
# but, we don't know exactly what they are yet. 
# However, we do have every OTU's taxonomy information stored in our phyloseq object (ie., the "t0.nutrients" phyloseq object from step 1 in this section)
##
t0.nutrients 
    #> t0.nutrients                     
    #phyloseq-class experiment-level object
    #otu_table()   OTU Table:         [ 4029 taxa and 14 samples ]
    #sample_data() Sample Data:       [ 14 samples by 12 sample variables ]
    #tax_table()   Taxonomy Table:    [ 4029 taxa by 8 taxonomic ranks ]

# get the associated taxonomy table:
t0.nutrients.tax <- data.frame(tax_table(t0.nutrients)) 
# and it looks like this:
head(t0.nutrients.tax)
    #> head(t0.nutrients.tax)           
    #             repseq     OTUS   domain                phylum
    #OTU_1    lp_0|26174    OTU_1 Bacteria        Proteobacteria
    #OTU_100   lp_0|5210  OTU_100 Bacteria unclassified_Bacteria
    #OTU_1003  lp_3|7409 OTU_1003  Archaea        Woesearchaeota
    #OTU_1004  lp_3|8356 OTU_1004 Bacteria          Spirochaetes
    #OTU_1005  lp_3|1731 OTU_1005  Archaea        Woesearchaeota
    #OTU_1007 lp_3|22037 OTU_1007  Archaea         Pacearchaeota
    #                               class                       order
    #OTU_1            Alphaproteobacteria             Rhodobacterales
    #OTU_100        unclassified_Bacteria       unclassified_Bacteria
    #OTU_1003 unclassified_Woesearchaeota unclassified_Woesearchaeota
    #OTU_1004                Spirochaetia              Spirochaetales
    #OTU_1005 unclassified_Woesearchaeota unclassified_Woesearchaeota
    #OTU_1007  unclassified_Pacearchaeota  unclassified_Pacearchaeota
    #                              family                              genus
    #OTU_1               Rhodobacteraceae                        Roseovarius
    #OTU_100        unclassified_Bacteria              unclassified_Bacteria
    #OTU_1003 unclassified_Woesearchaeota Woesearchaeota Incertae Sedis AR16
    #OTU_1004             Spirochaetaceae                        Spirochaeta
    #OTU_1005 unclassified_Woesearchaeota        unclassified_Woesearchaeota
    #OTU_1007  unclassified_Pacearchaeota  Pacearchaeota Incertae Sedis AR13

# add the taxonomy information to the OTUs that are significantly different in abundances using function `merge`
# first, let's check the dimension of "t0.nutrients.res.sig"
dim(t0.nutrients.res.sig)
    #> dim(t0.nutrients.res.sig)        
    #[1] 680   6
    ##
	# this means that our final merged data table should have 680 rows. If it doesn't, something went really wrong. 
# now, let's merge in the taxonomy information
t0.nutrients.res.sig <- merge(t0.nutrients.res.sig, t0.nutrients.tax, by = "row.names")
# check the dimension again:
dim(t0.nutrients.res.sig)
    #> dim(t0.nutrients.res.sig)        
    #[1] 680  15
	##
	# and it has the same number of rows before merging
	##
# and the table looks like this now:
head(t0.nutrients.res.sig)
    #> head(t0.nutrients.res.sig)       
    #  Row.names   baseMean log2FoldChange     lfcSE      stat       pvalue
    #1     OTU_1 25.0781417       3.106299 0.8492941  3.657507 2.546800e-04
    #2  OTU_1004 10.6669415      -6.394586 1.3484221 -4.742273 2.113336e-06
    #3  OTU_1017  8.1745704       3.504332 1.0595986  3.307226 9.422475e-04
    #4  OTU_1047  0.9562083       3.728892 1.4264269  2.614148 8.945020e-03
    #5  OTU_1052 15.3830749       3.017087 0.9068023  3.327171 8.773244e-04
    #6  OTU_1055 14.4562508       2.845759 1.0439347  2.725994 6.410824e-03
    #          padj     repseq     OTUS   domain         phylum               class
    #1 2.612289e-03 lp_0|26174    OTU_1 Bacteria Proteobacteria Alphaproteobacteria
    #2 5.578751e-05  lp_3|8356 OTU_1004 Bacteria   Spirochaetes        Spirochaetia
    #3 6.720534e-03 lp_3|17051 OTU_1017 Bacteria  Bacteroidetes    Sphingobacteriia
    #4 3.304901e-02 lp_4|13277 OTU_1047 Bacteria Actinobacteria      Actinobacteria
    #5 6.405955e-03  lp_4|3724 OTU_1052 Bacteria   Spirochaetes        Spirochaetia
    #6 2.610381e-02  lp_4|8991 OTU_1055 Bacteria  Bacteroidetes    Sphingobacteriia
    #               order           family       genus
    #1    Rhodobacterales Rhodobacteraceae Roseovarius
    #2     Spirochaetales  Spirochaetaceae Spirochaeta
    #3 Sphingobacteriales   Saprospiraceae   Lewinella
    #4   Acidimicrobiales        Iamiaceae       Iamia
    #5     Spirochaetales  Spirochaetaceae Salinispira
    #6 Sphingobacteriales   Saprospiraceae   Lewinella

# 5. save this results as a table for later use (e.g., plotting)
write.table(t0.nutrients.res.sig, "nutrients_vs_t0_significant_OTUs_w_taxa.txt", sep="\t", quote=F, row.names=F) 
 
