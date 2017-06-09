##												##
## Data wrangling and visulization on 			##
## OTUs with significantly different 			##
## abundances 									##
##												##
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
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(phyloseq)

# you can also use package `dplyr`. You can find the tutorial here: http://tracykteal.github.io/R-genomics/04-dplyr.html
# `dplyr` is the newer version of `plyr` and the two overlap quite a bit. 
# load either `dplyr` or `plyr`, not both. 

#################################
# import previously saved table #
#################################
# for this part of the tutorial, I will use the table with OTUs differed significantly in abundance in "Nutrients" and "T0" 
# from step 7 (DESeq2) as an example on how to wrangle the data for visualization. 
##
# import the table we saved at the end of the step 7:
t0.nutrients.res.sig <- read.delim("nutrients_vs_t0_significant_OTUs_w_taxa.txt") 
# check the table and make sure it has the same dimension as before (see step 7). 
dim(t0.nutrients.res.sig) 
    #> dim(t0.nutrients.res.sig)        
    #[1] 680  15

#################################
# simple plot                          #
#################################
# 1. we need to have a column defines all positive numbers (more abundant in "Nutrients") and negative numbers (more abundant in "T0")
t0.nutrients.res.sig$pos_neg <- ifelse(t0.nutrients.res.sig$log2FoldChange > 0, "positive", "negative")
# and we can see the table now looks like this:
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
    #               order           family       genus  pos_neg
    #1    Rhodobacterales Rhodobacteraceae Roseovarius positive
    #2     Spirochaetales  Spirochaetaceae Spirochaeta negative
    #3 Sphingobacteriales   Saprospiraceae   Lewinella positive
    #4   Acidimicrobiales        Iamiaceae       Iamia positive
    #5     Spirochaetales  Spirochaetaceae Salinispira positive
    #6 Sphingobacteriales   Saprospiraceae   Lewinella positive

# 2. plotting
p1 <- ggplot(data=t0.nutrients.res.sig, aes(x=phylum, y=log2FoldChange)) +
    geom_point(aes(color=pos_neg))+
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5)) 
p1
    # viola!
    # now the make it even better by adding things to "p1"
# see example figure here: https://github.com/germs-lab/shrimp_workshop/blob/master/example_figures/08_step2_1.pdf

# 3. fine-tuning
# 3.1. we can add a line to emphasize "0"
p3.1 <- p1 + geom_hline(yintercept=0)
p3.1
# see example figure here: https://github.com/germs-lab/shrimp_workshop/blob/master/example_figures/08_step2_3.1.pdf

# 3.2. we can change the color of dots
p3.2 <- p3.1 + scale_color_brewer(palette="Dark2")
p3.2
# see example figure here: https://github.com/germs-lab/shrimp_workshop/blob/master/example_figures/08_step2_3.2.pdf

# 3.3. we can relabel the x and y axes
p3.3 <- p3.2 + xlab("Bacterial Phyla") + ylab("Differences in Abundance (log2)")
p3.3
# see example figure here: https://github.com/germs-lab/shrimp_workshop/blob/master/example_figures/08_step2_3.3.pdf

# 3.4. we can change the position of the legend
p3.4 <- p3.3 + theme(legend.position = "top")
p3.4
# see example figure here: https://github.com/germs-lab/shrimp_workshop/blob/master/example_figures/08_step2_3.4.pdf

# 3.5. change the size of the axis text and labels
p3.5 <- p3.4 + theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size =12), axis.title.x = element_text(face="bold", size =14), axis.title.y = element_text(face="bold", size = 14))
p3.5
# see example figure here: https://github.com/germs-lab/shrimp_workshop/blob/master/example_figures/08_step2_3.5.pdf

# 3.6. we can also remove x-axis label
p3.6 <- p3.5 + theme(axis.title.x = element_blank())
p3.6
# see example figure here: https://github.com/germs-lab/shrimp_workshop/blob/master/example_figures/08_step2_3.6.pdf


#################################################
# use DESeq2 reults to guide data visualization #
#################################################
# just because an OTU were significantly more abundant in one condition than the other,
# it does not mean it's abundant. 
## 
# I like to use the DESeq2 results as a statistical test result (like, anova, t-test, etc) to indicate what's significant
# but plotting standardized abundance (ie, relative abundance or rarefied abundance) to actually show how abundant these OTUs are.
# I will use relative abundance for the demo purpose.

# 1. load relative abundance phyloseq object saved from step 3. 
data.mintax5.excluded.rela.phy <- readRDS("otu_sum_min_5_excluded_4samples_relative_abundance_phyloseq.RDS")

# 2. subset significantly different OTUs from "T0" and "Nutrients" in "Variables"
# subset "T0" and "Nutrients" first
t0.nutrients.sig.rela.phy <- subset_samples(data.mintax5.excluded.rela.phy, grepl("T0|Nutrients", Variable))
t0.nutrients.sig.rela.phy <- prune_taxa(taxa_sums(t0.nutrient.sig.rela.phy) > 0,t0.nutrient.sig.rela.phy) 
t0.nutrients.sig.rela.phy
# subset significant OTUs
t0.nutrients.sig.rela.phy <- subset_taxa(t0.nutrient.sig.rela.phy, taxa_names(t0.nutrient.sig.rela.phy) %in% t0.nutrients.res.sig$Row.names)
# check
t0.nutrients.sig.rela.phy
    #> t0.nutrients.sig.rela.phy         
    #phyloseq-class experiment-level object
    #otu_table()   OTU Table:         [ 680 taxa and 14 samples ]
    #sample_data() Sample Data:       [ 14 samples by 12 sample variables ]
    #tax_table()   Taxonomy Table:    [ 680 taxa by 8 taxonomic ranks ]

# 3. convert the phyloseq object into a 2 dimensional long table
t0.nutrients.sig.rela.psmelt <- psmelt(t0.nutrients.sig.rela.phy)
# and the table looks like this:
head(t0.nutrients.sig.rela.psmelt)
    #> head(t0.nutrients.sig.rela.psmelt)
    #         OTU Sample Abundance                      CIBNOR_id Sample_id
    #124 OTU_1140   lp_7 0.4499120 Sediment_nutrients_72h_E1_rep3     S_300
    #123 OTU_1140   lp_8 0.4306097 Sediment_nutrients_72h_E1_rep4     S_301
    #126 OTU_1140  lp_28 0.4108014    Water_nutrients_72h_E1_rep4     S_413
    #113 OTU_1140   lp_5 0.4070750 Sediment_nutrients_72h_E1_rep1     S_298
    #117 OTU_1140   lp_6 0.3816208 Sediment_nutrients_72h_E1_rep2     S_299
    #118 OTU_1140  lp_29 0.3243486    Water_nutrients_72h_E1_rep5     S_414
    #    Culture.media  Variable Replicate           plate well DNA.concentration
    #124  Sediment 72h Nutrients         3 Magallon_Plate2   E4              11.5
    #123  Sediment 72h Nutrients         4 Magallon_Plate2   E5              12.3
    #126    Water  72h Nutrients         4 Magallon_Plate2   F3               6.6
    #113  Sediment 72h Nutrients         1 Magallon_Plate2   E2              16.7
    #117  Sediment 72h Nutrients         2 Magallon_Plate2   E3              15.2
    #118    Water  72h Nutrients         5 Magallon_Plate2   F4              11.9
    #    X..of.reads X..of.singletons Experiment SAMPLES    repseq     OTUS   domain
    #124       22794              667        RK1    lp_7 lp_5|3426 OTU_1140 Bacteria
    #123       20330              545        RK1    lp_8 lp_5|3426 OTU_1140 Bacteria
    #126       19078              149        RK1   lp_28 lp_5|3426 OTU_1140 Bacteria
    #113       19106              703        RK1    lp_5 lp_5|3426 OTU_1140 Bacteria
    #117       19881              633        RK1    lp_6 lp_5|3426 OTU_1140 Bacteria
    #118       20053              205        RK1   lp_29 lp_5|3426 OTU_1140 Bacteria
    #            phylum               class       order       family  genus
    #124 Proteobacteria Gammaproteobacteria Vibrionales Vibrionaceae Vibrio
    #123 Proteobacteria Gammaproteobacteria Vibrionales Vibrionaceae Vibrio
    #126 Proteobacteria Gammaproteobacteria Vibrionales Vibrionaceae Vibrio
    #113 Proteobacteria Gammaproteobacteria Vibrionales Vibrionaceae Vibrio
    #117 Proteobacteria Gammaproteobacteria Vibrionales Vibrionaceae Vibrio
    #118 Proteobacteria Gammaproteobacteria Vibrionales Vibrionaceae Vibrio
	##
	# note: the column "Abundance" are the relative abundance of each OTU

# 4. with the above table, we know all OTUs in the table are significantly different between "T0" and "Nutrients". We can do several plots based on this.
###################################
# 4.1 top 10 most abundant OTUs   #
###################################
	# 4.1.1. identify the total relative abundance of every OTU
	otu.total <- ddply(t0.nutrients.sig.rela.psmelt, .(OTU), summarise, total = sum(Abundance))
	head(otu.total)
            #> head(otu.total)                    
            #       OTU       total
            #1    OTU_1 0.025315704
            #2 OTU_1004 0.006212894
            #3 OTU_1017 0.007012947
            #4 OTU_1047 0.001117648
            #5 OTU_1052 0.013385506
            #6 OTU_1055 0.012644564
    	# sort by "total" in descending order 
	otu.total <- otu.total[order(-otu.total$total), ]
	head(otu.total)
            #> head(otu.total)                  
            #         OTU     total
            #9   OTU_1140 3.5361366
            #625 OTU_6780 0.8713007
            #53  OTU_1531 0.7770710
            #12  OTU_1165 0.5900284
            #558 OTU_6475 0.3707745
            #14  OTU_1189 0.3077775

	# find the top 10 most abundant ones
	otu.top10 <- otu.total[1:10, ]
	# check
	otu.top10
            #> otu.top10
            #         OTU     total
            #9   OTU_1140 3.5361366
            #625 OTU_6780 0.8713007
            #53  OTU_1531 0.7770710
            #12  OTU_1165 0.5900284
            #558 OTU_6475 0.3707745
            #14  OTU_1189 0.3077775
            #561 OTU_6478 0.2665647
            #407  OTU_336 0.2609854
            #653 OTU_9093 0.2044837
            #11  OTU_1164 0.1874009

	# 4.1.2. subset these top 10 OTUs out of table "t0.nutrients.sig.rela.psmelt"
	t0.nutrients.sig.rela.top10 <- subset(t0.nutrients.sig.rela.psmelt, OTU %in% otu.top10$OTU)
	# check the new table's OTU number
	length(unique(t0.nutrients.sig.rela.top10$OTU))
            #> length(unique(t0.nutrients.sig.rela.top10$OTU))
            #[1] 10

	# 4.1.3. now we can plot
	p4.1.3.1 <- ggplot(t0.nutrients.sig.rela.top10, aes(x = OTU, y = Abundance, fill=Variable)) + stat_summary(fun.y = mean, geom="bar", position="dodge") + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5))
	p4.1.3.1
		# we can see that ggplot2 can calculate average of "Abundance" for each "OTU" across all samples. 



	# let's plot the OTU in descending order (based on their abundance in "Nutrients")
		# we need a table with unique OTU ids to to this:
	otu.top10$OTU <- reorder(otu.top10$OTU, -otu.top10$total)
	# check the levels:
	str(otu.top10)
            #> str(otu.top10)                   
            #'data.frame':   10 obs. of  2 variables:
            # $ OTU  : Factor w/ 10 levels "OTU_1140","OTU_6780",..: 1 2 3 4 5 6 7 8 9 10
            #  ..- attr(*, "scores")= num [1:10(1d)] -3.536 -0.187 -0.59 -0.308 -0.777 ...
            #  .. ..- attr(*, "dimnames")=List of 1
            #  .. .. ..$ : chr  "OTU_1140" "OTU_1164" "OTU_1165" "OTU_1189" ...
            # $ total: num  3.536 0.871 0.777 0.59 0.371 ...
	# reorder the "OTU" levels in "t0.nutrients.sig.rela.top10" based on the "OTU" levels in "otu.top10"
		# right now, this is what is like
		str(t0.nutrients.sig.rela.top10$OTU)
            #> str(t0.nutrients.sig.rela.top10$OTU)
            # chr [1:140] "OTU_1140" "OTU_1140" "OTU_1140" "OTU_1140" ...
	t0.nutrients.sig.rela.top10$OTU <- factor(t0.nutrients.sig.rela.top10$OTU, levels=levels(otu.top10$OTU))
		# now it looks like this:
		str(t0.nutrients.sig.rela.top10$OTU)
            #> str(t0.nutrients.sig.rela.top10$OTU)
            # Factor w/ 10 levels "OTU_1140","OTU_6780",..: 1 1 1 1 1 1 1 1 2 1 ...
	# replot:
	p4.1.3.2 <- ggplot(t0.nutrients.sig.rela.top10, aes(x = OTU, y = Abundance, fill=Variable)) + stat_summary(fun.y = mean, geom="bar", position="dodge") + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5))
	p4.1.3.2
		# again, average "Abundance" was shown

	# we should add errobars
	p4.1.3.3 <- p4.1.3.2 + stat_summary(fun.data = mean_se, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)	
	p4.1.3.3
	
	# we can see that some T0 OTUs are really small, we can plot them separately
	p4.1.3.4 <- p4.1.3.3 + facet_wrap(~Variable, scale="free")
	p4.1.3.4

###################################
# 4.2 grouping by genus           # 
###################################
# because "t0.nutrients.sig.rela.psmelt" has all OTUs that were significantly more abundant in one condition than the other, we can group OTUs by genus.

    # 4.2.1. sum all OTUs from the same genus in the same sample. we will also keep the phylum information for plotting purpose
    t0.nutrients.sig.genus <- ddply(t0.nutrients.sig.rela.psmelt, .(genus, SAMPLES, Variable, phylum), summarise, total_per_sample = sum(Abundance))
    # and it looks like this
    head(t0.nutrients.sig.genus)
        #> head(t0.nutrients.sig.genus)
        #         genus SAMPLES  Variable      phylum total_per_sample
        #1 Acholeplasma   lp_15        T0 Tenericutes     0.0000000000
        #2 Acholeplasma   lp_16        T0 Tenericutes     0.0000000000
        #3 Acholeplasma   lp_17        T0 Tenericutes     0.0000000000
        #4 Acholeplasma   lp_18        T0 Tenericutes     0.0000000000
        #5 Acholeplasma   lp_25 Nutrients Tenericutes     0.0007586870
        #6 Acholeplasma   lp_26 Nutrients Tenericutes     0.0008892653
    
    # 4.2.2. now let's plot it
    p4.2.2.1 <- ggplot(t0.nutrients.sig.genus, aes(x = genus, y = total_per_sample, fill=Variable)) + stat_summary(fun.y = mean, geom="bar", position="dodge") + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5)) + stat_summary(fun.data = mean_se, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)	
    p4.2.2.1
    # clearly, it's not pretty. There are too many genera and many of them are really low in relative abundance. 
    # we can subset the top 10 most abundant genera.

    # 4.2.3. subset the top 10 most abundant genera. Similar to what we did in step 4.1.1-4.1.2.
    # first, we find the average abundance of each genus across all samples 
