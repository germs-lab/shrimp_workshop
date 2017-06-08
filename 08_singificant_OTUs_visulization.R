##################################################
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
# data reduction                #
#################################
# as we can see, there are 680 different OTUs that were significantly different in abundance in "Nutrients" and "T0". 
# There are too many of them to plot all together.
# so we will have to do some grouping
##
# in this example, I'll group everything by "phylum" by taking the average of OTUs "log2FoldChange" by positive numbers and negative numbers

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

# 2. group by "phylum" and "pos_neg" then take average of "log2FoldChange"
res.sig.phylum <- ddply(t0.nutrients.res.sig, .(phylum, pos_neg), summarise, AVG = mean(log2FoldChange))
# the size of the newly created table is:
dim(res.sig.phylum)
    #> dim(res.sig.phylum)              
    #[1] 28  3
# and it looks like this:
head(res.sig.phylum)
    #> head(res.sig.phylum)             
    #          phylum  pos_neg       AVG
    #1  Acidobacteria positive  4.578075
    #2 Actinobacteria negative -5.292752
    #3 Actinobacteria positive  5.091601
    #4  Bacteroidetes negative -8.289403
    #5  Bacteroidetes positive  3.622159
    #6    Chloroflexi negative -5.088074

# now everything is much smaller and we can visualize the differences quickly
p1 <- ggplot(data=res.sig.phylum, aes(x=phylum, y=AVG, fill=pos_neg)) +
    geom_bar(stat="identity") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5)) 
p1
    # viola!
    # now the make it even better, we can simply add things to "p1"

# 3. fine-tuning the plot
# right now, the phyla are plotted in alphabetical order.
# let's sort the phyla to be plotted from the most changes to the least
# so we know that there are 28 rows in "res.sig.phylum"
# However:
length(unique(res.sig.phylum$phylum))
    #> length(unique(res.sig.phylum$phylum))
    #[1] 20
    ##
    # we can see that some phyla have both postive and negative numbers. 
    # 
