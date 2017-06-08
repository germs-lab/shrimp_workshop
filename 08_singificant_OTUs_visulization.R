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
library(RColorBrewer):39

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
# plot                          #
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



