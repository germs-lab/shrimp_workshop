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

