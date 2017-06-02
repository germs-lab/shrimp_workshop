##################################################
##						##
## Customized NMDS plots using ggplot2          ##
##						##
##################################################
## Author: Fan Yang

## As we showed you in the workshop, you can make a lot of cool plots using `phyloseq`.
## such as these: 
## http://joey711.github.io/phyloseq/plot_richness-examples.html
## http://joey711.github.io/phyloseq/plot_ordination-examples.html
## http://joey711.github.io/phyloseq/plot_bar-examples.html
## However, a lot of times, these may not be enough. 
## For this part of the scripts, I'm going to show you what I do to generate some of the customized plots in R 

## Before start, navigate to where your phyloseq objects were saved at. 
## It should be the same as where your "otu count table", "taxa table", and "sample meta data" are at.
## set it as your working directory. 
## For example, in my case (and my operating system is Mac OSX):
setwd("~/Box sync/Mexico/paola/data_for_R/rk1") 
## this means from now on, all of the files I import are coming directly out of the `rk1` folder

#####################
# load libraries    #
#####################
library(phyloseq)
library(vegan)
library(ggplot2)
library(plyr)
library(RColorBrewer)

######################################################################
# NMDS plot with ellipses and arrows                                 #
######################################################################
# frequently, showing the significant experimental factor groupoing is not enough. 
# we would also like to know how do some measurements (e.g., nitrogen content) is influencing our communities. 
# we can map the measurements (vector fitting)to the community variations (NMDS) using the function `envfit` in package `vegan`. 
# here is a good tutorial on how to use package `vegan` for community analysis: http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf
# but to plot the results using package `ggplot2`, we need to do some customization. 
    
# 1. if you are not continuing from previous step, load previously saved metaMDS and metadata object:    
data.mintax5.excluded.rela.mds <- readRDS("otu_sum_min_5_excluded_4samples_relative_abundance_metaMDS.RDS")
data.mintax5.excluded.rela.si <- readRDS("otu_sum_min_5_excluded_4samples_relative_abundance_metadata.RDS")

# 2. determine the measurements and how they map to the NMDS plot
# import the customized function:
source("~/Documents/repos/shrimp_workshop/misc_codes/mds.envfit.arrows.R")
# I need columns of numbers for vector fitting. Here, I'm going to use columns "DNA.concentration", "X..of.reads",and "X..of.singletons" as examples.  

data.mintax5.excluded.rela.envfit <- mds.envfit.arrows(data.mintax5.excluded.rela.mds, data.mintax5.excluded.rela.si[, c("DNA.concentration", "X..of.reads", "X..of.singletons", "SAMPLES")], "SAMPLES")
    # the customized function is `mds.envfit.arrows`
    # the function takes the metaMDS object we created earlier (see the above NMDS plot procedures): "data.mintax5.excluded.rela.mds". It can also be imported from step 1. 
    # we also need a data.frame that includes all of the measurements we would like to evaluate: "data.mintax5.excluded.rela.si, c("DNA.con    centration", "X..of.reads", "X..of.singletons", "SAMPLES")]". 
        # note: the measure data.frame should also include a column that can shares the same sample names as the "totu" table we used to generate the metaMDS object (see the above NMDS plot procedures). In this case, the sample name column is column "SAMPLES". 
    # then we need to specify the name of the column that contains the same sample names as they are in the "totu" table. In this case, it's "SAMPLES". 

# the result is a table showing vector fitting: 
data.mintax5.excluded.rela.envfit
    #> data.mintax5.excluded.rela.envfit
    #                        MDS1        MDS2        r2  pval
    #DNA.concentration -0.4351548  0.51170117 0.4511978 0.001
    #X..of.reads       -0.5421700 -0.08051237 0.3004305 0.007
    #X..of.singletons  -0.7269523  0.59152838 0.8783655 0.001
    ##
    # the result shows:
    # row names: the measurements we evaluated
    # MDS1 and MDS2: these are the cordinates where the measurements show extend to from the center of the NMDS plot
    # r2: how much does each measurements explain the sample community variations. eg. DNA.concentration explains 45% of the sample variations we observed in NMDS plot. 
    # pval: significance. 
    ##
    # not: we only plot the measurements that are significantly fitted. so based on this results, you can further subset for for significant measurements. In this demonstration, we don't need to worry about it. 

# 3. generate color schemes:
# again, I will use the column "Culture.media" as the experimental group as the example.  
# create a color pallete generating function using RColorBrewer pallete "Dark2"
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
# determine the number of colors needed:
colorCount = length(unique(data.mintax5.excluded.rela.si$Culture.media))
# get a set of colors:
colors = getPalette(colorCount)

# 4. plot the NMDS with ellipses and arrows
# import the customized plotting function:
source("~/Documents/repos/shrimp_workshop/misc_codes/ggplot.NMDS.ellipse.arrow.R")
# plot it!
ggplot.NMDS.ellipse.arrow(data.mintax5.excluded.rela.mds, data.mintax5.excluded.rela.envfit, data.mintax5.excluded.rela.si$Culture.media, colors)
    # `ggplot.NMDS.ellipse.arrow` is the function
    # "data.mintax5.excluded.rela.mds" is the metaMDS object we created earlier during NMDS with ellipse plot 
    # "data.mintax5.excluded.rela.envfit" is the vector fitting table from step 2.
    # "data.mintax5.excluded.rela.si$Culture.media" is the column with experimental grouping information that we would like to draw ellipses for.
    # "colors" are the colors scheme we created for the experiment groups in "data.mintax5.excluded.rela.si$Culture.media"
##
# you can see the plot generated from this example here:
# https://github.com/germs-lab/shrimp_workshop/blob/master/example_figures/nmds_ellipses_arrows_rk1.pdf

