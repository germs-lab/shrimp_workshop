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
# load the phyloseq object with relative abundance                   #
######################################################################
# if you recall, in step 3, we save a phyloseq object with relative abundance to the working directory using command: saveRDS(data.mintax5.excluded.rela.phy, "otu_sum_min_5_excluded_4samples_relative_abundance_phyloseq.RDS")
# again, you can also use the saved rarefied phyloseq object. 
# but for the demonstration purpose, I'm going to use the relative abundance phyloseq object.

# now we are going to start by read in the same phyloseq object
data.mintax5.excluded.rela.phy <- readRDS("otu_sum_min_5_excluded_4samples_relative_abundance_phyloseq.RDS")

# do a quick check of the phyloseq object:
data.mintax5.excluded.rela.phy
    #> data.mintax5.excluded.rela.phy
    #phyloseq-class experiment-level object
    #otu_table()   OTU Table:         [ 4942 taxa and 30 samples ]
    #sample_data() Sample Data:       [ 30 samples by 12 sample variables ]
    #tax_table()   Taxonomy Table:    [ 4942 taxa by 8 taxonomic ranks ]

######################################################################
# NMDS plot with ellipses                                            #
######################################################################
# Although `phyloseq` package allows you to plot NMDS plots, it could not plot ellipses (95% confidence level) around your experimental groups for a really long time. 
# Very recently, this could be done. See here: https://github.com/joey711/phyloseq/issues/323
# Please feel free to try the above method. 
# I have a similar function that adds the ellipses to each experimental group. The advantage of this function is that it automatically removes any NA groups in your experimental design. 
# I also likes this customized function, because it uses the orginal functions from package `vegan` and I know exactly what is happening. 
# to use it, we will have to take the phyloseq object apart:

# 1. get the standardized OTU count table. The OTU needs to be in columns, with sample names as row.names. 
totu<-t(data.frame(otu_table(data.mintax5.excluded.rela.phy)))
# check the dimension of the table:
dim(totu)
    #> dim(totu)
    #[1]   30 4942
    # 
    # and we can see that the number of samples are now in rows, and number of OTUs are in columns. 

# 2. create a metaMDS object:
data.mintax5.excluded.rela.mds <- metaMDS(totu, autotransform = F, k = 3, trymax = 100)
# function `metaMDS` takes the otutable, calculate "bray-curtis" distances for each pair of samples, then project them on a multi dimention space. 
# because we have already standardized our datat, I would turn the option "autotransfrom" off.
# "k" represents the number of dimention you would like to calculate to ensure a low stress score.
# "trymax" sets how many iteration of calculation will be performed. When your communities are really different, the default iteration of 20 does not ensure a solution. Therefore, I set it as 100. 
# Here is the site that explains NMDS really nicely: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
    #here is what the output looks like:
    #> data.mintax5.excluded.rela.mds <- metaMDS(totu, autotransform = F, k = 3, trymax = 100)
    #Run 0 stress 0.05309471 
    #Run 1 stress 0.05531997 
    #Run 2 stress 0.05255858 
    #... New best solution
    #... Procrustes: rmse 0.01513865  max resid 0.05879683 
    #Run 3 stress 0.05514236 
    #Run 4 stress 0.05788142 
    #Run 5 stress 0.05555345 
    #Run 6 stress 0.05913575 
    #Run 7 stress 0.05880336 
    #Run 8 stress 0.053255 
    #Run 9 stress 0.05325478 
    #Run 10 stress 0.05255946 
    #... Procrustes: rmse 0.0004554904  max resid 0.00145098 
    #... Similar to previous best
    #Run 11 stress 0.0546978 
    #Run 12 stress 0.05781932 
    #Run 13 stress 0.05555074 
    #Run 14 stress 0.05354062 
    #Run 15 stress 0.05732792 
    #Run 16 stress 0.05904354 
    #Run 17 stress 0.05732729 
    #Run 18 stress 0.05309078 
    #Run 19 stress 0.05256082 
    #... Procrustes: rmse 0.0004653287  max resid 0.0009364643 
    #... Similar to previous best
    #Run 20 stress 0.05255847 
    #... New best solution
    #... Procrustes: rmse 0.0002134823  max resid 0.000844039 
    #... Similar to previous best
    #*** Solution reached

# also the resulted object is:
data.mintax5.excluded.rela.mds
    #> data.mintax5.excluded.rela.mds
    #
    #Call:
    #metaMDS(comm = totu, k = 3, trymax = 100, autotransform = F)
    #
    #global Multidimensional Scaling using monoMDS
    #
    #Data:     totu
    #Distance: bray
    #
    #Dimensions: 3
    #Stress:     0.05255847
    #Stress type 1, weak ties
    #Two convergent solutions found after 20 tries
    #Scaling: centring, PC rotation, halfchange scaling
    #Species: expanded scores based on ‘totu’

# because every time you rerun the metaMDS, it takes time and your samples will be projected a little bit differently. 
# therefore, let's save the metaMDS object in case it's needed later as well.
saveRDS(data.mintax5.excluded.rela.mds, "otu_sum_min_5_excluded_4samples_relative_abundance_metaMDS.RDS")

# 3. get the sample metadata from the same phyloseq object:
data.mintax5.excluded.rela.si <- data.frame(sample_data(data.mintax5.excluded.rela.phy))
# also save this for later use:
saveRDS(data.mintax5.excluded.rela.si, "otu_sum_min_5_excluded_4samples_relative_abundance_metadata.RDS")

# 4. load the plot function:
source("~/Documents/repos/shrimp_workshop/misc_codes/ggplot.NMDS.ellipse.R")

# 5. determine what colors you would like to use for your experimental groups. 
## I very much dislike the default colors `phyloseq` and `ggplot2` have to offer. 
## I use package `RColorBrewer` to generate a set of subtle and easy to distinguish colors.
## I will use experimental groups defined by column "Culture.media" as an example:
    # create a color pallete generating function using RColorBrewer pallete "Dark2"
    getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
    # determine the number of colors needed:
    colorCount = length(unique(data.mintax5.excluded.rela.si$Culture.media))
    # get a set of colors:
    colors = getPalette(colorCount)

# 6. plot it!
ggplot.NMDS.ellipse(data.mintax5.excluded.rela.mds, data.mintax5.excluded.rela.si$Culture.media, colors)

# and you can repeat step 5 and 6 for other significant experimental factors (see step 3)

