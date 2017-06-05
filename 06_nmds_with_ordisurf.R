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
    
###############################################################################    
# NMDS plot showing two different experimental grouping and a surface fitting #    
###############################################################################    
# Measurement fitting showing the arrows on NMDS is great when you want to show multiple measurements.     
# But it suffers from suggesting linear changes only, and it's not ture a lot of times.     
# When we suspect that a measurement changes non-linearly according to the community variations, we can use "surface fitting"    
#     
# Because we are adding an entire surface contour lines to the background, this plot tens to get busy really quick.     
# therefore, the NMDS plot will not use ellipses to show the experimental groups. It also only plots one surface.    
# instead, it will use different colors and shapes to describe your experimental groups.     
# and it can show two different experimental factors (e.g., "Culture.media" and "Variable" in this example).    
#    
# in this example, I will plot a fitted surface contour using "X..of.singletons" as the measurement. Use different colors to represent different "Culture.media" groups. Use different shapes to represent different "Variable" groups.     
#    
# also note, in `ggplot2`, there are only 5 shapes that you can use different combinations of fill and color. So if your second experimental factor contains more than 5 groups, this will not work.     
##    
##    
# 1. if you are not continuing from previous step, load previously saved metaMDS and metadata object:        
data.mintax5.excluded.rela.mds <- readRDS("otu_sum_min_5_excluded_4samples_relative_abundance_metaMDS.RDS")    
data.mintax5.excluded.rela.si <- readRDS("otu_sum_min_5_excluded_4samples_relative_abundance_metadata.RDS")    
    
# 2. map a measurement (in this demo, it will be the "X..of.singletons" column in "data.mintax5.excluded.rela.si") to the NMDS ("data.mintax5.excluded.rela.mds")    
# import the function for mapping:     
source("~/Documents/repos/shrimp_workshop/misc_codes/ordi.sf.R")    
# mapping "X..of.singletons" to the NMDS    
data.mintax5.excluded.rela.ordisurf <- ordi.sf(data.mintax5.excluded.rela.mds, data.mintax5.excluded.rela.si[, c("X..of.singletons", "SAMPLES")], "SAMPLES")    
    # `ordi.sf`: is the customized function    
    # "data.mintax5.excluded.rela.mds": is the metaMDS object we created earlier during NMDS with ellipses plot or loaded in step 1.     
    # "data.mintax5.excluded.rela.si[, c("X..of.singletons", "SAMPLES")]" is the data.frame including the column "X..of.singletons" and the sample name column "SAMPLES", which is the same as the "totu" row.names (see NMDS with ellipses section)    
    # "SAMPLES" is the column name of column that contains the same sample name as those in "totu" table.     
    ##    
    # the output looks like below. It also returns the significance of the fitting.     
    # For statistical values, report the "p-value" in "Appromiate significance of smooth terms" and the goodness of fit from "Deviance explained"    
    #> data.mintax5.excluded.rela.ordisurf <- ordi.sf(data.mintax5.excluded.rela.mds, data.mintax5.excluded.rela.si[, c("X..of.singletons", "SAMPLES")], "SAMPLES")    
    #    
    #Family: quasipoisson     
    #Link function: log     
    #    
    #Formula:    
    #y ~ s(x1, x2, k = 10, bs = "tp", fx = FALSE)    
    #    
    #Parametric coefficients:    
    #            Estimate Std. Error t value Pr(>|t|)        
    #(Intercept)   6.2206     0.0324     192   <2e-16 ***    
    #---    
    #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1    
    #    
    #Approximate significance of smooth terms:    
    #           edf Ref.df     F p-value        
    #s(x1,x2) 7.577      9 57.57  <2e-16 ***    
    #---    
    #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1    
    #    
    #R-sq.(adj) =   0.95   Deviance explained = 96.9%    
    #-REML = 67.502  Scale est. = 12.092    n = 30    
    
# we can also take a look to see what "data.mintax5.excluded.rela.ordisurf" is like:    
head(data.mintax5.excluded.rela.ordisurf)    
    #> head(data.mintax5.excluded.rela.ordisurf)    
    #            x         y        z    
    #7  -0.5064395 -1.219561 264.8708    
    #8  -0.4374926 -1.219561 258.2193    
    #9  -0.3685457 -1.219561 251.8137    
    #10 -0.2995988 -1.219561 245.0261    
    #11 -0.2306519 -1.219561 237.4731    
    #12 -0.1617050 -1.219561 228.9169    
    ##    
    # as you can see, the result is a data.frame with 3 dimensions. These are the cordinates and modeled fitting results.     
    
# 3. generate color schemes for the experimental groups in "Culture.media":    
# create a color pallete generating function using RColorBrewer pallete "Dark2"    
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))    
# determine the number of colors needed:    
colorCount = length(unique(data.mintax5.excluded.rela.si$Culture.media))    
# get a set of colors:    
colors = getPalette(colorCount)    
    
# 4. plot the NMDS with 2 different experimental factors and the fitted measurement contour    
# import the customized plot function:    
source("~/Documents/repos/shrimp_workshop/misc_codes/ggplot.NMDS.ordisurf.2f.R")    
# plot it!    
ggplot.NMDS.ordisurf.2f(data.mintax5.excluded.rela.mds, data.mintax5.excluded.rela.si[, c("Culture.media","Variable")], colors, data.mintax5.excluded.rela.ordisurf, "X..of.singletons")    
    # `ggplot.NMDS.ordisurf.2f` is the customized function for plotting    
    # "data.mintax5.excluded.rela.mds" is the metaMDS object loaded in step 1 or from previous steps.    
    # "data.mintax5.excluded.rela.si[, c("Culture.media","Variable")]" is the data.frame with 2 different experimental factors. The first column ("Culture.media" in this case) will be plotted using different colors, and the second column ("Variable" in this case) will be plotted in different shapes.     
    # "colors" is the color scheme generated in step 2    
    # "data.mintax5.excluded.rela.ordisurf" is the fitted measurement from step 2.     
    # "X..of.singletons" tells the plot what the legend title is for the fitted surface contour    
# you can see the plot generated from this example here:
# https://github.com/germs-lab/shrimp_workshop/blob/master/example_figures/nmds_2factor_ordisurface_rk1.pdf      
