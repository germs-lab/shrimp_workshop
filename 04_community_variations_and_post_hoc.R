##################################################
##						##
## evaluate how microbial communities differ    ##
## based on different experimental factors      ##
## and associated post-hoc analysis             ##
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
library(vegan)

######################################################################
# load the phyloseq object with potentially errorness OTUs removed   #
######################################################################
# if you recall, in step 3, we made phyloseq objects with standardized OTU abundances (ie. rarefied and relative abundance)
# both would work, but
# for the purpose of this demo, we will use the relative abundance phyloseq object
data.mintax5.excluded.rela.phy <- readRDS("otu_sum_min_5_excluded_4samples_relative_abundance_phyloseq.RDS")

data.mintax5.excluded.rela.phy
    #> data.mintax5.excluded.rela.phy   
    #phyloseq-class experiment-level object
    #otu_table()   OTU Table:         [ 4942 taxa and 30 samples ]
    #sample_data() Sample Data:       [ 30 samples by 12 sample variables ]
    #tax_table()   Taxonomy Table:    [ 4942 taxa by 8 taxonomic ranks ]

######################################################################
# calculating community distance and determine experimental factor   #
# effects using PERMANOVA                                            #
######################################################################
# The following analysis uses standardized OTU counts. 
# either rarefied or relative abundance works. 
# For the example, I'll use the relative abundance phyloseq object we generated.
# first, calculate how different the samples are based on their community structures
rela.dist = phyloseq::distance(data.mintax5.excluded.rela.phy, "bray")
    # "bray" means we are using bray-curtis distance, which weights in individual OTU abundance for distance calculation.

# second, determine which experimental factor you would like to look into:
head(sample_data(data.mintax5.excluded.rela.phy))
    # and we can see that columns "Culture.media", "Variable", and "Replicate" are the experimental factors that could influence the community distance.

# third, calculate the experimental factor effects using PERMANOVA.
adonis(rela.dist ~ Culture.media + Variable + as.factor(Replicate), perm=9999,as(sample_data(data.mintax5.excluded.rela.phy),"data.frame"))
    # we have to use "as.factor" for column "Replicate", becase column "Replicate" only has numbers. R automatically take a column with all numbers as numeric vector. 
    ###########
    #########and we can see the results are:
    #> adonis(rela.dist ~ Culture.media + Variable + as.factor(Replicate), perm=9999,as(sample_data(data.mintax5.excluded.rela.phy,),"data.frame"))
    #
    #Call:
    #adonis(formula = rela.dist ~ Culture.media + Variable + as.factor(Replicate),      data = as(sample_data(data.mintax5.excluded.rela.phy, ),          "data.frame"), permutations = 9999) 
    #
    #Permutation: free
    #Number of permutations: 9999
    #
    #Terms added sequentially (first to last)
    #
    #                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    #Culture.media         2    2.7225 1.36125 12.1448 0.29864 0.0001 ***
    #Variable              2    3.7229 1.86144 16.6074 0.40838 0.0001 ***
    #as.factor(Replicate)  4    0.3171 0.07928  0.7073 0.03478 0.8353    
    #Residuals            21    2.3538 0.11208         0.25820           
    #Total                29    9.1163                 1.00000           
    #---
    #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
        # the result says:
        # Different sample grouns in both "Culture.media" and "Variable" have significantly different communities (Pr (>F)).
        # "Variable" describes 40.838% of the community variance. "Culture.media" describes 29.864% of the community variance (R2)
        # Communities among replicates are not significantly different, which is what one would like to see.

######################################################################
# post-hoc analysis on significant experimental factors 	    ##
# using pairwised `adonis`					    ##
######################################################################
# as we can see that "Variable" contributes to the most community variations we are observing (40%). 
# however, there are 4 different factor levels in "Variable" and we don't know if they are all different from each other or acted differently. 

# 1. get the overall metadata for the samples
si<-data.frame(sample_data(data.mintax5.excluded.rela.phy))

# 2. make sure the experimental factor we would like to evaluate is factor with levels (eg, "Variable" in our case)
str(si$Variable)
    #it should look something like this:
    #> str(si$Variable)
    # Factor w/ 4 levels "Control","Nutrients",..: 1 3 3 3 3 4 4 4 4 1 ...

# 3. create a dataframe with all pairwised comparisons
pairs <- t(combn(unique(si$Variable), 2))

# 4. use a for loop:
df <- data.frame() #start an empty data.frame
for (i in 1:nrow(pairs)){
	temp.rowname <- paste(pairs[i, 1], pairs[i, 2], sep="::")
	temp.phy <- subset_samples(data.mintax5.excluded.rela.phy, Variable %in% pairs[i, ])
	temp.phy <- prune_taxa(taxa_sums(temp.phy) > 0, temp.phy)
	temp.dist <- phyloseq::distance(temp.phy, "bray")
	temp.result <- adonis(temp.dist ~ Variable, perm=9999, as(sample_data(temp.phy), "data.frame"))
	temp.df <- data.frame(temp.rowname, temp.result$aov.tab[4][1, ], temp.result$aov.tab[5][1, ], temp.result$aov.tab[6][1, ])
	df <- rbind(df, temp.df)
}
names(df)<-c("factor", "rela.F.model", "rela.adonis_R2", "rela.Pr(>F)") #add meaningful columnnames. 
    #the result looks something like this:
    #> df
    #                factor rela.F.model rela.adonis_R2 rela.Pr(>F)
    #1   Control::Prebiotic     4.293772      0.2347122      0.0041
    #2          Control::T0     4.613875      0.3157188      0.0046
    #3   Control::Nutrients    16.057036      0.5008896      0.0001
    #4        Prebiotic::T0    10.175646      0.5043529      0.0020
    #5 Prebiotic::Nutrients    15.442680      0.4911375      0.0002
    #6        T0::Nutrients    29.198579      0.7087278      0.0006
    ##
    # the results suggest that all pairs were contributing to the community variations observed significantly. 
    # and "T0 VS Nutrients" was describing 71% of the community variations. 
