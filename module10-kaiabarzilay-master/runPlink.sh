#!/usr/bin/env bash

# get basic summary stats about the file
plink --file hapmap1

# make a binary PED file
plink --file hapmap1 --make-bed --out hapmap1

# specify that data being used is in binary format
plink --bfile hapmap1

# generate some summary stats on rates of missing data
plink --bfile hapmap1 --missing --out miss_stat

# look at specific file
# see the number of missing individuals the proportion
more miss_stat.lmiss

# look at specific file
# see genoytyping rate
more miss_stat.imiss

# use plink to analyze data for chromosome 1
plink --bfile hapmap1 --chr 1 --out res1 --missing

# use plink to analyze data for chromosome 2
plink --bfile hapmap1 --chr 2 --out res1 --missing

# generates a file called freq_stat.frq which contains the minor allele frequency and allele codes for each SNP
plink --bfile hapmap1 --freq --out freq_stat

# perform a stratified analysis
plink --bfile hapmap1 --freq --within pop.phe --out freq_stat

# view file
more freq_stat.frq.strat

# frequency in the two populations for a specific SNP
plink --bfile hapmap1 --snp rs1891905 --freq --within pop.phe --out snp1_frq_stat

# basic association analysis on the disease trait for all single SNPs
plink --bfile hapmap1 --assoc --out as1

# sort the list of association statistics and print out the top ten
sort --key=7 -nr as1.assoc | head

# get a sorted list of association results
plink --bfile hapmap1 --assoc --adjust --out as2

# look at most significant associations
more as2.assoc.adjusted

# directly look at the inflation factor that results from having population membership as the phenotype in a case/control analysis
# use alternate phenotype option
plink --bfile hapmap1 --pheno pop.phe --assoc --adjust --out as3

# calculate association stats
#  the 2-by-3 genotype table as well as the standard allelic test
# calculate tests that assume dominant or recessive action of the minor allele
# perform the Cochran-Armitage trend test instead of the basic allelic test
plink --bfile hapmap1 --model --snp rs2222162 --out mod1

# to force the genotypic tests for this particular SNP just for illustrative purposes
plink --bfile hapmap1 --model --cell 0 --snp rs2222162 --out mod2

# perform a cluster analysis that pairs up individuals on the basis of genetic identity.
plink --bfile hapmap1 --cluster --mc 2 --ppc 0.05 --out str1

# read file that contains the results of clustering
more str1.cluster1

# use the Cochran-Mantel-Haenszel (CMH) association statistic
# get a sorted list of CMH association results
plink --bfile hapmap1 --mh --within str1.cluster2 --adjust --out aac1

# look at adjusted results file
more aac1.cmh.adjusted

# perform clustering
# do not impose a maximum cluster size: rather we request that each cluster contains at least 1 case and 1 control 
plink --bfile hapmap1 --cluster --cc --ppc 0.01 --out version2

# repeat association analysis
plink --bfile hapmap1 --mh --within version2.cluster2 --adjust --out aac2

# peform stratification analysis
# specify the number of clusters one wants in the final solution
plink --bfile hapmap1 --cluster --K 2 --out version3

# use external clustering in analysis
plink --bfile hapmap1 --mh --within pop.phe --adjust --out aac3

# generate a visualisation of the substructure in the sample by creating a matrix of pairwsie IBS distances
# use a statistical package such as R to generate a multidimensional scaling plot
plink --bfile hapmap1 --cluster --matrix --out ibd_view

# use R to perform following commands to generate plot
# color coded plot
m <- as.matrix(read.table("ibd_view.mdist"))
mds <- cmdscale(as.dist(1-m))
k <- c( rep("green",45) , rep("blue",44) )
plot(mds,pch=20,col=k)

# tell plink to use the quantitative trait
plink --bfile hapmap1 --assoc --pheno qt.phe --out quant1

# request clustered permutation
plink --bfile hapmap1 --assoc --pheno qt.phe --perm --within str1.cluster2 --out quant2

# compare each observed test statistic against the maximum of all permuted statistics in each replicate
plink --bfile hapmap1 --assoc --pheno qt.phe --mperm 1000 --within str1.cluster2 --out quant3

# test whether the association with the continuous phenotype differs between the two populations
# perform this analysis for the main SNP of interest rather than all SNPs
plink --bfile hapmap1 --pheno qt.phe --gxe --covar pop.phe --snp rs2222162 --out quant3

# extract single SNP
plink --bfile hapmap1 --snp rs2222162 --recodeAD --out rec_snp1

# to repeat the main analysis as a simple logistic regression using the R package
d <- read.table("rec_snp1.recode.raw" , header=T)
summary(glm(PHENOTYPE-1 ~ rs2222162_A, data=d, family="binomial"))

