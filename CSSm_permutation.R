#! /usr/bin/env Rscript

# (c) David A. Marques, 2016/2020/2021
# Script to generate permutations for the cluster-separation score (CSS) estimated with the script CSSm.R
# Group membership of individuals is permuted (optionally within strata) and the CSS recomputed for those permutations of the distance matrix

# Version history
# v1.1: index bugfix, change of group file format to more standard format
# v1.2: added optional stratified permutation, bugfix if individual order in VCF file and grouplist file are different

usage="PCACSSm_permutation.R file.vcf file.CSSm.dmat.gz file.CSSm.txt file.grouplist npermutations"
detail=paste0("\n\tfile.vcf : VCF file input (uncompressed or gzip / bgzip compressed)",
              "\n\tfile.CSSm.dmat.gz : Outfile from the CSSm.R script containing distance matrices for each window",
              "\n\tfile.CSSm.txt : Outfile from the CSSm.R script containing window coordinates and CSS values",
              "\n\tfile.grouplist : Group assignment file for all individuals in the VCF file.",
              "\n\t               Format:individual<TAB>group(Optional:<TAB>strata) (one individual per line)",
              "\n\t               (If the optional column strata is given, permutations are preformed within each stratum)",
              "\n\tnpermutations : Number of permutations to perform to estimate the empirical quantile of the observed CSS value")

# Reads input parameters
args<-commandArgs(trailingOnly=T)
vcf=args[1]
dmatn=args[2]
pmatn=args[3]
grp=args[4]
nperm=as.integer(args[5])
if(length(args)<5){
  stop(paste0("Aborted, not enough arguments given.\nUsage: ",usage,detail))
}

print("reading files")
dmat=read.delim(dmatn,h=T)
pmat=read.table(pmatn,h=T)
print("files read")

# Read grouping file
# Format of grouping file: individual\tgroup(optional:\tstrata)
grpfile<-read.delim(grp,h=F,col.names=c("ind","grp","strata")[1:dim(read.delim(grp))[2]],stringsAsFactors = F)
# Read individual order from VCF file
vcfin=readLines(vcf)
sample.id=unlist(strsplit(grep("#CHROM",vcfin,v=T),"\t"))[-c(1:9)]
# Reorder grpfile based on VCF file
grpfile=grpfile[match(sample.id,grpfile$ind),]

#######################################
# (Optionally stratified) permutation #
#######################################

# Randomly permute group assignments: either among all individuals, or within strata
# For unstratified permutation, assign a single stratum to all individuals
if(length(grpfile$stratum)==0){grpfile=data.frame(grpfile,strata=1)}

# Get list of indices for each stratum
cats=list()
for(i in unique(grpfile$strata)){
  cats[[i]]=which(grpfile$strata==i)
}

# Create permutation matrix of group assignments by resampling within strata
print("Creating permutation matrix.")
grpmat=NULL
for(i in 1:nperm){
  y=grpfile$grp
  for(j in unique(grpfile$strata)){
    y[cats[[j]]]=sample(y[cats[[j]]])
  }
  grpmat=rbind(grpmat,y)
}
# In case of identical assignments, re-run permutation and overwrite one of the duplicates
while(length(which(duplicated(grpmat)))>0){
  grpmat=grpmat[-which(duplicated(grpmat)),]
  for(i in 1:(nperm-dim(grpmat)[1])){
    y=grpfile$grp
    for(j in unique(grpfile$strata)){
      y[cats[[j]]]=sample(y[cats[[j]]])
    }
    grpmat=rbind(grpmat,y)
  }
}

# Compute CSS from permuted group assignments
print("Computing CSS scores for permutations.")

# The Euclidean distance matrices above have been converted into vectors
# with indices (1,2), (1,3), (1,4), ..., (1,n), (2,3), (2,4), ..., (2,n), (3,4),  ..., (n-1,n)

# 0) First create a matrix with these index pairs
# with first partner on first, second partner on second row (always larger than first partner / row)
pairidx=combn(dim(grpfile)[1],2)

# Permute group memership and compute CSS
cmat=sapply(1:nperm,function(x){
  
  # Draw permuted group membership
  grpfile$grp=grpmat[x,]

  ###############
  # Compute css #
  ###############

  # 1) Get indices for between group and within group comparisons
  # 1.1) Identify group assignment of individuals
  #      Thereby use sorting of vcf/genofile in "samples"
  grp1=which(grpfile$grp==unique(grpfile$grp)[1])
  grp2=which(grpfile$grp==unique(grpfile$grp)[2])
  
  # 1.2) Get indices of pairwise Euclidean distances
  # within group 1 (i/i), between groups (i/j) and within group 2(j/j)
  idxii=which(apply(pairidx,2,function(x){1*(x[1]%in%grp2)+1*(x[2]%in%grp2)})==0)
  idxij=which(apply(pairidx,2,function(x){1*(x[1]%in%grp2)+1*(x[2]%in%grp2)})==1)
  idxjj=which(apply(pairidx,2,function(x){1*(x[1]%in%grp2)+1*(x[2]%in%grp2)})==2)
  
  # 2) compute css across distance matrices and combine with position information
  # According to formula by Miller et al. 2019 Curr Biol, but without weighting by the number of sequenced bases
  s=length(grp1)
  n=length(grp2)
  unlist(apply(dmat,1,function(x){
    sum(x[idxij])/(s*n)-(1/(s+n))*(sum(x[idxii])/((s-1)/2)+sum(x[idxjj])/((n-1)/2))}))
  
})
print("CSS score compuatation done, starting q-value computation.")

# Calculate empirical quantile value for each window, add to pmat
pmat=cbind(pmat,q=apply(cbind(pmat$css,cmat),1,function(x){if(sum(is.na(x))==0){m=ecdf(x[-1]);m(x[1])}else{NA}}))
print("Q-values computed, writing output.")

# Write position/css/qval to file
write.table(pmat,paste0(sub(".dmat.gz","",dmatn),".",nperm,"perm.txt"),quote = F,sep="\t",
            row.names=F,col.names=c("chr","pos","sta","end","css","qval"))