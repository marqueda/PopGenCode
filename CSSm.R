#! /usr/bin/env Rscript

# (c) David A. Marques, 2016/2020
# Script to compute the cluster-separation score used by Jones et al. 2012 Nature
# with corrections to the formula published by Miller et al. 2019 Curr Biol
# from a VCF file and and grouping file (format: individual\tgroup, one individual per line)
# using either principal component analysis or multidimensional scaling

# Version history
# v1.1: index bugfix, change of group file format to more standard format
# v1.2: included pca or mds option
# v1.3: replaced snpgdsSlidingWindow function with loop over chromosomes and apply function due to bug in SNPRelate
# v1.4: fixed bugs in output (number of SNPs per window, new: tab separated values not in scientific notation)
#       and added optional specifying minor allele frequency
# v1.5: cleaning up script

usage="CSSm.R file.vcf file.grouplist windowsize stepsize minsnpperwindow [locus|basepair] [pca|mds] minorallelefrequency"
detail=paste0("\n\tfile.vcf : VCF file input (uncompressed or gzip / bgzip compressed)",
              "\n\tfile.grouplist : Group assignment file for all individuals in the VCF file.",
              "\n\t               Format:individual<TAB>group(Optional:<TAB>strata) (one individual per line)",
              "\n\t               (Optional column strata relevant for permutations only -> see script CSSm_permutation.R)",
              "\n\twindowsize : Size of window",
              "\n\tstepsize : Size between window starts (stepsize<windowsize -> sliding windows)",
              "\n\tminsnpperwindow : Minimum number of SNPs for CSS to be computed",
              "\n\t[locus|basepair] : Unit of window- and stepsizes as number of SNPs (locus) or base pairs (basepair)",
              "\n\t[pca|mds] : Method to be used to compute CSS, either",
              "\n\t            Principal Component Analysis (pca), using the first two PC axes as PCA space, or",
              "\n\t            Multidimensional Scaling (mds) with two factors",
              "\n\tminorallelefrequency : Minor allele frequency filter for SNPs to be used, e.g. 0.05 means MAF>=5%")

# Load libraries
library(SNPRelate)
library(intervals)

# Reads input parameters
args<-commandArgs(trailingOnly=T)
vcf=args[1]
grp=args[2]
win=as.integer(args[3])
step=as.integer(args[4])
minsnp=as.integer(args[5])
uni=args[6]
method=args[7]
maf=as.numeric(args[8])
if(length(args)<8){
  stop(paste0("Aborted, not enough arguments given.\nUsage: ",usage,detail))
}

# Load datafiles
# Convert VCF to GDS file
gdsfile=paste0(sub(".vcf.*","",vcf),".gds")
if(file.exists(gdsfile)==FALSE){
  snpgdsVCF2GDS(vcf,gdsfile, method="biallelic.only")
}else{
  print(paste0("File ",gdsfile," already exists -> not overwritten."))
}

# Read gds file
genofile<-snpgdsOpen(gdsfile)

# Read grouping file
# Format of group file: individual\tgroup
# individual names must be equivalent to VCF individual names
grpfile<-read.delim(grp,h=F,col.names=c("ind","grp","strata")[1:dim(read.delim(grp))[2]],stringsAsFactors = F)
sample.id=read.gdsn(index.gdsn(genofile, "sample.id"))

# Create sample and SNP subsets based on group file
# Order of the group file does not matter - the relevant order (for later) is the VCF/genofile
samples=sample.id[sort(match(grpfile$ind,sample.id))]
snps=snpgdsSelectSNP(genofile,maf=maf,autosome.only=F,sample.id=samples)

# Get SNP coordinates
snpdf=snpgdsSNPList(genofile,sample.id=samples)[snps,2:3]

# Check whether all individuals in grouping file are present in VCF
if(sum(is.na(match(grpfile$ind,sample.id))>0)){
  stop("Aborted, because sample names in VCF file and grouping file not identical.")
}

# Define function to compute Euclidean distance in PCA space for a list of SNPs
pcadist=function(snpid){
  if(length(snpid)>minsnp & sum(is.na(snpid))==0){
    # Perform principal components analysis
    m=snpgdsPCA(genofile,snp.id=snpid,sample.id=samples,num.thread=1,
                autosome.only=F,verbose=F,eigen.cnt=2)$eigenvect
    # Compute pairwise Euclidean distances in PC1 / PC2 space
    c(dist(m,method="euclidean"))
  }else{
    NA
  }
}

# Define function to compute Euclidean distance in MDS space for a list of SNPs
mdsdist=function(snpid){
  if(length(snpid)>minsnp & sum(is.na(snpid))==0){
    # Compute IBS pairwise distances for a subset of SNPs
    ibs=snpgdsIBS(genofile,snp.id=snpid,sample.id=samples,num.thread=1,
                  autosome.only=F,verbose=F)
    # Perform multidimensional scaling into two components
    mds=cmdscale(1-ibs$ibs,k=2)
    # Compute pairwise Euclidean distances in MDS space
    c(dist(mds,method="euclidean"))
  }else{
    NA
  }
}

# Define function to compute cluster separation score (CSs)
{
  # The Euclidean distance matrices above have been converted into vectors
  # with indices (1,2), (1,3), (1,4), ..., (1,n), (2,3), (2,4), ..., (2,n), (3,4),  ..., (n-1,n)
  
  # 0) First create a matrix with these index pairs
  # with first partner on first, second partner on second row (always larger than first partner / row)
  pairidx=combn(dim(grpfile)[1],2)
  
  # 1) Get indices for between group and within group comparisons
  # 1.1) Identify group assignment of individuals
  #      Thereby use sorting of vcf/genofile in "samples"
  grp1=which(grpfile$grp[match(samples,grpfile$ind)]==unique(grpfile$grp)[1])
  grp2=which(grpfile$grp[match(samples,grpfile$ind)]==unique(grpfile$grp)[2])
  
  # 1.2) Get indices of pairwise Euclidean distances
  # within group 1 (i/i), between groups (i/j) and within group 2(j/j)
  idxii=which(apply(pairidx,2,function(x){1*(x[1]%in%grp2)+1*(x[2]%in%grp2)})==0)
  idxij=which(apply(pairidx,2,function(x){1*(x[1]%in%grp2)+1*(x[2]%in%grp2)})==1)
  idxjj=which(apply(pairidx,2,function(x){1*(x[1]%in%grp2)+1*(x[2]%in%grp2)})==2)
  
  # 2) compute css across distance matrices and combine with position information
  # According to formula by Miller et al. 2019 Curr Biol, but without weighting by the number of sequenced bases
  s=length(grp1)
  n=length(grp2)
  
  # 3) CSS function
  css=function(distmattable){
    return(unlist(apply(distmattable,1,function(x){
      sum(x[idxij])/(s*n)-(1/(s+n))*(sum(x[idxii])/((s-1)/2)+sum(x[idxjj])/((n-1)/2))})))
  }
}

# Loop over available chromsomes
chr=unique(read.gdsn(index.gdsn(genofile, "snp.chromosome")))
for(i in chr){
  # Define sliding windows for chromosome
  snpsub=snps[snpdf$chromosome==i]
  # Get coordinates of SNPs
  snpcoord=snpdf$position[snpdf$chromosome==i]
  if(uni=="basepair"){
    # Create windows and overlap with SNPs
    winsta=seq(1,max(snpcoord),by=step)
    winend=winsta+win-1
    # Overlap coordinates
    wlist=lapply(interval_overlap(Intervals(cbind(winsta,winend)),snpcoord),function(x){snpsub[x]})
  }else if(uni=="locus"){
    # Group loci into windows of length step
    wlist=split(snpsub,ceiling(seq_along(snpsub)/step))
    # Get start and end positions from grouping
    winsta=unlist(lapply(wlist,function(x){snpcoord[x[1]]}))
    winend=unlist(lapply(wlist,function(x){snpcoord[rev(x)[1]]}))
  }
  
  # Compute distance matrices in windows
  if(method=="pca"){
    dlist=lapply(wlist,pcadist)
  }else if(method=="mds"){
    dlist=lapply(wlist,mdsdist)
  }
  
  # Combine window distance matrices and coordinates into matrix
  dmat=do.call(rbind,c(dlist))
  pmat=data.frame(chr=i,sta=winsta,end=winend,nsnps=unlist(lapply(wlist,length)))
  
  # Remove empty rows in both matrices
  dmat=dmat[!is.na(pmat[,2]),]
  pmat=pmat[!is.na(pmat[,2]),]
  pmat=pmat[!is.na(dmat[,1]),]
  dmat=dmat[!is.na(dmat[,1]),]

  # Compute CSS
  pmat=cbind(pmat,css=css(dmat))

  # Save to file
  if(i==chr[1]){
    # Write chromosome, window coordinates, no. of SNPs and CSS estimate to file *css.txt
    write.table(format(pmat,scientific=F),paste0(sub(".vcf.*","",vcf),".",format(win,scientific=F),uni,format(step,scientific=F),
                            "step.window.",method,".",
                            tools::file_path_sans_ext(grp),".CSSm.txt"),quote = F,sep="\t",
                row.names=F,col.names=c("chr","sta","end","nsnps","css"))
    # Write Euclidean distances for each window to file *CSSm.dmat.gz
    gz1=gzfile(paste0(sub(".vcf.*","",vcf),".",format(win,scientific=F),uni,format(step,scientific=F),
                      "step.window.",method,".",
                      tools::file_path_sans_ext(grp),".CSSm.dmat.gz"),"w")
    write.table(dmat,gz1,quote = F,sep="\t",
                row.names=F,col.names=apply(pairidx,2,function(x){paste("d",x[2],"_",x[1],sep="")}))
  }else{
    # Append chromosome, mean SNP position of the window and CSS estimate to file *css.txt
    write.table(format(pmat,scientific=F),paste0(sub(".vcf.*","",vcf),".",format(win,scientific=F),uni,format(step,scientific=F),
                            "step.window.",method,".",
                            tools::file_path_sans_ext(grp),".CSSm.txt"),quote = F,sep="\t",
                row.names=F,col.names=F,append=T)
    # Write Euclidean distances for each window to file *CSSm.dmat.gz
    write.table(dmat,gz1,quote = F,append=T,sep="\t",
                row.names=F,col.names=F)
  }
}

# Close infile and outfile
snpgdsClose(genofile)
close(gz1)
