#! /usr/bin/env Rscript

# (c) David A. Marques, 2018
# Usage: wsfs_pi.R -w sfsfile -b bedfile -o outfile
# Computes nucleotide diversity pi from multidimensional SFS

# Load libaries
library(optparse)

# Read input arguments
option_list = list(
  make_option(c("-w", "--wsfs"), type="character", default=F, 
              help="filename of multidimensional site-frequency spectrum file with mulitple SFS as replicates (one per line)", metavar="character"),
  make_option(c("-b", "--bedfile"), type="character", default=F, 
              help="BED file name", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=F, 
              help="output filename", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#opt$wsfs="GRA.BOH.all.2500bpwin.k5_DSFS.obs"
#opt$bedfile="winstats.all.GLA_12dmin2500loci.filt.bed"
#opt$outfile="GRA.BOH.all.2500bpwin.k5.pi"

# Reading infile into vector
infile<-scan(file=opt$wsfs,character(0), sep ="\t")
# Determine: minor or derived SFS
if(grepl("MSFS",opt$wsfs)){sfstype<-"M"}else{sfstype<-"D"}
# Determine: observed or expected SFS
if(grepl("obs",opt$wsfs)){ending<-".obs"}else{ending<-".txt"}
# Get the number of SFS stored in the file
nsfs<-as.integer(gsub(" obs.*","",infile[1]))
# Get the number of populations from the file
npop<-as.integer(infile[2])
# Read population sizes from the file
popsiz<-as.integer(infile[2+1:npop])

# Create list of multidimensional arrays and parse data into them
arr=list()
# Parse data into multidmensional array
# Note: index = allele count + 1, i.e. arr[1,1,1] = entry (0,0,0)
# IMPORTANT: fastsimcoal2 and dadi store multi-dimensional SFS in the format:
#            (0,0,0),(0,0,1),(0,0,2),(0,1,0),(0,1,1),(0,1,2),...
#            first iterating the last population!
#            thus,
#            a) population sizes are reversed while parsing the data into an array
#            and then,
#            b) the array is reversed to restore the population order as in the file
for(j in 1:nsfs){
  arr[[j]]<-array(as.numeric(infile[-c(1,2+0:npop)][((j-1)*prod(popsiz+1)+1)+0:(prod(popsiz+1)-1)]),rev(popsiz+1))
  arr[[j]]<-aperm(arr[[j]])
}

# Option 1-A: convert multidimensional to marginal 1D-SFS and fold 1D-SFS
# Create list to contain matrix
marSFS=list()
# Create matrix, one row per population and columns as in the population with most individuals
for(j in 1:nsfs){
  marSFS[[j]]<-matrix(NA,nrow=npop,ncol=max(popsiz)+1)
}
# Function to fold 2D-SFS
derived1maf=function(d1dsfs){
  n1=length(d1dsfs)
  m1dsfs=matrix(0,ncol=n1)
  m1dsfs[1:ceiling(n1/2)]=unlist(d1dsfs[1:ceiling(n1/2)])
  for(k in 0:(ceiling(n1/2)-1)){
    m1dsfs[1,k+1]=m1dsfs[1,k+1]+d1dsfs[n1-k]
  }
  m1dsfs
}
# Fill the matrix with the marginal SFS entries for each population
# i.e. the sum of all entries across the mSFS for this 1D-category
# e.g. all entries / SNPs where population 1 (size=10) has 9 major and 1 minor allele (-> category MAC=1)
for(j in 1:nsfs){
  for(i in 1:npop){
    marSFS[[j]][i,1:dim(arr[[j]])[i]]<-derived1maf(apply(arr[[j]],i,function(x) sum(x)))
  }
  colnames(marSFS[[j]])<-paste("d_",0:max(popsiz),sep="")
}

# Compute Pi for each population and SFS
# Function to compute Pi
sfs=marSFS[[1]][1,]
nucpi=function(sfs){
  n=sum(sfs)
  m=length(sfs)-1
  sum(sapply(0:(m/2),function(x){(2*x*(m-x)*sfs[x+1])/(m*(m-1))}))/n
}

# Read BED file to get coordinates of windows
bed=read.delim(opt$bedfile,col.names=c("scaffold","start","end","nsites"))
# Append columns for Pi estimates
for(i in 1:npop){
  bed=cbind(bed,pi=NA)
  names(bed)[dim(bed)[2]]=paste0("pi.",i)
}
# Fill values in matrix with pi estimates
for(j in 1:nsfs){
  for(i in 1:npop){
    bed[j,paste0("pi.",i)]=nucpi(marSFS[[j]][i,1+0:(popsiz[i])])
  }
}

# Write output pi-file
write.table(bed,opt$outfile,sep="\t",quote=F,row.names=F)
