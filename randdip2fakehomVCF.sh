#! /bin/bash

# Usage: randdip2fakehomVCF.sh filename.vcf randseed
# (c) 2015-10-11 David Marques

grep "^#" $1 > ${1%.vcf}.randhapfakehom.seed"$2".vcf
cat $1 | grep -v "^#" |\
awk -v s=$2 'BEGIN{srand(s);OFS="\t"}{l=$1;for(i=2;i<=9;i++){l=l"\t"$i};for(i=10;i<=NF;i++){
split($i,a,":");split(a[1],b,"/");c=b[1+(int(rand()*2))];sub(a[1],c"/"c,$i)
if(c==1 && length(a)>2){sub(a[length(a)],"50,50,0",$i)}else if(c==0 && length(a)>2){sub(a[length(a)],"0,50,50",$i)};l=l"\t"$i};print l}'\
>> ${1%.vcf}".randhapfakehom.seed"$2".vcf"
