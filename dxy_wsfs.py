#! /usr/bin/env python

# Author: (c) David Marques, June 1, 2016, Victoria BC, Canada
# Written for Python 3.4.3
# Changelog:
# 2016-09-21: adapted to loop across array of window SFS

import argparse, os
import numpy as np

parser=argparse.ArgumentParser(description='dxy from two angsd .saf.idx files and one bed file specifying window coordinates')

parser.add_argument('-w', '--wsfs', dest='w', help='window sfs file [required]', required=True)
parser.add_argument('-m', '--npop1', dest='m', help='no. of individuals in population 1 [required]', required=True)
parser.add_argument('-n', '--npop2', dest='n', help='no. of individuals in population 2 [required]', required=True)
parser.add_argument('-o', '--out', dest='o', help='name of outfile [required]', required=True)

args=parser.parse_args()

# Calculate number of alleles
m=int(args.m)*2
n=int(args.n)*2

# function to get all pairwise pairs of two arrays
def cartesian(arrays, out=None):
	arrays = [np.asarray(x) for x in arrays]
	dtype = arrays[0].dtype

	n = np.prod([x.size for x in arrays])
	if out is None:
		out = np.zeros([n, len(arrays)], dtype=dtype)

	m = n / arrays[0].size
	out[:,0] = np.repeat(arrays[0], m)
	if arrays[1:]:
		cartesian(arrays[1:], out=out[0:m,1:])
		for j in range(1, arrays[0].size):
			out[j*m:(j+1)*m,1:] = out[0:m,1:]
	return out

mat=np.zeros((m+1,n+1))

# Calculate weighting vector for SFS entries
idp=cartesian((np.arange(m+1),np.arange(n+1)))
wei=[]
for j in range(0,np.size(idp,0)):
	wei=np.append(wei,np.mean(np.absolute(np.diff(cartesian((np.append(np.repeat([0],idp[j][0]),np.repeat([1],m-idp[j][0])),np.append(np.repeat([0],idp[j][1]),np.repeat([1],n-idp[j][1]))))))))

# Loop through bed file and calculate dxy
wsfs=open(args.w,'r')
out=open(args.o,'w')
out.write("scaffold\tstart\tend\tdxy\n")
for line in wsfs:
	cols=line.strip("\n").split(" ")
	if len(cols)>4:
		sfs=cols[3:(len(cols)-1)]
		sfs=[float(x) for x in sfs]
		dxy=np.divide(np.sum(np.multiply(wei,np.asarray(sfs))),np.sum(np.asarray(sfs)))
		out.write('\t'.join(cols[0:3])+"\t"+str(dxy)+'\n')
wsfs.close()
out.close()
