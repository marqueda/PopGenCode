#! /usr/bin/env python

# Author: (c) David Marques, Dec 14, 2018, Bern
# Written for Python 3.4.3
# Modified from dxy_wsfs.py

import argparse, os, re, itertools
import numpy as np

parser=argparse.ArgumentParser(description='Computes Dxy from window-sfs file.')

parser.add_argument('-w', '--wsfs', dest='w', help='Multidimensional SFS file (MSFS / DSFS) with 2 populations and windows stored as replicate SFS (e.g. multiple lines) [required]', required=True)
parser.add_argument('-b', '--bed', dest='b', help='BED file with window coordinates, needs to have the same no. of lines as no. of replicate SFS in the SFS file [required]', required=True)
parser.add_argument('-o', '--out', dest='o', help='name of outfile [required]', required=True)

args=parser.parse_args()

# Get number of windows and population sizes
wsfs=open(args.w)
for i, line in enumerate(wsfs):
	if i == 0:
		# Read number of SFS windows in the file
		nsfs=int(re.sub(r' obs.*\n','',line))
	elif i == 1:
		# Read population sizes into list
		m=[int(x) for x in line.strip('\n').split('\t')][1]
		n=[int(x) for x in line.strip('\n').split('\t')][2]
	elif i > 1:
		break
wsfs.close()

# Read window coordinates into lists
bed=open(args.b,'r')
wincoord=[]
for line in bed:
	if not (re.match('^track',line) or re.match('^browser',line) or re.match('^#',line)):
		wincoord.append((line.strip("\n").split("\t"))[0:3])
bed.close()

# Calculate weighting vector for SFS entries
idp=[list(x) for x in itertools.product(*[list(range(0,m+1)),list(range(0,n+1))])]
wei=[]
for j in range(0,np.size(idp,0)):
	wei=np.append(wei,np.mean(np.absolute(np.diff([list(x) for x in itertools.product(*[list(np.append(np.repeat([0],idp[j][0]),np.repeat([1],m-idp[j][0]))),list(np.append(np.repeat([0],idp[j][1]),np.repeat([1],n-idp[j][1])))])]))))

# Loop through SFS file and calculate pi, dxy
wsfs=open(args.w,'r')
outfile=open(args.o,'w')
outfile.write("scaffold\tstart\tend\tnsites\tnsnps\tdxy\n")
for i, line in enumerate(wsfs):
	if i > 1:
		cols=line.strip("\n").split("\t")
		sfs=cols[0:len(cols)]
		sfs=[float(x) for x in sfs]
		if sum(sfs) > 0:
			dxy=np.divide(np.sum(np.multiply(wei,np.asarray(sfs))),np.sum(np.asarray(sfs)))
		else:
			dxy=0
		outfile.write('\t'.join(wincoord[i-2])+"\t"+str(int(sum(sfs)))+"\t"+str(int(sum(sfs[1:-1])))+"\t"+str(dxy)+'\n')
wsfs.close()
outfile.close()
