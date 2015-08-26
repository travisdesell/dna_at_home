#!/usr/bin/python

import sys

from pylab import *
from matplotlib import rc, rcParams

#rc('text',usetex=True)
#rc('font',**{'family':'serif','serif':['Computer Modern']})


#  Import the data from a text file and save as a 2D matrix.
dataMatrix1 = genfromtxt(sys.argv[1])

# Slice out required abscissae and ordinate vectors.

x = dataMatrix1[:,0]
minval = dataMatrix1[:,1]
avgval = dataMatrix1[:,2]
medval = dataMatrix1[:,3]
maxval = dataMatrix1[:,4]
stddev = dataMatrix1[:,5]

plot(x,minval,label='min')
plot(x,maxval,label='max')
plot(x,avgval,label='avg')
plot(x,medval,label='med')

legend(loc='upper right')

xlabel('step',fontsize=16)
plt.ylabel('distance',fontsize=16)

plt.savefig(sys.argv[2])
