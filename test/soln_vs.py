#!/usr/bin/env python
import sys
import os
import h5py
import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

# used for pfm testcase only, you can modify it for your problem
# compare the difference between two result files, especially used for checking the numerical parameters.

def get_color():
    for item in ['r','g','b','c','m','y','k']:
        yield item

def main(argv):
    if(len(argv)<4):
        print 'usage   : python tsvs.py f1 f2 label'
        print 'example : python tsvs.py result1.h5 result2.h5 1vs2'
        return 1

    f1 = argv[1]
    f2 = argv[2]
    lname = argv[3]
    print 'input files are : {0} .vs. {1}\n'.format(f1, f2)

    xlen = 1.0
    ylen = 1.0
    
    pf1 = h5py.File(f1, 'r')
    ds1 = pf1.get('solution')

    pf2 = h5py.File(f2, 'r')
    ds2 = pf2.get('solution')
    
    (totalcount,)= ds1.shape
    pos = np.arange(totalcount) 

    d1 = ds1[...]
    d2 = ds2[...]
    grid = (d1 - d2)/d2  # when d2 is very small, the relative error is not accurate 
    grid = np.absolute(grid)
    grid.sort()
    gmean = np.mean(grid)
    gstd = np.std(grid)
    pdf = stats.norm.pdf(grid, gmean, gstd)

    weights = np.ones_like(grid)/float(len(grid))

    print 'the size of the result is  {0}, mean is {1},  std is {2}.\n'.format(totalcount, gmean, gstd)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plt.scatter(posx[indxy], posy[indxy], c=xygrid, edgecolor='none', cmap=cm.binary)  # x-y plane projection
    #plt.plot(pos, grid)         # x direction plot
    #xs = np.linspace(grid.min(), grid.max(), 20)
    #density = stats.gaussian_kde(grid)
    #density.covariance_factor = lambda : .25
    #density._compute_covariance()
    #plt.plot(xs, density(xs)/1000.0, '-o')
    #plt.plot(grid, pdf, '-o') # including h here is crucial
    plt.hist(grid,50, weights=weights, facecolor='g', alpha=0.75)  
    #plt.colorbar()
    #fig.colorbar(p)
    ax = plt.gca()
    ax.set_xlabel('Relative error')
    ax.set_ylabel('Distribution')
    xvals = ax.get_xticks()
    yvals = ax.get_yticks()
    #ax.set_xticklabels(['{:0.3f}%'.format(x*100) for x in xvals])
    #ax.set_yticklabels(['{:0.3f}%'.format(y*100/totalcount) for y in yvals])

    plt.title(lname)
    #axes = plt.gca()
    #ax.set_xlim([0.00001, 0.0001])
    #ax.set_xlim([-0.00001,0.00026])
    #ax.set_ylim([0, 0.00006*totalcount])
    #ax.set_ylim([0, 0.004*totalcount])
    #plt.savefig('pfm3d.png',bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main(sys.argv)
