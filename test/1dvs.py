#!/usr/bin/env python
import sys
import os
import h5py
import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt 
import matplotlib.cm as cm

# used for pfm testcase only, you can modify it for your problem
# compare the absolute difference between two result files, especially used for checking the numerical parameters.

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
    
    ds = pf1.get('coords')
    (totalcount, dims) = ds.shape
    print 'coords shape is ({0},{1})\n'.format(totalcount, dims)
    posx = ds[..., 0]

    d1 = ds1[...]
    d2 = ds2[...]
    grid = (d1 - d2)

    color = get_color()
    acolor = next(color)  
    plt.plot(posx, grid, color=acolor, linewidth=2.0)

    ax = plt.gca()
    ax.set_xlabel('X position')
    ax.set_ylabel('absolute error')
    xvals = ax.get_xticks()
    yvals = ax.get_yticks()
    ax.set_xticks([0.0, 0.1, 0.2, 0.25, 0.3, 0.5, 0.6, 0.8, 1.0])
    #ax.set_xticklabels(['{:0.3f}%'.format(x*100) for x in xvals])
    #ax.set_yticklabels(['{:0.3f}%'.format(y*100/totalcount) for y in yvals])

    plt.title(lname)
    #axes = plt.gca()
    #ax.set_xlim([-0.001, 0.018])
    #ax.set_xlim([-0.00001,0.00026])
    #ax.set_ylim([0, 0.00006*totalcount])
    #ax.set_ylim([0, 0.004*totalcount])
    #plt.savefig('pfm3d.png',bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main(sys.argv)
