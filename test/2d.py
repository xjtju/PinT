#!/usr/bin/env python
import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm

def get_color():
    for item in ['r','g','b','c','m','y','k']:
        yield item

def main(argv):
    if(len(argv)<3):
        print 'usage   : python pfm.py input_file label'
        print 'example : python pfm.py 1.txt bd4'
        return 1

    fname = argv[1]
    lname = argv[2]
    print 'input file is {0} \n'.format(fname)

    #labels = ['0.001', '0.01', '0.02', '0.1']

    color = get_color()
    ind = 0
    pf = h5py.File(fname, 'r')
    ds1 = pf.get('solution')
    
    ds2 = pf.get('coords')
    (totalcount, dims) = ds2.shape
    print 'coords shape is ({0},{1})\n'.format(totalcount, dims)

    posx = ds2[..., 0]
    posy = ds2[..., 1]

    grid = ds1[...]
    
    plt.scatter(posx, posy, c=grid, edgecolor='none')
    #p = ax.scatter3D(posx, posy, posz, c=grid, s=1, edgecolor='none', cmap=cm.binary)
    plt.colorbar()
    plt.title(lname)

    axes = plt.gca()
    axes.set_ylim([0.0, 1.0])
    axes.set_xlim([0.0, 1.0])
    axes.set_xlabel(r"$X$")
    axes.set_ylabel(r"$Y$")
    #p.savefig('pfm3d.png',bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main(sys.argv)

