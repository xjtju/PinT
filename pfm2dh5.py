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
    if(len(argv)<2):
        print 'usage   : python pfm.py input_file xpoints'
        print 'example : python pfm.py 1.txt 1000'
        return 1

    fname = argv[1]
    cnt   = int(argv[2])
    print 'input file is {0} with {1} points \n'.format(fname,cnt)

    xlen = 1.0
    ylen = 1.0
    #posx  = np.zeros(cnt)
    #posy  = np.zeros(cnt)
    
    labels = ['0.001', '0.01', '0.02', '0.1']
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
   # grid = tmp.reshape(cnt, cnt) 
    
    #plt.imshow(grid, interpolation='nearest')
    plt.scatter(posx, posy, c=grid)
    #plt.colorbar()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(fname)
    #axes = plt.gca()
    #axes.set_ylim([-0.1,1.1])
    #axes.set_xlim([-0.02, xlen+0.02])
    #P.savefig('pfm.eps',bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main(sys.argv)
