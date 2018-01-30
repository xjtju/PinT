#!/usr/bin/env python
import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib

xlen = 1.0
ylen = 1.0
def get_color():
    for item in ['r','g','b','c','m','y','k']:
        yield item

def main(argv):
    if(len(argv)<2):
        print 'usage   : python pfm.py input_file label'
        print 'example : python pfm.py dbbug.h5 bd4'
        return 1

    fname = argv[1]
    lname = argv[2]
    print 'input file is {0}. \n'.format(fname)

    color = get_color()
    ind = 0
    
    pf = h5py.File(fname, 'r')

    ds = pf.get('coords')
    (totalcount, dims) = ds.shape
    print 'coords shape is ({0},{1})\n'.format(totalcount, dims)
    posx = ds[..., 0]

    ds1 = pf.get('solution')
    

    #with open(fname) as fp:
    #    for line in fp:
        #line = fp.readline();
        #ds = [ float(x) for x in line.split()]
    acolor = next(color)  
    plt.plot(posx, ds1, color=acolor, linewidth=2.0, label=lname)

    plt.legend(loc='upper right', shadow=True)
   #plt.tick_params(axis='both', which='major', labelsize=20,length=9)
   #plt.tick_params(axis='both', which='minor', labelsize=12,length=5)

    #plt.gca().xaxis.set_major_locator(MultipleLocator(1.0))
    #plt.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
    #plt.gca().yaxis.set_major_locator(MultipleLocator(0.5))
    #plt.gca().yaxis.set_minor_locator(MultipleLocator(0.1))

    #P.xlabel('x',fontsize=24)
    #P.xlim(-0.02, 1.02)
    #P.ylabel('U',fontsize=24)
    #P.ylim(-0.1, 1.1)

    plt.title(fname)
    axes = plt.gca()
    axes.set_xlim([-0.02, xlen+0.02])
    axes.set_ylim([-0.02, ylen+0.02])
    axes.set_xlabel(r'$X$')
    axes.set_ylabel(r'$\Phi$', fontsize=18)
    #P.savefig('pfm.eps',bbox_inches='tight')
    plt.grid()
    plt.xticks(np.arange(0, 1.0, 0.05))
    plt.show()

if __name__ == '__main__':
    main(sys.argv)
