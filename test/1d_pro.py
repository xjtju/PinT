#!/usr/bin/env python
import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rcParams['text.latex.unicode']=True
matplotlib.rc('axes', linewidth=2)

xlen = 1.0
ylen = 1.0

def get_color():
    for item in ['r','g', 'm', 'b', 'c', 'y','k']:
        yield item

# plot from multi files for comparison  

def main(argv):
    if(len(argv)<3):
        print 'usage   : python pfm.py number f1 label1 f2 lable2 ...'
        print 'example : python pfm.py 3 1.h5 l1 2.h5 l2 3.h5 l3'
        return 1

    fnum = int(argv[1])
    idx = 1

    print 'the numbers of input file is {0}. \n'.format(fnum)

    color = get_color()
    while (idx/2 < fnum) : 
        print 'processing the {0}th file.\n'.format(idx/2+1)

        idx = idx+1
        fname = argv[idx]
        idx = idx+1
        lname = argv[idx]
        pf = h5py.File(fname, 'r')

        ds = pf.get('coords')
        (totalcount, dims) = ds.shape
        print 'coords shape is ({0},{1})\n'.format(totalcount, dims)

        posx = ds[..., 0]
        ds1 = pf.get('solution')

        acolor = next(color)  
        plt.plot(posx, ds1, color=acolor, label=lname)

    plt.legend(loc='lower center', shadow=True)
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

    #plt.title(r'Heat equation', fontsize=18)
    plt.title(r'G1: ${k=1600~}$,  ${\beta}=0.0$', fontsize=18)
    #plt.title(r'G2: ${k=16000}$, ${\beta}=-0.128$', fontsize=18)
    #plt.title(r'G3: ${k=1600~}$,  ${\beta}=-0.128$', fontsize=18)

    axes = plt.gca()
    axes.set_xlim([-0.02, xlen+0.02])
    axes.set_ylim([-0.02, ylen+0.02])
    axes.set_xlabel(r'x position', fontsize=16)
    axes.set_ylabel(r'$\Phi$', fontsize=16)
    plt.grid()
    plt.xticks(np.arange(0, 1.02, 0.1))

    #plt.savefig('pfm.png',bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main(sys.argv)
