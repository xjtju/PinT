#!/usr/bin/env python
import sys
import os
import numpy as np
#import pylab as P
import matplotlib.pyplot as plt 
import matplotlib
from matplotlib.ticker import MultipleLocator

#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rcParams['text.latex.unicode']=True
#matplotlib.rc('axes', linewidth=2)
#matplotlib.rc('legend', fontsize=22)

xlen = 1.0
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

    posx  = np.zeros(cnt)
    for i in range(cnt):
        posx[i] = i*xlen/cnt
    labels = ['0.001', '0.01', '0.02', '0.1']
    color = get_color()
    ind = 0
    with open(fname) as fp:
        for line in fp:
        #line = fp.readline();
            ds = [ float(x) for x in line.split()]
            acolor = next(color)  
            plt.plot(posx, ds, color=acolor, linewidth=2.0, label=labels[ind])
            ind = ind + 1

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

    plt.xlabel('x')
    plt.ylabel('U')
    plt.title(fname)
    axes = plt.gca()
    axes.set_ylim([-0.1,1.1])
    axes.set_xlim([-0.02, xlen+0.02])
    #P.savefig('pfm.eps',bbox_inches='tight')
    plt.grid()
    plt.xticks(np.arange(0, 1.0, 0.05))
    plt.show()

if __name__ == '__main__':
    main(sys.argv)
