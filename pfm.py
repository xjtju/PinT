#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt 

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
    index = 0
    print 'input file is {0} with {1} points \n'.format(fname,cnt)

    with open(fname) as fp:
        #for line in fp:
        line = fp.readline();
        ds = [ float(x) for x in line.split()]

    posx  = np.zeros(cnt)
    for i in range(cnt):
        posx[i] = i*1.0/cnt

    acolor = get_color()

    plt.plot(posx,ds)

    plt.xlabel('x')
    plt.ylabel('T')
    plt.title(fname)
    axes = plt.gca()
    axes.set_ylim([-0.1,1.1])
    axes.set_xlim([-0.02,1.02])
    plt.grid()
    plt.xticks(np.arange(0, 1.0, 0.05))
    plt.show()

if __name__ == '__main__':
    main(sys.argv)
