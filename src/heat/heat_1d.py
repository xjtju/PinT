#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt 

def get_color():
    for item in ['r','g','b','c','m','y','k', 'r','g','b','c','m','y','k', 'r','g','b','c','m','y','k', 'r','g','b','c','m','y','k']:
        yield item

def main(argv):
    if(len(argv)<2):
        print 'usage   : python heat.py input_file count'
        print 'example : python heat.py 1.txt 100'
        return 1

    fname = argv[1]
    cnt   = int(argv[2])
    ds    = np.zeros(cnt)
    index = 0
    print 'input file is {0} with {1} lines \n'.format(fname,cnt)
    with open(fname) as fp:
        for line in fp:
            ds[index] = float(line)
            index += 1

    posx  = np.zeros(cnt)
    for i in range(cnt):
        posx[i] = i*1.0/cnt

    acolor = get_color()

    plt.plot(posx,ds)

    plt.xlabel('x')
    plt.ylabel('T')
    plt.title(fname)
    axes = plt.gca()
    axes.set_ylim([0,1.0])

    plt.show()

if __name__ == '__main__':
    main(sys.argv)
