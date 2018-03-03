#!/usr/bin/env python
import sys
import os
import h5py
import numpy as np
import matplotlib
import pylab as p 
import matplotlib.pyplot as plt 
import matplotlib.ticker as mtick 

matplotlib.rc('text', usetex=True)
plt.rc('font', family='serif')
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rcParams['text.latex.unicode']=True
matplotlib.rc('axes', linewidth=2)

def get_color():
    for item in ['r','g','b','c','m','y','k']:
        yield item

def main(argv):
    if(len(argv) < 2):
        print "usage  : python plot.py var"
        return 1
    name = argv[1]

    print 'plot object is {0}\n'.format(name)

    #with open(fname) as fp:
    #    for line in fp:
        #line = fp.readline();
        #ds = [ float(x) for x in line.split()]
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    if( name == 'cn-eu'):
        x05  = np.arange(5) +1  
        x10  = np.arange(10)+1  
        x20  = np.arange(20)+1  
        v05  = np.array([1.03e-2, 2.89e-5, 4.53e-6, 2.26e-06, 4.53e-7])
        v10  = np.array([1.04e-2, 7.86e-4, 1.33e-4, 9.50e-5,1.14e-4,9.50e-5,5.43e-5,2.04e-5,4.53e-6,4.53e-7])
        v20  = np.array([1.19E-02,1.13E-02 ,3.14E-03 ,2.19E-03 ,7.02E-03 ,1.75E-02 ,3.51E-02 ,5.70E-02 ,7.60E-02 ,8.36E-02 ,7.60E-02 ,5.70E-02 ,3.51E-02 ,1.75E-02 ,7.02E-03 ,2.19E-03 ,5.16E-04 ,8.60E-05 ,9.05E-06 ,4.53E-07])
        plt.plot(x05, v05, color='royalblue', linestyle='--', linewidth=2, marker='d', label='5')
        plt.plot(x10, v10, color='firebrick', linestyle='-',  linewidth=2, marker='s', label='10')
        plt.plot(x20, v20, color='seagreen',  linestyle=':',  linewidth=2, marker='p', label='20')
        plt.legend(loc='lower right', title='time slices', frameon=False, columnspacing=None)
        plt.title(r'CN-EU : ${\delta}T/{\delta}t=100$, ${\delta}t=2.5e-6s$')

    elif( name == 'bd4-eu'):
        x05  = np.arange(5) +1  
        x10  = np.arange(7)+1  
        x20  = np.arange(8)+1  
        x50  = np.arange(12)+1  
        x100  = np.arange(27)+1  
        x100n  = np.arange(21)+1  
        v05 = np.array([7.36E-01, 1.46E-02, 5.67E-05, 9.27E-08, 3.63E-11])
        v10 = np.array([1.75E+00 ,2.87E-02 ,2.29E-04 ,7.30E-07 ,9.20E-09 ,1.13E-10 ,1.17E-12])
        v20 = np.array([4.13E+00 ,7.75E-02 ,3.11E-03 ,6.28E-05 ,1.42E-06 ,3.95E-08 ,1.20E-09 ,3.37E-11])
        v50 = np.array([1.38E+01 ,4.63E-01 ,5.43E-02 ,4.54E-03 ,3.95E-04 ,3.41E-05 ,2.96E-06 ,2.51E-07 ,2.04E-08 ,1.56E-09 ,1.11E-10 ,7.57E-12])
        v100 = np.array([8.32E+01 ,5.34E+00 ,4.17E+00 ,1.16E+00 ,5.17E-01 ,1.94E-01 ,7.68E-02 ,2.96E-02 ,1.15E-02 ,4.41E-03 ,1.69E-03 ,6.47E-04 ,2.46E-04 ,9.29E-05 ,3.46E-05 ,1.27E-05 ,4.58E-06 ,1.62E-06 ,5.61E-07 ,1.90E-07 ,6.32E-08 ,2.06E-08 ,6.57E-09 ,2.06E-09 ,6.31E-10 ,1.90E-10 ,5.56E-11])
        v100n = np.array([1.05E+00 ,3.17E-01 ,2.43E-02 ,5.90E-03 ,1.89E-03 ,6.04E-04 ,1.99E-04 ,6.98E-05 ,2.45E-05 ,8.57E-06 ,3.00E-06 ,1.05E-06 ,3.69E-07 ,1.30E-07 ,4.56E-08 ,1.58E-08 ,5.41E-09 ,1.82E-09 ,6.06E-10 ,1.98E-10 ,6.36E-11])
        
        plt.plot(x05, v05, color='royalblue',  linewidth=2, marker='d', label='5')
        plt.plot(x10, v10, color='firebrick',  linewidth=2, marker='s', label='10')
        plt.plot(x20, v20, color='seagreen',   linewidth=2, marker='p', label='20')
        plt.plot(x50, v50, color='blue',       linewidth=2, marker='v', label='50')
        plt.plot(x100,v100, color='green',       linewidth=2, marker='h', label='100')
        plt.plot(x100n,v100n, color='red',   linewidth=2, marker='*', label='100n')
        plt.legend(loc='lower right', title='time slices', frameon=False, columnspacing=None)
        plt.title(r'BD4-EU : ${\delta}T/{\delta}t=100$, ${\delta}t=2.5e-6s$')

    elif( name == 'cn-cn'):
        x10  = np.arange(10)+1  
        x20  = np.arange(20)+1  
        x20n = np.arange(20)+1  
        v10  = np.array([ 1.04E-02 ,7.86E-04 ,1.33E-04 ,9.50E-05 ,1.14E-04 ,9.50E-05 ,5.43E-05 ,2.04E-05 ,4.52E-06 ,4.52E-07])
        v20  = np.array([ 1.18E-02 ,1.13E-02 ,3.14E-03 ,2.19E-03 ,7.02E-03 ,1.75E-02 ,3.51E-02 ,5.70E-02 ,7.60E-02 ,8.36E-02 ,7.60E-02 ,5.70E-02 ,3.51E-02 ,1.75E-02 ,7.02E-03 ,2.19E-03 ,5.16E-04 ,8.60E-05 ,9.05E-06 ,4.53E-07])
        v20n = np.array([ 5.27E-03 ,5.10E-03 ,9.23E-04 ,3.86E-04 ,1.94E-04 ,1.10E-04 ,7.50E-05 ,6.28E-05 ,5.76E-05 ,5.53E-05 ,5.03E-05 ,3.77E-05 ,2.32E-05 ,1.16E-05 ,4.64E-06 ,1.45E-06 ,3.41E-07 ,5.69E-08 ,5.99E-09 ,2.99E-10])

        plt.plot(x10, v10, color='firebrick', linewidth=2, marker='s', label='10')
        plt.plot(x20, v20, color='seagreen',  linewidth=2, marker='p', label='20')
        plt.plot(x20, v20n, color='blue',     linewidth=2, marker='v', label='20n')
        plt.legend(loc='lower right', title='time slices', frameon=False, columnspacing=None)
        plt.title(r'CN-CN : ${\delta}T/{\delta}t=100$, ${\delta}t=2.5e-6s$')

    ax = plt.gca()
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    ax.set_xlim([1, 100])
    ax.set_xticks([1, 5, 10, 20, 50, 100])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    ax.set_ylim([1e-15, 1e2])
    ax.set_yticks([1.0e-15, 1.0e-12, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-3, 1.0])
    ax.set_xlabel(r'$Iteration number K^{par}$')
    ax.set_ylabel(r'$Residual$')
    ax.yaxis.grid(True)

    plt.savefig(name+'.eps',bbox_inches='tight')
    #plt.show()

if __name__ == '__main__':
    main(sys.argv)
