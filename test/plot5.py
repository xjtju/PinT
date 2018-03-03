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

    if( name == 'cn-eu-5'):
        x10  = np.arange(10)+1  
        x20  = np.arange(20)+1  
        x50  = np.arange(20)+1  
        x100 = np.arange(20)+1  
        x125 = np.arange(20)+1  

        v10  = np.array([ 1.93E-03 ,7.83E-08 ,6.20E-11 ,4.36E-11 ,3.21E-11 ,7.91E-11 ,5.61E-11 ,1.52E-11 ,7.62E-12 ,1.29E-15])
        v20  = np.array([ 1.89E-03 ,2.88E-07 ,2.40E-09 ,1.22E-10 ,7.54E-11 ,1.56E-10 ,1.53E-10 ,6.01E-11 ,5.56E-11 ,8.04E-11 ,3.84E-11 ,5.45E-11 ,6.56E-11 ,6.77E-11 ,1.63E-11 ,1.31E-11 ,1.79E-11 ,2.33E-14 ,1.26E-15 ,7.26E-16])
        v50  = np.array([ 1.80E-03 ,1.50E-04 ,1.22E-05 ,1.26E-06 ,1.77E-07 ,3.83E-08 ,8.44E-09 ,1.87E-09 ,4.57E-10 ,2.66E-10 ,1.80E-10 ,4.92E-11 ,1.25E-10 ,8.49E-11 ,9.60E-11 ,6.18E-11 ,7.38E-11 ,1.20E-10 ,9.60E-11 ,9.17E-11])
        v100 = np.array([ 3.72E-03 ,3.48E-03 ,5.99E-04 ,2.08E-04 ,7.93E-05 ,3.78E-05 ,2.33E-05 ,1.64E-05 ,1.21E-05 ,9.08E-06 ,6.98E-06 ,5.43E-06 ,4.28E-06 ,3.40E-06 ,2.72E-06 ,2.20E-06 ,1.79E-06 ,1.46E-06 ,1.21E-06 ,1.00E-06])
        v125 = np.array([ 8.56E-03 ,8.39E-03 ,1.79E-03 ,8.02E-04 ,4.24E-04 ,2.64E-04 ,2.06E-04 ,1.92E-04 ,1.91E-04 ,1.95E-04 ,2.03E-04 ,2.15E-04 ,2.29E-04 ,2.47E-04 ,2.69E-04 ,2.95E-04 ,3.26E-04 ,3.63E-04 ,4.05E-04 ,4.56E-04])
        plt.plot(x10, v10, color='firebrick', linewidth=2, marker='s', label='10')
        plt.plot(x20, v20, color='seagreen',  linewidth=2, marker='p', label='20')
        plt.plot(x50, v50, color='royalblue', linewidth=2, marker='d', label='50')
        plt.plot(x100,v100, color='red',      linewidth=2, marker='v', label='100')
        plt.plot(x125,v125, color='black',    linewidth=2, marker='h', label='125')

        plt.legend(loc='lower right', title='time slices', frameon=False, columnspacing=None)
        plt.title(r'CN-EU : ${\delta}T/{\delta}t=100$, ${\delta}t=1.0e-6s$')

    elif( name == 'bd4-eu-5'):
        x10  = np.arange(10)+1  
        x20  = np.arange(20)+1  
        x50  = np.arange(17)+1  
        x100 = np.arange(11)+1  
        x125 = np.arange(13)+1  

        v10 = np.array([6.51E-01 ,4.37E-03 ,1.43E-05 ,1.67E-08 ,8.63E-11 ,1.04E-10 ,7.67E-11 ,5.88E-11 ,5.89E-11 ,6.54E-13])
        v20 = np.array([1.44E+00 ,1.09E-02 ,1.87E-04 ,1.57E-06 ,1.43E-08 ,1.46E-10 ,1.03E-10 ,5.45E-11 ,5.04E-11 ,1.79E-10 ,1.45E-10 ,1.34E-10 ,1.84E-10 ,1.23E-10 ,1.09E-10 ,2.01E-11 ,3.98E-11 ,2.07E-12 ,2.38E-13 ,2.17E-13])
        v50 = np.array([ 4.23E+00 ,6.32E-02 ,3.58E-03 ,1.46E-04 ,6.17E-06 ,2.53E-07 ,1.01E-08 ,3.71E-10 ,9.33E-11 ,9.31E-11 ,8.76E-11 ,8.99E-11 ,9.37E-11 ,1.04E-10 ,2.63E-10 ,1.79E-10 ,1.19E-10])


        v100 = np.array([ 9.92E+00 ,2.69E-01 ,2.99E-02 ,2.71E-03 ,2.58E-04 ,2.40E-05 ,2.19E-06 ,1.95E-07 ,1.71E-08 ,1.52E-09 ,2.79E-10])
        v125 = np.array([ 1.37E+01 ,4.53E-01 ,6.32E-02 ,7.40E-03 ,9.08E-04 ,1.09E-04 ,1.29E-05 ,1.49E-06 ,1.69E-07 ,1.89E-08 ,2.13E-09 ,3.09E-10 ,8.31E-11])
        
        plt.plot(x10, v10, color='firebrick',  linewidth=2, marker='s', label='10')
        plt.plot(x20, v20, color='seagreen',   linewidth=2, marker='p', label='20')
        plt.plot(x50, v50, color='blue',       linewidth=2, marker='d', label='50')
        plt.plot(x100,v100, color='red',     linewidth=2, marker='v', label='100')
        plt.plot(x125,v125, color='black',   linewidth=2, marker='h', label='125')
        plt.legend(loc='lower right', title='time slices', frameon=False, columnspacing=None)
        plt.title(r'BD4-EU : ${\delta}T/{\delta}t=100$, ${\delta}t=1.0e-6s$')

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
