#! /usr/bin/env python

## Required Software:
# Python: https://www.python.org/downloads/
# NumPy and SciPy: http://www.scipy.org/scipylib/download.html
# matplotlib: http://matplotlib.org/downloads.html

import math
import scipy
import scipy.stats as stats
import numpy as np
import sys
import matplotlib.pyplot as plt
import csv
from argparse import ArgumentParser
from scipy.optimize import curve_fit
## from scipy.optimize import *

print ("============================================================")
print ("PyReweighting-3D: Python script used to reweight accelerated molecular dynamics simulations")
print ("                  and calculate 3D potential of mean force (PMF) profiles and ligand binding free energies.")
print ("  ")
print ("Author: Yinglong Miao <miao@ku.edu>, Copyright 2019-2020.")
print ("\n\
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"PyReweighting\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following citation: \n\
\n\
Miao Y, Sinko W, Pierce L, Bucher D, Walker RC, McCammon JA (2014) Improved reweighting of accelerated molecular dynamics simulations for free energy calculation. J Chemical Theory and Computation. 10(7): 2677-2689.")
print (" ")

###########MAIN
def main():
## Set control parameters
    plt_figs = 0
    lig_dG = 1

    args = cmdlineparse()   
    
    data=loadfiletoarray(args.input)
    
    rows = len(data[:,0])
    weights,dV = weightparse(rows, args)

    if args.discX:
        discX=float(args.discX)
    else :
        discX = 6

    if args.discY:
        discY=float(args.discY)
    else :
        discY = 6

    if args.discZ:
        discZ=float(args.discZ)
    else :
        discZ = 6

    Vcell = discX * discY * discZ

    if args.Xdim:
        binsX= assignbins(args.Xdim, discX)
    else:
        max_data = discX * (int(np.amax(data[:,0])/discX) + 1)
        min_data = discX * (int(np.amin(data[:,0])/discX) - 1)
        binsX= assignbins([min_data,max_data], discX)  ## Default bin size

    if args.Ydim:
        binsY= assignbins(args.Ydim, args)
    else:
        max_data = discY * (int(np.amax(data[:,1])/discY) + 1)
        min_data = discY * (int(np.amin(data[:,1])/discY) - 1)
        binsY= assignbins([min_data,max_data], discY)  ## Default bin size

    if args.Zdim:
        binsZ= assignbins(args.Zdim, args)
    else:
        max_data = discZ * (int(np.amax(data[:,2])/discZ) + 1)
        min_data = discZ * (int(np.amin(data[:,2])/discZ) - 1)
        binsZ= assignbins([min_data,max_data], discZ)  ## Default bin size

##  SET MAX ENERGY FOR ALL INFINITY VALUES
    if args.Emax:
        cb_max=float(args.Emax)
    else :
        cb_max = 8

##  SET HISTOGRAM CUTOFF
    if args.cutoff:
        hist_min=int(args.cutoff)
    else :
        hist_min = 10    # minimum number of configurations in one bin

##  SET ORDER of McLaurin series expansion
    if args.order:
        order=int(args.order)
    else :
        order = 10    # default

##  SET TEMPERATURE
    if args.T:
        T=float(args.T)
    else :
        T = 300    # simulation temperature
    beta = 1.0/(0.001987*T)

## SET flag for calculating ligand binding free energy
    if args.lig_dG:
        lig_dG=args.lig_dG
    else :
        lig_dG=False

##  SET LIGAND BOUND STATE CUTOFF DISTANCE
    if args.rb:
        rb=float(args.rb)
    else :
        rb = 7.5

##  SET LIGAND UNBOUND STATE CUTOFF DISTANCE
    if args.ru:
        ru=float(args.ru)
    else :
        ru = 7.5

##REWEIGHTING
##  SET flag for Gaussian fitting of deltaV
    if args.fit:
        fit=args.fit
    else :
        fit=False    # simulation temperature
##    print "gaussian fitting:", fit

##REWEIGHTING
    if args.job == "amdweight_CE":
        hist3,newedgesX,newedgesY,newedgesZ,c1,c2,c3 = reweight_CE(data,hist_min,binsX,discX,binsY,discY,binsZ,discZ,dV,T,fit)
        pmf = hist3pmf3D(hist3,hist_min,T)
        c1 = -np.multiply(1.0/beta,c1)
        c2 = -np.multiply(1.0/beta,c2)
        c3 = -np.multiply(1.0/beta,c3)
        
        c12 = np.add(c1,c2)
        c123 = np.add(c12,c3)
        pmf_c1 = np.add(pmf, c1)
        print ("pmf_min-c1 = ", np.min(pmf_c1))
        pmf_c1 = normalize3D(pmf_c1,cb_max)
        pmf_c2 = np.add(pmf, c12)
        print ("pmf_min-c2 = ", np.min(pmf_c2))
        pmf_c2 = normalize3D(pmf_c2,cb_max)
        pmf_c3 = np.add(pmf, c123)
        print ("pmf_min-c3 = ", np.min(pmf_c3))
        pmf_c3 = normalize3D(pmf_c3,cb_max)
    elif args.job == "amdweight_MC":
        n=order
        MCweight=np.zeros(len(dV))
        beta_dV=np.multiply(dV,beta)
        for x in range(0,n+1):
          MCweight=np.add(MCweight,(np.divide(np.power(beta_dV, x), float(scipy.misc.factorial(x)))))
        weights=MCweight
        hist3,(newedgesX,newedgesY,newedgesZ) = np.histogramdd(data, bins = (binsX, binsY, binsZ), weights=weights)
        hist3=prephist(hist3,T,cb_max)
    elif args.job == "amdweight":
        hist3,(newedgesX,newedgesY,newedgesZ) = np.histogramdd(data, bins = (binsX, binsY, binsZ), weights=weights)
        hist3=prephist(hist3,T,cb_max)
    else :
        hist3,(newedgesX,newedgesY,newedgesZ) = np.histogramdd(data, bins = (binsX, binsY, binsZ), weights=None)
        hist3=prephist(hist3,T,cb_max)

##SAVE FREE ENERGY DATA INTO A FILE
    if args.job == "amdweight_MC" or args.job == "amdweight" or args.job == "noweight" :
        pmffile = 'pmf-'+str(args.input)+'.xvg'
        output_pmf3D(pmffile,hist3,binsX,binsY,binsZ)
    if args.job == "amdweight_CE" :
        hist3 = pmf_c1
        pmffile = 'pmf-c1-'+str(args.input)+'.xvg'
        output_pmf3D(pmffile,hist3,binsX,binsY,binsZ)
        
        hist3 = pmf_c3
        pmffile = 'pmf-c3-'+str(args.input)+'.xvg'
        output_pmf3D(pmffile,hist3,binsX,binsY,binsZ)
        
        hist3 = pmf_c2
        pmffile = 'pmf-c2-'+str(args.input)+'.xvg'
        output_pmf3D(pmffile,hist3,binsX,binsY,binsZ)

    if args.job == "histo" :
        hist3,(newedgesX,newedgesY,newedgesZ) = np.histogramdd(data, bins = (binsX, binsY, binsZ), weights=None)
        pmffile = 'histo-'+str(args.input)+'.xvg'
        output_dV_anharm3D(pmffile,binsX,binsY,binsZ,hist3)

    if args.job == "amd_dV":
        plt_figs = 0
        hist3,newedgesX,newedgesY,newedgesZ,binfX,binfY,binfZ,dV_avg,dV_std,dV_anharm,dV_mat = reweight_dV(data,hist_min,binsX,discX,binsY,discY,binsZ,discZ,dV,T)
        
        pmffile = 'dV-hist-3D-'+str(args.input) + '.xvg'
        output_dV(pmffile,dV)
        
        alpha = anharm(dV)
        print ("Anharmonicity of all dV = " + str(alpha))
        
        pmffile = 'dV-anharm-3D-'+str(args.input)+'.xvg'
        output_dV_anharm3D(pmffile,binsX,binsY,binsZ,dV_anharm)
        
        pmffile = 'dV-stat-3D-'+str(args.input)+'.xvg'
        output_dV_stat3D(pmffile,binsX,binsY,binsZ,dV_avg,dV_std,dV_anharm)
        
        pmffile = 'dV-mat-3D-'+str(args.input)+'.xvg'
        output_dV_mat3D(pmffile,binsX,binsY,binsZ,hist3,dV_avg,dV_std,dV_anharm,dV_mat)

    if lig_dG :
        V0 = 1661
        Vb, Vb0 = calc_Vb(hist3,binsX,discX,binsY,discY,binsZ,discZ,rb,T)
        dW, Vu, Vu0 = calc_dW(hist3,binsX,discX,binsY,discY,binsZ,discZ,ru,T)
        dG = - (0.001987*T)*np.log(Vb/V0) - dW 
        print ("\nLigand binding free energy calculation: ")
        print ("Input (rb, ru)   = %8.3f, %8.3f" % (rb, ru))
        print ("dG               = %8.3f kcal/mol with %s reweighting" % (dG, args.job))
        print ("(Vb, Vb0)        = %8.3f, %8.3f" % (Vb, Vb0))
        print ("(dW, Vu, Vu0)    = %8.3f, %8.3f, %8.3f" % (dW, Vu, Vu0))

###PLOTTING FUNCTION FOR FREE ENERGY FIGURE
    if plt_figs :
        cbar_ticks=[0, cb_max*.25, cb_max*.5, cb_max*.75, 8.0]
        plt.figure(2, figsize=(11,8.5))
        extent = [newedgesX[0], newedgesX[-1], newedgesY[-1], newedgesY[0], newedgesZ[-1], newedgesZ[0]]
        print (extent)
        plt.imshow(hist3.transpose(), extent=extent, interpolation='gaussian')
        cb = plt.colorbar(ticks=cbar_ticks, format=('% .1f'), aspect=10) # grab the Colorbar instance
        imaxes = plt.gca()
        plt.sca(cb.ax)
        plt.clim(vmin=0,vmax=8.0)
        plt.yticks(fontsize=18)
        plt.sca(imaxes)
        axis=(min(binsX), max(binsX), min(binsY), max(binsY), min(binsZ), max(binsZ))
        plt.axis(axis)
        plt.xticks(size='18')
        plt.yticks(size='18')
        plt.xlabel('RC1',fontsize=18)
        plt.ylabel('RC2',fontsize=18)
##        plt.xlabel(r'$\phi$',fontsize=18)
##        plt.ylabel(r'$\psi$',fontsize=18)
##        plt.xlabel(r'$\chi$1',fontsize=18)
##        plt.ylabel(r'$\chi$2',fontsize=18)
        plt.savefig('3D_Free_energy_surface.png',bbox_inches=0)
        print ("FIGURE SAVED 3D_Free_energy_surface.png")
    
###PLOTTING FUNCTION FOR WEIGHTS histogram
        [hist, edges] = np.histogram(weights, bins=100)
        width=np.absolute(np.subtract(edges[0], edges[1]))
        plt.figure(1, figsize=(11,8.5))
        plt.bar(edges[:100], hist, width=width, log=True)
        plt.yscale('log')   ###if typerror is thrown delete .matplotlib/fontList.cache  file
        plt.xticks(fontsize='18')
        plt.yticks(fontsize='18')
        plt.savefig('weights.png',bbox_inches=0)
        print ("FIGURE SAVED weights.png")

    print (" ")
    print ("END")

def cmdlineparse():
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument("-input", dest="input", required=True, help="3D input file", metavar="<3D input file>")
    parser.add_argument("-job", dest="job", required=True, help="Reweighting method to use: <noweight>, <weighthist>, <amd_time>, <amd_dV>, <amdweight>, <amdweight_MC>, <amdweight_CE>, <histo>", metavar="<Job type reweighting method>")
    parser.add_argument("-weight", dest="weight", required=False, help="weight file", metavar="<weight file>")
    parser.add_argument("-Xdim", dest="Xdim", required=False, nargs="+", help="Xdimensions", metavar="<Xmin Xmax >")
    parser.add_argument("-Ydim", dest="Ydim", required=False, nargs="+", help="Ydimension", metavar="<Ymin Ymax >")
    parser.add_argument("-Zdim", dest="Zdim", required=False, nargs="+", help="Zdimension", metavar="<Zmin Zmax >")
    parser.add_argument("-discX", dest="discX", required=False,  help="Discretization size in X dimension", metavar="<discretization-X>")
    parser.add_argument("-discY", dest="discY", required=False,  help="Discretization size in Y dimension", metavar="<discretization-Y>")
    parser.add_argument("-discZ", dest="discZ", required=False,  help="Discretization size in Z dimension", metavar="<discretization-Z>")
    parser.add_argument("-cutoff", dest="cutoff", required=False,  help="histogram cutoff", metavar="<cutoff>")
    parser.add_argument("-T", dest="T", required=False,  help="Temperature", metavar="<Temperature>")
    parser.add_argument("-Emax", dest="Emax", required=False,  help="Maximum free energy", metavar="<Emax>")
    parser.add_argument("-fit", dest="fit", required=False, help="Fit deltaV distribution", metavar="<fit>")
    parser.add_argument("-order", dest="order", required=False, help="Order of Maclaurin series", metavar="<order>")
    parser.add_argument("-lig_dG", dest="lig_dG", required=False,  help="Flag for calculating ligand binding free energy", metavar="<lig_dG>")
    parser.add_argument("-rb", dest="rb", required=False,  help="Ligand bound cutoff distance", metavar="<rb>")
    parser.add_argument("-ru", dest="ru", required=False,  help="Ligand unbound cutoff distance", metavar="<ru>")
    args=parser.parse_args()
    return args
    
def loadfiletoarray(file):
    loaded=np.loadtxt(file, usecols=[0,1,2])
    print ("DATA LOADED:    "+file)
    return loaded

def weightparse(rows, args):
    if args.job == "weighthist":
        data=np.loadtxt(args.weight)
        weights=data[:,0]
        dV = np.zeros(rows)
    elif args.job == "amd_time" or args.job == "amd_dV" or args.job == "amdweight" or args.job == "amdweight_MC" or args.job == "amdweight_CE" :
        data=np.loadtxt(args.weight)
        weights = np.exp(data[:,0])
        dV = data[:,2]
    elif args.job == "noweight" or args.job == "histo":
        weights = np.zeros(rows)
        weights = weights + 1
        dV = np.zeros(rows)
    else:
        print ("ERROR: JOBTYPE"+ args.job+ " NOT RECOGNIZED")
        del data
        del weights
        del dV
    return weights,dV

def assignbins(dim, disc):
    minimum=float(dim[0])
    maximum=float(dim[1])
    bins =np.arange(minimum,(maximum+disc),disc)
    return bins

def normalize3D(pmf,cb_max):
    pmf=pmf-np.min(pmf)  ## zero value to lowest energy state
    temphist=pmf
    #set infinity free energy values to is cb_max
    for jz in range(len(temphist[0,0,:])):
      for jy in range(len(temphist[0,:,0])):
        for jx in range(len(temphist[:,0,0])):
          if np.isinf(temphist[jx,jy,jz]):
                temphist[jx,jy,jz]=cb_max
    return temphist

def prephist(hist3,T,cb_max):
    hist3=np.add(hist3,0.000000000000000001)  ###so that distrib
    hist3=(0.001987*T)*np.log(hist3) ####Convert to free energy in Kcal/mol
    hist3=np.max(hist3)-hist3  ## zero value to lowest energy state
##    print np.max(hist3)
    temphist3=hist3
    #set infinity free energy values to is cb_max
    for jz in range(len(temphist3[0,0,:])):
      for jy in range(len(temphist3[0,:,0])):
        for jx in range(len(temphist3[:,0,0])):
            if np.isinf(temphist3[jx,jy,jz]):
                temphist3[jx,jy,jz]=cb_max
    return temphist3

# memory usage is much reduced with multidimensional list for dV_mat; pretty fast ~ O(N)
def reweight_CE(data,hist_min,binsX,discX,binsY,discY,binsZ,discZ,dV,T,fit):
    hist3, (newedgesX, newedgesY, newedgesZ) = np.histogramdd(data, bins = (binsX, binsY, binsZ), weights=None)

    beta = 1.0/(0.001987*T)
    nf = len(data[:,0])
    nbinsX = len(hist3[:,0,0])
    nbinsY = len(hist3[0,:,0])
    nbinsZ = len(hist3[0,0,:])

    c1 = np.zeros((nbinsX,nbinsY,nbinsZ)) 
    c2 = np.zeros((nbinsX,nbinsY,nbinsZ)) 
    c3 = np.zeros((nbinsX,nbinsY,nbinsZ)) 

    binfX = np.zeros(nf) # array for storing assigned bin of each frame
    binfY = np.zeros(nf) # array for storing assigned bin of each frame
    binfZ = np.zeros(nf) # array for storing assigned bin of each frame
    nA = np.zeros((nbinsX,nbinsY,nbinsZ),dtype=np.int) # nA is equivalent to hist here
    dV_avg = np.zeros((nbinsX,nbinsY,nbinsZ)) 
    dV_avg2 = np.zeros((nbinsX,nbinsY,nbinsZ)) 
    dV_avg3 = np.zeros((nbinsX,nbinsY,nbinsZ)) 
    dV_std = np.zeros((nbinsX,nbinsY,nbinsZ)) 

    dV_avg_all=np.average(dV)
    dV_std_all=np.std(dV)
    print ('dV all: avg = ', dV_avg_all, 'std = ', dV_std_all)

    dV_mat = [[[[[] for i in range(1)] for i in range(nbinsZ)] for i in range(nbinsY)] for i in range(nbinsX)]
    for i in range(len(data[:,0])):
        jx = int((data[i,0]-binsX[0])/discX)
        jy = int((data[i,1]-binsY[0])/discY)
        jz = int((data[i,2]-binsZ[0])/discZ)
        if jx < nbinsX and jy < nbinsY and jz < nbinsZ :
            binfX[i] = jx
            binfY[i] = jy
            binfZ[i] = jz
            dV_mat[jx][jy][jz].append(dV[i])
            nA[jx,jy,jz] = nA[jx,jy,jz]+1

    for jx in range(nbinsX):
      for jy in range(nbinsY):
        for jz in range(nbinsZ):
          if nA[jx,jy,jz]>=hist_min :
            num = int(nA[jx,jy,jz])
            atemp = np.asarray(dV_mat[jx][jy][jz][1:num+1])
            atemp2 = np.power(atemp,2)
            atemp3 = np.power(atemp,3)
            dV_avg[jx,jy,jz] = np.average(atemp)
            dV_std[jx,jy,jz]=np.std(atemp)
            dV_avg2[jx,jy,jz]=np.average(atemp2)
            dV_avg3[jx,jy,jz]=np.average(atemp3)
            del atemp
            del atemp2
            del atemp3
            c1[jx,jy,jz] = beta*dV_avg[jx,jy,jz]
            c2[jx,jy,jz] = 0.5*beta**2*dV_std[jx,jy,jz]**2
            c3[jx,jy,jz] = (1.0/6.0)*beta**3*(dV_avg3[jx,jy,jz]-3.0*dV_avg2[jx,jy,jz]*dV_avg[jx,jy,jz]+2.0*dV_avg[jx,jy,jz]**3)

    del dV_mat
    del dV_avg 
    del dV_avg2
    del dV_avg3
    del dV_std 
    return hist3,newedgesX,newedgesY,newedgesZ,c1,c2,c3

def reweight_dV(data,hist_min,binsX,discX,binsY,discY,binsZ,discZ,dV,T):
    hist3, (newedgesX, newedgesY, newedgesZ) = np.histogramdd(data, bins = (binsX, binsY, binsZ), weights=None)

    nf = len(data[:,0])
    nbinsX = len(hist3[:,0,0])
    nbinsY = len(hist3[0,:,0])
    nbinsZ = len(hist3[0,0,:])

    binfX = np.zeros(nf) # array for storing assigned bin of each frame
    binfY = np.zeros(nf) # array for storing assigned bin of each frame
    binfZ = np.zeros(nf) # array for storing assigned bin of each frame
    nA = np.zeros((nbinsX,nbinsY,nbinsZ),dtype=np.int) # nA is equivalent to hist here
    dV_avg = np.zeros((nbinsX,nbinsY,nbinsZ)) 
    dV_avg2 = np.zeros((nbinsX,nbinsY,nbinsZ)) 
    dV_avg3 = np.zeros((nbinsX,nbinsY,nbinsZ)) 
    dV_std = np.zeros((nbinsX,nbinsY,nbinsZ)) 
    dV_anharm = np.zeros((nbinsX,nbinsY,nbinsZ)) 

    dV_avg_all=np.average(dV)
    dV_std_all=np.std(dV)
    print ('dV all: avg = ', dV_avg_all, 'std = ', dV_std_all)

    dV_mat = [[[[[] for i in range(1)] for i in range(nbinsZ)] for i in range(nbinsY)] for i in range(nbinsX)]
    for i in range(len(data[:,0])):
        jx = int((data[i,0]-binsX[0])/discX)
        jy = int((data[i,1]-binsY[0])/discY)
        jz = int((data[i,2]-binsZ[0])/discZ)
        if jx < nbinsX and jy < nbinsY and jz < nbinsZ :
            binfX[i] = jx
            binfY[i] = jy
            binfZ[i] = jz
            dV_mat[jx][jy][jz].append(dV[i])
            nA[jx,jy,jz] = nA[jx,jy,jz]+1

    for jx in range(nbinsX):
      for jy in range(nbinsY):
        for jz in range(nbinsZ):
          if nA[jx,jy,jz]>=hist_min :
            num = int(nA[jx,jy,jz])
            atemp = np.asarray(dV_mat[jx][jy][jz][1:num+1])
            dV_avg[jx,jy,jz] = np.average(atemp)
            dV_std[jx,jy,jz] = np.std(atemp)
            dV_anharm[jx,jy,jz] = anharm(atemp)
            del atemp
    return hist3,newedgesX,newedgesY,newedgesZ,binfX,binfY,binfZ,dV_avg,dV_std,dV_anharm,dV_mat

##  Convert histogram to free energy in Kcal/mol
def hist3pmf3D(hist,hist_min,T):
        nbinsX = len(hist[:,0,0])
        nbinsY = len(hist[0,:,0])
        nbinsZ = len(hist[0,0,:])
        pmf = np.zeros((nbinsX,nbinsY,nbinsZ))
        pmf_min = 100
        for jx in range(len(hist[:,0,0])):
          for jy in range(len(hist[0,:,0])):
            for jz in range(len(hist[0,0,:])):
                if hist[jx,jy,jz]>=hist_min :
                    pmf[jx,jy,jz]=-(0.001987*T)*np.log(hist[jx,jy,jz])
                    if pmf_min > pmf[jx,jy,jz] :
                       pmf_min=pmf[jx,jy,jz]
##        pmf=pmf-pmf_min  ## zero value to lowest energy state
        return pmf

def output_pmf3D(pmffile,hist,binsX,binsY,binsZ):
        fpmf = open(pmffile, 'w')
        strpmf='#RC1\tRC2\tRC3\tPMF(kcal/mol)\n\n@    xaxis  label \"RC1\"\n@    yaxis  label \"RC2\"\n@    zaxis  label \"RC3\"\n@TYPE xy\n'
        fpmf.write(strpmf)
        for jx in range(len(hist[:,0,0])):
          for jy in range(len(hist[0,:,0])):
            for jz in range(len(hist[0,0,:])):
                strpmf=str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(binsZ[jz]) + ' \t' + str(hist[jx,jy,jz]) + '\n'
                fpmf.write(strpmf)
        fpmf.closed
        return fpmf

def output_pmf3D_dx(pmffile,hist,binsX,binsY,binsZ):
        fpmf = open(pmffile, 'w')
        strpmf='#RC1\tRC2\tRC3\tPMF(kcal/mol)\n\n@    xaxis  label \"RC1\"\n@    yaxis  label \"RC2\"\n@    zaxis  label \"RC3\"\n@TYPE xy\n'
        fpmf.write(strpmf)
        for jx in range(len(hist[:,0,0])):
          for jy in range(len(hist[0,:,0])):
            for jz in range(len(hist[0,0,:])):
                strpmf=str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(binsZ[jz]) + ' \t' + str(hist[jx,jy,jz]) + '\n'
                fpmf.write(strpmf)
        fpmf.closed
        return fpmf

def calc_dW(pmf,binsX,discX,binsY,discY,binsZ,discZ,ru,T):
        beta = 1.0/(0.001987*T)
        Vu=0.0
        Vu0=0.0
        ncells=0
        for jx in range(len(pmf[:,0,0])):
          for jy in range(len(pmf[0,:,0])):
            for jz in range(len(pmf[0,0,:])):
                r2=binsX[jx]**2 + binsY[jy]**2 + binsZ[jz]**2
                if r2 > ru**2 :
                    ncells = ncells + 1
                    Vu = Vu + np.exp(- beta * pmf[jx,jy,jz])
        Vu0 = ncells * discX * discY * discZ
        Vu = Vu * discX * discY * discZ
        dW = - 0.001987*T*np.log(Vu/Vu0)
        return dW,Vu,Vu0

def calc_Vb(pmf,binsX,discX,binsY,discY,binsZ,discZ,rb,T):
        beta = 1.0/(0.001987*T)
        Vb=0.0
        Vb0=0.0
        ncells=0
        for jx in range(len(pmf[:,0,0])):
          for jy in range(len(pmf[0,:,0])):
            for jz in range(len(pmf[0,0,:])):
                r2=binsX[jx]**2 + binsY[jy]**2 + binsZ[jz]**2
                if r2 <= rb**2 :
                    ncells = ncells + 1
                    Vb = Vb + np.exp(- beta * pmf[jx,jy,jz])
        Vb0 = ncells * discX * discY * discZ
        Vb = Vb * discX * discY * discZ
        return Vb,Vb0

def output_dV(pmffile,dV):
        fpmf = open(pmffile, 'w')
        strpmf='#dV \tp(dV) \n\n@    xaxis  label \"dV\"\n@    yaxis  label \"p(dV)\"\n@TYPE xy\n'
        hist_dV, bin_dV = np.histogram(dV, bins=50)
        for k in range(len(hist_dV)):
            strpmf=strpmf + str(bin_dV[k]) + ' \t' + str(hist_dV[k]) + ' \n'
        fpmf.write(strpmf)
        fpmf.closed
        return fpmf

def output_dV_anharm3D(pmffile,binsX,binsY,binsZ,dV_anharm):
        fpmf = open(pmffile, 'w')
        strpmf='#RC \tdV_anharm \tError\n\n@    xaxis  label \"RC\"\n@    yaxis  label \"dV_anmarm\"\n@TYPE xy\n'
        fpmf.write(strpmf)
        for jx in range(len(dV_anharm[:,0,0])):
          for jy in range(len(dV_anharm[0,:,0])):
            for jz in range(len(dV_anharm[0,0,:])):
                strpmf=str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(binsZ[jz]) + ' \t' + str(dV_anharm[jx,jy,jz]) + '\n'
                fpmf.write(strpmf)
        fpmf.closed
        return fpmf

def output_dV_stat3D(pmffile,binsX,binsY,binsZ,dV_avg,dV_std,dV_anharm):
        fpmf = open(pmffile, 'w')
        strpmf='#RC \tdV_avg(kcal/mol) \tError\n\n@    xaxis  label \"RC\"\n@    yaxis  label \"dV(kcal/mol)\"\n@TYPE xydy\n'
        fpmf.write(strpmf)
        for jx in range(len(dV_anharm[:,0,0])):
          for jy in range(len(dV_anharm[0,:,0])):
            for jz in range(len(dV_anharm[0,0,:])):
                strpmf=str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(binsZ[jz]) + ' \t' + str(dV_avg[jx,jy,jz]) + ' \t' + str(dV_std[jx,jy,jz]) + ' \t' + str(dV_anharm[jx,jy,jz]) + '\n'
                fpmf.write(strpmf)
        fpmf.closed
        return fpmf

def output_dV_mat3D(pmffile,binsX,binsY,binsZ,hist,dV_avg,dV_std,dV_anharm,dV_mat):
        fpmf = open(pmffile, 'w')
        strpmf='#RC \tNf \tdV_avg \tdV_std \tdV_ij \n\n@    xaxis  label \"RC\"\n@    yaxis  label \n@    zaxis  label \"dV(kcal/mol)\"\n@TYPE xy\n'
        fpmf.write(strpmf)
        for jx in range(len(hist[:,0,0])):
          for jy in range(len(hist[0,:,0])):
            for jz in range(len(hist[0,0,:])):
                nf_j = int(hist[jx,jy,jz])
                strpmf=str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(binsZ[jz]) + ' \t' + str(hist[jx,jy,jz]) + ' \t' + str(dV_avg[jx,jy,jz]) + ' \t' + str(dV_std[jx,jy,jz]) + ' \t' + str(dV_anharm[jx,jy,jz])
                strpmf=strpmf + ' \t' + str(dV_mat[jx][jy][jz][1:nf_j+1])
                strpmf=strpmf + '\n'
                fpmf.write(strpmf)
        fpmf.closed
        return fpmf

def anharm(data):
#    print "Compute anharmonicity"
    var=np.var(data)
    hist, edges=np.histogram(data, 50, normed=True)
    hist=np.add(hist,0.000000000000000001)  ###so that distrib
    dx=edges[1]-edges[0]
    S1=-1*np.trapz(np.multiply(hist, np.log(hist)),dx=dx)
    S2=0.5*np.log(2.00*np.pi*np.exp(1.0)*var+0.000000000000000001)
    alpha=S2-S1
    if np.isinf(alpha):
       alpha = 100
    return alpha
 
if __name__ == '__main__':
    main()
    
