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
print ("PyReweighting: Python scripts used to reweight accelerated molecular dynamics simulations.")
print ("  ")
print ("Authors: Yinglong Miao <yinglong.miao@gmail.com>")
print ("         Bill Sinko <wsinko@gmail.com>")
print ("\n\
Copyright <2014-2019> <Yinglong Miao and William Sinko> \n\
\n\
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"PyReweighting\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following citation: \n\
\n\
Miao Y, Sinko W, Pierce L, Bucher D, Walker RC, McCammon JA (2014) Improved reweighting of accelerated molecular dynamics simulations for free energy calculation. J Chemical Theory and Computation. 10(7): 2677-2689.")
print (" ")

###########MAIN
def main():
## Set control parameters
    plt_figs = 0

    args = cmdlineparse()   
    
    data=loadfiletoarray(args.input)
    
    rows = len(data[:])
    weights,dV = weightparse(rows, args)

    if args.disc:
        discX=float(args.disc)
    else :
        discX = 6

    if args.Xdim:
        binsX= assignbins(args.Xdim, discX)
    else:
        max_data = discX * (int(np.amax(data)/discX) + 1)
        min_data = discX * (int(np.amin(data)/discX) - 1)
        binsX= assignbins([min_data,max_data], discX)  ## Default bin size

##  SET MAX ENERGY FOR ALL INFINITY VALUES
    if args.Emax:
        cb_max=float(args.Emax)
    else :
        cb_max = 8

##  SET HISTOGRAM CUTOFF
    if args.cutoff:
        hist_min=int(args.cutoff)
    else :
        hist_min = 10	# minimum number of configurations in one bin

##  SET ORDER of McLaurin series expansion
    if args.order:
        order=int(args.order)
    else :
        order = 10	# default

##  SET TEMPERATURE
    if args.T:
        T=float(args.T)
    else :
        T = 300	# simulation temperature
    beta = 1.0/(0.001987*T)

##  SET flag for Gaussian fitting of deltaV
    if args.fit:
        fit=args.fit
    else :
        fit=False	# simulation temperature

##REWEIGHTING
    if args.job == "amdweight_CE":
        ##CALCULATE effective acceleration factor
        accel = accel_amd(dV,T)
        print ("Effective acceleration factor with Gaussian approximation: ", accel)
        
        hist,newedgesX,c1,c2,c3 = reweight_CE(data,hist_min,binsX,discX,dV,T,fit)
        pmf=prephist(hist,T,cb_max)
        c1 = -np.multiply(1.0/beta,c1)
        c2 = -np.multiply(1.0/beta,c2)
        c3 = -np.multiply(1.0/beta,c3)
        c12 = np.add(c1,c2)
        c123 = np.add(c12,c3)
        pmf_c1 = np.add(pmf, c1)
        print ("pmf_min-c1 = ", np.min(pmf_c1))
        pmf_c1 = normalize(pmf_c1,cb_max)
        pmf_c2 = np.add(pmf, c12)
        print ("pmf_min-c2 = ", np.min(pmf_c2))
        pmf_c2 = normalize(pmf_c2,cb_max)
        pmf_c3 = np.add(pmf, c123)
        print ("pmf_min-c3 = ", np.min(pmf_c3))
        pmf_c3 = normalize(pmf_c3,cb_max)
    elif args.job == "amd_time" :
        ##CALCULATE effective acceleration factor
        accel = accel_amd(dV,T)
        print ("Effective acceleration factor with Gaussian approximation: ", accel)
        
        hist,newedgesX,c1,c2,c3 = reweight_CE(data,hist_min,binsX,discX,dV,T,fit)
        c1 = np.multiply(1.0/beta,c1)
        c2 = np.multiply(1.0/beta,c2)
        c3 = np.multiply(1.0/beta,c3)
        c12 = np.add(c1,c2)
        c123 = np.add(c12,c3)
        
        c1=prepdV(c1,cb_max)
        print ("accel_min-c1 = ", np.min(c1))
        c1 = c1 - np.min(c1)
        c12=prepdV(c12,cb_max)
        print ("accel_min-c2 = ", np.min(c12))
        c12 = c12 - np.min(c12)
        c123=prepdV(c123,cb_max)
        print ("accel_min-c3 = ", np.min(c123))
        c123 = c123 - np.min(c123)
        
        itb=10
        ite=len(c12)-5
        dV=c12[itb:ite]
        accel = accel_amd(dV,T)
        print ("Time-averaged acceleration factor with Gaussian approximation: ", accel)

    elif args.job == "amdweight_MC":
        n=order
        MCweight=np.zeros(len(dV))
        beta_dV=np.multiply(dV,beta)
        for x in range(0,n+1):
          MCweight=np.add(MCweight,(np.divide(np.power(beta_dV, x), float(scipy.misc.factorial(x)))))
        weights=MCweight
        hist, newedgesX = np.histogram(data, bins = binsX, weights=weights)
        hist=prephist(hist,T,cb_max)
    elif args.job == "amdweight":
        hist, newedgesX = np.histogram(data, bins = binsX, weights=weights)
        hist=prephist(hist,T,cb_max)
    else :
        hist, newedgesX = np.histogram(data, bins = binsX, weights=None)
        hist=prephist(hist,T,cb_max)

##SAVE FREE ENERGY DATA INTO A FILE
    if args.job == "amdweight_MC" or args.job == "amdweight" or args.job == "noweight" :
        pmffile = 'pmf-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX)

    if args.job == "amdweight_CE" :
        hist = pmf_c1
        pmffile = 'pmf-c1-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX)

        hist = pmf_c2
        pmffile = 'pmf-c2-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX)

        hist = pmf_c3
        pmffile = 'pmf-c3-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX)

    if args.job == "amd_time" :
        hist = c1
        pmffile = 'accel-c1-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX)
        
        hist = c12
        pmffile = 'accel-c2-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX)
        
        hist = c123
        pmffile = 'accel-c3-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX)

##SAVE WEIGHTS
    if args.job == "amdweight_MC" or args.job == "amdweight" :
        pmffile = 'weights-'+str(args.input)+'.xvg'
        output_pmf(pmffile,weights,data)
    if args.job == "amdweight_CE" :
        hist = np.exp(c1)
        pmffile = 'weights-c1-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX)

        hist = np.exp(c12)
        pmffile = 'weights-c2-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX)

        hist = np.exp(c123)
        pmffile = 'weights-c3-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX)

    if args.job == "histo" :
        hist,newedgesX = histo(data,hist_min,binsX,discX)
        pmffile = 'histo-'+str(args.input)+'.xvg'
        output_dV_anharm(pmffile,binsX,hist)

    if args.job == "amd_dV" :
        hist, newedgesX, binf, dV_avg, dV_std, dV_anharm, dV_mat = reweight_dV(data,hist_min,binsX,discX,dV,T)

        pmffile = 'dV-hist-'+str(args.input) + '.xvg'
        output_dV(pmffile,dV)

        dV_avg_all=np.average(dV)
        dV_std_all=np.std(dV)
        print ('dV all: avg = ', dV_avg_all, 'std = ', dV_std_all)
        
        alpha = anharm(dV)
        print ("Anharmonicity of all dV = " + str(alpha))
        
        pmffile = 'dV-anharm-'+str(args.input)+'.xvg'
        output_dV_anharm(pmffile,binsX,dV_anharm)
        
        pmffile = 'dV-stat-'+str(args.input)+'.xvg'
        output_dV_stat(pmffile,binsX,dV_avg,dV_std,dV_anharm)
        
        pmffile = 'dV-mat-'+str(args.input)+'.xvg'
        output_dV_mat(pmffile,binsX,hist,dV_avg,dV_std,dV_anharm,dV_mat)

###PLOTTING FUNCTION FOR WEIGHTS histogram
    if plt_figs :
        [hist, edges] = np.histogram(weights, bins=100)
        width=np.absolute(np.subtract(edges[0], edges[1]))
        plt.figure(1, figsize=(11,8.5))
        plt.bar(edges[:100], hist, width=width, log=True)
        plt.yscale('log')   ###if typerror is thrown delete .matplotlib/fontList.cache  file
        plt.xticks(fontsize='18')
        plt.yticks(fontsize='18')
        plt.savefig('weights.png',bbox_inches=0)
        print ("FIGURE SAVED weights.png")
        plt.show()

    print (" ")
    print ("END")
 
########READ datafiles and print weights

def cmdlineparse():
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument("-input", dest="input", required=True, help="input file", metavar="<input file>")
    parser.add_argument("-job", dest="job", required=True, help="Reweighting method to use: <noweight>, <weighthist>, <amd_time>, <amd_dV>, <histo>, <amdweight>, <amdweight_MC>, <amdweight_CE>", metavar="<Job type reweighting method>")
    parser.add_argument("-weight", dest="weight", required=False, help="weight file", metavar="<weight file>")
    parser.add_argument("-Xdim", dest="Xdim", required=False, nargs="+", help="Xdimensions", metavar="<Xmin Xmax>")
    parser.add_argument("-Ydim", dest="Ydim", required=False, nargs="+", help="Ydimension", metavar="<Ymin Ymax>")
    parser.add_argument("-disc", dest="disc", required=False,  help="Discretization size", metavar="<discretization>")
    parser.add_argument("-cutoff", dest="cutoff", required=False,  help="histogram cutoff", metavar="<cutoff>")
    parser.add_argument("-T", dest="T", required=False,  help="Temperature", metavar="<Temperature>")
    parser.add_argument("-Emax", dest="Emax", required=False,  help="Maximum free energy", metavar="<Emax>")
    parser.add_argument("-fit", dest="fit", required=False, help="Fit deltaV distribution", metavar="<fit>")
    parser.add_argument("-order", dest="order", required=False, help="Order of Maclaurin series", metavar="<order>")
    args=parser.parse_args()
    return args
    
def loadfiletoarray(file):
    loaded=np.loadtxt(file, usecols=[0])
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
        print ("ERROR JOBTYPE"+ args.job+ " NOT RECOGNIZED")
        del data
        del weights
        del dV
    return weights,dV

def reweight_CE(data,hist_min,binsX,discX,dV,T,fit):
    hist, newedgesX = np.histogram(data, bins = binsX, weights=None)

    beta = 1.0/(0.001987*T)
    nf = len(data)
    nbins = len(hist)

    c1 = np.zeros(nbins)
    c2 = np.zeros(nbins)
    c3 = np.zeros(nbins)

    binf = np.zeros(nf) # array for storing assigned bin of each frame
    nA = np.zeros(nbins,dtype=np.int) # nA is equivalent to hist here
    dV_avg = np.zeros(nbins)
    dV_avg2 = np.zeros(nbins)
    dV_avg3 = np.zeros(nbins)
    dV_std = np.zeros(nbins)

    dV_avg_all=np.average(dV)
    dV_std_all=np.std(dV)
    print ('dV all: avg = ', dV_avg_all, 'std = ', dV_std_all)

    diff_tol_avg = dV_std_all * 3
    diff_tol_std = dV_std_all
    dV_binsize = 50

    dV_mat = [[[] for i in range(1)] for i in range(nbins)]
    for i in range(len(data)):
        j = int((data[i]-binsX[0])/discX)
        if j == nbins :
          j = nbins-1
        binf[i] = j
        dV_mat[j].append(dV[i])
        nA[j] = nA[j]+1

    for j in range(nbins):
        if nA[j]>=hist_min :
          num = int(nA[j])
          atemp = np.asarray(dV_mat[j][1:num+1])
          atemp2 = np.power(atemp,2)
          atemp3 = np.power(atemp,3)
          if fit :
            ## calculate average/std through gaussian fitting
            hist_temp, bin_edges_temp = np.histogram(atemp, bins=dV_binsize)
            bin_centres_temp = (bin_edges_temp[:-1] + bin_edges_temp[1:])/2
            ## output original histograms
            pmffile = 'dV-hist-forFit-RC'+str('%#08.2f' % binsX[j]) + '.xvg'
            output_pmf(pmffile,hist_temp,bin_centres_temp)
            # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
            mean=np.average(atemp)
            std=np.std(atemp)
            p0 = [0., 1., mean, std]
            ## coeff, var_matrix = curve_fit(gauss, bin_centres_temp, hist_temp, p0=p0)
            coeff, var_matrix = curve_fit(gauss, bin_centres_temp, hist_temp, p0=p0)
            # Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
            print (binsX[j], ': mean = ', coeff[2], 'sigma = ', coeff[3])
            dV_avg[j]=coeff[2]
            dV_std[j]=coeff[3]
            # Get the fitted curve
            hist_fit = gauss(bin_centres_temp, *coeff)
            ## output fitted histograms
            pmffile = 'dV-hist-gaussFit-RC'+str('%#08.2f' % binsX[j]) + '.xvg'
            output_pmf(pmffile,hist_fit,bin_centres_temp)
          else :
            ## calculate average/std directly
            dV_avg[j]=np.average(atemp)
            dV_std[j]=np.std(atemp)

##	  if np.absolute(dV_avg[j]-dV_avg_all)>diff_tol_avg or np.absolute(dV_std[j]-dV_std_all)>diff_tol_std :
##	     dV_avg[j]=0
##	     dV_std[j]=0

          dV_avg2[j]=np.average(atemp2)
          dV_avg3[j]=np.average(atemp3)
          del atemp
          del atemp2
          del atemp3
          c1[j] = beta*dV_avg[j]
          c2[j] = 0.5*beta**2*dV_std[j]**2
          c3[j] = (1.0/6.0)*beta**3*(dV_avg3[j]-3.0*dV_avg2[j]*dV_avg[j]+2.0*dV_avg[j]**3)

    del dV_mat
    del dV_avg 
    del dV_avg2
    del dV_avg3
    del dV_std 
    return hist,newedgesX,c1,c2,c3

def reweight_dV(data,hist_min,binsX,discX,dV,T):
    hist, newedgesX = np.histogram(data, bins = binsX, weights=None)

    nf = len(data)
    nbins = len(hist)

    binf = np.zeros(nf) # array for storing assigned bin of each frame
    nA = np.zeros(nbins,dtype=np.int) # nA is equivalent to hist here
    dV_avg = np.zeros(nbins)
    dV_std = np.zeros(nbins)
    dV_anharm = np.zeros(nbins)

    dV_mat = [[[] for i in range(1)] for i in range(nbins)]
    for i in range(len(data)):
        j = int((data[i]-binsX[0])/discX)
        if j >= nbins :
          j = nbins-1
        binf[i] = j
        dV_mat[j].append(dV[i])
        nA[j] = nA[j]+1

    for j in range(nbins):
        dV_anharm[j] = 100
        if nA[j]>0 :
          num = int(nA[j])
          atemp = np.asarray(dV_mat[j][1:num+1])
          dV_avg[j] = np.average(atemp)
          dV_std[j] = np.std(atemp)
          dV_anharm[j] = anharm(atemp)
          del atemp
    return hist,newedgesX,binf,dV_avg,dV_std,dV_anharm,dV_mat

def histo(data,hist_min,binsX,discX):
    hist, newedgesX = np.histogram(data, bins = binsX, weights=None)
    return hist,newedgesX

def assignbins(dim, disc):
    minimum=float(dim[0])
    maximum=float(dim[1])
    bins =np.arange(minimum,(maximum+disc),disc)
    return bins

def normalize(pmf,cb_max):
    pmf=pmf-np.min(pmf)  ## zero value to lowest energy state
    temphist=pmf
    #set infinity free energy values to is cb_max
    for x in range(len(temphist[:])):
      if np.isinf(temphist[x]):
                temphist[x]=cb_max
    return temphist

def prephist(hist,T,cb_max):
    hist=np.add(hist,0.000000000000000001)  ###so that distrib
    hist=(0.001987*T)*np.log(hist) ####Convert to free energy in Kcal/mol
    print ("PMF_min = ", -np.max(hist))
    hist=np.max(hist)-hist  ## zero value to lowest energy state
    temphist=hist
    #set infinity free energy values to is cb_max
    for x in range(len(temphist[:])):
      if np.isinf(temphist[x]):
                temphist[x]=cb_max
    return temphist

def prepdV(hist,cb_max):
    for x in np.where(hist == 0)[0]:
                hist[x]=cb_max
    return hist

##  Convert histogram to free energy in Kcal/mol
def hist2pmf(hist,hist_min,T):
    nbins = len(hist)
    pmf = np.zeros(nbins)
    pmf_min = 100
    for j in range(len(hist)):
       if hist[j]>=hist_min :
          pmf[j]=-(0.001987*T)*np.log(hist[j])
          if pmf_min > pmf[j] :
             pmf_min=pmf[j]
    return pmf

def output_pmf(pmffile,hist,binsX):
    fpmf = open(pmffile, 'w')
    strpmf='#RC \tPMF(kcal/mol)\n\n@    xaxis  label \"RC\"\n@    yaxis  label \"PMF(kcal/mol)\"\n@TYPE xy\n'
    fpmf.write(strpmf)
    for j in range(len(hist[:])):
        strpmf=str(binsX[j]) + ' \t' + str(hist[j]) + '\n'
        fpmf.write(strpmf)
    fpmf.closed
    return fpmf

def output_dV(pmffile,dV):
    fpmf = open(pmffile, 'w')
    strpmf='#dV \tp(dV) \n\n@    xaxis  label \"dV\"\n@    yaxis  label \"p(dV)\"\n@TYPE xy\n'
    hist_dV, bin_dV = np.histogram(dV, bins=50)
    for k in range(len(hist_dV)):
        strpmf=strpmf + str(bin_dV[k]) + ' \t' + str(hist_dV[k]) + ' \n'
    fpmf.write(strpmf)
    fpmf.closed
    return fpmf

def output_dV_anharm(pmffile,binsX,dV_anharm):
    fpmf = open(pmffile, 'w')
    strpmf='#RC \tdV_anharm \tError\n\n@    xaxis  label \"RC\"\n@    yaxis  label \"dV_anmarm\"\n@TYPE xy\n'
    fpmf.write(strpmf)
    for j in range(len(dV_anharm[:])):
        strpmf=str(binsX[j]) + ' \t' + str(dV_anharm[j]) + '\n'
        fpmf.write(strpmf)
    fpmf.closed
    return fpmf

def output_dV_stat(pmffile,binsX,dV_avg,dV_std,dV_anharm):
    fpmf = open(pmffile, 'w')
    strpmf='#RC \tdV_avg(kcal/mol) \tError\n\n@    xaxis  label \"RC\"\n@    yaxis  label \"dV(kcal/mol)\"\n@TYPE xydy\n'
    fpmf.write(strpmf)
    for j in range(len(dV_avg[:])):
        strpmf=str(binsX[j]) + ' \t' + str(dV_avg[j]) + ' \t' + str(dV_std[j]) + ' \t' + str(dV_anharm[j]) + '\n'
        fpmf.write(strpmf)
    fpmf.closed
    return fpmf

def output_dV_mat(pmffile,binsX,hist,dV_avg,dV_std,dV_anharm,dV_mat):
    fpmf = open(pmffile, 'w')
    strpmf='#RC \tNf \tdV_avg \tdV_std \tdV_anharm \tdV_ij \n\n@    xaxis  label \"RC\"\n@    yaxis  label \"dV(kcal/mol)\"\n@TYPE xy\n'
    fpmf.write(strpmf)
    for j in range(len(dV_avg[:])):
        nf_j = int(hist[j])
        strpmf=str(binsX[j]) + ' \t' + str(hist[j]) + ' \t' + str(dV_avg[j]) + ' \t' + str(dV_std[j]) + ' \t' + str(dV_anharm[j])
        strpmf=strpmf + ' \t' + str(dV_mat[j][1:nf_j+1])
        strpmf=strpmf + '\n'
        fpmf.write(strpmf)
    fpmf.closed
    return fpmf

def accel_amd_window(c12,itb,ite):
    avg = np.average(dV)
    std = np.std(dV)
    accel = np.exp(c12)
    return accel

def accel_amd(dV,T):
    avg = np.average(dV)
    std = np.std(dV)
    accel = np.exp(avg/(0.001987*T)+(std/(0.001987*T))**2/2.)
    return accel

def gauss(x, *p):
    y0, A, mu, sigma = p
    return y0+A*np.exp(-(x-mu)**2/(2.*sigma**2))

def anharm(data):
    var=np.var(data)
    # hist, edges=np.histogram(data, 50, normed=True)
    hist, edges=np.histogram(data, 50, density=True)
    hist=np.add(hist,0.000000000000000001)  ###so that distrib
    dx=edges[1]-edges[0]
    S1=-1*np.trapz(np.multiply(hist, np.log(hist)),dx=dx)
    S2=0.5*np.log(np.add(2.00*np.pi*np.exp(1)*var,0.000000000000000001))
    alpha=S2-S1
    if np.isinf(alpha):
       alpha = 100
    return alpha

def anharmND(datafull):
    print ("Performing error estimation")
    width=datafull[0,:]
    lendata=len(datafull[:,0])
    for x in range(len(width)):
        var=np.var(datafull[:,x])
        std=np.std(datafull[:,x])
        print ("variance of "+str(x)+" is : " +str(var)+" Standard Deviation:  "+str(std))
        hist, edges=np.histogram(datafull[:,x], 100, normed=True)
        hist=np.add(hist,0.000000000000000001)  ###so that distrib
        dx=edges[1]-edges[0]
        S1=-1*np.trapz(np.multiply(hist, np.log(hist)),dx=dx)
        S2=0.5*np.log(np.add(2.00*np.pi*np.exp(1)*var,0.000000000000000001))
        alpha=S2-S1
        print (str(x)+"dimension   S1 = "+str(S1)+"  S2 = "+str(S2)+" Alpha = "+str(alpha))
    return var, std, alpha
  
if __name__ == '__main__':
    main()
    
