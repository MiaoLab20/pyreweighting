#!~/amber-git-install/miniconda/bin/python
#!/proj/ymiaolab/software/amber22//miniconda/bin/python

import os
import numpy as np
import argparse
from pathlib import Path
from scipy.special import factorial

print ("============================================================")
print ("PyReweighting-GaMD-TI: Python scripts used to reweight GaMD-TI simulations.")
print ("  ")
print ("Authors: Yinglong Miao <Yinglong_Miao@med.unc.edu>")
print ("         Jinan Wang <Jinan_Wang@med.unc.edu>")
print ("\n\
Copyright <2024-2026> <Yinglong Miao and Jinan Wang> \n\
\n\
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"PyReweighting-GaMD-TI\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so.")
print (" ")

def read_from_gamd_file(filename):
    """Reads dvdl and boost values from a given file."""
    time = []
    dvdl = []
    boost = []

    with open(filename, 'r') as file:
        for line in file:
            values = line.split()
            if len(values) >= 3:
                try:
                    time.append(float(values[0]))  # First column
                    dvdl.append(float(values[1]))  # Second column
                    boost.append(float(values[2]))  # Third column
                except ValueError:
                    continue  # Skip lines with invalid floats

    time = np.array(time, dtype=np.float64)
    dvdl = np.array(dvdl, dtype=np.float64)
    boost = np.array(boost, dtype=np.float64)

    return time, dvdl, boost

def extract_dvdl_mdout(filename):
    start_collecting = False
    stop_collecting = False
    dvdl = []

    with open(filename) as f:
        for line in f:
            if "End of dvdl summary" in line:
                stop_collecting = True

            if start_collecting and not stop_collecting:
                dvdl.append(float(line.strip()))

            if "Summary of dvdl values" in line:
                start_collecting = True

    dvdl = np.array(dvdl, dtype=np.float64)
    return dvdl

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

def calculate_average_noweight(dvdl):
    """Calculates the average value of dv/dl without taking boost into account (un-reweighted)."""

    average = np.sum(dvdl) / len(dvdl)
    se = np.std(dvdl)

    return average, se

def reweight_Jarzynski(dvdl, boost, b):
    """Reweight dvdl with Jarzynski's equality approximation."""
    dvdl = dvdl-boost

    # exp_cutoff = 709
    # for dV in dvdl :
        # if dV > exp_cutoff / b :
            # print(f"Warning: modified potential is out of range for calculating exponential: {dV}")

    dfdl = -np.log(np.sum(np.exp(-b*dvdl))/len(dvdl))/b
    dfdl_std = np.std(dvdl)
    # dfdl_std = np.sqrt(dvdl_std**2+boost_std**2)

    return dfdl, dfdl_std

def reweight_exp_average(dvdl, boost, b):
    """Reweight dvdl by average of boost potential values."""

    exp_cutoff = 709    # dV ~ 1200 kcal/mol
    for dV in boost :
        if dV > exp_cutoff / b :
            print(f"Warning: modified potential is out of range for calculating exponential: {dV}")

    exp_b_boost = np.exp(b * boost)
    weighted_exp = dvdl * exp_b_boost
    ave_exp_b_boost = np.sum(exp_b_boost) / len(boost)
    ave_weighted_exp = np.sum(weighted_exp) / len(weighted_exp)
    dfdl = ave_weighted_exp / ave_exp_b_boost

    dvdl_std = np.std(dvdl)
    boost_std = np.std(boost)
    dfdl_std = np.sqrt(dvdl_std**2+boost_std**2)

    # print(f"len_boost = \t{len(boost)}\tlen_weighted_exp = \t{len(weighted_exp)}\tave_weighted_exp = \t{ave_weighted_exp:.2f}\tave_exp_b_boost = {ave_exp_b_boost:.2f}")
    return dfdl, dfdl_std

def reweight_average(dvdl, boost, b):
    """Reweight dvdl by average of boost potential values OR Cumulant Expansion to the 1st order."""

    dvdl_avg = np.average(dvdl)
    dvdl_std = np.std(dvdl)
    boost_avg = np.average(boost)
    boost_std = np.std(boost)

    dfdl = dvdl_avg - boost_avg
    dfdl_std = np.sqrt(dvdl_std**2+boost_std**2)

    return dfdl, dfdl_std, dvdl_avg, dvdl_std, boost_avg, boost_std

def reweight_CE12(dvdl, boost, beta):
    """Reweight dvdl with boost potential using cumulant expansion on different orders."""

    boost_avg = np.average(boost)
    boost_std = np.std(boost)

    # cov_xy = np.mean((x - np.mean(x)) * (y - np.mean(y)))  [np calculate covariance of two vectors]
    cov_dvdl_boost = np.mean((dvdl - np.mean(dvdl)) * (boost - np.mean(boost)))

    ## Full cross-correlation (2nd order correlation function)
    # corr = np.correlate(x - np.mean(x), y - np.mean(y), mode='full')
    ## Normalize (so correlation at lag=0 is 1 if signals are identical)
    # corr = corr / (np.std(x) * np.std(y) * len(x))
    # corr_dvdl_boost = np.correlate(dvdl - np.mean(dvdl), boost - np.mean(boost), mode='full')
    corr_dvdl_boost_1_2 = third_order_corr_value(dvdl, boost, p=1, q=2, normalize=False)

    c1 = beta*cov_dvdl_boost
    # c2 = 0.5*beta**2*corr_dvdl_boost
    c2 = 0.5*beta**2*corr_dvdl_boost_1_2

    return c1, c2

def reweight_CE12_21(dvdl, boost, beta):
    """Reweight dvdl with boost potential using cumulant expansion on different orders."""

    boost_avg = np.average(boost)
    boost_std = np.std(boost)

    # cov_xy = np.mean((x - np.mean(x)) * (y - np.mean(y)))  [np calculate covariance of two vectors]
    cov_dvdl_boost = np.mean((dvdl - np.mean(dvdl)) * (boost - np.mean(boost)))

    ## Full cross-correlation (2nd order correlation function)
    # corr = np.correlate(x - np.mean(x), y - np.mean(y), mode='full')
    ## Normalize (so correlation at lag=0 is 1 if signals are identical)
    # corr = corr / (np.std(x) * np.std(y) * len(x))
    # corr_dvdl_boost = np.correlate(dvdl - np.mean(dvdl), boost - np.mean(boost), mode='full')
    corr_dvdl_boost_1_2 = third_order_corr_value(dvdl, boost, p=1, q=2, normalize=False)
    corr_dvdl_boost_2_1 = third_order_corr_value(dvdl, boost, p=2, q=1, normalize=False)

    c1 = beta*cov_dvdl_boost
    # c2 = 0.5*beta**2*corr_dvdl_boost
    c2 = 0.5*beta**2*(corr_dvdl_boost_1_2 + corr_dvdl_boost_2_1)

    return c1, c2

def reweight_CE(dvdl, boost, beta):
    """Reweight dvdl with boost potential using cumulant expansion on different orders."""

    boost_avg = np.average(boost)
    boost_std = np.std(boost)

    # cov_xy = np.mean((x - np.mean(x)) * (y - np.mean(y)))  [np calculate covariance of two vectors]
    cov_dvdl_boost = np.mean((dvdl - np.mean(dvdl)) * (boost - np.mean(boost)))

    ## Full cross-correlation (2nd order correlation function)
    # corr = np.correlate(x - np.mean(x), y - np.mean(y), mode='full')
    ## Normalize (so correlation at lag=0 is 1 if signals are identical)
    # corr = corr / (np.std(x) * np.std(y) * len(x))
    # corr_dvdl_boost = np.correlate(dvdl - np.mean(dvdl), boost - np.mean(boost), mode='full')
    corr_dvdl_boost_1_2 = third_order_corr_value(dvdl, boost, p=1, q=2, normalize=False)
    corr_dvdl_boost_2_1 = third_order_corr_value(dvdl, boost, p=2, q=1, normalize=False)

    c1 = beta*cov_dvdl_boost
    # c2 = 0.5*beta**2*corr_dvdl_boost
    c2 = 0.5*beta**2*(corr_dvdl_boost_1_2 + corr_dvdl_boost_2_1)

    return c1, c2

def third_order_corr_value(x, y, p=2, q=1, normalize=True):
    """
    Calculate the 3rd-order correlation value between two arrays.

    Parameters:
        x, y : array_like
            Input arrays of same length.
        p, q : int
            Exponents such that p+q=3 (default = 2,1).
        normalize : bool
            If True, normalize by std deviations raised to (p+q).

    Returns:
        float : 3rd-order correlation value.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    # Center the data
    x_c = x - np.mean(x)
    y_c = y - np.mean(y)

    # Compute raw moment
    val = np.mean((x_c**p) * (y_c**q))

    if normalize:
        val /= (np.std(x)**p * np.std(y)**q)

    return val

def integrate_dfdl_std(g, se, lambdas):
    # given lambdas = [0.0, 0.1, ... ,1.0]
    # and g = [<dH/dlambda>_lambda_i], se = [std error per window]
    # trapezoid integration
    dg = 0.0
    for i in range(len(lambdas)-1):
        dl = lambdas[i+1]-lambdas[i]
        dg += 0.5*(g[i]+g[i+1])*dl
    # for error propagation (independent approx)
    var = 0.0
    for i in range(len(lambdas)-1):
        dl = lambdas[i+1]-lambdas[i]
        var += 0.25*(se[i]**2 + se[i+1]**2)*(dl**2)
    se_total = np.sqrt(var)
    # print("ΔG =", dg, "+/-", se_total)
    return dg, se_total

def integrate_dfdl(g, lambdas):
    # given lambdas = [0.0, 0.1, ... ,1.0]
    # and g = [<dH/dlambda>_lambda_i]
    # trapezoid integration
    dg = 0.0
    for i in range(len(lambdas)-1):
        dl = lambdas[i+1]-lambdas[i]
        dg += 0.5*(g[i]+g[i+1])*dl
    return dg

def output_dV_histo(pmffile,binsX,dV_anharm):
    fpmf = open(pmffile, 'w')
    strpmf='#RC \tdV_histo \tError\n\n@    xaxis  label \"dV (kcal/mol)\"\n@    yaxis  label \"Probability\"\n@TYPE xy\n'
    fpmf.write(strpmf)
    for j in range(len(dV_anharm[:])):
        strpmf=str(binsX[j]) + ' \t' + str(dV_anharm[j]) + '\n'
        fpmf.write(strpmf)
    fpmf.closed
    return fpmf

def output_dG_mat(dGfile,dG_mat,time):
    fdG = open(dGfile, 'w')
    strdG='#Time \tdG_runs \n\n@    xaxis  label \"Time (ns)\"\n@    yaxis  label \"\\f{Symbol}DD\\f{}G (kcal/mol)\"\n@TYPE xy\n'
    fdG.write(strdG)
    for i in range(len(dG_mat[:,0])):
        strdG=str(time[i]) + ' \t'
        for j in range(len(dG_mat[0,:])):
            strdG=strdG + str(dG_mat[i,j]) + ' \t'
        strdG=strdG + '\n'
        fdG.write(strdG)
    fdG.closed
    return fdG

def output_dG_avg(dGfile,dG_avg,dG_std,time):
    fdG = open(dGfile, 'w')
    strdG='#Time \tdG_avg \tdG_std \n\n@    xaxis  label \"Time (ns)\"\n@    yaxis  label \"\\f{Symbol}DD\\f{}G (kcal/mol)\"\n@TYPE xydy\n'
    fdG.write(strdG)
    for j in range(len(dG_avg)):
        strdG=str(time[j]) + ' \t' + str(dG_avg[j]) + ' \t' + str(dG_std[j]) + ' \n'
        fdG.write(strdG)
    fdG.closed
    return fdG

def main(method):
    """Main function to calculate binding free energy using specified method."""
    lambdas = [0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738,
               0.56262, 0.68392, 0.79366, 0.88495, 0.95206, 0.99078]

##  SET INPUT FILENAME
    if args.input:
        filename = Path(args.input).stem
    else :
        filename = 'gamd-ti-energy'

##  SET TEMPERATURE
    if args.T:
        T=float(args.T)
    else :
        T = 300	# simulation temperature
    beta = 1.0 / (0.001987 * T)

##  SET NUMBER of SIMULATION RUNS
    if args.nruns:
        nruns=int(args.nruns)
    else :
        nruns = 3

##  SET NUMBER of SIMULATION LENGTHS USED TO CALCULATE dG
    if args.nlen:
        nlen=int(args.nlen)
    else :
        nlen = 20

    ltime = np.zeros(nlen)
    time = []
    dvdl = []
    boost = []

    # calculate statistics
    if method == 'stats':
        print(f"Statistics:")
        for l in range(1,nlen+1):
            print(f"\nSimulation Length: {l/nlen:.2f}")

            for run in range(1,nruns+1):
                for y in range(len(lambdas)):
                    gamd_file = f'run_{run}/{y}_lambda/{filename}.dat'
                    mdout_file = f'run_{run}/{y}_lambda/{filename}.out'
                    if os.path.exists(gamd_file) :
                        time, dvdl, boost = read_from_gamd_file(gamd_file)
                        if run==1 and y==0 :
                            print(f"Run\tWindow\tLength\tdvdl_avg\tstd\tmin\tmax\tboost_avg\tstd\tmin\tmax\tanharm")
                    elif os.path.exists(mdout_file) and not os.path.exists(gamd_file):
                        dvdl = extract_dvdl_mdout(mdout_file)
                        if run==1 and y==0 :
                            print(f"Run\tWindow\tLength\tdvdl_avg\tstd\tmin\tmax")
                    else :
                        print("Error with reading input file.")
                        return

                    # calculate simulation length
                    if run == 1 and y == 0 :
                        # tl=int(time[-1]*l/nlen)
                        tl=l/nlen
                        ltime[l-1]=tl

                    # dvdl
                    if len(dvdl)>0 :
                        dvdl=dvdl[0:int(len(dvdl)*l/nlen)]
                        dvdl_avg = np.average(dvdl)
                        dvdl_std = np.std(dvdl)
                        dvdl_min = np.min(dvdl)
                        dvdl_max = np.max(dvdl)
                        # print(f"{y}\t{dvdl_avg:.2f}\t\t{dvdl_std:.2f}\t{dvdl_min:.2f}\t{dvdl_max:.2f}")

                    # boost
                    if len(boost)>0 :
                        boost=boost[0:int(len(boost)*l/nlen)]
                        boost_avg = np.average(boost)
                        boost_std = np.std(boost)
                        boost_min = np.min(boost)
                        boost_max = np.max(boost)
                        boost_anharm = anharm(boost)
                        # print(f"{boost_avg:.2f}\t\t{boost_std:.2f}\t{boost_min:.2f}\t{boost_max:.2f}")
                        # calculate histogram of boost
                        hist, edges = np.histogram(boost, bins = "auto", weights=None)
                        if nlen == 5 :
                            histo_file = f'run_{run}/{y}_lambda/{filename}-len{l}-dV-histo.xvg'
                        else :
                            histo_file = f'run_{run}/{y}_lambda/{filename}-dV-histo.xvg'
                        output_dV_histo(histo_file,edges,hist)

                    if len(dvdl)>0 and len(boost)>0 :
                        print(f"{run}\t{y}\t{l/nlen:.2f}\t{dvdl_avg:.2f}\t\t{dvdl_std:.2f}\t{dvdl_min:.2f}\t{dvdl_max:.2f}\t{boost_avg:.2f}\t\t{boost_std:.2f}\t{boost_min:.2f}\t{boost_max:.2f}\t{boost_anharm:.2e}")
                    else :
                        print(f"{run}\t{y}\t{l/nlen:.2f}\t{dvdl_avg:.2f}\t\t{dvdl_std:.2f}\t{dvdl_min:.2f}\t{dvdl_max:.2f}")
        return

    # calculate free energy
    ltime = np.zeros(nlen)
    dG_avg = np.zeros(nlen)
    dG_std = np.zeros(nlen)
    dG_mat = np.zeros((nlen, nruns))
    for l in range(1,nlen+1):
        print(f"\nSimulation Length: {l/nlen:.2f}")
        dG_list = np.zeros(nruns)
        for run in range(1, nruns+1):
            # calculate dfdl
            dfdl_list = []
            dfdl_std_list = []

            if method == 'Exp':
                print(f"Window\tlambda\tdfdl\tboost_avg")
            elif method == 'ExpAverage' or method == 'JE':
                print(f"Window\tlambda\tdfdl\tstd")
            elif method == 'Average':
                print(f"Window\tlambda\tdfdl\tstd\tdvdl\tstd\tboost\tstd")
            elif method == 'Gaussian' or method == 'CE12' or method == 'CE12_21':
                print(f"Window\tlambda\tdfdl\tstd\tc1\tc2")

            for y in range(len(lambdas)):
                gamd_file = f'run_{run}/{y}_lambda/{filename}.dat'
                mdout_file = f'run_{run}/{y}_lambda/{filename}.out'
                if os.path.exists(gamd_file) :
                    time, dvdl, boost = read_from_gamd_file(gamd_file)
                elif os.path.exists(mdout_file) and not os.path.exists(gamd_file):
                    dvdl = extract_dvdl_mdout(mdout_file)
                else :
                    print("Error with reading input file.")
                    return

                # calculate simulation length
                if run == 1 and y == 0 :
                    # tl=int(len(dvdl)*l/nlen)
                    tl=l/nlen
                    ltime[l-1]=tl
                    print(f"Simulation Length: {tl:.2f}")

                if len(dvdl)>0 :
                    dvdl=dvdl[0:int(len(dvdl)*l/nlen)]
                    dvdl_avg = np.average(dvdl)
                    dvdl_std = np.std(dvdl)

                if len(boost)>0 :
                    boost=boost[0:int(len(boost)*l/nlen)]

                if method == 'noweight':
                    dfdl, dfdl_std = calculate_average_noweight(dvdl)
                    dfdl_std_list.append(dfdl_std)
                elif method == 'ExpAverage':
                    dfdl, dfdl_std = reweight_exp_average(dvdl, boost, beta)
                    dfdl_std_list.append(dfdl_std)
                    print(f"{y}\t{lambdas[y]}\t{dfdl:.2f}\t{dfdl_std:.2f}")
                elif method == 'JE':
                    dfdl, dfdl_std = reweight_Jarzynski(dvdl, boost, beta)
                    dfdl_std_list.append(dfdl_std)
                    print(f"{y}\t{lambdas[y]}\t{dfdl:.2f}\t{dfdl_std:.2f}")
                elif method == 'Average':
                    dfdl, dfdl_std, dvdl_avg, dvdl_std, boost_avg, boost_std = reweight_average(dvdl, boost, beta)
                    dfdl_std_list.append(dfdl_std)
                    print(f"{y}\t{lambdas[y]}\t{dfdl:.2f}\t{dfdl_std:.2f}\t{dvdl_avg:.2f}\t{dvdl_std:.2f}\t{boost_avg:.2f}\t{boost_std:.2f}")
                elif method == 'Gaussian':
                    c1, c2 = reweight_CE12(dvdl, boost, beta)
                    c12 = np.add(c1,c2)
                    dfdl = dvdl_avg+c12
                    dfdl_std = np.std(dvdl) # NEED correction
                    dfdl_std_list.append(dfdl_std)
                    print(f"{y}\t{lambdas[y]}\t{dfdl:.2f}\t{dfdl_std:.2f}\t{c1:.2f}\t{c2:.2f}")
                elif method == 'CE12':
                    c1, c2 = reweight_CE12(dvdl, boost, beta)
                    c12 = np.add(c1,c2)
                    dfdl = dvdl_avg+c12
                    dfdl_std = np.std(dvdl) # NEED correction
                    dfdl_std_list.append(dfdl_std)
                    print(f"{y}\t{lambdas[y]}\t{dfdl:.2f}\t{dfdl_std:.2f}\t{c1:.2f}\t{c2:.2f}")
                elif method == 'CE12_21':
                    c1, c2 = reweight_CE12_21(dvdl, boost, beta)
                    c12 = np.add(c1,c2)
                    dfdl = dvdl_avg+c12
                    dfdl_std = np.std(dvdl) # NEED correction
                    dfdl_std_list.append(dfdl_std)
                    print(f"{y}\t{lambdas[y]}\t{dfdl:.2f}\t{dfdl_std:.2f}\t{c1:.2f}\t{c2:.2f}")

                dfdl_list.append(dfdl)
                # print(f"Window lambda dfdl: {y} {lambdas[y]} {dfdl:.2f}")

            # calculate free energy
            if method == 'noweight':
                dG, dG_std_dvdl = integrate_dfdl_std(dfdl_list,dfdl_std_list,lambdas)
                print(f"Run {run} ΔG = {dG:.2f} +/- {dG_std_dvdl:.4f}")
            # elif method == 'ExpAverage' or method == 'JE' or method == 'Average' or method == 'Gaussian' or method == 'CE12' or method == 'CE12_21':
            else :
                dG = integrate_dfdl(dfdl_list,lambdas)
                print(f"Run {run} ΔG = {dG:.2f}")
            dG_list[run-1] = dG
            dG_mat[l-1,run-1] = dG

        dG_avg[l-1] = np.mean(dG_list)
        dG_std[l-1] = np.std(dG_list)
        print(f"Average of all runs: ΔG = {np.mean(dG_list):.2f} +/- {np.std(dG_list):.4f}")
        # print(f"Average of all runs: ΔG = {np.mean(dG_list):.2f} +/- {np.std(dG_list, ddof=1) / np.sqrt(len(dG_list)):.4f}")

    dGfile = 'dG-avg-'+str(filename)+'-'+str(method)+'.xvg'
    output_dG_avg(dGfile,dG_avg,dG_std,ltime)

    dGfile = 'dG-mat-'+str(filename)+'-'+str(method)+'.xvg'
    output_dG_mat(dGfile,dG_mat,ltime)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate binding free energy using different methods.')
    parser.add_argument('-method', type=str, choices=['stats', 'noweight', 'ExpAverage', 'JE', 'Average', 'Gaussian', 'CE12', 'CE12_21'], default='Gaussian',
                        help='Method to use for calculation (stats, noweight, ExpAverage, JE, Gaussian).')
    parser.add_argument("-input", dest="input", required=True, help="input file", metavar="<input file>")
    parser.add_argument("-T", dest="T", required=False,  help="Temperature", metavar="<Temperature>")
    parser.add_argument("-nruns", type=int, default='3', dest="nruns", required=False,  help="Number of runs", metavar="<nruns>")
    parser.add_argument("-nlen", type=int, default='2', dest="nlen", required=False,  help="Number of simulation lengths used to calculate free energy", metavar="<nlen>")
    args = parser.parse_args()
    main(args.method)

