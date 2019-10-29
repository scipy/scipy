#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 15:34:53 2019
@author: Solomon Vimal (solomon.vimal@gmail.com)
Based on R code (from Prof. Rong Fu, Oct, 2019 AOS2019 Fall 2019, UCLA)
& Gauss program by Vogelsang (https://msu.edu/~tjv/working.html) from 8/13/96
"""

from __future__ import division, print_function, absolute_import
import sys
import numpy as np
import pandas as pd
from numpy.linalg import inv

def vogelsang_trend(y):
    """
    Computes the individual tests using the t-PS test of Vogelsang 
    (1998, Econometrica)
    
    Test the null hypothesis that a time series with serial correlation has a 
    slope with statistically significant trend.
    
    The test is evaluated at 90%, 95%, 97.5%, and 99.0% confidence levels.
    
    Parameters
    ----------
    y : array-like 
        The time series for which the trend is to be calculated and tested for 
        statistical significance.
        
    Returns
    -------
    slope: float
        The slope of the trend function
    pvalue: binary 
        Binary result (1 or 0) for pass or fail for 4 levels of statistical 
        significance at: 90%, 95%, 97.5%, and 99.0% confidence levels.
        
    See Also
    --------
    https://www.jstor.org/stable/2998543
    
    Notes
    -----
    This test is applicable when the time series has a serial auto correlation. 
    which in mann-kendall and modified mann-kendall test will have a tendency 
    for over rejection.
    
    This procedure computes the nonparametric estimator of the long run variance 
    using Andrews (1991) AR(1) plug in optimal truncation lag.  
    v is the vector of residuals and M is the trucation lag.  
    M<0 uses the automatic bandwidth.  
    Kernel is an integer representing the kernel that is used 1 = bartlett, 
    2 = parzen, 3 = Tukey-Hanning, 4 = quadratic spectral
    if prewhite == 1 then prewhitening using an AR(1) model is used
    if prewhite == 0 no prewhitening is used.

    References
    ----------
    Vogelsang, Timothy J. 1998. “Trend Function Hypothesis Testing in the 
    Presence of Serial Correlation.” Econometrica 66 (1): 123–48. 
    https://doi.org/10.2307/2998543.
    
    """
    
    print("\n###############################################\n")
    print("Trend Function Hypothesis Test (Vogelsang, 1998)")
    print("https://www.jstor.org/stable/2998543")
    print("\n###############################################\n")

    # Critical values for right tailed 5%, 2.5%, 1% test and one-sided
    # Confidence intervals

    cvtps = [1.331, 1.720, 2.152, 2.647]
    btps = [0.494, 0.716, 0.995, 1.501]
    
    z = np.cumsum(y)    
    T = len(y)
    t1 = range(1,T+1)
    X1 = pd.DataFrame({"T":[1]*T, "t1":t1})

    # Compute the partial sums of X1
    nr = X1.shape[0] #nc = X1.shape[1]
    X2 = np.cumsum(X1)
    
    # Compute inverse matrices
    X1inv = inv(np.matmul(X1.T, X1))
    X2inv = inv(np.matmul(X2.T, X2))
    
    # Compute orthonormal trends for the J statistics
    Xj = np.array(X1["T"] / np.sqrt(nr))
    Xj = pd.DataFrame(Xj)
    i = 1;
    for i in range(1,10):
        t1 = np.arange(1,len(y)+1)
        ehat = tt - np.matmul(Xj, np.matmul(inv(np.matmul(Xj.T, Xj)), 
                                            np.matmul(Xj.T,tt)))
                                            np.matmul(Xj.T,tt)))
        ehat = ehat / np.sqrt(np.matmul(ehat.T, ehat));
        Xj = pd.concat([Xj, pd.DataFrame(ehat)], axis=1)
        
    Xjinv = inv(np.matmul(Xj.T, Xj))
    
    # Compute stats in standard regressions
    bhat = np.matmul( np.matmul(X1inv, X1.T) , np.array(y))
    uhat = y - np.matmul(X1,bhat)
    rss1 = np.matmul(uhat.T, uhat)
    
    #s2dan = sig2np(v, M, my_kernel=1, prewhite=0)
    btild = np.matmul(np.matmul(X2inv, X2.T), z)
    jbeta = np.matmul(np.matmul(Xjinv, Xj.T), y)
    rssj = np.matmul((y-np.matmul(Xj,jbeta)).T, (y - np.matmul( Xj, jbeta))) 
    s2z = np.matmul((z-np.matmul(X2,btild)).T, 
                    (z - np.matmul( X2, btild))) / (X2.shape[0]-X2.shape[1]) 
    J = (rss1-rssj) / rssj
    
    tps = (btild[1]/np.sqrt(T*s2z*X2inv[1,1]))*np.exp(-np.array(btps)*J)[0]
    tps_100 = [1 if ((tps[0][0]>cvtps[0]) | (tps[0][0]<-1*cvtps[0])) else 0][0]
    tps_050 = [1 if ((tps[0][1]>cvtps[1]) | (tps[0][1]<-1*cvtps[1])) else 0][0]
    tps_025 = [1 if ((tps[0][2]>cvtps[2]) | (tps[0][2]<-1*cvtps[2])) else 0][0]
    tps_010 = [1 if ((tps[0][3]>cvtps[3]) | (tps[0][3]<-1*cvtps[3])) else 0][0]
    
    results = [bhat[1], tps_100, tps_050, tps_025, tps_010]
    
    #results = vogelsang_trend(df)
    trend = results[0]
    # pValue equals to 0 means the trend is not significant at the 
    #corresponding confidence level, i.e. cannot pass the test
    # pValue equals to 1 means the trend is significant at the corresponding 
    #cofindence level, i.e. pass the test
    alpha_values = [0.1, 0.05, 0.025, 0.01]
    confidence_levels = [100*(1-x) for x in alpha_values]
    
    print("The trend slope is " + str("{:.3f}".format(trend[0])) + "\n")
    for i, p in enumerate(results[1:]):
        p
        if p == 0:
            print("Not significant at " + str(confidence_levels[i]) 
            + "% confidence (alpha=" + str(alpha_values[i]) + ")")
        else:
            print("Significant at " + str(confidence_levels[i]) 
            + "% confidence (alpha=" + str(alpha_values[i]) + ")")
    print("\n")
    
    return results

if __name__ == "__main__":
    filename = sys.argv[1:][0]
    #filename = "/home/svimal/Desktop/Fall_2019_courses/AOS_Climate_Stats_Viz/week1/DJFLArainfall.txt"
    y = pd.read_csv(filename,header=None)
    vogelsang_trend(y)
