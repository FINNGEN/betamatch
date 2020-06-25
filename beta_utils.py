#! /usr/bin/env python3

from scipy.stats import pearsonr, norm
import statsmodels.regression.linear_model as lm
import numpy as np #type: ignore
import pandas as pd #type: ignore
from typing import List, Dict, Tuple, Optional
from collections import namedtuple

def weighted_cov(x: np.array, y: np.array, w: np.array) -> float:
    """Weighted covariance between vectors x and y, with weights w
    Args:
        x (np.array): numerical vector x
        y (np.array): numerical vector y
        w (np.array): weights for calculating weighted average of x
    Returns:
        (float): weighted covariance
    """
    return np.average( ( (x-np.average(x, weights=w) ) * (y - np.average(y, weights=w) ) ) , weights=w )

def weighted_pearsonr(x: np.array, y: np.array, w: np.array) -> float:
    return weighted_cov(x, y, w) / np.sqrt( weighted_cov(x, x, w) * weighted_cov(y, y, w) )

def calculate_r2(dataset: pd.DataFrame, x_label: str, y_label: str, stderr_label: str) -> Tuple[float, float, int, int]:
    """Calculate r2 values for dataset
    Args:
        dataset (pd.DataFrame): dataset
        x_label (str): x values column label in dataset
        y_label (str): y values column label in dataset
        stderr_label (str): standard error column label in dataset
    Returns:
        (Tuple[float,float,int,int]): Tuple of values: r^2 without weights, r^2 weighted by inverse stderr, N of r^2 calculation, N of r^2 weighted calculation
    """
    data_w=dataset[[x_label,y_label,stderr_label]].copy()
    data_r2 = data_w[[x_label,y_label]].copy()
    data_r2=data_r2.dropna(how="any")
    data_w=data_w.dropna(how="any")
    r_2=np.nan
    r_w=np.nan
    N_r=np.nan
    N_w=np.nan
    if data_r2.shape[0]>= 2:
        x_array = data_r2[x_label].values
        y_array = data_r2[y_label].values
        r,_=pearsonr(x_array,y_array)
        r_2=r**2
        N_r=data_r2.shape[0]
    if data_w.shape[0]>=2:
        x_array = data_w[x_label].values
        y_array = data_w[y_label].values
        stderr = data_w[stderr_label].values
        weight_array= 1/(stderr**2 + 1e-9) #weights as inverse of variance 
        r_w=weighted_pearsonr(x_array,y_array,weight_array)
        r_w=r_w**2
        N_w=data_w.shape[0]
    return (r_2,r_w,N_r,N_w)


RegressionResults = namedtuple('RegressionResults',['slope', 'stderr', 'tstat', 'pval', 'rsquared'])

def calculate_regression(x: np.array,
                        y: np.array,
                        weights: Optional[np.array] = None)-> RegressionResults :
    """Calculate regression y=bx coefficients from data.
    Model uses statsmodels.regression.linear_model.WLS.
    Args:
        x (np.array): numpy array of x-coordinates
        y (np.array): numpy array of y-coordinates
        weights (Optional[np.array]): numpy array of point weights. These are passed directly to WSL.
    Returns:
        (RegressionResults): Named tuple with variables slope, stderr, t-statistic, pvalue, rsquared
    """
    weights = weights if type(weights) != type(None) else 1.0
    model = lm.WLS(y,x,weights)
    results = model.fit()
    slope=results.params[0]
    stderr = results.bse[0]
    tstat = results.tvalues[0]
    pval = results.pvalues[0]
    rsquared = results.rsquared

    return RegressionResults(slope, stderr, tstat, pval,rsquared)