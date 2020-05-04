#! /usr/bin/env python3

from scipy.stats import pearsonr, norm
import numpy as np #type: ignore
import pandas as pd #type: ignore
from typing import List, Dict, Tuple, Optional

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

def standard_error_slope(x: np.array, y: np.array, yhat: np.array, N: int) -> float:
    """Calculate standard error of slope estimate
    Args:
        x (np.array): numpy array of observed estimators
        y (np.array): numpy array of observed values
        yhat (np.array): numpy array of estimated values
        N (int): number of observations, >= 3
    Returns:
        (float): standard error of slope estimate
    """
    x_avg = np.average(x)
    return np.sqrt( np.sum( (y-yhat)**2 ) / (N-2)) / np.sqrt(np.sum( (x-x_avg)**2 ) )

def calculate_regression(x: np.array, y: np.array, weights: Optional[np.array] = None) -> List[float] :
    """Calculate regression coefficients from data
    Use np.linalg.polyfit for linear regression. std.err is std.err of slope estimate.
    Args:
        x (np.array): numpy array of x-coordinates
        y (np.array): numpy array of y-coordinates
        weights (Optional[np.array]): numpy array of point weights
    Returns:
        (List[float]): List with [intercept, slope, stderr]
    """
    coeff = np.polynomial.polynomial.polyfit(x,y,deg=[1],w=weights)
    intercept=coeff[0]
    slope = coeff[1]
    N=x.shape[0]
    yhat = intercept + slope*x
    stderr = standard_error_slope(x,y,yhat,N)
    return [ intercept, slope, stderr]