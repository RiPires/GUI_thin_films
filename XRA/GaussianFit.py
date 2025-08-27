from XRA.thickness import film_thickness, uncertainty
import numpy as np
from sympy import sympify
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from math import erf, sqrt, pi
import csv


# Define the Gaussian function
def gaussian(x, A, x0, sigma, B):
    return A * np.exp(-(x - x0)**2 / (2 * sigma**2)) + B

def gaussian_fit(func, x_range, y_data):
    """
    Fit a Gaussian function to the data.
    
    Parameters:
    - func: The Gaussian function to fit.
    - x_range: The range of x values.
    - y_data: The y data to fit.
    - p0: Initial guess for the parameters.
    
    Returns:
    - popt: Optimal parameters for the fit.
    """
    # Initial guess
    init = [max(y_data), x_range[np.argmax(y_data)], 1, 0]
    ## Fit data to gaussian profile  
    popt, _ = curve_fit(func, x_range, y_data, p0=init)
    return popt

def gaussian_integral_erf(A, x0, sigma, B, x1, x2):
    erf_part = 0.5 * (erf((x2 - x0) / (sqrt(2) * sigma)) - erf((x1 - x0) / (sqrt(2) * sigma)))
    area_gaussian = A * sigma * sqrt(2 * pi) * erf_part
    area_below= B * (x2 - x1)
    return area_gaussian + area_below

