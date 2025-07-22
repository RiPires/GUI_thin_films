from new import film_thickness
import numpy as np
from sympy import sympify
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from math import erf, sqrt, pi

def load_mca_data(filename):
    with open(filename, 'r', encoding='latin1') as file:
        lines = file.readlines()
    # Finding where the data starts and ends
    start = 12
    end = [i for i, line in enumerate(lines) if '<<END>>' in line][0]
    # Reading data lines
    data = []
    for line in lines[start:end]:
        parts = line.strip().split()
        if parts:
            data.append(float(parts[0]))
    return data

source_data = load_mca_data('Data/NoFilm_for_C12_4h_2.mca')
source_film_data = load_mca_data('Data/SampleC12_4h.mca')

# Setting ROI and attenuation coefficient
x1 = int(input("Enter the start channel: "))  # start channel
x2 = int(input("Enter the end channel: "))  # end channel 
mu = float(sympify(input("Enter the attenuation coefficient in nm^-1: "))) # attenuation coefficient (nm^-1)

x_range = list(range(x1, x2 + 1))
y_source = source_data[x1:x2 + 1]
y_film = source_film_data[x1:x2 + 1]

# Define the Gaussian function
def gaussian(x, A, x0, sigma, B):
    return A * np.exp(-(x - x0)**2 / (2 * sigma**2)) + B

# Fit source data
popt_source, _ = curve_fit(gaussian, x_range, y_source, p0=[max(y_source), x_range[np.argmax(y_source)], 1, 0])
A0, x0_0, sigma0, B0 = popt_source

# Fit film data
popt_film, _ = curve_fit(gaussian, x_range, y_film, p0=[max(y_film), x_range[np.argmax(y_film)], 1, 0])
A1, x0_1, sigma1, B1 = popt_film

def gaussian_integral_erf(A, x0, sigma, B, x1, x2):
    erf_part = 0.5 * (erf((x2 - x0) / (sqrt(2) * sigma)) - erf((x1 - x0) / (sqrt(2) * sigma)))
    area_gaussian = A * sigma * sqrt(2 * pi) * erf_part
    area_below= B * (x2 - x1)
    return area_gaussian + area_below

# Calculate Gaussian integrals
I0 = gaussian_integral_erf(*popt_source, x1, x2)
I = gaussian_integral_erf(*popt_film, x1, x2)

## Using the previous method
I0_ = A0 * sigma0
I_ = A1 * sigma1

# Film thickness calculation
thickness = film_thickness(I, I0, mu)
thickness_ = film_thickness(I_, I0_, mu) ## previous method
print(f"Film thickness: {round(thickness)} nm")
print(f"Film thickness (prev method): {round(thickness_)} nm")


# Plot fits
# Plot 1: Source (no film)
plt.figure()
plt.plot(x_range, y_source, '*', label='Source', color='blue')
plt.plot(x_range, gaussian(x_range, *popt_source), '--', label='Gaussian Fit (Source)', color='cyan')
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.title('Gaussian Fit - Source (No Film)')
plt.legend()
plt.grid(True)

# Plot 2: Film
plt.figure()
plt.plot(x_range, y_film,'x', label='Source + Film', color = 'green')
plt.plot(x_range, gaussian(x_range, *popt_film), '--', label='Gaussian Fit (Source + Film)', color='lime')
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.title('Gaussian Fit - Source + Film')
plt.legend()
plt.grid(True)

plt.show()