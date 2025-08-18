from thickness import film_thickness, uncertainty
import numpy as np
from sympy import sympify
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from math import erf, sqrt, pi
import csv


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


def get_measurement_time(filename):

    with open(filename, 'r', encoding='latin1') as file:
        reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
        data = list(reader)
        time_str = data[7][0]  #8th line
        time = float(time_str[time_str.find('- ') + 2:])
        return time


time_source = get_measurement_time('Data/NoFilm_2h_17Jul_CuSource.mca')
time_film = get_measurement_time('Data/Al_foil_CuSource_01Aug.mca')
time_bkg = get_measurement_time('Data/Combined_bkg_Internship.txt')

print(f"Measurement time (source): {time_source} s")
print(f"Measurement time (film): {time_film} s")
print(f"Measurement time (background): {time_bkg} s")


source_data = load_mca_data('Data/NoFilm_2h_17Jul_CuSource.mca')
source_film_data = load_mca_data('Data/Al_foil_CuSource_01Aug.mca')
bkg_data = load_mca_data('Data/Combined_bkg_Internship.txt')

# Setting ROI and attenuation coefficient
x1 = int(input("Enter the start channel: "))  # start channel
x2 = int(input("Enter the end channel: "))  # end channel 
mu = float(sympify(input("Enter the attenuation coefficient in nm^-1: "))) # attenuation coefficient (nm^-1)

x_range = list(range(x1, x2 + 1))
y_source = source_data[x1:x2 + 1]
y_film = source_film_data[x1:x2 + 1]
y_bkg = bkg_data[x1:x2 + 1]  

m = 0.030826941169
b = 0.092673711109335
energy_range = [m * ch + b for ch in x_range] #converting channels to energy
x_smooth = np.linspace(x1, x2, 500)
energy_smooth = [m * ch + b for ch in x_smooth]

# Define the Gaussian function
def gaussian(x, A, x0, sigma, B):
    return A * np.exp(-(x - x0)**2 / (2 * sigma**2)) + B

# Fit source data
popt_source, _ = curve_fit(gaussian, x_range, y_source, p0=[max(y_source), x_range[np.argmax(y_source)], 1, 0])
A0, x0_0, sigma0, B0 = popt_source

print(f"Source centroid: {x0_0:.1f}")

# Fit film data
popt_film, _ = curve_fit(gaussian, x_range, y_film, p0=[max(y_film), x_range[np.argmax(y_film)], 1, 0])
A1, x0_1, sigma1, B1 = popt_film

# Fit background data
popt_bkg, _ = curve_fit(gaussian, x_range, y_bkg, p0=[max(y_bkg), x_range[np.argmax(y_bkg)], 1, 0])
A2, x0_2, sigma2, B2 = popt_bkg

def gaussian_integral_erf(A, x0, sigma, B, x1, x2):
    erf_part = 0.5 * (erf((x2 - x0) / (sqrt(2) * sigma)) - erf((x1 - x0) / (sqrt(2) * sigma)))
    area_gaussian = A * sigma * sqrt(2 * pi) * erf_part
    area_below= B * (x2 - x1)
    return area_gaussian + area_below

# Calculate Gaussian integrals
N0 = gaussian_integral_erf(*popt_source, x1, x2)
N = gaussian_integral_erf(*popt_film, x1, x2)
Nb = gaussian_integral_erf(*popt_bkg, x1, x2)


# Film thickness calculation
thickness = film_thickness(N, N0, Nb, mu)
print(f"Film thickness: {round(thickness)} nm")
uncertainty_value = uncertainty(N, N0, Nb, mu, time_source, time_film, time_bkg)
print(f"Uncertainty: {uncertainty_value:.2f}")
print(f"Full answer: {thickness:.2f} \u00B1 {uncertainty_value:.2f} nm")

# Plot fits
# Plot 1: Source (no film)
plt.figure()
plt.plot(x_range, y_source, '*', label='Source', color='blue')
plt.plot(x_smooth, gaussian(x_smooth, *popt_source), '--', label='Gaussian Fit (Source)', color='cyan')
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.title('Gaussian Fit - Source (No Film)')
plt.legend()
plt.grid(True)

# Plot 2: Film
plt.figure()
plt.plot(x_range, y_film,'x', label='Source + Film', color = 'green')
plt.plot(x_smooth, gaussian(x_smooth, *popt_film), '--', label='Gaussian Fit (Source + Film)', color='lime')
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.title('Gaussian Fit - Source + Film')
plt.legend()
plt.grid(True)

#Plot 3: Source - counts vs energy
plt.figure()
plt.plot(energy_range, y_source, '*', label='Source (Energy)', color='navy')
plt.plot(energy_smooth, gaussian(x_smooth, *popt_source), '--', label='Gaussian Fit', color='dodgerblue')
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Gaussian Fit - Source (No Film) [vs Energy]')
plt.legend()
plt.grid(True)

# Plot 4: Film - counts vs energy
plt.figure()
plt.plot(energy_range, y_film, 'x', label='Source + Film (Energy)', color='darkgreen')
plt.plot(energy_smooth, gaussian(x_smooth, *popt_film), '--', label='Gaussian Fit', color='limegreen')
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Gaussian Fit - Source + Film [vs Energy]')
plt.legend()
plt.grid(True)

plt.show()