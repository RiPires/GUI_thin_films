from new import peakSum, film_thickness
import numpy as np
from sympy import sympify

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
background_data = load_mca_data('Data/Fundo_acqTest_1hour_01.mca')

# Setting ROI and attenuation coefficient
x1 = int(input("Enter the start channel: "))  # start channel
x2 = int(input("Enter the end channel: "))  # end channel 
mu = float(sympify(input("Enter the attenuation coefficient in nm^-1: "))) # attenuation coefficient (nm^-1)


# Net peak areas
I0 = peakSum(x1, x2, source_data)
I  = peakSum(x1, x2, source_film_data)
# Calculate film thickness
thickness = film_thickness(I, I0, mu)
print(f"Film thickness: {round(thickness)} nm")