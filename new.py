import numpy as np

####################################
def peakSum(x1, x2, list):
    peaksum = 0.0

    if x1 < 0 or x2 < 0 or x2 < x1:
        peaksum = float('nan')
    else:
        for i in range(int(x1), int(x2)):
            peaksum += list[i]

    return peaksum
####################################


def background_file(x1, x2, background_data):
    """Return background values from file for ROI [x1, x2] (inclusive)."""
    return background_data[int(x1)-1:int(x2)]

####################################

def calc_bckg(x1, x2, values_list):
    """Calculate linear background for ROI [x1, x2] (inclusive)."""
    yL = values_list[int(x1) - 1]
    yU = values_list[int(x2) - 1]
    m = (yU - yL) / (x2 - x1) if x2 - x1 != 0 else 0
    return [m * (i - x1) + yL for i in range(int(x1), int(x2) + 1)]

###################################

def peakNet(x1, x2, values_list, background_data=None):
    if x1 <= 0 or x2 <= 0 or x2 < x1:
        return "N/A"
    if x1 > len(values_list) or x2 > len(values_list):
        return "N/A"
    if background_data is not None and (x1 > len(background_data) or x2 > len(background_data)):
        return "N/A"

    soma = 0
    for i in range(int(x1), int(x2) + 1):
        myCellValue = values_list[i - 1]
        if background_data is not None:
            bgd = background[i - 1]
        else:
            yL = values_list[int(x1) - 1]
            yU = values_list[int(x2) - 1]
            m = (yU - yL) / (x2 - x1) if x2 - x1 != 0 else 0
            bgd = m * (i - x1) + yL
        soma += myCellValue - bgd

    return soma

def film_thickness(I, I0, mu):
    """
    Calculate film thickness using t = -1/mu * ln(I/I0)
    I: net peak area with film (source+film - background_data)
    I0: net peak area without film (source - background)
    mu: attenuation coefficient
    """
    if I <= 0 or I0 <= 0 or mu == 0 or I < I0:
        return float('nan')
    return -1.0 / mu * np.log(I / I0)


print("text")
## More code
print("AGAIN")