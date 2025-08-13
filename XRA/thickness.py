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
    
    #both backgrounds

    bckg_file = background_file(x1, x2, background_data) if background_data is not None else None
    bckg_calc = calc_bckg(x1, x2, values_list)

    soma = 0
    for idx, i in enumerate(range(int(x1), int(x2) + 1)):
        myCellValue = values_list[i - 1]
        if bckg_file is not None:
            bgd = max(bckg_file[idx], bckg_calc[idx])
        else:
            bgd = bckg_calc[idx]
        soma += myCellValue - bgd

    return soma

"""
def film_thickness(I, I0, mu):  
    
    
    Calculate film thickness using t = -1/mu * ln(I/I0)
    I: net peak area with film (source+film - background_data)
    I0: net peak area without film (source - background)
    mu: attenuation coefficient
    

    if I <= 0 or I0 <= 0 or mu == 0:
        return float('nan')
    return -1.0 / mu * np.log(I / I0)

"""
def film_thickness(N, N0, Nb, mu):
    N_net = N - Nb
    N0_net = N0 - Nb


    if N_net <= 0 or N0_net <= 0 or mu == 0:
        print("i am here")
        return float('nan')
    
    return -1.0 / mu * np.log(N_net / N0_net)

def uncertainty(N, N0, Nb, mu, time_source, time_film, time_bkg):
    I_net = N/time_film - Nb/time_bkg
    comp1 = Nb/time_bkg - N/time_film
    comp2 = N0/time_source - Nb/time_bkg
    comp3 = N0/time_source - N/time_film

    if N <= 0 or N0 <= 0 or mu == 0:
        return float('nan')
    

    num = N0 / (time_source * time_source)*comp1*comp1 + N / (time_film *time_film)*comp2*comp2 + Nb / (time_bkg * time_bkg) *comp3*comp3
    den = (comp2*comp2 * I_net*I_net)*mu*mu

    return np.sqrt( num/den)  #((I0 / time_source)*comp1**2 + (I / time_film)*comp2**2 + (Ib / time_bkg)*comp3**2) / ((comp2**2)*(I_net**2)))