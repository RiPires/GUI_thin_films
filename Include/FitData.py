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

####################################
def peakNet(x1, x2, values_list):
    if x1 <= 0 or x2 <= 0 or x2 < x1:
        return "N/A"

    if x1 > len(values_list) or x2 > len(values_list):
        return "N/A"

    yL = values_list[int(x1) - 1]
    yU = values_list[int(x2) - 1]

    if x2 - x1 != 0:
        m = (yU - yL) / (x2 - x1)
    else:
        m = 0

    soma = 0
    for i in range(int(x1), int(x2) + 1):
        myCellValue = values_list[i - 1]
        bgd = m * (i - x1) + yL  # background
        soma += myCellValue - bgd

    return soma
####################################

####################################
def peakCentroid(x1, x2, values):
    if x1 <= 0 or x2 <= 0 or x2 < x1:
        return "N/A"

    yL = values[int(x1) - 1]
    yU = values[int(x2) - 1]

    if x2 - x1 != 0:
        m = (yU - yL) / (x2 - x1)
    else:
        m = 0

    soma = 0
    xsoma = 0
    for i in range(int(x1), int(x2) + 1):
        myCellValue = values[i - 1]

        if m == 0:
            bgd = 0
        else:
            bgd = m * (i - x1) + yL  # background

        soma += myCellValue
        xsoma += i * myCellValue

    print("here" + str(bgd))
    if soma > 0:
        return xsoma / soma
    else:
        return 0
####################################

####################################
def peakSigma(x1, x2, values_list):
    if x1 < 0 or x2 < 0 or x2 < x1:
        return "N/A"

    yL = values_list[int(x1) - 1] if x1 > 0 else 0
    yU = values_list[int(x2) - 1] if x2 > 0 else 0

    if x2 - x1 != 0:
        m = (yU - yL) / (x2 - x1)
    else:
        m = 0

    soma = 0
    xsoma = 0
    xxsoma = 0
    for i in range(int(x1), int(x2) + 1):
        if i > 0 and i <= len(values_list):
            myCellValue = values_list[i - 1]
        else:
            myCellValue = 0

        if m == 0:
            bgd = 0
        else:
            bgd = m * (i - x1) + yL  # background

        deviation = myCellValue - bgd
        soma += deviation
        xsoma += i * deviation
        xxsoma += i * i * deviation

    if soma > 0:
        variance = xxsoma / soma - xsoma * xsoma / (soma * soma)
        if variance > 0:
            peakSigma = variance ** 0.5
        else:
            peakSigma = 0
    else:
        peakSigma = 0

    return peakSigma
####################################

####################################
def peakFWHM(x1, x2, values_list):
    if x1 <= 0 or x2 <= 0 or x2 < x1:
        return "N/A"

    yL = values_list[int(x1) - 1] if x1 > 0 else 0
    yU = values_list[int(x2) - 1] if x2 > 0 else 0

    if x2 - x1 != 0:
        m = (yU - yL) / (x2 - x1)
    else:
        m = 0

    soma = 0
    xsoma = 0
    xxsoma = 0
    for i in range(int(x1), int(x2) + 1):
        if i > 0 and i <= len(values_list):
            myCellValue = values_list[i - 1]
        else:
            myCellValue = 0

        if m == 0:
            bgd = 0
        else:
            bgd = m * (i - x1) + yL  # background

        deviation = myCellValue - bgd
        soma += deviation
        xsoma += i * deviation
        xxsoma += i * i * deviation

    if soma > 0:
        variance = xxsoma / soma - xsoma * xsoma / (soma * soma)
        if variance > 0:
            peakFWHM = 2.35482 * (variance ** 0.5)
        else:
            peakFWHM = 0
    else:
        peakFWHM = 0

    return peakFWHM
####################################

####################################
def peakMAX(x1, x2, values_list):
    if x1 <= 0 or x2 <= 0 or x2 < x1:
        return "N/A"

    maxval = 0  # reset max. value
    xmax = 0.0
    for i in range(int(x1), int(x2) + 1):
        if i > 0 and i <= len(values_list):
            myCellValue = values_list[i - 1]
        else:
            myCellValue = 0

        if myCellValue > maxval:
            maxval = myCellValue
            xmax = i

    return xmax
####################################