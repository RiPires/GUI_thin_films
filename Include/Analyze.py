from Include.FitData import*
import numpy as np

def Analyze(y, ROId, ROIu):
    centroids = []
    sigmas = []
    error = []

    for r in range(len(ROId)):
        try:
            x1 = float(ROId[r])
            x2 = float(ROIu[r])
        except (ValueError, TypeError):
            continue  # Skip invalid input

        # Skip if ROI is zero-width, negative, or too narrow to be useful
        if x2 <= x1 or (x2 - x1) < 2:
            continue

        Cent = peakCentroid(x1, x2, y)
        Sigma = peakSigma(x1, x2, y)
        Net = peakNet(x1, x2, y)

        centroids.append(Cent)
        sigmas.append(Sigma)
        error.append(Sigma / np.sqrt(Net) if Net > 0 else 0)

    return centroids, error, sigmas
