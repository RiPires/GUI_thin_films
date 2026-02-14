from Include.FitData import*
import numpy as np
import warnings

def Analyze(y, ROId, ROIu):
    centroids = []
    sigmas = []
    error = []

    for r in range(len(ROId)):
        try:
            x1 = float(ROId[r])
            x2 = float(ROIu[r])
        except (ValueError, TypeError):
            warnings.warn(f"Invalid ROI input at index {r}: ROId = {ROId[r]}, ROIu = {ROIu[r]}", UserWarning)
            continue  # Skip invalid input

        # Skip if ROI is zero-width, negative, or too narrow to be useful
        if x2 <= x1:
            warnings.warn(f"ROI upper bound must be greater than lower bound at index {r}: x1 = {x1}, x2 = {x2}", UserWarning)
            continue
        elif (x2 - x1) < 2:
            warnings.warn(f"ROI too narrow at index {r} (width = {x2 - x1}). Minimum width of 2 required.", UserWarning)
            continue

        Cent = peakCentroid(x1, x2, y)
        Sigma = peakSigma(x1, x2, y)
        Net = peakNet(x1, x2, y)

        centroids.append(Cent)
        sigmas.append(Sigma)
        error.append(Sigma / np.sqrt(Net) if Net > 0 else 0)

    return centroids, error, sigmas
