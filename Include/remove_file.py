####################################################################################################
##                                                                                                ##       
## ARC-TF stands for Alpha particles' energy loss and Rutherford backscattering                   ##
## spectrometry methods for Characterization of Thin Films.                                       ##
##                                                                                                ##
## ARC-TF is a GUI (Guided User Interface) intended to expedite the process of                    ##
## characterizing thin films, via an interface with ease of use, and backend                      ##
## algorithms that accelerate the data analysis.                                                  ##
##                                                                                                ##
## This project is the result of a LIP Summer Internship, within the NUC-RIA group.               ##
## A publication resulting from the internship, resuming the work taken to develop                ##
## the first version of ARC-TF, can be searched for by the reference LIP-STUDENTS-23-15.          ##
## Directly available at https://www.lip.pt/files/training/papers/2023/pdf/2023-PAPER-179-15.pdf  ##
##                                                                                                ## 
####################################################################################################
import os

def remove_file(filepath):
    """Remove file safely if it exists, ignoring FileNotFoundError."""
    try:
        os.remove(filepath)
    except FileNotFoundError:
        pass