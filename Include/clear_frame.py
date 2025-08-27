##########################################################################################################
##                                                                                                      ##       
## ARC-TF stands for Alpha particles' energy loss and Rutherford backscattering                         ##
## spectrometry methods for Characterization of Thin Films.                                             ##
##                                                                                                      ##
## ARC-TF is a Python based GUI (Guided User Interface) intended to expedite the process of             ##
## characterizing thin films, via an interface with ease of use, and backend algorithms that            ##
## accelerate the data analysis.                                                                        ##
##                                                                                                      ##
## This project is the result of a LIP Summer Internship, within the NUC-RIA group.                     ##
## A publication resulting from the internship, resuming the work taken to develop                      ##
## the first version of ARC-TF, can be searched for by the reference LIP-STUDENTS-23-15.                ##
## Directly available at https://www.lip.pt/files/training/papers/2023/pdf/2023-PAPER-179-15.pdf.       ##
##                                                                                                      ##
## Authors: Alexandre Gusmão (1),           alex_vng@hotmail.com,   https://github.com/AlexVnGit        ##
##          Ricardo Matoza Pires (1,2),     rpires@lip.pt,          https://github.com/RiPires          ##
##          Tomás Campante Tavares (1,2),   tmctavares@lip.pt,      https://github.com/TomasCampante    ##
##                                                                                                      ##
##  (1) Faculty of Sciences of the Universty of Lisbon (FCUL),                                          ##
##      Rua Ernesto de Vasconcelos, 1749-016 Lisboa, Portugal,                                          ##
##      Building C8, 5th floor, room 8.5.15                                                             ##
##                                                                                                      ##
##  (2) Laboratory of Instrumentation and Experimental Particle Physics (LIP)                           ##
##      NUC-RIA group (Nuclear Reactions, Instrumentation and Astrophysics)                             ##
##      Av. Prof. Gama Pinto, 2, 1649-003 Lisboa, Portugal                                              ##
##                                                                                                      ##
##########################################################################################################

def clear_frame(frame, remove_grid=True):
    """Destroy all widgets in the given frame and optionally remove it from the grid."""
    for widget in frame.winfo_children():
        widget.destroy()
    if remove_grid:
        frame.grid_remove()