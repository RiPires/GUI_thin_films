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

## ------------------------------- Import necessary librarires ---------------------------------- ##
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
from PIL import ImageTk, Image
import webbrowser
from matplotlib.pyplot import axhline
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import os
import math
from shutil import copy2
import numpy as np
import ctypes
import sys

from Include.Analyze import*
from Include.Calibration import*
from Include.FitData import*
from Include.Eloss import*
from Include.Thick import*
from Include.remove_file import*
from Include.clear_frame import*
## ---------------------------------------------------------------------------------------------- ##

#####################################################
# Handles display scaling                           #
# works on Windows 10/11 and Linux/Kubuntu 20.04    #
#####################################################

# Adjust DPI awareness on Windows only
if sys.platform == "win32":
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception as e:
        print(f"Warning: DPI awareness setting failed: {e}")

################################################
# Returns the index of the Tab the user is on  #
################################################
def Current_Tab():
    """
    Returns the ID of the tab where the user is on,
    by getting the index of the selected tab in the notebook.
    """
    ## Get index of the current tab in the notebook
    tabID = tab_manager.notebook.index(tab_manager.notebook.select()) - 1

    ## Check if the user is in the "final results" tab (has negative index)
    if tabID < 0:
        pass
    else:
        return tabID
  
#####################################################################
# Clears and resets specific UI frames and associated data files    #
# depending on the type of frame specified                          #
#####################################################################
def ClearWidget(Frame, parameter):
    """
     Function: ClearWidget
     --------------------------------------------
     Purpose:
       Clears and resets specific UI frames and 
       associated data files depending on the 
       type of frame specified.
    
     Parameters:
       Frame     (str) : The name of the frame to clear.
                         Options include: 'Graphic', 'Algorithm',
                         'Results', 'Source', 'Popup', 'Image',
                         'Linear', 'Thickness', 'Final', 'Everything'
       parameter (int) : Used to control conditional deletion of 
                         files or reset of variables.
                         (1 = full reset with file deletion)
    
     Dependencies:
       - Uses the global TabList and TabTracker structures.
       - Assumes each frame is a Tkinter Frame with children.
       - Depends on external objects like `tab_manager`, `warnings_manager`.
    """
    num = Current_Tab()
    tab = TabList[num][1]

    if Frame == 'Graphic':
        # Clear all widgets from graphic-related frames
        clear_frame(tab.GraphicFrame)
        clear_frame(tab.Extra_Frame)

    elif Frame == 'Algorithm':
        # Clear previous algorithm UI and reset selection
        clear_frame(tab.AlgFrame)
        if parameter == 1:
            tab.Algorithm_Method.set('Select Algorithm to Run')
            tab.Algorithm.set(0)

    elif Frame == 'Results':
        # Clear result display and delete results file if required
        clear_frame(tab.ResultFrame)
        if parameter == 1:
            remove_file(TabList[num][3])  # Output file

    elif Frame == 'Source':
        # Clear source options UI and reset selection
        clear_frame(tab.SourceOptionsFrame)
        if parameter == 1:
            tab.Source.set('Radiation Sources')

    elif Frame == 'Popup':
        # Destroy warning popup
        for widget in warnings_manager.warning.winfo_children():
            widget.destroy()
        warnings_manager.warning.destroy()

    elif Frame == 'Image':
        # Destroy decay image popup
        for widget in warnings_manager.decay.winfo_children():
            widget.destroy()
        warnings_manager.decay.destroy()

    elif Frame == 'Linear':
        # Clear linear regression result display and optionally delete file
        clear_frame(tab.LinearRegressionFrame)
        clear_frame(tab_manager.Calib_Result2)
        if parameter == 1:
            remove_file(TabList[num][4])  # Calibration file

    elif Frame == 'Thickness':
        # Clear thickness calculation results and optionally delete file
        clear_frame(tab.ThicknessFrame)
        clear_frame(tab_manager.Mat_Result2)
        if parameter == 1:
            tab.Mat.set('Select Material')
            remove_file(TabList[num][4])  # Thickness file

    elif Frame == 'Final':
        # Clear final results display
        clear_frame(tab_manager.Mat_Result2)
        clear_frame(tab_manager.Calib_Result2)

    elif Frame == 'Everything':
        # Full reset of all UI elements and files in the current tab
        clear_frame(tab.GraphicFrame)
        clear_frame(tab.Extra_Frame)
        clear_frame(tab.AlgFrame)
        tab.Algorithm_Method.set('Select Algorithm to Run')
        tab.Algorithm.set(0)
        clear_frame(tab.ResultFrame)

        # Conditional clearing depending on tab type
        if TabTracker[num] < 0:
            clear_frame(tab.SourceOptionsFrame)
            tab.Source.set('Radiation Sources')
            clear_frame(tab.LinearRegressionFrame)

        elif TabTracker[num] > 0:
            clear_frame(tab.ThicknessFrame)
            tab.Mat.set('Select Material')

        # Delete associated output files if requested
        if parameter == 1:
            remove_file(TabList[num][3])
            remove_file(TabList[num][4])

############################################################################################
# Reads data from a text file and returns it as a 1D or 2D list, depending on the format.  #
# Handles different column separators, number formats (int/float/string), and optional     #
# real-time data extraction for GUI updates.                                               #
############################################################################################
def File_Reader(Document, Separator, Decimal, Upload):
    """
     Function: File_Reader
     ------------------------------------------------------------------------------------------
     Purpose:
       Reads data from a text file and returns it as a 1D or 2D list, depending on the format.
       Handles different column separators, number formats (int/float/string), and optional
       real-time data extraction for GUI updates.
    
     Parameters:
       Document  (str)  : Path to the text file to be read.
       Separator (str)  : Character used to separate columns (e.g. '\t', ',', ' ').
                          Use '0' to indicate a 1D list (no column separation).
       Decimal   (str)  : Controls data type conversion:
                          'Yes'  -> Convert values to float
                          'No'   -> Convert values to int
                          'String' -> Return raw strings (1D only)
       Upload    (str)  : If set to 'Yes', updates the GUI with a real-time value from the file.
    
     Returns:
       List of values:
         - 2D list if Separator != '0'
         - 1D list if Separator == '0'
         - Type of values (string/float/int) based on Decimal input
    
     Notes:
       - Assumes that if Upload is 'Yes', line 9 (index 8) contains real-time info (in seconds).
       - Assumes valid format and no missing separators for 2D inputs.
    """
    num = Current_Tab()

    with open(Document, 'r', encoding='latin1') as OpenFile:
        lines = OpenFile.read().splitlines()
    
    # Split the content by lines
    #lines = lines.splitlines()

    # If needed, update GUI with the "measurement live-time" value from line 9
    if Upload == 'Yes':
        TabList[num][1].Real_Time.set(lines[8] + ' s')

        if TabList[num][1].tab_kind == 5:
            try:
                TabList[num][1].Real_Time_2.set(lines[8] + ' s')
            except AttributeError:
                TabList[num][1].Real_Time_2 = tk.StringVar(value=lines[8] + ' s')

    start = 12
    try:
        end = next(i for i, line in enumerate(lines) if '<<END>>' in line)
    except StopIteration:
        end = len(lines)  # fallback if <<END>> not found


    # Slice the lines to desired range
    lines = lines[:end]


    # If it's a 2D file (e.g. table with multiple columns)
    if Separator != '0':
        Results = [[0] for _ in range(len(lines))]  # Create a placeholder matrix
        for i, line in enumerate(lines):
            Results[i] = line.split(Separator)  # Split each line by separator


        # Convert each value to float or int
        for j in range(len(Results)):
            for i in range(len(Results[0])):  # Assumes all lines have equal number of columns
                if Results[j][i] == '':
                    pass  # Skip empty entries
                elif Decimal == 'Yes':
                    Results[j][i] = float(Results[j][i])
                elif Decimal == 'No':
                    Results[j][i] = int(Results[j][i])
        return Results

    else:
        # Handle 1D vector (no separator between values)
        if Decimal == 'String':
            return lines  # Return raw strings
        else:
            for i in range(len(lines)):
                if Decimal == 'Yes':
                    lines[i] = float(lines[i])
                elif Decimal == 'No':
                    lines[i] = int(lines[i])
            return lines

###########################################################################
# Returns the number of significant figures based on an uncertainty value #
###########################################################################
def Precision(value):
    """
    Calculate the number of decimal places required to represent the first significant digit
    of based on an uncertainty value.

    Parameters:
        value (float or str): The uncertainty value.

    Returns:
        int: Number of decimal places needed for the first significant digit.
             Returns 0 if the value is zero or invalid.    
    """
    try:
        number = abs(float(value))
        if number == 0:
            return 0
        # Get exponent of the first significant digit
        exponent = int(math.floor(math.log10(number)))
        # If exponent >= 0, no decimal places needed
        if exponent >= 0:
            return 0
        else:
            return abs(exponent)
    except Exception:
        return 0

####################################################################
# Reads and displays the final results on the first tab of the GUI #
####################################################################
def Final_Results(tracker):
    """
     Function: Final_Results
     ------------------------
     Reads and displays the final results on the first tab of the GUI.

     The behavior depends on:
       - The current method (e.g., 'ROI Select') chosen for the analysis.
       - The tracker value, which indicates whether calibration (-), material thickness (+),
         or idle (0) operations are being requested.

     If calibration: calls Linearize or LinearizeWithErrors depending on the method.
     If material thickness: calls Calculate_Thickness (handles both standard and ROI Select).

     For each tab (loop over TabTracker):
       - Reads result files and displays values using tkinter labels.
       - If calibration: shows slope/intersect with uncertainties.
       - If material: shows thickness per peak and average with uncertainty.
    """
    num = Current_Tab()
    method = TabList[num][1].Algorithm_Method.get()

    # Decide which analysis to apply
    if tracker < 0:
        LinearRegression()
    elif tracker > 0:
        Calculate_Thickness()
    elif tracker == 0:
        pass

    # Loop through all tabs to update results
    for i in range(len(TabTracker)):
        # CALIBRATION RESULTS DISPLAY
        if (tracker <= 0) and TabTracker[i] < 0:
            if os.path.isfile(TabList[i][4]):
                # Determine energy units
                energy_val = TabList[i][1].energy.get()
                unit_energy = 'MeV' if energy_val == 1000 else 'keV'

                # Read calibration results
                Results = File_Reader(TabList[i][4], '0', 'String', 'No')

                # Display calibration trial info
                tk.Label(tab_manager.Calib_Result2, text=f'Calibration Trial {-TabTracker[i]} - {TabList[i][1].Source.get()}').grid(row=4*i+i, columnspan=3)
                tk.Label(tab_manager.Calib_Result2, text=f'({unit_energy})').grid(row=4*i+i+1, column=0)
                tk.Label(tab_manager.Calib_Result2, text='Values').grid(row=4*i+i+1, column=1)
                tk.Label(tab_manager.Calib_Result2, text='Uncertainty').grid(row=4*i+i+1, column=2)
                tk.Label(tab_manager.Calib_Result2, text='Slope').grid(row=4*i+i+2, column=0)
                tk.Label(tab_manager.Calib_Result2, text='Intersect').grid(row=4*i+i+3, column=0)
                tk.Label(tab_manager.Calib_Result2, text='').grid(row=4*i+i+4, column=0)

                # Display slope and intercept with uncertainty
                tk.Label(tab_manager.Calib_Result2, text='%.*f' % (int(Results[5]), float(Results[1]))).grid(row=4*i+i+2, column=1)
                tk.Label(tab_manager.Calib_Result2, text='%.*f' % (int(Results[5]), float(Results[2]))).grid(row=4*i+i+2, column=2)
                tk.Label(tab_manager.Calib_Result2, text='%.*f' % (int(Results[6]), float(Results[3]))).grid(row=4*i+i+3, column=1)
                tk.Label(tab_manager.Calib_Result2, text='%.*f' % (int(Results[6]), float(Results[4]))).grid(row=4*i+i+3, column=2)

            else:
                # Update canvas if no result file
                tab_manager.calib_canvas.update_idletasks()
                tab_manager.calib_canvas.config(scrollregion=tab_manager.Calib_Result2.bbox())
                tab_manager.Calib_Result2.bind('<Configure>',
                    lambda e: tab_manager.calib_canvas.configure(
                        scrollregion=tab_manager.calib_canvas.bbox('all'), width=e.width))

        # MATERIAL RESULTS DISPLAY
        if (tracker >= 0) and TabTracker[i] > 0:
            if os.path.isfile(TabList[i][4]) and os.path.isfile(TabList[i][3]):
                # Read result and peak files
                Results = File_Reader(TabList[i][4], '0', 'Yes', 'No')
                Peaks = File_Reader(TabList[i][3], ',', 'Yes', 'No')
                Peaks.sort()

                # Units and conversion setup
                units_list = ['nm', 'μm', 'μg cm⁻²', '10¹⁵ Atoms cm⁻³']
                units_values = [1e9, 1e6, 0.0, -1.0]
                index = units_values.index(TabList[i][1].units.get())

                size = len(Peaks)

                for j in range(size):
                    if j == 0:
                        # Trial header and column labels
                        tk.Label(tab_manager.Mat_Result2, text=f'Material Trial {TabTracker[i]} - {TabList[i][1].Mat.get()}').grid(row=(4+size)*i+i+j, columnspan=2)
                        tk.Label(tab_manager.Mat_Result2, text='Peak Centroid').grid(row=(4+size)*i+i+j+1, column=0)
                        tk.Label(tab_manager.Mat_Result2, text='Thickness').grid(row=(4+size)*i+i+j+1, column=1)

                    # Peak energy and thickness
                    tk.Label(tab_manager.Mat_Result2, text=f'{Peaks[j][0]:.1f}').grid(row=(4+size)*i+i+j+2, column=0)
                    tk.Label(tab_manager.Mat_Result2, text='%.*f %s' % (int(Results[-1]), Results[j], units_list[index])).grid(row=(4+size)*i+i+j+2, column=1)

                # Display average and uncertainty
                tk.Label(tab_manager.Mat_Result2, text='\nAverage').grid(row=(4+size)*i+i+j+3, column=0)
                tk.Label(tab_manager.Mat_Result2, text='Uncertainty').grid(row=(4+size)*i+i+j+4, column=0)
                tk.Label(tab_manager.Mat_Result2, text='\n%.*f %s' % (int(Results[-1]), Results[j+1], units_list[index])).grid(row=(4+size)*i+i+j+3, column=1)
                tk.Label(tab_manager.Mat_Result2, text='%.*f %s' % (int(Results[-1]), Results[j+2], units_list[index])).grid(row=(4+size)*i+i+j+4, column=1)
                tk.Label(tab_manager.Mat_Result2, text='').grid(row=(4+size)*i+i+j+5, columnspan=2)

                # Update scrollable canvas
                tab_manager.mat_canvas.update_idletasks()
                tab_manager.mat_canvas.config(scrollregion=tab_manager.Mat_Result2.bbox())
                tab_manager.Mat_Result2.bind('<Configure>',
                    lambda e: tab_manager.mat_canvas.configure(
                        scrollregion=tab_manager.mat_canvas.bbox('all'), width=e.width))
    return


###############################################################################


###############################################################################
# Allows the user to choose which linear regression(s) to use for calibration #
###############################################################################
def Calib_Choice():
    """
    Opens a popup window allowing the user to select which linear regression(s) to use for calibration.

    This function scans all calibration tabs (where TabTracker < 0) and checks for the existence of valid
    regression result files. For each valid regression, a checkbox is presented to the user. The user can
    select one or more regressions to be used in subsequent calculations (e.g., thickness determination).
    If no valid regressions are found, a warning popup is shown.

    Behavior:
        - Resets all previous regression selections for the current tab.
        - Lists all calibration tabs with valid regression files.
        - If none are found, displays a warning and instructions.
        - Otherwise, displays a selection menu with checkboxes for each regression.
        - User selections are stored in Regression_List variables for the current tab.

    Dependencies:
        - Uses global TabList and TabTracker structures.
        - Relies on tkinter for GUI elements and warnings_manager for popup management.

    Returns:
        None
    """
    num = Current_Tab()    # Get the index of the currently active tab
    validCalib = []        # List of calibration tabs with valid regressions
    validCalib_index = []  # Corresponding index in TabList for each valid regression
    regression_value_map = {}

    # Reset all previous selections (checkboxes set to -1 = unchecked)
    for i in range(len(TabList[num][1].Regression_List)):
        TabList[num][1].Regression_List[i].set(-1)

    # Identify all calibration tabs (TabTracker < 0) with valid regression result files
    for i in range(len(TabTracker)):
        if TabTracker[i] < 0 and os.path.isfile(TabList[i][4]):
            validCalib.append(TabTracker[i])    # e.g., -1, -2, ...
            validCalib_index.append(i)          # Store actual tab index for reference

    # If no valid calibration regressions found, show warning popup
    if not validCalib:
        warnings_manager.popup('No Linear Regressions detected')
        tk.Label(warnings_manager.warning, text=(
            'No linear regressions were detected.\n\n'
            'Please perform a Calibration Trial before calculating the film\'s thickness.\n\n')).pack()
        tk.Button(warnings_manager.warning, text='Return', command=lambda: warnings_manager.warning.destroy()).pack()

    # Otherwise, show selection popup for user to choose regressions
    else:
        warnings_manager.popup('Linear Regression Selection Menu')
        tk.Label(warnings_manager.warning, text=(
            'Please select one or more calibration trials.\n'
            'Choosing multiple calibrations will average the slopes and intercepts.\n\n')).pack()

        for i in range(len(validCalib)):
            tab_idx = validCalib_index[i]  # Actual index in TabList and TabTracker
            if i >= len(TabList[num][1].Regression_List):
                print(f"Warning: More regressions than Regression_List slots (i={i})")
                break
            # Create a checkbox for each valid calibration regression
            # Use a unique positive value for onvalue
            on_value = i + 1
            regression_value_map[on_value] = TabTracker[tab_idx]
            button_Choice = tk.Checkbutton(
                            warnings_manager.warning,
                            text=f'Linear Regression of Calibration Trial {-validCalib[i]}',
                            variable=TabList[num][1].Regression_List[i],  # Correct indexing
                            onvalue=on_value,
                            offvalue=-1)
            TabList[num][1].Regression_List[i].set(-1)  # Ensure unchecked by default
            button_Choice.pack()
        
        # Store the mapping for later use (e.g., as an attribute of the tab or globally)
        TabList[num][1].regression_value_map = regression_value_map

        # Button to simply close the popup (could be extended with a Confirm button)
        tk.Button(warnings_manager.warning, text='Return', command=lambda: ClearWidget('Popup', 0)).pack()

#####################################################################################
# Calculate the material thickness for each detected peak and the average thickness #
# based on selected calibration regressions and stopping power data                 #
#####################################################################################
def Calculate_Thickness():
    """
    Calculates the material thickness for each detected peak and the average thickness
    for the current material trial tab, using the selected calibration regression(s)
    and stopping power data.

    This function supports both the standard and 'ROI Select' analysis methods:
      - For 'ROI Select', it uses a single selected calibration and calculates thickness
        based on energy loss between calibration and film peaks.
      - For the standard method, it allows averaging over multiple selected calibration
        regressions, computes the average slope/intercept, and determines thickness for
        each detected peak.

    The results (thickness per peak, average, and uncertainty) are written to a results file
    and displayed in the GUI.

    Steps:
        1. Clears previous thickness results from the GUI.
        2. Determines the analysis method and retrieves relevant calibration(s).
        3. Loads material stopping power data and peak information.
        4. Calculates thickness for each peak and the average/uncertainty.
        5. Updates the results file and displays results in the GUI.

    Dependencies:
        - Uses global TabList, TabTracker, and tkinter for GUI elements.
        - Relies on helper functions: File_Reader, Eloss, Thickness, Precision, ClearWidget.

    Returns:
        None
    """
    num = Current_Tab()
    method = TabList[num][1].Algorithm_Method.get()
    ClearWidget('Thickness', 0)
    TabList[num][1].ThicknessFrame.grid(row=5, columnspan=3, pady=5)

    # Units and conversion setup
    units_list = ['nm', '\u03bcm',
                  '\u03bcg cm\u207B\u00B2',
                  '10\u00B9\u2075 Atoms cm\u207B\u00B3']
    units_values = [1e9, 1e6, 0.0, -1.0]
    index = units_values.index(TabList[num][1].units.get())

    if method == 'ROI Select':
        # --- ROI Select logic ---
        calibration = TabList[num][1].Regression_List[0].get()  # vamos selecionar apenas uma calibração

        # Get the selected material and stop. pow.
        Material_choice = TabList[num][1].Mat.get()
        Material_choice = 'Files\\Materials\\' + Material_choice + '.txt'
        material_data = File_Reader(Material_choice, '|', 'Yes', 'No')

        # Get energies of selected source
        energies = [TabList[0][1].DecayList[k].get() for k in range(len(TabList[0][1].DecayList))]
        energies = [e for e in energies if e != -1.0]  # Remove unselected energies

        # Get slope and intercept of the selected calib
        calib_params = File_Reader(TabList[0][4], '0', 'String', 'No')
        # Check energy units to match the stop. pow. ones
        if calib_params[0] == 'keV':
            m = str(float(calib_params[1]) * 1000)  # slope
            dm = str(float(calib_params[2]) * 1000)  # slope uncertainty
        else:
            m = str(float(calib_params[1]) * 1000)
            dm = str(float(calib_params[2]) * 1000)

        # Get calibration centroids and errors
        calibCents = [File_Reader(TabList[0][3], ',', 'Yes', 'No')[k][0] for k in range(len(File_Reader(TabList[0][3], ',', 'Yes', 'No')))]
        calibErr = [File_Reader(TabList[0][3], ',', 'Yes', 'No')[k][1] for k in range(len(File_Reader(TabList[0][3], ',', 'Yes', 'No')))]
        filmCents = [File_Reader(TabList[num][3], ',', 'Yes', 'No')[k][0] for k in range(len(File_Reader(TabList[num][3], ',', 'Yes', 'No')))]
        filmErr = [File_Reader(TabList[num][3], ',', 'Yes', 'No')[k][1] for k in range(len(File_Reader(TabList[num][3], ',', 'Yes', 'No')))]

        # Calculate energy loss and uncertainty, returning min and max energy of alphas after crossing the film
        Emin, Emax, eloss = Eloss(energies, calibCents, filmCents, calibErr, filmErr, m, dm)

        # Get selected material stop. pow.
        Material_choice = TabList[num][1].Mat.get()
        Material_choice = 'Files\\Materials\\' + Material_choice + '.txt'
        material_data = File_Reader(Material_choice, '|', 'Yes', 'No')

        # Calculate thickness from energy loss for each peak
        thickPeak = Thickness(energies, Emin, Emax, material_data)
        meanThick = np.mean([thick for thick in thickPeak])
        stdDevThick = np.std(thickPeak)

        # Write results to file
        with open(TabList[num][4], 'w') as result_file:
            for val in thickPeak:
                result_file.write(str("{:.0f}".format(val)) + '\n')
            result_file.write(str(meanThick) + '\n')
            result_file.write(str(stdDevThick) + '\n')
            result_file.write(str(0))

        # Display results in the frame
        # header
        tk.Label(TabList[num][1].ThicknessFrame, text='Peak Energy\n(MeV)').grid(row=0, column=0)
        tk.Label(TabList[num][1].ThicknessFrame, text=' ').grid(row=0, column=1)  # spacer
        tk.Label(TabList[num][1].ThicknessFrame, text='Channel\n').grid(row=0, column=2)
        tk.Label(TabList[num][1].ThicknessFrame, text=' ').grid(row=0, column=3)  # spacer
        tk.Label(TabList[num][1].ThicknessFrame, text='Eloss\n(keV)').grid(row=0, column=4)
        tk.Label(TabList[num][1].ThicknessFrame, text=' ').grid(row=0, column=5)  # spacer
        tk.Label(TabList[num][1].ThicknessFrame, text='Thickness\n(' + units_list[index] + ')').grid(row=0, column=6)
        # peak info
        for k in range(len(eloss)):
            tk.Label(TabList[num][1].ThicknessFrame, text='%.3f' % (energies[k])).grid(row=k + 1, column=0)
            tk.Label(TabList[num][1].ThicknessFrame, text=' ').grid(row=0, column=1)  # spacer
            tk.Label(TabList[num][1].ThicknessFrame, text='%.1f' % (filmCents[k])).grid(row=k + 1, column=2)
            tk.Label(TabList[num][1].ThicknessFrame, text=' ').grid(row=0, column=3)  # spacer
            tk.Label(TabList[num][1].ThicknessFrame, text='%.0f' % (eloss[k])).grid(row=k + 1, column=4)
            tk.Label(TabList[num][1].ThicknessFrame, text=' ').grid(row=0, column=5)  # spacer
            tk.Label(TabList[num][1].ThicknessFrame, text='%.0f' % (thickPeak[k])).grid(row=k + 1, column=6)
        # thickness final result
        tk.Label(TabList[num][1].ThicknessFrame,
                 text='Average Thickness (' + units_list[index] + ')').grid(row=k + 2, column=0)
        tk.Label(TabList[num][1].ThicknessFrame,
                 text='Uncertainty (' + units_list[index] + ')').grid(row=k + 2, column=4)
        tk.Label(TabList[num][1].ThicknessFrame,
                 text="{:.0f}".format(meanThick)).grid(row=k + 3, column=0)
        tk.Label(TabList[num][1].ThicknessFrame,
                 text="{:.0f}".format(stdDevThick)).grid(row=k + 3, column=4)
        tk.Button(TabList[num][1].ThicknessFrame,
                  command=lambda: ClearWidget('Thickness', 1),
                  text='Reset Results').grid(row=k + 4, columnspan=2)
    else:
        # --- Standard (non-ROI) logic ---
        # Get the index of the currently selected units from the units_values list
        try:
            index = units_values.index(TabList[num][1].units.get())
        except ValueError:
            print("Error: Selected unit not found in units_values list.")
            return

        # Determine the material file path based on user selection
        material_filename = 'Files/Materials/' + TabList[num][1].Mat.get() + '.txt'
        material_data = File_Reader(material_filename, '|', 'Yes', 'No')

        # Initialize accumulators for slope and intercept from selected regressions
        slope_sum = 0.0
        intercept_sum = 0.0
        selectedCalibs = []
        regression_value_map = getattr(TabList[num][1], 'regression_value_map', {})

        # Build the list of indices from TabTracker corresponding to selected regressions
        for i, reg_var in enumerate(TabList[num][1].Regression_List):
            val = reg_var.get()
            if val != -1 and val in regression_value_map:
                try:
                    idx = TabTracker.index(regression_value_map[val])
                    selectedCalibs.append(idx)
                except ValueError:
                    print(f"Warning: Regression selection value {regression_value_map[val]} not found in TabTracker")

        if not selectedCalibs:
            print("No valid regression selected. Aborting calculation.")
            return

        # Assuming all selected calibrations share the same DecayList source, get channels used for stopping power range
        selectedEnergies = []
        first_reg_idx = selectedCalibs[0]
        if first_reg_idx >= len(TabList):
            print(f"Error: First regression index {first_reg_idx} out of range for TabList.")
            return

        for j in range(len(TabList[first_reg_idx][1].DecayList)):
            decays = TabList[first_reg_idx][1].DecayList[j].get()
            if decays != -1:
                selectedEnergies.append(decays)

        selectedEnergies.sort()
        nrPeaks = len(selectedEnergies)

        # Sum slopes and intercepts from each selected regression file to compute average later
        for reg_idx in selectedCalibs:
            if reg_idx >= len(TabList):
                print(f"Warning: Regression index {reg_idx} out of range for TabList, skipping.")
                continue
            reg_file = TabList[reg_idx][4]
            reg_data = File_Reader(reg_file, '0', 'String', 'No')
            if reg_data[0] == 'keV':
                reg_data[1] = str(float(reg_data[1]) * 1000)
                reg_data[3] = str(float(reg_data[3]) * 1000)
            slope_sum += float(reg_data[1])
            intercept_sum += float(reg_data[3])

        if len(selectedCalibs) == 0:
            print("No regressions to average.")
            return

        slope_avg = slope_sum / len(selectedCalibs)
        intercept_avg = intercept_sum / len(selectedCalibs)

        # Read the peaks data (material points), sort ascending by channel or energy
        points = File_Reader(TabList[num][3], ',', 'No', 'No')
        points.sort()
        thickness_list = []

        # Setup labels for GUI
        tk.Label(TabList[num][1].ThicknessFrame, text=f'Thickness ({units_list[index]})').grid(row=0, column=0)
        tk.Label(TabList[num][1].ThicknessFrame, text='Channel').grid(row=0, column=1)

        for j in range(nrPeaks):
            # Calibrate each point with average regression slope and intercept
            calibrated_energy = (slope_avg * points[j][0]) + intercept_avg
            summed_stopping_power = 0.0

            # Sum inverse stopping power over energy interval [calibrated_energy, selectedEnergies[j]]
            for i in range(1, len(material_data)):
                energy = material_data[i][0]
                stopping = material_data[i][1]
                if calibrated_energy <= energy <= selectedEnergies[j]:
                    if stopping != 0:
                        summed_stopping_power += (1 / stopping)
                    else:
                        print(f"Warning: Stopping power zero at energy {energy}")

            # Unit conversions according to selected unit
            if index == 0:
                summed_stopping_power = summed_stopping_power / material_data[0][1] * 10000
            elif index == 1:
                summed_stopping_power = summed_stopping_power / material_data[0][1] * 10
            elif index == 2:
                summed_stopping_power *= 0.001
            elif index == 3:
                summed_stopping_power *= 1000 * (6.02214076e-23 / material_data[0][0])

            thickness_list.append(summed_stopping_power)

            # Display thickness and channel values in GUI
            tk.Label(TabList[num][1].ThicknessFrame, text=f'{summed_stopping_power:.2f}').grid(row=j+1, column=0)
            tk.Label(TabList[num][1].ThicknessFrame, text=str(points[j][0])).grid(row=j+1, column=1)

        # Calculate average thickness and uncertainty
        avg_thickness = sum(thickness_list) / len(thickness_list)
        uncertainty_sum = sum((x - avg_thickness) ** 2 for x in thickness_list) / (len(thickness_list) - 1)
        uncertainty = math.sqrt(uncertainty_sum)
        sig_fig = Precision(f'{uncertainty:.2g}')

        # Save results to file
        with open(TabList[num][4], 'w') as my_file:
            for val in thickness_list:
                my_file.write(f'{val:.3f}\n')
            my_file.write(f'{avg_thickness}\n')
            my_file.write(f'{uncertainty}\n')
            my_file.write(str(sig_fig))

        # Display average thickness and uncertainty in GUI
        tk.Label(TabList[num][1].ThicknessFrame, text=f'Average Thickness ({units_list[index]})').grid(row=nrPeaks+2, column=0)
        tk.Label(TabList[num][1].ThicknessFrame, text=f'Uncertainty ({units_list[index]})').grid(row=nrPeaks+2, column=1)
        tk.Label(TabList[num][1].ThicknessFrame, text=f'{avg_thickness:.{sig_fig}f}').grid(row=nrPeaks+3, column=0)
        tk.Label(TabList[num][1].ThicknessFrame, text=f'{uncertainty:.{sig_fig}f}').grid(row=nrPeaks+3, column=1)
        tk.Button(TabList[num][1].ThicknessFrame, command=lambda: ClearWidget('Thickness', 1),
                  text='Reset Results').grid(row=nrPeaks+4, columnspan=2)
    
    return

#########################################
# Displays algorithm results in the GUI #
#########################################
def ResultManager():
    """
    Displays the results of the selected analysis algorithm in the Results frame of the current tab.

    This function handles both standard and 'ROI Select' analysis methods:
      - For 'ROI Select', it displays a list of detected peak centroids, their uncertainties (σ), and σ/√N,
        each with a selectable checkbox.
      - For standard methods, it displays a list of detected channels and their counts, each with a selectable checkbox.

    For each result, a Tkinter Checkbutton is created to allow the user to select or deselect individual results.
    The function also resets the Results frame before displaying new results.

    Dependencies:
        - Uses global TabList and tkinter for GUI elements.
        - Relies on the File_Reader helper function to load results from file.

    Returns:
        None
    """
    num = Current_Tab()

    ClearWidget('Results', 0)  # Reset previous results
    TabList[num][1].ResultFrame.grid(row=3, columnspan=2, pady=5)

    method = TabList[num][1].Algorithm_Method.get()

    # ROI Select: expects centroids, sigma/sqrt(N), sigma
    if method == 'ROI Select':
        values = File_Reader(TabList[num][3], ',', 'Yes', 'No')
        for j in range(len(values)):
            Result_Button = tk.Checkbutton(
                TabList[num][1].ResultFrame,
                variable=TabList[num][1].Var_Data[j],
                onvalue=1, offvalue=-1,
                text='Centroid: ' + str("{:.1f}".format(values[j][0]))
            )
            Result_Button.grid(row=j, column=0)
            Result_Button.select()
            tk.Label(
                TabList[num][1].ResultFrame,
                text='\t \u03C3 =  ' + str("{:.1f}".format(values[j][2]))
            ).grid(row=j, column=1)
            tk.Label(
                TabList[num][1].ResultFrame,
                text='\t \u03C3/\u221aN =  ' + str("{:.3f}".format(values[j][1]))
            ).grid(row=j, column=2)
    # Standard: expects channel and counts
    else:
        values = File_Reader(TabList[num][3], ',', 'No', 'No')
        for j in range(len(values)):
            Result_Button = tk.Checkbutton(
                TabList[num][1].ResultFrame,
                variable=TabList[num][1].Var_Data[j],
                onvalue=1, offvalue=-1,
                text='Channel: ' + str(values[j][0])
            )
            Result_Button.grid(row=j, column=0)
            Result_Button.select()
            tk.Label(
                TabList[num][1].ResultFrame,
                text='\t Counts: ' + str(values[j][1])
            ).grid(row=j, column=1)
    return
        
#####################################################################
# Remove unchecked results from the GUI and update the results file #
#####################################################################       
def Unchecked_Results():
    """
    Unchecked_Results: Remove unchecked results from the GUI and update the results file
    
    This function:
    - Scans all result widgets (checkbuttons and labels)
    - Removes GUI elements corresponding to unchecked channels
    - Updates the results file with only the checked channels and counts
    """
    num = Current_Tab()  # Get the current tab index

    Eraser = []  # List to store widgets (checkbuttons and labels) to be potentially destroyed
    Aux = []     # List to store checked lines for rewriting the results file

    # Collect all widgets (both checkbuttons and labels) from the ResultFrame
    for widget in TabList[num][1].ResultFrame.winfo_children():
        Eraser.append(widget)

    # Read the results file as raw lines (strings) containing channel and counts
    values = File_Reader(TabList[num][3], '0', 'String', 'No')

    j = 0  # Index for Eraser list (widgets)
    k = 0  # Index for values list (file lines)

    for i in range(len(TabList[num][1].Var_Data)):

        var_value = TabList[num][1].Var_Data[i].get()

        if var_value == 1:
            # If checked, keep this line for rewriting the file
            Aux.append(values[k])

        elif var_value == -1:
            # If unchecked (-1), remove widgets and reset variable to 0
            TabList[num][1].Var_Data[i].set(0)

            # Destroy checkbutton and label widgets (assumes pairs)
            Eraser[j].destroy()
            Eraser[j + 1].destroy()

        elif var_value == 0:
            # If variable is already 0, adjust indices to keep synchronization
            k -= 1
            j -= 2

        # Increment widget index by 2 (checkbutton + label)
        j += 2
        # Increment values line index by 1
        k += 1

    # Move all variables with value 0 to the end of Var_Data list to keep selected variables grouped
    # Do this after main loop to avoid index shifting issues during iteration
    for i in range(len(TabList[num][1].Var_Data)):
        if TabList[num][1].Var_Data[i].get() == 0:
            TabList[num][1].Var_Data.append(TabList[num][1].Var_Data[i])
            TabList[num][1].Var_Data.pop(i)

    # Rewrite the results file with only checked (selected) entries
    with open(TabList[num][3], "w") as file:
        for line in Aux:
            file.write(line + '\n')
    
    return

#####################################################################################
# Performs linear regression (with or without errors) on selected data and displays #
# results in the GUI.                                                               #
#####################################################################################
def LinearRegression():
    """
    Performs linear regression on selected data and displays results in the GUI.

    This function automatically detects whether to use errors (uncertainties) in the regression
    based on the current analysis method ('ROI Select' uses errors, others do not).

    - For standard methods:
        * Reads channel (x-axis) data from results file.
        * Collects selected decay values (y-axis) from GUI.
        * Checks if x and y data match in length.
        * Calculates slope (m), intercept (b), and their uncertainties (sigma_m, sigma_b).
        * Converts units if needed (MeV or keV).
        * Saves regression results to a file.
        * Displays the results in the GUI with options to clear.

    - For 'ROI Select':
        * Reads centroids and errors from results file.
        * Collects selected decay energies from GUI.
        * Checks if centroids and energies match in length.
        * Performs weighted linear regression using uncertainties.
        * Converts units if needed (MeV or keV).
        * Saves regression results to a file.
        * Displays the results in the GUI with options to clear.

    Dependencies:
        - Uses global TabList, TabTracker, and tkinter for GUI elements.
        - Relies on helper functions: File_Reader, Calib, Precision, ClearWidget.

    Returns:
        None
    """
    num = Current_Tab()
    method = TabList[num][1].Algorithm_Method.get()

    # --- ROI Select: Use errors in regression ---
    if method == 'ROI Select':
        centroids = []
        errors = []
        energies = []
        # Read centroids and errors from results file
        values = File_Reader(TabList[num][3], ',', 'Yes', 'No')
        for i in range(len(values)):
            centroids.append(values[i][0])
            errors.append(values[i][1])
        # Collect selected decay energies from GUI
        for i in range(len(TabList[num][1].DecayList)):
            if TabList[num][1].DecayList[i].get() != -1:
                energies.append(TabList[num][1].DecayList[i].get())
        # Sort for consistency
        centroids = sorted(centroids)
        errors = sorted(errors)
        energies = sorted(energies)
        # Check for matching lengths
        if len(centroids) != len(energies):
            warnings_manager.popup('Invalid Linear Regression Configuration')
            tk.Label(warnings_manager.warning, 
                     text="Number of Radiation Decay does not match the number of Peaks detected.\n").pack()
            tk.Label(warnings_manager.warning, 
                     text="Please adjust the Searching Algorithms or the number of Decay Energy.\n\n").pack()
            tk.Button(warnings_manager.warning, text='Return', command=lambda: warnings_manager.warning.destroy()).pack()
            return

        # Clear previous regression output and show regression frame
        ClearWidget('Linear', 0)
        TabList[num][1].LinearRegressionFrame.grid(row=3, columnspan=2, pady=5)

        # Perform weighted linear regression using Calib()
        m, b, sigma_m, sigma_b = Calib(energies, centroids, errors)

        # Handle units and significant digits
        if TabList[num][1].energy.get() == 1000:
            unit_string = 'MeV'
        elif TabList[num][1].energy.get() == 1:
            unit_string = 'keV'
            m = m * 1000
            sigma_m = sigma_m * 1000
            b = b * 1000
            sigma_b = sigma_b * 1000

        # Save regression results to file
        with open(TabList[num][4], 'w') as my_file:
            my_file.write(unit_string + '\n')
            my_file.write(str(m) + '\n')
            my_file.write(str(sigma_m) + '\n')
            my_file.write(str(b) + '\n')
            my_file.write(str(sigma_b) + '\n')
            my_file.write(str(6) + '\n')
            my_file.write(str(6))

        # Display results in GUI
        tk.Label(TabList[num][1].LinearRegressionFrame, text='(' + unit_string + ')').grid(row=0, column=0)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='Values').grid(row=0, column=1)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='Uncertainty').grid(row=0, column=2)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='Slope').grid(row=1, column=0)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='Intersect').grid(row=2, column=0)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='%.*f' % (6, m)).grid(row=1, column=1)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='%.*f' % (6, sigma_m)).grid(row=1, column=2)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='%.*f' % (6, b)).grid(row=2, column=1)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='%.*f' % (6, sigma_b)).grid(row=2, column=2)
        tk.Button(TabList[num][1].LinearRegressionFrame, text='Clear Regression', 
                  command=lambda: ClearWidget('Linear', 1)).grid(row=3, column=0, columnspan=3)

    # --- Standard: No errors in regression ---
    else:
        xaxis = []
        yaxis = []
        # Read peak data from results file: returns list of [channel, counts]
        values = File_Reader(TabList[num][3], ',', 'No', 'No')
        # Extract channel numbers as x-axis values
        for val in values:
            xaxis.append(val[0])
        # Extract selected decay energies as y-axis values
        for var in TabList[num][1].DecayList:
            if var.get() != -1:
                yaxis.append(var.get())
        # Sort both datasets for consistency
        xvalues = sorted(xaxis)
        yvalues = sorted(yaxis)
        # Ensure both datasets have the same length (necessary for linear regression)
        if len(xvalues) != len(yvalues):
            warnings_manager.popup('Invalid Linear Regression Configuration')
            tk.Label(warnings_manager.warning, 
                     text="Number of Radiation Decay does not match the number of Peaks detected.\n").pack()
            tk.Label(warnings_manager.warning, 
                     text="Please adjust the Searching Algorithms or the number of Decay Energy.\n\n").pack()
            tk.Button(warnings_manager.warning, text='Return', command=lambda: warnings_manager.warning.destroy()).pack()
            return

        # Clear previous regression output and show regression frame
        ClearWidget('Linear', 0)
        TabList[num][1].LinearRegressionFrame.grid(row=3, columnspan=2, pady=5)

        # Calculate averages of x and y
        avgx = np.sum(xvalues) / len(xvalues)
        avgy = np.sum(yvalues) / len(yvalues)

        # Calculate slope (m) and intercept (b)
        numerator = np.sum((x - avgx) * (y - avgy) for x, y in zip(xvalues, yvalues))
        denominator = np.sum((x - avgx) ** 2 for x in xvalues)
        m = numerator / denominator
        b = avgy - m * avgx

        # Calculate uncertainties (sigma_m and sigma_b)
        residual_sum = np.sum((y - m * x - b) ** 2 for x, y in zip(xvalues, yvalues))
        residual_variance = residual_sum / (len(xvalues) - 2)
        sum_x2 = np.sum(x ** 2 for x in xvalues)
        sum_x = np.sum(xvalues)
        n = len(xvalues)
        sigma_m = math.sqrt(residual_variance / (sum_x2 - (sum_x ** 2) / n))
        sigma_b = math.sqrt(residual_variance * sum_x2 / (n * sum_x2 - sum_x ** 2))

        # Handle units and significant digits depending on energy unit selected
        energy_unit = TabList[num][1].energy.get()
        if energy_unit == 1000:
            unit_string = 'MeV'
            significant_digits_m = Precision('%.2g' % sigma_m)
            significant_digits_b = Precision('%.2g' % sigma_b)
        elif energy_unit == 1:
            unit_string = 'keV'
            # Convert slope and intercept and their uncertainties to keV
            m *= 1000
            sigma_m *= 1000
            b *= 1000
            sigma_b *= 1000
            significant_digits_m = 0 if sigma_m > 1 else Precision('%.2g' % sigma_m)
            significant_digits_b = 0 if sigma_b > 1 else Precision('%.2g' % sigma_b)

        # Write regression results to file for use in other functions
        with open(TabList[num][4], 'w') as my_file:
            my_file.write(unit_string + '\n')
            my_file.write(str(m) + '\n')
            my_file.write(str(sigma_m) + '\n')
            my_file.write(str(b) + '\n')
            my_file.write(str(sigma_b) + '\n')
            my_file.write(str(significant_digits_m) + '\n')
            my_file.write(str(significant_digits_b))

        # Display results in GUI with proper formatting
        tk.Label(TabList[num][1].LinearRegressionFrame, text='(' + unit_string + ')').grid(row=0, column=0)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='Values').grid(row=0, column=1)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='Uncertainty').grid(row=0, column=2)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='Slope').grid(row=1, column=0)
        tk.Label(TabList[num][1].LinearRegressionFrame, text='Intercept').grid(row=2, column=0)
        tk.Label(TabList[num][1].LinearRegressionFrame, text=f'{m:.6f}').grid(row=1, column=1)
        tk.Label(TabList[num][1].LinearRegressionFrame, text=f'{sigma_m:.6f}').grid(row=1, column=2)
        tk.Label(TabList[num][1].LinearRegressionFrame, text=f'{b:.6f}').grid(row=2, column=1)
        tk.Label(TabList[num][1].LinearRegressionFrame, text=f'{sigma_b:.6f}').grid(row=2, column=2)
        tk.Button(TabList[num][1].LinearRegressionFrame, text='Clear Regression', 
                  command=lambda: ClearWidget('Linear', 1)).grid(row=3, column=0, columnspan=3)

    return

##########################################################################
# Manual aelection algorithm, finds the data point closest to the given  # 
# input point (Valuex, Valuey) by calculating the minimum Euclidean      # 
# distance between the input point and the points in the dataset         #
##########################################################################
def ManSelec_Alg(Valuex, Valuey):
    """
    Finds the data point closest to the given input point (Valuex, Valuey)
    by calculating the minimum Euclidean distance between the input point
    and the points in the dataset (where x = channel index, y = counts).
    
    Inputs:
      Valuex - the x-coordinate (channel) of the input point
      Valuey - the y-coordinate (count) of the input point
    
    The function reads the counts data from the current tab file, finds the
    closest data point, writes the result (channel,count) to an output file
    (appending or creating it), and then calls ResultManager() to update results.
    """
    num = Current_Tab()
    Counts = File_Reader(TabList[num][2], '0', 'No', 'No')

    # Function to calculate Euclidean distance between input point and data point at index i
    def dist(i):
        channel = i + 1  # Channels start at 1
        count = Counts[i]
        return math.sqrt((Valuex - channel)**2 + (Valuey - count)**2)

    # Find the index of the data point with the minimum distance to the input point
    min_index = min(range(len(Counts)), key=dist)

    # Update Valuex and Valuey with the closest data point found
    Valuex = min_index + 1
    Valuey = Counts[min_index]

    # Write the closest data point (channel,count) to the output file
    mode = 'a' if os.path.isfile(TabList[num][3]) else 'w'
    with open(TabList[num][3], mode) as results:
        results.write(f"{Valuex},{Valuey}\n")

    # Update results after selection
    ResultManager()

    return
    
################################################################
# Algorithm to detect and record peaks above a given threshold #
################################################################
def Threshold_Alg():
    """
    Algorithm to detect and record peaks above a given threshold.
    
    The function scans through the count data to identify peaks where
    counts exceed a specified threshold after a certain channel cut-off.
    For each peak found, it records the maximum count and the corresponding channel.
    It also merges peaks that are closer than a specified width, keeping only
    the highest peak in such cases.
    
    Results are saved to an output file (appending or writing new).
    """
    num = Current_Tab()
    Counts = File_Reader(TabList[num][2], '0', 'No', 'No')

    Threshold = TabList[num][1].Algorithm.get()  # Threshold value for peak detection
    yaxis = []               # Stores peak heights (max counts)
    xaxis = []               # Stores peak positions (channels)
    current_peak = []        # Temporary list to hold counts belonging to current peak
    current_peak_counts = 0  # Sum of counts in current peak (not used further here but can be for FWHM)
    counter = 0              # Tracks overall position in counts array
    i = 0                    # Index for looping over counts
    j = 0                    # Index for tracking peaks stored in xaxis/yaxis

    cut = TabList[num][1].channels_cut.get()  # Channel index cut-off, ignore counts below this channel
    width = TabList[num][1].peaks_widths.get() # Minimum width to consider two peaks separate

    # Iterate over all counts
    for count in Counts:
        # If count is above threshold and past the cut-off channel, accumulate peak data
        if count > Threshold and i > cut:
            current_peak_counts += count
            current_peak.append(count)
        else:
            counter += 1
            # If current_peak is not empty, we have reached the end of a peak
            if current_peak:
                # Record the maximum value of the current peak and its position
                max_peak = max(current_peak)
                peak_pos = counter + current_peak.index(max_peak)
                yaxis.append(max_peak)
                xaxis.append(peak_pos)

                # Merge peaks closer than 'width'
                if len(xaxis) > 1:
                    if xaxis[j] - xaxis[j-1] < width:
                        # Decide which peak to keep (the highest one)
                        decider = max(yaxis[j-1], yaxis[j])
                        if decider == yaxis[j]:
                            yaxis.pop(j-1)
                            xaxis.pop(j-1)
                        elif decider == yaxis[j-1]:
                            yaxis.pop(j)
                            xaxis.pop(j)
                        j -= 1

                # Reset for next peak
                counter = counter + len(current_peak)
                current_peak_counts = 0
                current_peak = []
                j += 1
        i += 1

    # Write detected peaks (channel, count) to output file
    mode = 'a' if os.path.isfile(TabList[num][3]) else 'w'
    with open(TabList[num][3], mode) as results:
        for i in range(len(xaxis)):
            results.write(f"{xaxis[i]},{yaxis[i]}\n")

    # Note: If the last peak extends to the end of the data, it may not be added.
    # This should be checked if needed.

    ResultManager()

################################################################################
# ROI_Select_Alg: Detects peaks within user-defined ROIs and records their     #
# centroids, uncertainties, and sigma/sqrt(N) for the ROI Select algorithm.    #
################################################################################
def ROI_Select_Alg():
    """
    Detects and analyzes peaks within user-defined Regions of Interest (ROIs)
    for the 'ROI Select' analysis method. For each ROI, calculates the centroid,
    uncertainty, and sigma/sqrt(N) using the Analyze() helper.

    Results are appended or written to the results file for the current tab,
    and the GUI is updated to display the detected peaks and their statistics.

    Steps:
        1. Reads the current tab's counts data.
        2. Retrieves lower and upper bounds for up to 6 ROIs from the GUI.
        3. Calls Analyze() to compute centroids, uncertainties, and sigma/sqrt(N).
        4. Writes results to the output file (appending if it exists, creating otherwise).
        5. Calls ResultManager() to update the results display in the GUI.

    Dependencies:
        - Uses global TabList and tkinter for GUI elements.
        - Relies on File_Reader and Analyze helper functions.

    Returns:
        None
    """
    num = Current_Tab()

    # Read counts data for the current tab
    counts = File_Reader(TabList[num][2], '0', 'Yes', 'No')

    # Retrieve ROI lower and upper bounds from the GUI (up to 6 ROIs)

    # Check if the current tab is XRA (TabTracker value == 5)

    if TabList[num][1].tab_kind == 5:

    # Only use ROI 1 for XRA
        roi_down = [TabList[num][1].ROIdown1.get()]
        roi_up = [TabList[num][1].ROIup1.get()]
        print("I am here")

    # Hide ROI 2–6 Entry widgets
        for i in range(1, min(6, len(TabList[num][1].ROIdown_entries))):
            try:
                TabList[num][1].ROIdown_entries[i].grid_remove()
                TabList[num][1].ROIup_entries[i].grid_remove()
                if hasattr(TabList[num][1], 'ROIlabels'):
                    TabList[num][1].ROIlabels[i].grid_remove()
                TabList[num][1].ROIdown_entries[i].delete(0, tk.END)
                TabList[num][1].ROIup_entries[i].delete(0, tk.END)
                getattr(TabList[num][1], f'ROIdown{i+1}').set(0)
                getattr(TabList[num][1], f'ROIup{i+1}').set(0)
            except IndexError:
                print(f"Skipping ROI entry {i} — not defined for this tab.")

    else:
        # Use all 6 if not XRA
        roi_down = [TabList[num][1].ROIdown1.get(),
                    TabList[num][1].ROIdown2.get(),
                    TabList[num][1].ROIdown3.get(),
                    TabList[num][1].ROIdown4.get(),
                    TabList[num][1].ROIdown5.get(),
                    TabList[num][1].ROIdown6.get()]
        roi_up = [TabList[num][1].ROIup1.get(),
                TabList[num][1].ROIup2.get(),
                TabList[num][1].ROIup3.get(),
                TabList[num][1].ROIup4.get(),
                TabList[num][1].ROIup5.get(),
                TabList[num][1].ROIup6.get()]

    print("TabTracker:", TabTracker[num])
    print("roi_down:", roi_down)
    print("roi_up:", roi_up)
    print("Number of ROIs passed to Analyze:", len(roi_down))

    for idx, (d, u) in enumerate(zip(roi_down, roi_up)):
        try:
            x1 = float(d)
            x2 = float(u)
        except (ValueError, TypeError):
            tk.messagebox.showerror("Error", f"Invalid ROI input at peak {idx+1}: \nROId = {d}, ROIu = {u}")
            return

        if x2 <= x1 or (x2 - x1) < 2:
            tk.messagebox.showerror("Error", f"Invalid ROI range at peak {idx+1}: \nx1 = {x1}, x2 = {x2} "
                                f"\nUpper bound must be greater than lower bound and difference must be at least 2.")
            return

    # Analyze the counts within each ROI to get centroids, uncertainties, and sigma/sqrt(N)
    cents, errs, sigmas = Analyze(counts, roi_down, roi_up)

    # Write results to file: append if file exists, otherwise create new
    if os.path.isfile(TabList[num][3]):
        with open(TabList[num][3], 'a') as results:
            for i in range(len(cents)):
                results.write(f"{cents[i]},{errs[i]},{sigmas[i]}\n")
    else:
        with open(TabList[num][3], 'w') as results:
            for i in range(len(cents)):
                results.write(f"{cents[i]},{errs[i]},{sigmas[i]}\n")

    # Update the results display in the GUI
    ResultManager()

    return
  
################################################################
# Displays the decay chain image for the selected alpha source #
################################################################
def showimage():
    """
    Displays the decay chain image for the currently selected alpha source.

    This function scans the 'Files/Sources/Images' directory for an image file
    whose name matches the selected alpha source. If a match is found, it constructs
    the image path and, for specific isotopes, assigns a relevant URL for more information.
    The function then attempts to clear any previous image widgets and display the
    selected image in a popup window, optionally with a hyperlink.

    Steps:
        1. Get the currently selected alpha source from the GUI.
        2. Scan the images directory for a matching file.
        3. If found, construct the image path and assign a URL if applicable.
        4. Clear any previous image widgets.
        5. Display the image and optional link in a popup window.

    Dependencies:
        - Uses global TabList and tkinter for GUI elements.
        - Relies on the warnings_manager.Images helper for displaying images and links.

    Returns:
        None
    """
    num = Current_Tab()
    alphas = TabList[num][1].Source.get()  # Get the selected alpha source name
    Dir = os.scandir('Files\Sources\Images')  # Scan the images directory

    for entry in Dir:
        if entry.is_file():
            temp = os.path.splitext(entry.name)  # Split filename and extension
            if alphas == temp[0]:  # Check if file name matches the selected source
                picture = 'Files\Sources\Images\\' + alphas + temp[1]  # Construct image path

                # Assign relevant URL for specific isotopes
                if alphas == '226Ra':
                    domain = 'https://www.nist.gov/image-23773'
                elif alphas == '232U':
                    domain = ('https://www.nuclear-power.com/nuclear-power-plant/'
                              'nuclear-fuel/uranium/uranium-232/decay-half-life-uranium-232/')
                else:
                    domain = ''  # No URL for other isotopes

                try:
                    ClearWidget('Image', 0)  # Clear previous image widget(s)
                    warnings_manager.Images('Decay Chain of ' + alphas, picture, domain)  # Display image and link
                except:
                    # If clearing fails, still try to display the image
                    warnings_manager.Images('Decay Chain of ' + alphas, picture, domain)

    return

############################################################################
# Event handler for capturing points directly from the plot on mouse click #
############################################################################
def onclick(event):
    """
    This function is triggered when the user clicks on the graph. It checks if
    the current algorithm method is set to 'Manual Selection' and the right mouse
    button (button 3) is clicked. If both conditions are met, it captures the
    x and y coordinates of the click event and passes them to the manual selection
    algorithm function (ManSelec_Alg).
    """
    num = Current_Tab()

    decider = TabList[num][1].Algorithm_Method.get()

    # Only process right-click events in 'Manual Selection' mode
    if decider == 'Manual Selection' and event.button == 3:
        xpoint = event.xdata  # Get x-coordinate from the click event
        ypoint = event.ydata  # Get y-coordinate from the click event
        ManSelec_Alg(xpoint, ypoint)  # Call function to process selected point
    
    return

##############################################################
# Opens a file dialog to upload a data file and processes it #
##############################################################
def DataUploader(): 
    """
    It restricts the selectable files to '.mca' text files and all files.
    After the user selects a file, it reads the data, processes the structure,
    plots the graph, and applies the selected algorithm (e.g., thresholding).
    The tab label is also updated with the filename.
    """
    num = Current_Tab()

    # Define the allowed file types for upload
    domain = (('Text Files', '*.mca'), ('All Files', '*.*'))

    tab_type = getattr(TabList[num][1], 'tab_kind', None)

    # XRA tab
    if tab_type == 5:
        # Ask XRA user: Source only or Source + Film?
        mode_window = tk.Toplevel()
        mode_window.title("XRA Plot Mode")

        tk.Label(mode_window, text="XRA tab detected.\nDo you want to plot:").pack(padx=10, pady=10)

        def upload_source_only():
            mode_window.destroy()
            filename = fd.askopenfilename(title='Open Source File', filetypes=[('Text Files', '*.mca'), ('All Files', '*.*')])
            if not filename:
                return

            file = File_Reader(filename, '0', 'string', 'Yes')
            TabList[num][5].source_data = file
            TabList[num][5].source_file = filename
            TabList[num][5].Structure(file, filename)
            TabList[num][5].subplots()
            Tabs.RenameTab(filename.split('/')[-1])

        def upload_film():
            mode_window.destroy()

            try:
                source_data = TabList[num][5].source_data
                source_file = TabList[num][5].source_file
            except AttributeError:
                tk.messagebox.showerror("Missing Source File", "Please upload a Source file first using 'Source Only'.")
                return

            film_file = fd.askopenfilename(title='Open Film File', filetypes=[('Text Files', '*.mca'), ('All Files', '*.*')])
            if not film_file:
                return

            film_data = File_Reader(film_file, '0', 'string', 'No')

            TabList[num][5].Structure(source_data, source_file, film_data)
            TabList[num][5].subplots()
            Tabs.RenameTab(film_file.split('/')[-1])

        tk.Button(mode_window, text="Source Only", command=upload_source_only).pack(padx=10, pady=5)
        tk.Button(mode_window, text="Source + Film", command=upload_film).pack(padx=10, pady=5)

    # All other tabs
    else: 

        # Open a file dialog window for user to select the data file
        filename = fd.askopenfilename(title='Open a file', initialdir=',', filetypes=domain) 

        # If no file was selected, do nothing and exit
        if not filename:
            pass

        else:
            # Read the uploaded file's contents using File_Reader function
            file = File_Reader(filename, '0', 'string', 'Yes')

            # Process the file structure and generate the initial plot
            TabList[num][5].Structure(file, filename)
            TabList[num][5].subplots()

            # Retrieve the current algorithm selection and apply it if needed
            value = TabList[num][1].Algorithm.get()
            if value != 0:
                TabList[num][5].threshold(value)

        # Rename the tab to the filename (only the file's name, not the full path)
        try:
            Tabs.RenameTab(filename.split('/')[-1])
        except:
            pass  # Ignore errors in renaming tab

        return

######################################################
# Handles the event of changing tabs in the notebook #
######################################################
def handleTabChange(event):
    """
    If the user selects the last tab (the '+' tab) and the number of tabs is less than 15,
    a popup menu appears to allow the user to select a new type of tab to open.
    
    If the number of tabs reaches 15 or more, the '+' tab is hidden to prevent adding more tabs.
    """

    # Check if the currently selected tab is the last one (the '+' tab)
    # and if the total number of tabs is less than 15
    if tab_manager.notebook.select() == tab_manager.notebook.tabs()[-1] and len(tab_manager.notebook.tabs()) < 15:
        # Show popup menu for tab selection
        warnings_manager.popup('Tab Selector Menu')

        tk.Label(warnings_manager.warning, text='Please Select New Type of Tab to Open \n').pack()
        tk.Button(warnings_manager.warning, command=lambda: Tabs.tab_change(1),
                  text='Calibration Trial').pack()
        tk.Button(warnings_manager.warning, command=lambda: Tabs.tab_change(2), 
                  text='Material Trial').pack()
        tk.Button(warnings_manager.warning, command=lambda: Tabs.tab_change(5), text='XRA').pack()
        tk.Label(warnings_manager.warning, text='\n').pack()
        tk.Button(warnings_manager.warning, text='Return',
                  command=lambda: Tabs.tab_change(3)).pack()
        # The numbers passed to tab_change() indicate which type of tab to add

    # If there are already 15 or more tabs, hide the '+' tab to prevent adding more
    elif len(tab_manager.notebook.tabs()) >= 15:
        tab_manager.notebook.hide(14)

    return

#############################################################
# Handles the setup and display of available alpha sources, #
# decay energies, and related controls for the current tab  #
#############################################################
def SourceReader(*args):
    """
    Sets up the source options for the current tab, including available alpha decay energies,
    checkbuttons for user selection, and controls for displaying the decay chain image and
    performing linear regression.

    Steps:
        1. Clears and resets the source options frame for the current tab.
        2. Loads the available alpha decay energies from the corresponding file.
        3. Creates checkbuttons for each decay energy, all selected by default.
        4. Adds a button to display the decay chain image for the selected source.
        5. Adds a button to perform linear regression and display results.

    Path Handling:
        - Uses os.path.join for all file and directory paths to ensure compatibility
          across different operating systems.

    Dependencies:
        - Uses global TabList, TabTracker, and tkinter for GUI elements.
        - Relies on the showimage and Final_Results helper functions.

    Returns:
        None
    """
    num = Current_Tab()
    value = TabTracker[num]

    ClearWidget('Source', 0)  # Reset widgets and layout geometry
    TabList[num][1].SourceOptionsFrame.grid(row=2, columnspan=2)

    # Reset all decay energies in DecayList to -1 to avoid interference with linear regression
    for i in range(len(TabList[num][1].DecayList)):
        TabList[num][1].DecayList[i].set(-1)

    # Construct path to the file containing alpha decay energies based on user selection
    Alpha_name = TabList[num][1].Source.get()
    Alpha_path = os.path.join('Files', 'Sources', 'Values', f'{Alpha_name}.txt')

    # Read all decay energies from the file
    with open(Alpha_path, 'r') as file:
        Decay = [float(line) for line in file]

    # Create checkbuttons for each decay energy for user selection
    for i in range(len(Decay)):
        checkbutton = tk.Checkbutton(TabList[num][1].SourceOptionsFrame,
                                     text=str(Decay[i]) + ' MeV',
                                     variable=TabList[num][1].DecayList[i],
                                     onvalue=Decay[i],
                                     offvalue=-1)
        checkbutton.grid(row=i, columnspan=2)
        checkbutton.select()  # By default, select all options

    # Button to show the decay chain image corresponding to the selected source
    tk.Button(TabList[num][1].SourceOptionsFrame,
              text='Show Decay Chain',
              command=showimage).grid(row=i + 1, column=0)

    # Button to perform linear regression and display results on the first tab
    tk.Button(TabList[num][1].SourceOptionsFrame,
              text='Linear Regression',
              command=lambda: Final_Results(value)).grid(row=i + 1, column=1)
    return

###################################################################
# Updates the algorithm input interface for the selected analysis #
# method in the current tab                                       #
###################################################################
def Method(*args):
    """
    Updates the algorithm input interface for the selected analysis method in the current tab.

    This function dynamically configures the controls and widgets in the algorithm frame
    based on the user's choice of analysis method:
      - Manual Selection: Provides instructions and buttons for manual peak selection.
      - Threshold Input: Provides entry for threshold value and buttons for peak detection.
      - ROI Select: Provides entry fields for up to 6 Regions of Interest (ROIs) and buttons
        for peak detection within those regions.

    Steps:
        1. Clears the previous algorithm UI and resets the algorithm frame.
        2. Determines the selected analysis method from the GUI.
        3. Sets up the appropriate controls and widgets for the chosen method.

    Notes:
        - Avoids unnecessary repetition by grouping similar widget creation.
        - Ensures only relevant controls are shown for the selected method.
        - Calls helper functions to clear results or trigger algorithm execution.

    Dependencies:
        - Uses global TabList and tkinter for GUI elements.
        - Relies on helper functions: ClearWidget, Unchecked_Results, ROI_Select_Alg, Threshold_Alg.

    Returns:
        None
    """
    num = Current_Tab()

    # Clear previous algorithm UI and reset the algorithm frame
    ClearWidget('Algorithm', 0)
    TabList[num][1].AlgFrame.grid(row=2, columnspan=2)

    decider = TabList[num][1].Algorithm_Method.get()  # Get selected algorithm method

    if decider == 'Manual Selection':
        # Manual Selection: Instructions and buttons for manual peak selection
        tk.Label(TabList[num][1].AlgFrame, 
                 text='Right-Click on/near the Peaks in the Graphic ').grid(row=2, columnspan=2)
        tk.Label(TabList[num][1].AlgFrame, 
                 text='For an automatic point detection: ').grid(row=3, columnspan=2)
        tk.Button(TabList[num][1].AlgFrame, text='Remove Unchecked', 
                  command=Unchecked_Results).grid(row=4, column=0)
        tk.Button(TabList[num][1].AlgFrame, text='Remove All',
                  command=lambda: ClearWidget('Results', 1)).grid(row=4, column=1)

        TabList[num][1].Algorithm.set(0)
        TabList[num][5].destroyer()

    elif decider == 'Threshold Input':
        # Threshold Input: Entry for threshold value and buttons for peak detection
        tk.Label(TabList[num][1].AlgFrame, 
                 text='Please input Threshold: ').grid(row=2, columnspan=3)
        tk.Entry(TabList[num][1].AlgFrame, textvariable=TabList[num][1].Algorithm, relief='sunken',
                 borderwidth=2).grid(row=3, columnspan=3)
        tk.Button(TabList[num][1].AlgFrame, text='Search', 
                  command=Threshold_Alg).grid(row=4, column=0)
        tk.Button(TabList[num][1].AlgFrame, text='Remove Unchecked',
                  command=Unchecked_Results).grid(row=4, column=1)
        tk.Button(TabList[num][1].AlgFrame, text='Remove All',
                  command=lambda: ClearWidget('Results', 1)).grid(row=4, column=2)

    elif decider == 'ROI Select':
        # ROI Select: Entry fields for up to 6 ROIs and buttons for peak detection
        tk.Label(TabList[num][1].AlgFrame, text='ROI Down: ').grid(row=2, column=0)
        tk.Label(TabList[num][1].AlgFrame, text='ROI Up: ').grid(row=2, column=1)

        # Create labels for each peak
        ## Check if we are in a XRA tab
        
        if TabList[num][1].tab_kind == 5:
            range_limit = 1
        else:
            range_limit = 6

        
        # Initialize storage
        TabList[num][1].ROIdown_entries = []
        TabList[num][1].ROIup_entries = []
        TabList[num][1].ROIlabels = []
        for idx in range(range_limit):
            down_entry = tk.Entry(TabList[num][1].AlgFrame, textvariable=getattr(TabList[num][1], f'ROIdown{idx+1}'))
            down_entry.grid(row=3+idx, column=0)
            TabList[num][1].ROIdown_entries.append(down_entry)

            up_entry = tk.Entry(TabList[num][1].AlgFrame, textvariable=getattr(TabList[num][1], f'ROIup{idx+1}'))
            up_entry.grid(row=3+idx, column=1)
            TabList[num][1].ROIup_entries.append(up_entry)

            label = tk.Label(TabList[num][1].AlgFrame, text=f'Peak {idx+1}')
            label.grid(row=3+idx, column=2)
            TabList[num][1].ROIlabels.append(label)

        # Create entry fields for each ROI lower and upper bound
        #for idx in range(6):
            """
            tk.Entry(TabList[num][1].AlgFrame, textvariable=getattr(TabList[num][1], f'ROIdown{idx+1}'),
                     relief='sunken', borderwidth=2).grid(row=3+idx, column=0)
            tk.Entry(TabList[num][1].AlgFrame, textvariable=getattr(TabList[num][1], f'ROIup{idx+1}'),
                     relief='sunken', borderwidth=2).grid(row=3+idx, column=1)
            """
        
        tk.Button(TabList[num][1].AlgFrame, text='Search', 
                  command=ROI_Select_Alg).grid(row=9, column=0)
        tk.Button(TabList[num][1].AlgFrame, text='Remove Unchecked',
                  command=Unchecked_Results).grid(row=9, column=1)
        tk.Button(TabList[num][1].AlgFrame, text='Remove All',
                  command=lambda: ClearWidget('Results', 1)).grid(row=9, column=2)
         
    return



#############################################################################
# Esta funcao muda uma linha de threshold, caso o utilizador escreva um numero
#############################################################################
def on_entry_change(*args):

    num = Current_Tab()

    try:
        value = TabList[num][1].Algorithm.get()
        TabList[num][5].threshold(value)
    except:
        pass

##############################################################################
# Apaga os ficheiros listados e recebidos, faz parte do File_Manager
#############################################################################
def Delete(deletes, names, directory, last):

    dir = os.getcwd()

    directory = dir + '\\' + directory + '\\'

    for i in range(0, len(deletes)):
        if deletes[i].get() == 1:
            os.remove(directory + names[i] + last)

    warnings_manager.warning.destroy()

##########################################################################
# Handles uploading and deleting of permanent data files used            #
# by the application, including alpha source value files, source images, #
# and material files                                                     #
##########################################################################
def File_Manager(Choice, Nature, Action):
    """
    Manages the upload and deletion of permanent data files used by the application.

    Parameters:
        Choice (str): Type of file to manage ('Source' or 'Material').
        Nature (int): For 'Source', 1 = value file (.txt), 0 = image file (.jpg/.jpeg/.png).
                      For 'Material', always 0.
        Action (int): 1 = Upload file, 0 = Delete file.

    Behavior:
        - For uploading, opens a file dialog for the user to select a file and copies it to the
          appropriate directory.
        - For deleting, shows a popup with checkboxes for each file in the relevant directory,
          allowing the user to select files to delete.
        - After any change, updates the source or material dropdown menus in all relevant tabs.

    Notes:
        - Uses os.path and os.scandir for directory/file handling.
        - Uses tkinter.filedialog for file selection and tkinter for popups.
        - Uses shutil.copy2 for file copying.
        - Updates global lists (source_list, materials_list) after changes.

    Returns:
        None
    """
    # Handle Source files (value files or images)
    if Choice == 'Source':
        # Value files (.txt)
        if Nature == 1:
            if Action == 1:
                # Upload: Open file dialog and copy selected file to source values directory
                filename = fd.askopenfilename(filetypes=(('Text Files', '*.txt'), ('All Files', '*.*')),
                                              title='Add Alpha Source Energy File')
                target_dir = os.path.join(os.getcwd(), 'Files', 'Sources', 'Values')
                if filename:
                    copy2(filename, target_dir, follow_symlinks=True)
            else:
                # Delete: Show popup with checkboxes for each file in source values directory
                domain = os.path.join('Files', 'Sources', 'Values')
                dir_entries = os.scandir(domain)
                warnings_manager.popup('Delete Alpha Source Files')
                name_list = []
                delete_vars = []
                for entry in dir_entries:
                    if entry.is_file():
                        name, _ = os.path.splitext(entry.name)
                        name_list.append(name)
                        delete_vars.append(tk.IntVar())
                tk.Label(warnings_manager.warning, text='Files available for deletion\n').pack()
                for i, name in enumerate(name_list):
                    tk.Checkbutton(warnings_manager.warning, text=name, variable=delete_vars[i],
                                   onvalue=1, offvalue=0).pack()
                tk.Button(warnings_manager.warning, command=lambda: Delete(delete_vars, name_list, domain, '.txt'),
                          text='Delete Files').pack()
                tk.Button(warnings_manager.warning, command=lambda: warnings_manager.warning.destroy(),
                          text='Return').pack()
            # Update source_list after changes
            Dir = os.scandir(os.path.join('Files', 'Sources', 'Values'))
            source_list.clear()
            for entry in Dir:
                if entry.is_file():
                    name, _ = os.path.splitext(entry.name)
                    source_list.append(name)
        # Image files (.jpg, .jpeg, .png)
        elif Nature == 0:
            if Action == 1:
                # Upload: Open file dialog and copy selected image to source images directory
                filename = fd.askopenfilename(filetypes=(('Image Files', '.jpg .jpeg .png'),
                                                         ('All Files', '*.*')),
                                              title='Add Alpha Source Energy Decay Image')
                target_dir = os.path.join(os.getcwd(), 'Files', 'Sources', 'Images')
                if filename:
                    copy2(filename, target_dir, follow_symlinks=True)
            else:
                # Delete: Show popup with checkboxes for each file in source images directory
                domain = os.path.join('Files', 'Sources', 'Images')
                dir_entries = os.scandir(domain)
                warnings_manager.popup('Delete Alpha Source Chain Images')
                name_list = []
                delete_vars = []
                for entry in dir_entries:
                    if entry.is_file():
                        name, _ = os.path.splitext(entry.name)
                        name_list.append(name)
                        delete_vars.append(tk.IntVar())
                tk.Label(warnings_manager.warning, text='Files available for deletion\n').pack()
                for i, name in enumerate(name_list):
                    tk.Checkbutton(warnings_manager.warning, text=name, variable=delete_vars[i],
                                   onvalue=1, offvalue=0).pack()
                tk.Button(warnings_manager.warning, command=lambda: Delete(delete_vars, name_list, domain, '.txt'),
                          text='Delete Files').pack()
                tk.Button(warnings_manager.warning, command=lambda: warnings_manager.warning.destroy(),
                          text='Return').pack()

    # Handle Material files (.txt)
    elif Choice == 'Material':
        if Action == 1:
            # Upload: Open file dialog and copy selected file to materials directory
            filename = fd.askopenfilename(filetypes=(('Text Files', '*.txt'), ('All Files', '*.*')),
                                          title='Add Material File')
            target_dir = os.path.join(os.getcwd(), 'Files', 'Materials')
            if filename:
                copy2(filename, target_dir, follow_symlinks=True)
        else:
            # Delete: Show popup with checkboxes for each file in materials directory
            domain = os.path.join('Files', 'Materials')
            dir_entries = os.scandir(domain)
            warnings_manager.popup('Delete Material Files')
            name_list = []
            delete_vars = []
            for entry in dir_entries:
                if entry.is_file():
                    name, _ = os.path.splitext(entry.name)
                    name_list.append(name)
                    delete_vars.append(tk.IntVar())
            tk.Label(warnings_manager.warning, text='Files available for deletion\n').pack()
            for i, name in enumerate(name_list):
                tk.Checkbutton(warnings_manager.warning, text=name, variable=delete_vars[i],
                               onvalue=1, offvalue=0).pack()
            tk.Button(warnings_manager.warning, command=lambda: Delete(delete_vars, name_list, domain, '.txt'),
                      text='Delete Files').pack()
            tk.Button(warnings_manager.warning, command=lambda: warnings_manager.warning.destroy(),
                      text='Return').pack()
        # Update materials_list after changes
        Dir = os.scandir(os.path.join('Files', 'Materials'))
        materials_list.clear()
        for entry in Dir:
            if entry.is_file():
                name, _ = os.path.splitext(entry.name)
                materials_list.append(name)

    # Update dropdown menus in all tabs after file changes
    for i in range(len(TabTracker)):
        if TabTracker[i] < 0:
            TabList[i][1].Source_Menu.destroy()
            TabList[i][1].Source_Menu = tk.OptionMenu(
                TabList[i][1].SourceFrame, TabList[i][1].Source,
                *source_list, command=SourceReader)
            TabList[i][1].Source_Menu.grid(row=1, columnspan=2)
        elif TabTracker[i] > 0:
            TabList[i][1].Mat_Menu.destroy()
            TabList[i][1].Mat_Menu = tk.OptionMenu(
                TabList[i][1].SourceFrame, TabList[i][1].Mat, *materials_list)
            TabList[i][1].Mat_Menu.grid(row=1, columnspan=2)
    return

##########################################
# Allows the user to save all results    #
# obtained in the program to a text file #
##########################################
def Save_Results():
    """
    Saves all results obtained in the program to a user-specified text file.

    This function opens a file dialog for the user to choose a save location and filename.
    It then writes a summary of all calibration and material trials, including:
      - Detected peaks and counts for each calibration trial.
      - Radiation source used and linear regression results (slope, intercept, uncertainties).
      - Detected peaks, counts, and calculated thickness for each material trial.
      - Average thickness and uncertainty for each material.

    Steps:
        1. Opens a save file dialog for the user to specify the output file.
        2. Iterates through all tabs (TabTracker) to collect calibration and material results.
        3. For calibration trials:
            - Writes detected peaks and counts.
            - Writes source and regression results.
        4. For material trials:
            - Writes detected peaks, counts, and thickness.
            - Writes average thickness and uncertainty.
        5. Formats all numerical results with appropriate units and significant digits.

    Notes:
        - Uses File_Reader to load peak and regression/thickness data.
        - Handles both calibration (TabTracker < 0) and material (TabTracker > 0) tabs.
        - Units and formatting are handled according to user settings.

    Returns:
        None
    """
    # Define file types for save dialog
    filetypes = (('Text Files', '*.txt'), ('All Files', '*.*'))
    file = fd.asksaveasfile(
        title='Save Results',
        initialdir=".",
        filetypes=filetypes,
        defaultextension=".txt")

    if file:
        file.write("The Results calculated by NUC-RIA's ARC-TF were the following:\n\n")
        # Loop through all tabs for calibration trials
        for i in range(len(TabTracker)):
            if TabTracker[i] < 0:
                file.write(f'Calibration Trial {-TabTracker[i]}\n\n')
                file.write('The detected Peaks and Counts were:\n\n')
                peaks = File_Reader(TabList[i][3], ',', 'Yes', 'No')
                for peak in peaks:
                    file.write(f'Channel: {peak[0]:.1f}\tCounts: {peak[1]:.1f}\n')
                file.write('\nThe Radiation source used for this trial was: ' +
                           TabList[i][1].Source.get() + '\n\n')
                regression = File_Reader(TabList[i][4], '0', 'String', 'No')
                file.write('The Linear Regression calculated was the following, using ' +
                           regression[0] + ' units.\n\n')
                file.write('The Slope: ' +
                           '%.*f' % (int(regression[5]), float(regression[1])) + ' ± ' +
                           '%.*f' % (int(regression[5]), float(regression[2])) + '\n')
                file.write('The Intersect: ' +
                           '%.*f' % (int(regression[6]), float(regression[3])) + ' ± ' +
                           '%.*f' % (int(regression[6]), float(regression[4])) + '\n')
                file.write('_______________________________________________________________\n\n')

        # Loop through all tabs for material trials
        for i in range(len(TabTracker)):
            if TabTracker[i] > 0:
                file.write(f'Material Trial {TabTracker[i]}\n\n')
                file.write('The detected Channels, Counts and respective Thickness approximation were:\n\n')
                peaks = File_Reader(TabList[i][3], ',', 'Yes', 'No')
                peaks.sort()
                thickness = File_Reader(TabList[i][4], '0', 'String', 'No')

                units_list = [
                    'nm', '\u03bcm',
                    '\u03bcg cm\u207B\u00B2',
                    '10\u00B9\u2075 Atoms cm\u207B\u00B3'
                ]
                units_values = [1e9, 1e6, 0.0, -1.0]
                index = units_values.index(TabList[i][1].units.get())

                for j, peak in enumerate(peaks):
                    file.write(
                        f'Channel: {peak[0]:.1f}\tCounts: {peak[1]:.1f}\tThickness: '
                        f'%.*f {units_list[index]}\n' % (int(thickness[-1]), float(thickness[j]))
                    )

                file.write('\nThe average thickness, of material ' +
                           TabList[i][1].Mat.get() + ' was calculated to be: (' +
                           '%.*f' % (int(thickness[-1]), float(thickness[j + 1])) +
                           ' ± ' +
                           '%.*f' % (int(thickness[-1]), float(thickness[j + 2])) + ') ' +
                           units_list[index] + '\n\n')
                file.write('_______________________________________________________________\n\n')
    return
        
#############################################################################
# Skeleton: Main application window and menu bar class for ARC-TF GUI       #
#                                                                           #
# This class creates the main Tkinter window, sets up the menu bar, and     #
# provides access to all major program functions (file operations, tab      #
# management, settings, help, etc.).                                        #
#############################################################################
class Skeleton:
    """
    Main application window and menu bar class for the ARC-TF GUI.

    Responsibilities:
        - Initializes the main Tkinter window with title, size, and background.
        - Sets up the main menu bar with submenus for:
            * File operations (plot, save, clear, exit)
            * Settings (opens settings dialog)
            * Tab management (add/remove calibration/material tabs)
            * Data file management (add/remove source files, images, materials)
            * Help and About (opens help dialog or README link)
        - Each menu item is linked to the appropriate callback function.

    Notes:
        - Variable and method names are in English for clarity.
        - Menu and submenu variables use double underscores to indicate internal use.
        - All menu commands use lambda for parameterized callbacks.
        - The run() method starts the Tkinter main event loop.

    Attributes:
        main (tk.Tk): The main application window.
        menu (tk.Menu): The main menu bar.

    Methods:
        __init__(): Initializes the window and menu bar.
        run(): Starts the Tkinter main event loop.
    """
    def __init__(self):
        # Main application window setup
        self.main = tk.Tk()
        self.main.title('ARC-TF')
        try:
            if sys.platform == 'win32':
                self.main.state('zoomed')
            else:
                self.main.attributes('-zoomed', True)  # May work on Linux (e.g., GNOME)
        except Exception as e:
            print(f"Window zoom/maximize not supported: {e}")
        self.main.configure(background = 'dark grey')

        # Main menu bar
        self.menu = tk.Menu(self.main)
        self.main.config(menu=self.menu)

        # File menu: plot, save, clear, exit
        __file_menu = tk.Menu(self.menu, tearoff=False)
        self.menu.add_cascade(label='File', menu=__file_menu)
        __file_menu.add_command(label='Plot Data', command=DataUploader)
        __file_menu.add_command(label="Save Results", command=Save_Results)
        __file_menu.add_separator()
        __file_menu.add_command(label='Remove Current Plot', command=lambda: ClearWidget('Graphic', 0))
        __file_menu.add_command(label='Remove Algorithm Results', command=lambda: ClearWidget('Results', 0))
        __file_menu.add_command(label="Reset All Data from Current Tab", command=lambda: ClearWidget('Everything', 1))
        __file_menu.add_separator()
        __file_menu.add_command(label='Exit', command=self.main.quit)

        # Settings menu
        self.menu.add_command(label='Settings', command=lambda: warnings_manager.Settings())

        # Tab management menu
        __tabs_menu = tk.Menu(self.menu, tearoff=False)
        self.menu.add_cascade(label='Manage Tabs', menu=__tabs_menu)
        __tabs_menu.add_command(label='Add Calibration Tab', command=lambda: Tabs.tab_change(1))
        __tabs_menu.add_command(label='Add Material Tab', command=lambda: Tabs.tab_change(2))
        __tabs_menu.add_command(label='Add XRA Tab', command=lambda: Tabs.tab_change(5))        
        __tabs_menu.add_separator()
        __tabs_menu.add_command(label='Remove Current Tab', command=lambda: Tabs.tab_change(4))

        # Data file management menu
        __files_data = tk.Menu(self.menu, tearoff=False)
        self.menu.add_cascade(label='Manage Data Files', menu=__files_data)
        __files_data.add_command(label='Add Alpha Source File', command=lambda: File_Manager('Source', 1, 1))
        __files_data.add_command(label='Remove Alpha Source File', command=lambda: File_Manager('Source', 1, 0))
        __files_data.add_separator()
        __files_data.add_command(label='Add Alpha Source Image', command=lambda: File_Manager('Source', 0, 1))
        __files_data.add_command(label='Remove Alpha Source Image', command=lambda: File_Manager('Source', 0, 0))
        __files_data.add_separator()
        __files_data.add_command(label='Add Material File', command=lambda: File_Manager('Material', 0, 1))
        __files_data.add_command(label='Remove Material File', command=lambda: File_Manager('Material', 0, 0))

        # Help menu
        self.menu.add_command(label='Help', command=lambda: warnings_manager.Help())

        # About menu (opens README in browser)
        self.menu.add_command(
            label="About",
            command=lambda: webbrowser.open(
                'https://github.com/AlexVnGit/GUI_thin_films/blob/master/README.md', new=1
            )
        )

    def run(self):
        """
        Starts the Tkinter main event loop to keep the application running.
        """
        self.main.mainloop()

##############################################################################
# Warnings: Versatile class for creating popups and dialogs in the ARC-TF    #
# GUI. Handles warning messages, image displays, settings dialogs, and help  #
##############################################################################
class Warnings:
    """
    Provides popup dialogs and utility windows for the ARC-TF GUI.

    Responsibilities:
        - Display warning and information popups that must be closed before returning to the main window.
        - Show decay chain images with optional source links.
        - Present a settings dialog for adjusting general and algorithm-specific parameters.
        - Display the help window with scrollable content.

    Methods:
        popup(name): Opens a modal popup window with the given title.
        Images(name, picture, site): Shows a decay chain image and optional source link.
        Settings(): Opens the settings dialog for energy units, thickness units, and algorithm parameters.
        Help(): Opens a scrollable help window with content from the help file.

    Notes:
        - Uses self.warning, self.decay, self.configuration, and self.helping as Toplevel windows.
        - Uses global window.main as the parent for all popups.
        - Variable and method names are in English for clarity.
        - Settings dialog supports applying changes to the current tab or all tabs.

    Attributes:
        warning (tk.Toplevel): Popup window for warnings and selection menus.
        decay (tk.Toplevel): Popup window for displaying decay chain images.
        configuration (tk.Toplevel): Popup window for settings dialog.
        helping (tk.Toplevel): Popup window for help content.

    Returns:
        None
    """
    def popup(self, name):
        # Create a modal popup window with the given title
        self.warning = tk.Toplevel(main_window.main)
        self.warning.title(name)
        self.warning.wait_visibility()
        self.warning.geometry('700x300')
        self.warning.grab_set()
    
    def show(message):
        tk.messagebox.showwarning("Warning", message)    

    def Images(self, name, picture, site):
        # Display a decay chain image in a popup window, with optional source link
        self.decay = tk.Toplevel(main_window.main)
        self.decay.title(name)
        self.decay.geometry('1000x600')

        load = Image.open(picture)
        load_resize = load.resize((700, 500))
        img = ImageTk.PhotoImage(load_resize)
        label = tk.Label(self.decay, image=img)
        label.image = img
        label.pack()
        site_label = tk.Label(self.decay, text='From:  ' + site)
        site_label.pack()

    def Settings(self):
        # Open the settings dialog for energy units, thickness units, and algorithm parameters
        num = Current_Tab()

        self.configuration = tk.Toplevel(main_window.main)
        self.configuration.title('Settings')
        self.configuration.geometry('700x400')
        self.configuration.grab_set()

        self.parameter = ttk.Notebook(self.configuration)
        self.parameter.pack(expand=True, fill='both')
        self.parameter.enable_traversal()
        self.general_tab = tk.Frame(self.parameter)
        self.algorithm_tab = tk.Frame(self.parameter)
        self.parameter.add(self.general_tab, text='General Settings')
        self.parameter.add(self.algorithm_tab, text='Algorithm Settings')

        self.general_tab.columnconfigure(0, weight=1)
        self.general_tab.columnconfigure(1, weight=1)
        self.general_tab.columnconfigure(2, weight=1)

        self.algorithm_tab.columnconfigure(0, weight=1)
        self.algorithm_tab.columnconfigure(1, weight=1)
        self.algorithm_tab.columnconfigure(2, weight=1)
        self.algorithm_tab.columnconfigure(3, weight=1)

        # Store current settings for cancel/restore
        self.energy_value = TabList[num][1].energy.get()
        self.unit_value = TabList[num][1].units.get()
        self.peak_interval = TabList[num][1].peaks_widths.get()
        self.cut_low_energy = TabList[num][1].channels_cut.get()

        self.buttons_frame = tk.Frame(self.configuration)
        self.buttons_frame.pack()

        def close(choice):
            # Handle closing the settings dialog and applying/cancelling changes
            num = Current_Tab()
            if choice == 0:
                # Cancel and restore previous values
                TabList[num][1].units.set(warnings_manager.unit_value)
                TabList[num][1].energy.set(warnings_manager.energy_value)
                TabList[num][1].peaks_widths.set(warnings_manager.peak_interval)
                TabList[num][1].channels_cut.set(warnings_manager.cut_low_energy)
                warnings_manager.configuration.destroy()
            elif choice == 1:
                # Apply to current tab only
                try:
                    TabList[num][1].peaks_widths.set(warnings_manager.entry2.get())
                    TabList[num][1].channels_cut.set(warnings_manager.entry1.get())
                except ValueError:
                    TabList[num][1].peaks_widths.set(warnings_manager.peak_interval)
                    TabList[num][1].channels_cut.set(warnings_manager.cut_low_energy)
                warnings_manager.configuration.destroy()
            elif choice == 2:
                # Apply to all tabs and update global settings
                try:
                    TabList[num][1].peaks_widths.set(warnings_manager.entry2.get())
                    TabList[num][1].channels_cut.set(warnings_manager.entry1.get())
                except ValueError:
                    TabList[num][1].peaks_widths.set(warnings_manager.peak_interval)
                    TabList[num][1].channels_cut.set(warnings_manager.cut_low_energy)
                Energy_settings.set(TabList[num][1].energy.get())
                Unit_settings.set(TabList[num][1].units.get())
                Channel_cut.set(TabList[num][1].channels_cut.get())
                Peak_Width.set(TabList[num][1].peaks_widths.get())
                warnings_manager.configuration.destroy()

        # General settings UI
        tk.Label(self.general_tab, text=' ').grid(row=0, column=0, pady=10, padx=10)
        tk.Label(self.general_tab, text='Specify Energy Units').grid(row=1, column=0, pady=10, padx=10)
        tk.Label(self.general_tab, text='Specify Thickness Units').grid(row=1, column=2, pady=10, padx=10)
        tk.Radiobutton(self.general_tab, text='kev', variable=TabList[num][1].energy, value=1).grid(row=2, column=0)
        tk.Radiobutton(self.general_tab, text='Mev', variable=TabList[num][1].energy, value=1000).grid(row=3, column=0)
        tk.Radiobutton(self.general_tab, text='nm', variable=TabList[num][1].units, value=10.0**9).grid(row=2, column=2)
        tk.Radiobutton(self.general_tab, text='\u03bcm', variable=TabList[num][1].units, value=10.0**6).grid(row=3, column=2)
        tk.Radiobutton(self.general_tab, text='\u03bcg cm\u207B\u00B2', variable=TabList[num][1].units, value=0.0).grid(row=4, column=2)
        tk.Radiobutton(self.general_tab, text='10\u00b9\u2075 Atoms cm\u207B\u00b3', variable=TabList[num][1].units, value=-1.0).grid(row=5, column=2)

        # Algorithm settings UI
        tk.Label(self.algorithm_tab, text='Threshold Input Algorithm Settings').grid(row=0, columnspan=3, pady=10)
        tk.Label(self.algorithm_tab, text='Low Energy Cut: ').grid(row=1, column=1, pady=10)
        self.entry1 = tk.Entry(self.algorithm_tab)
        self.entry1.grid(row=1, column=2, pady=10)
        self.entry1.insert(0, TabList[num][1].channels_cut.get())
        tk.Label(self.algorithm_tab, text='Approximate width of peaks: ').grid(row=2, column=1, pady=10)
        self.entry2 = tk.Entry(self.algorithm_tab)
        self.entry2.grid(row=2, column=2, pady=10)
        self.entry2.insert(0, TabList[num][1].peaks_widths.get())

        # Buttons for applying/cancelling settings
        tk.Button(self.buttons_frame, text='Cancel and Return', command=lambda: close(0)).grid(row=1, column=2, padx=5, pady=5)
        tk.Button(self.buttons_frame, text='Apply to current Tab', command=lambda: close(1)).grid(row=1, column=0, padx=5, pady=5)
        tk.Button(self.buttons_frame, text='Apply to all Tabs', command=lambda: close(2)).grid(row=1, column=1, padx=5, pady=5)

    def Help(self):
        # Open the help window with scrollable content from the help file
        try:
            warnings_manager.helping.destroy()
        except:
            pass

        with open('Files\\Help.txt', 'r') as OpenFile:
            lines = OpenFile.read()

        self.helping = tk.Toplevel(main_window.main)
        self.helping.title('Help')
        self.helping.resizable(0, 0)
        self.frame = tk.Frame(self.helping)
        self.frame.pack(expand=True, fill='both')
        self.canvas = tk.Canvas(self.frame)
        self.frame2 = tk.Frame(self.canvas)
        self.scrollbar = tk.Scrollbar(self.frame)
        self.canvas.config(yscrollcommand=self.scrollbar.set, highlightthickness=0)
        self.scrollbar.config(orient=tk.VERTICAL, command=self.canvas.yview)
        self.canvas.grid(row=1, column=0)
        self.scrollbar.grid(row=1, column=2, sticky='ns')
        self.canvas.create_window(0, 0, window=self.frame2, anchor=tk.NW)
        self.frame2.bind('<Configure>', lambda e: self.canvas.configure(
            scrollregion=self.canvas.bbox('all'), width=e.width))

        tk.Label(self.frame2, text=lines).grid()

#############################################################################
# Manages the tabbed interface (Notebook) for ARC-TF, including             #
# creation, configuration, and management of calibration and material tabs. #
# Handles the layout and variable setup for each tab, as well as tab events.#
#############################################################################
class Tabs:
    """
    Manages the tabbed interface (Notebook) for the ARC-TF GUI.

    Responsibilities:
        - Creates and configures the main Notebook widget for tabbed navigation.
        - Sets up the initial "Final Results" and "+" (add tab) tabs.
        - Handles creation and layout of calibration and material analysis tabs.
        - Manages tab-specific frames (graphics, data, source, algorithm, results, etc.).
        - Initializes and tracks all variables needed for each tab (e.g., algorithm settings, ROI, regression, etc.).
        - Handles tab events (adding, removing, renaming, switching).
        - Provides static methods for tab management and renaming.

    Notes:
        - Uses both English and Portuguese variable names for legacy and clarity.
        - Some variable names (e.g., value, index) could be more descriptive.
        - Some code repetition exists in tab creation (could be refactored for DRY).
        - All tab-specific variables are initialized per tab instance for isolation.
        - Comments clarify the structure and intent of each section.

    Attributes:
        material_tab_counter (int): Counter for material tabs.
        calibration_tab_counter (int): Counter for calibration tabs.
        notebook (ttk.Notebook): The main tabbed widget.
        value (int): Tracks the number of tabs.
        [Many per-tab attributes for frames, variables, and controls.]

    Methods:
        First_Tabs(): Initializes the main Notebook and result frames.
        AnalysisTab(choice): Sets up a calibration or material analysis tab.
        tab_change(num): Static method to add/remove/switch tabs.
        RenameTab(name): Static method to rename a tab.

    Returns:
        None
    """
    material_tab_counter = 0
    calibration_tab_counter = 0
    xra_tab_counter = 0

    def First_Tabs(self):

        # Create the main Notebook widget for tab navigation
        self.notebook = ttk.Notebook(main_window.main)   
        self.notebook.pack(expand = True, fill = 'both') # Expand the notebook to fill the window
        self.notebook.enable_traversal() # Enable Ctrl+Tab navigation
        
        # Main frames: results and add-tab frame
        self.CRFrame = tk.Frame(self.notebook, bg = 'dark grey')
        self.PlusFrame = tk.Frame(self.notebook, bg = 'dark grey')
        self.CRFrame.columnconfigure(0, weight = 3)
        self.CRFrame.columnconfigure(1, weight = 3)
        self.CRFrame.rowconfigure(0, weight = 3)
        self.CRFrame.rowconfigure(1, weight = 3)

        # Add the main frames to the notebook
        self.notebook.add(self.CRFrame, text = 'Final Results')
        self.notebook.add(self.PlusFrame, text = '+')
        
        # Frames for displaying final results and scrollbars
        self.Calib_Result = tk.Frame(self.CRFrame, borderwidth = 5, relief = 'ridge')
        self.Calib_Result.grid(row = 0, column = 0, pady = 10, padx = 30, sticky = 'nw', rowspan = 2)
        self.Calib_Result.columnconfigure(0, weight = 3)
        self.Calib_Result.columnconfigure(1, weight = 1)

        self.Calib_Result_Title = tk.Label(self.Calib_Result, 
                 text = '                       Calibration Trials Results                       \n')
        self.Calib_Result_Title.grid(row = 0, column = 0, columnspan = 3)

        self.calib_canvas = tk.Canvas(self.Calib_Result)
        self.Calib_Result2 = tk.Frame(self.calib_canvas)
        self.calib_scrollbar = tk.Scrollbar(self.Calib_Result)

        self.calib_canvas.config(yscrollcommand = self.calib_scrollbar.set, highlightthickness = 0 )
        self.calib_scrollbar.config(orient = tk.VERTICAL, command = self.calib_canvas.yview)

        self.calib_canvas.grid(row = 1, column = 0)
        self.calib_scrollbar.grid(row = 1, column = 2, sticky = 'ns')
        self.calib_scrollbar.place
        self.calib_canvas.create_window(0, 0, window = self.Calib_Result2, anchor = tk.NW)
        self.Calib_Result2.bind('<Configure>', 
                              lambda e: self.calib_canvas.configure(
                                  scrollregion = self.calib_canvas.bbox('all'), width = e.width)) 

        self.Mat_Result = tk.Frame(self.CRFrame, borderwidth = 5, relief = 'ridge')
        self.Mat_Result.grid(row = 0, column = 1, pady = 10, padx = 30, sticky = 'ne', rowspan = 2)
        self.Mat_Result.columnconfigure(0, weight = 3)
        self.Mat_Result.columnconfigure(1, weight = 1)

        self.Mat_Result_Title = tk.Label(self.Mat_Result, 
                 text = '                       Material Trials Results                       \n')
        self.Mat_Result_Title.grid(row = 0, column = 0, columnspan = 3)

        self.mat_canvas = tk.Canvas(self.Mat_Result)
        self.Mat_Result2 = tk.Frame(self.mat_canvas)
        self.mat_scrollbar = tk.Scrollbar(self.Mat_Result)

        self.mat_canvas.config(yscrollcommand = self.mat_scrollbar.set, highlightthickness = 0 )
        self.mat_scrollbar.config(orient = tk.VERTICAL, command = self.mat_canvas.yview)

        self.mat_canvas.grid(row = 1, column = 0)
        self.mat_scrollbar.grid(row = 1, column = 2, sticky = 'ns')
        self.mat_canvas.create_window(0, 0, window = self.Mat_Result2, anchor = tk.NW )
        self.Mat_Result2.bind('<Configure>', 
                              lambda e: self.mat_canvas.configure(
                                  scrollregion = self.mat_canvas.bbox('all'), width = e.width)) 
        
        # Variable to count the number of analysis tabs
        self.value = 0

        # Bind the tab change event to handle adding new tabs
        self.notebook.bind("<<NotebookTabChanged>>", handleTabChange)      

    def AnalysisTab(self, tab_kind):
        self.tab_kind = tab_kind

        ### Configuracao de geometria e Frames comuns a Calib e Material Tabs  ####
        TabList[tab_manager.value][0].columnconfigure(0, weight = 4)
        TabList[tab_manager.value][0].columnconfigure(1, weight = 1)
        TabList[tab_manager.value][0].rowconfigure(0, weight = 1)
        TabList[tab_manager.value][0].rowconfigure(1, weight = 1)

        self.GraphicFrame = tk.Frame(TabList[tab_manager.value][0], borderwidth = 5, relief = 'ridge')
        self.GraphicFrame.grid(column = 0, row = 0, sticky = "nw", pady = 5, columnspan = 2)

        self.DataFrame = tk.Frame(TabList[tab_manager.value][0], borderwidth = 5, relief = 'ridge')
        self.DataFrame.grid(column = 2, row = 0, sticky = "ne", pady = 5)

        self.SourceFrame = tk.Frame(TabList[tab_manager.value][0], borderwidth = 5, relief = 'ridge')
        self.SourceFrame.grid(column = 1, row = 0, sticky = "ne", pady = 5)

        self.Extra_Frame = tk.Frame(TabList[tab_manager.value][0], borderwidth = 5, relief = 'ridge')

        self.AlgFrame = tk.Frame(self.DataFrame, borderwidth = 0)
        self.AlgFrame.grid(row = 2, columnspan = 2, pady = 5)
       
        self.ResultFrame = tk.Frame(self.DataFrame)
        
        ##################################### Variaveis para cada instance ##########################

        self.Algorithm_Method = tk.StringVar()
        self.Algorithm_Method.set('Select Algorithm to Run')

        self.Algorithm = tk.IntVar()
        self.Algorithm.trace('w', on_entry_change)
        self.Algorithm.set(0)

        self.ROIdown1 = tk.IntVar()
        self.ROIdown1.set(1927)
        self.ROIup1 = tk.IntVar()
        self.ROIup1.set(1973)
        self.ROIdown2 = tk.IntVar()
        self.ROIdown2.set(1466)
        self.ROIup2 = tk.IntVar()
        self.ROIup2.set(1522)
        self.ROIdown3 = tk.IntVar()
        self.ROIdown3.set(1364)
        self.ROIup3 = tk.IntVar()
        self.ROIup3.set(1412)
        self.ROIdown4 = tk.IntVar()
        self.ROIdown4.set(1310)
        self.ROIup4 = tk.IntVar()
        self.ROIup4.set(1360)
        self.ROIdown5 = tk.IntVar()
        self.ROIdown5.set(1228)
        self.ROIup5 = tk.IntVar()
        self.ROIup5.set(1284)        
        self.ROIdown6 = tk.IntVar()
        self.ROIdown6.set(1138)
        self.ROIup6 = tk.IntVar()
        self.ROIup6.set(1183)

        self.units = tk.DoubleVar()
        unit = Unit_settings.get()
        self.units.set(unit)

        self.energy = tk.IntVar()
        energies = Energy_settings.get()
        self.energy.set(energies)

        self.channels_cut = tk.IntVar()
        self.channels_cut.set(Channel_cut.get())

        self.peaks_widths = tk.IntVar()
        self.peaks_widths.set(Peak_Width.get())

        self.Real_Time = tk.StringVar()
        self.Real_Time_2 = tk.StringVar()
        self.Real_Time_2.set('5')
    

        self.Total_Counts = tk.StringVar()

        self.variable1 = tk.IntVar()
        self.variable1.set(0)
        self.variable2 = tk.IntVar()
        self.variable2.set(0)
        self.variable3 = tk.IntVar()
        self.variable3.set(0)
        self.variable4 = tk.IntVar()
        self.variable4.set(0)
        self.variable5 = tk.IntVar()
        self.variable5.set(0)
        self.variable6 = tk.IntVar()
        self.variable6.set(0)
        self.variable7 = tk.IntVar()
        self.variable7.set(0)
        self.variable8 = tk.IntVar()
        self.variable8.set(0)
        self.variable9 = tk.IntVar()
        self.variable9.set(0)
        self.variable10 = tk.IntVar()
        self.variable10.set(0)
        self.variable11 = tk.IntVar()
        self.variable11.set(0)
        self.variable12 = tk.IntVar()
        self.variable12.set(0)
        self.variable13 = tk.IntVar()
        self.variable13.set(0)
        self.variable14 = tk.IntVar()
        self.variable14.set(0)
        self.variable15 = tk.IntVar()
        self.variable15.set(0)

        self.Var_Data = [
            self.variable1, self.variable2, self.variable3, self.variable4, self.variable5,
            self.variable6, self.variable7, self.variable8, self.variable9, self.variable10,
            self.variable11, self.variable12, self.variable13, self.variable14, self.variable15]
    
        ## Check if we are in an XRA analysis tab
        if tab_kind == 5:
            self.Algorithm_Method.set("ROI Select")
            Method()

            """ for i in range(1, 6):  # index 1 to 5 = ROI 2 to 6
                self.ROIdown_entries[i].grid_remove()
                self.ROIup_entries[i].grid_remove()
                self.ROIlabels[i].grid_remove() """


        else:
            tk.Label(self.DataFrame, text='Analysis Method Selected:').grid(row=0, columnspan=2)
            Algs = ["Manual Selection", "Threshold Input", "ROI Select"]
            tk.OptionMenu(self.DataFrame, self.Algorithm_Method, *Algs, command=Method).grid(row=1, columnspan=2)

        def MatTab(self):

            self.Mat = tk.StringVar()
            self.Mat.set('Select Material')

            self.Regression1 = tk.DoubleVar()
            self.Regression1.set(-1)
            self.Regression2 = tk.DoubleVar()
            self.Regression2.set(-1)
            self.Regression3 = tk.DoubleVar()
            self.Regression3.set(-1)
            self.Regression4 = tk.DoubleVar()
            self.Regression4.set(-1)
            self.Regression5 = tk.DoubleVar()
            self.Regression5.set(-1)
            self.Regression6 = tk.DoubleVar()
            self.Regression6.set(-1)
            self.Regression7 = tk.DoubleVar()
            self.Regression7.set(-1)
            self.Regression8 = tk.DoubleVar()
            self.Regression8.set(-1)
            self.Regression9 = tk.DoubleVar()
            self.Regression9.set(-1)
            self.Regression10 = tk.DoubleVar()
            self.Regression10.set(-1)

            self.Regression_List = [self.Regression1, self.Regression2, self.Regression3, self.Regression4, 
                                    self.Regression5, self.Regression6, self.Regression7, self.Regression8,
                                    self.Regression9, self.Regression10]
            
            self.Calib_Sel_Frame = tk.Frame(self.SourceFrame, borderwidth = 0)
            self.Calib_Sel_Frame.grid(row = 4, columnspan = 2)

            self.ThicknessFrame = tk.Frame(self.SourceFrame)

            tk.Label(self.SourceFrame, text = 'Material of Film Used: ').grid(row = 0, columnspan = 2)
            self.Mat_Menu = tk.OptionMenu(self.SourceFrame, self.Mat, *materials_list)
            self.Mat_Menu.grid(row = 1, columnspan = 2)

            num = Current_Tab()
            value = TabTracker[num]
            
            tk.Label(self.SourceFrame, text = 'Choose the Calibration Regression \n'+ 
                     'Make sure the Peaks in the Calibration ' + 
                     'match the Peaks found for this trial').grid(row = 2, columnspan = 2, pady = 5)
            tk.Button(self.SourceFrame, text = 'Select Calibration', 
                      command = Calib_Choice).grid(row = 3, column = 0)
            tk.Button(self.SourceFrame, text = 'Calculate Thickness', 
                      command = lambda: Final_Results(value)).grid(row = 3, column = 1)

        def CalibTab(self):

            self.Source = tk.StringVar()
            self.Source.set('Radiation Sources')

            self.decay1 = tk.DoubleVar()
            self.decay1.set(-1)
            self.decay2 = tk.DoubleVar()
            self.decay2.set(-1)
            self.decay3 = tk.DoubleVar()
            self.decay3.set(-1)           
            self.decay4 = tk.DoubleVar()
            self.decay4.set(-1) 
            self.decay5 = tk.DoubleVar()
            self.decay5.set(-1) 
            self.decay6 = tk.DoubleVar()
            self.decay6.set(-1) 
            self.decay7 = tk.DoubleVar()
            self.decay7.set(-1) 

            self.DecayList = [self.decay1, self.decay2, self.decay3, self.decay4, self.decay5,
                              self.decay6, self.decay7]

            tk.Label(self.SourceFrame, 
                     text = 'Radiation Source Selected: ').grid(row = 0, columnspan = 2)
            self.Source_Menu = tk.OptionMenu(self.SourceFrame, self.Source, 
                          *source_list, command = SourceReader)
            self.Source_Menu.grid(row = 1, columnspan = 2)
            self.SourceOptionsFrame = tk.Frame(self.SourceFrame, borderwidth = 0)
            self.SourceOptionsFrame.grid(row = 2, columnspan = 2)
            self.LinearRegressionFrame = tk.Frame(self.SourceFrame, borderwidth = 1)

        #def XraTab(self):

            
        if tab_kind == 1:
            CalibTab(self)
            
        elif tab_kind == 2:
            MatTab(self)

    @staticmethod
    def tab_change(tab_type):
        """
        Handles adding, removing, and switching tabs in the Notebook.
        tab_type: 1 = add calibration, 2 = add material, 3 = cancel, 4 = remove current
        """
        value = Current_Tab()
        index = len(tab_manager.notebook.tabs()) - 1             

        if tab_type == 1:
            Tabs.calibration_tab_counter -= 1
            Data = os.path.join('Temp', 'Data' + str(Tabs.calibration_tab_counter) + '.txt')
            Analysis = os.path.join('Temp', 'Analysis' + str(Tabs.calibration_tab_counter) + '.txt')
            Result =  os.path.join('Temp', 'Result' + str(Tabs.calibration_tab_counter) + '.txt')
            ROIs = os.path.join('Temp', 'ROIs' + str(Tabs.calibration_tab_counter) + '.txt')
            TabList.append([tk.Frame(tab_manager.notebook, bg = 'dark grey'), Tabs(),  Data,
                        Analysis,  Result, Plot(), ROIs]) 

            TabTracker.append(Tabs.calibration_tab_counter)
            TabList[tab_manager.value][1].AnalysisTab(1)
            tab_manager.notebook.insert(index, TabList[tab_manager.value][0], 
                                     text = "Calibration Trial " + str(-Tabs.calibration_tab_counter))
            tab_manager.notebook.select(index)
            tab_manager.value += 1
            try:
                warnings_manager.warning.destroy()
            except:
                ()

        elif tab_type == 2:
            Tabs.material_tab_counter += 1
            Data = os.path.join('Temp', 'Data' + str(Tabs.material_tab_counter) + '.txt')
            Analysis = os.path.join('Temp', 'Analysis' + str(Tabs.material_tab_counter) + '.txt')
            Result =  os.path.join('Temp', 'Result' + str(Tabs.material_tab_counter) + '.txt')
            ROIs = os.path.join('Temp', 'ROIs' + str(Tabs.material_tab_counter) + '.txt')
            TabList.append([tk.Frame(tab_manager.notebook, bg = 'dark grey'), Tabs(),  Data,
                        Analysis,  Result, Plot(), ROIs]) 

            TabTracker.append(Tabs.material_tab_counter)
            TabList[tab_manager.value][1].AnalysisTab(2)
            tab_manager.notebook.insert(index, TabList[tab_manager.value][0], 
                                     text = "Material Trial " + str(Tabs.material_tab_counter))
            tab_manager.notebook.select(index)
            tab_manager.value += 1
            try:
                warnings_manager.warning.destroy()
            except:
                ()
        elif tab_type == 5:
            Tabs.xra_tab_counter += 1
            Data = os.path.join('Temp', 'Data' + str(Tabs.xra_tab_counter) + '.txt')
            Analysis = os.path.join('Temp', 'Analysis' + str(Tabs.xra_tab_counter) + '.txt')
            Result =  os.path.join('Temp', 'Results' + str(Tabs.xra_tab_counter) + '.txt')
            ROIs = os.path.join('Temp', 'ROIs' + str(Tabs.xra_tab_counter) + '.txt')
            TabList.append([tk.Frame(tab_manager.notebook, bg = 'dark grey'), Tabs(),  Data,
                        Analysis,  Result, Plot(), ROIs]) 

            TabTracker.append(Tabs.xra_tab_counter)
            TabList[tab_manager.value][1].AnalysisTab(5)
            tab_manager.notebook.insert(index, TabList[tab_manager.value][0],text = " XRA " + str(Tabs.xra_tab_counter))
            tab_manager.notebook.select(index)
            tab_manager.value += 1
            try:
                warnings_manager.warning.destroy()
            except:
                ()

        
        
        elif tab_type == 3:
            tab_manager.notebook.select(index - 1)
            warnings_manager.warning.destroy()

        elif tab_type == 4:
            if tab_manager.notebook.select() == '.!notebook.!frame' or tab_manager.notebook.select() == '.!notebook.!frame2':
                warnings_manager.popup('Bad Tab Deletion')
                tk.Label(warnings_manager.warning, text = '\n This Tab cannot be deleted.\nPlease delete Analysis Tabs').pack()
                tk.Label(warnings_manager.warning, text = '\n\n').pack()
                tk.Button(warnings_manager.warning, text = 'Return', command = lambda: warnings_manager.warning.destroy()).pack()

            else:
                tab_manager.notebook.forget("current")
                tab_manager.notebook.select(index - 2)

                if os.path.isfile(TabList[value][2]) == True:
                    os.remove(TabList[value][2])
                if os.path.isfile(TabList[value][3]) == True:
                    os.remove(TabList[value][3])
                if os.path.isfile(TabList[value][4]) == True:
                    os.remove(TabList[value][4])

                TabList.pop(value)

                tab_manager.value -= 1
                TabTracker.pop(value)
                ClearWidget('Final', 0)
                Final_Results(0)

    def RenameTab(name):

        tab_num = Current_Tab()
        tab_manager.notebook.tab(tab_num+1, text = str(name))

        return
   
###########################################################################
# Handles data loading, plotting, and graphical updates for ARC-TF.       #
# This class reads external data files, inserts plots into the main GUI   #
# frame, and manages plot overlays (e.g., threshold lines).               #
###########################################################################
class Plot:
    """
    Handles data loading, plotting, and graphical updates for the ARC-TF GUI.

    Responsibilities:
        - Reads data from external files and prepares it for plotting.
        - Inserts the plot into the main graphic frame of the GUI.
        - Creates a temporary text file for other functions to access the data.
        - Displays a threshold line if the threshold input algorithm is selected.
        - Provides methods to update, clear, or overlay lines on the plot.

    Methods:
        Structure(File, Name): Loads data from file, processes channels and counts, updates GUI.
        subplots(): Plots the data (channels vs. counts) and sets up axes and event handlers.
        destroyer(): Removes the last overlay line from the plot (e.g., threshold line).
        threshold(height): Draws a horizontal threshold line at the specified height.

    Notes:
        - Uses matplotlib for plotting and FigureCanvasTkAgg for embedding in Tkinter.
        - Variable names are now in English for clarity.
        - Comments explain each step of the plotting and data handling process.
        - No unnecessary code detected; logic is clear and concise.

    Attributes:
        Channel (list): List of channel indices (x-axis).
        Counts (list): List of count values (y-axis).
        line (list): List of overlay lines (e.g., threshold lines).
        Title (str): Title for the plot, set based on tab type.
        figure (Figure): Matplotlib Figure object for the plot.
        figure_canvas (FigureCanvasTkAgg): Canvas for embedding the plot in Tkinter.
        axes (Axes): Matplotlib Axes object for plotting.
    """
    def Structure(self, File, Name, File2 = None):
        # Initialize lists for channel and counts
        self.Channel = []
        self.Counts = []
        self.Channel2 = []
        self.Counts2 = []
        self.line = []

        num = Current_Tab()
        total_sum = 0
        total_sum2 = 0
        j = 0

        # Clear previous plot and set up the graphic frame
        ClearWidget('Graphic', 0)
        TabList[num][1].GraphicFrame.grid(column=0, row=0, sticky="nw", pady=5, columnspan=2)

        # Open the data file for writing processed counts
        Data = open(TabList[num][2], "w")

        # If the file is an .mca file, process accordingly (specific to AEL machine format)
        if Name[-4:] == ".mca":
            for i in range(12, len(File) - 1):
                self.Counts.append(int(File[i]))
                total_sum += self.Counts[j]
                self.Channel.append(i - 11)
                Data.write(str(self.Counts[i - 12]) + "\n")
                j += 1
        

        Data.close()

        # Update total counts and display in the extra frame
        TabList[num][1].Total_Counts.set(total_sum)
        TabList[num][1].Extra_Frame.grid(column=0, row=1, sticky="nw")
        
        tk.Label(TabList[num][1].Extra_Frame, text="Source:").grid(row=0, column=0, sticky="w")
        tk.Label(TabList[num][1].Extra_Frame, text=TabList[num][1].Real_Time.get()).grid(row=0, column=1, sticky="w")
        tk.Label(TabList[num][1].Extra_Frame, text='Total: ' + str(TabList[num][1].Total_Counts.get())).grid(row=0, column=2, sticky="w")
        
        #For XRA tab, also show File 2 time
        if TabList[num][1].tab_kind == 5:
            tk.Label(TabList[num][1].Extra_Frame, text="Source+Film:").grid(row=1, column=0, sticky="w")
            tk.Label(TabList[num][1].Extra_Frame, text=TabList[num][1].Real_Time_2.get()).grid(row=1, column=1, sticky="w")

        # Set plot title based on tab type
        if TabTracker[num] < 0:
            self.Title = 'Calibration Trial ' + str(-TabTracker[num])
        elif TabTracker[num] > 0:
            self.Title = 'Material Trial ' + str(TabTracker[num])
        else:
            self.Title = 'XRA Analysis'

        # If this is an XRA tab and File2 is provided, parse second dataset
        if getattr(TabList[num][1], 'tab_kind', None) == 5 and File2:
            for i in range(12, len(File2) - 1):
                self.Counts2.append(int(File2[i]))
                self.Channel2.append(i - 11)

        # Create the matplotlib figure and embed it in the Tkinter frame
        self.figure = Figure(figsize=(6, 4), dpi=100)
        self.figure_canvas = FigureCanvasTkAgg(self.figure, TabList[num][1].GraphicFrame)
        NavigationToolbar2Tk(self.figure_canvas, TabList[num][1].GraphicFrame)

    def subplots(self):
        # Create the plot with channels on x-axis and counts on y-axis
        self.axes = self.figure.add_subplot()
        self.axes.plot(self.Channel, self.Counts, '*', label='Source')

        # If second dataset exists, plot it
        if self.Counts2:
            self.axes.plot(self.Channel2, self.Counts2, '+', label='Source + Film', color='red')

        self.axes.set_title(self.Title)
        self.axes.set_xlabel('Channel')
        self.axes.set_ylabel('Counts')
        self.axes.legend()

        # Connect mouse click event for manual selection
        self.figure.canvas.mpl_connect('button_press_event', onclick)

        # Pack the plot widget into the Tkinter frame
        self.figure_canvas.get_tk_widget().pack()

    def destroyer(self):
        # Remove the last overlay line (e.g., threshold line) from the plot
        if self.line:
            self.line.pop().remove()
            self.figure_canvas.draw()

    def threshold(self, height):
        # Draw a horizontal threshold line at the specified height
        if self.line:
            self.line.pop().remove()
        self.line.append(self.axes.axhline(y=height, color='r', linestyle='-'))
        self.figure_canvas.draw()

#############################################################################
# Initialize source and material lists by scanning the relevant directories #
#############################################################################

# Scan the directory for alpha source value files and populate source_list
source_values_dir = os.scandir(os.path.join('Files', 'Sources', 'Values'))
Dir = os.scandir(os.path.join('Files', 'Sources', 'Values'))
source_list = []
for entry in source_values_dir:
    if entry.is_file():
        name, _ = os.path.splitext(entry.name)
        source_list.append(name)

# Scan the directory for material files and populate materials_list
materials_dir = os.scandir(os.path.join('Files', 'Materials'))
Dir = os.scandir(os.path.join('Files', 'Materials'))
materials_list = []
for entry in materials_dir:
    if entry.is_file():
        name, _ = os.path.splitext(entry.name)
        materials_list.append(name)

#############################################################################
# Initialize main application components and global variables               #
#############################################################################

# Instantiate main utility classes
warnings_manager = Warnings()
main_window = Skeleton()
tab_manager = Tabs()
tab_manager.First_Tabs()

# Lists to keep track of tab objects and their types (calibration/material)
TabList = []
TabTracker = []

# Global settings for energy, units, channel cut, and peak width
Energy_settings = tk.IntVar(value=1000)
Unit_settings = tk.DoubleVar(value=1e9)
Channel_cut = tk.IntVar(value=100)
Peak_Width = tk.IntVar(value=35)

#############################################################################
# Start the application: create the first calibration tab and run the GUI   #
#############################################################################

Tabs.tab_change(1)  # Add initial calibration tab
tab_manager.notebook.select(1)  # Select the first analysis tab

# Create a temporary directory for storing intermediate files
if not os.path.exists('Temp'):
    os.mkdir('Temp')

main_window.run()  # Start the Tkinter main event loop

#############################################################################
# Cleanup: Remove temporary files and directory after the application exits #
#############################################################################

# Remove all temporary files created for each tab
for i in range(tab_manager.value):
    for file_index in [2, 3, 4]:  # Data, Algorithm Results, Regression Results
        temp_file = TabList[i][file_index]
        if os.path.isfile(temp_file):
            os.remove(temp_file)

# Remove the temporary directory if it exists
if os.path.isdir('Temp'):
    os.rmdir('Temp')
