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

from Include.Analyze import*
from Include.Calibration import*
from Include.FitData import*
from Include.Eloss import*
from Include.Thick import*
from Include.clear_frame import*
from Include.remove_file import*
from Include.update_file_list import*
from Include.show_deletion_popup import*
## ---------------------------------------------------------------------------------------------- ##

###############################################
# Handles display scaling                     #
# works on Windows 10/11 and Linux/Kubuntu    #
###############################################
ctypes.windll.shcore.SetProcessDpiAwareness(1) 

###############################################
# Returns the index of the Tab the user is on #
###############################################
def Current_Tab():
    """
    Returns the ID of the tab where the user is on,
    by getting the index of the selected tab in the notebook.
    """
    ## Get index of the current tab in the notebook
    tabID = Notebook.notebook.index(Notebook.notebook.select()) - 1

    ## Check if the user is in the "final results" tab (has negative index)
    if tabID < 0:
        pass
    else:
        return tabID
  
#########################################################################################
# Clears and resets specific UI frames and associated data files depending on the 
# type of frame specified.
#########################################################################################
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
       - Depends on external objects like `Notebook`, `wng`.
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
        for widget in wng.warning.winfo_children():
            widget.destroy()
        wng.warning.destroy()

    elif Frame == 'Image':
        # Destroy decay image popup
        for widget in wng.decay.winfo_children():
            widget.destroy()
        wng.decay.destroy()

    elif Frame == 'Linear':
        # Clear linear regression result display and optionally delete file
        clear_frame(tab.LinearRegressionFrame)
        clear_frame(Notebook.Calib_Result2)
        if parameter == 1:
            remove_file(TabList[num][4])  # Calibration file

    elif Frame == 'Thickness':
        # Clear thickness calculation results and optionally delete file
        clear_frame(tab.ThicknessFrame)
        clear_frame(Notebook.Mat_Result2)
        if parameter == 1:
            tab.Mat.set('Select Material')
            remove_file(TabList[num][4])  # Thickness file

    elif Frame == 'Final':
        # Clear final results display
        clear_frame(Notebook.Mat_Result2)
        clear_frame(Notebook.Calib_Result2)

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

    # Open file and read its entire content
    with open(Document, 'r') as OpenFile:
        lines = OpenFile.read()
    
    # Split the content by lines
    lines = lines.splitlines()

    # If needed, update GUI with the "measurement live-time" value from line 9
    if Upload == 'Yes':
        TabList[num][1].Real_Time.set(lines[8] + ' s')

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
                tk.Label(Notebook.Calib_Result2, text=f'Calibration Trial {-TabTracker[i]} - {TabList[i][1].Source.get()}').grid(row=4*i+i, columnspan=3)
                tk.Label(Notebook.Calib_Result2, text=f'({unit_energy})').grid(row=4*i+i+1, column=0)
                tk.Label(Notebook.Calib_Result2, text='Values').grid(row=4*i+i+1, column=1)
                tk.Label(Notebook.Calib_Result2, text='Uncertainty').grid(row=4*i+i+1, column=2)
                tk.Label(Notebook.Calib_Result2, text='Slope').grid(row=4*i+i+2, column=0)
                tk.Label(Notebook.Calib_Result2, text='Intersect').grid(row=4*i+i+3, column=0)
                tk.Label(Notebook.Calib_Result2, text='').grid(row=4*i+i+4, column=0)

                # Display slope and intercept with uncertainty
                tk.Label(Notebook.Calib_Result2, text='%.*f' % (int(Results[5]), float(Results[1]))).grid(row=4*i+i+2, column=1)
                tk.Label(Notebook.Calib_Result2, text='%.*f' % (int(Results[5]), float(Results[2]))).grid(row=4*i+i+2, column=2)
                tk.Label(Notebook.Calib_Result2, text='%.*f' % (int(Results[6]), float(Results[3]))).grid(row=4*i+i+3, column=1)
                tk.Label(Notebook.Calib_Result2, text='%.*f' % (int(Results[6]), float(Results[4]))).grid(row=4*i+i+3, column=2)

            else:
                # Update canvas if no result file
                Notebook.calib_canvas.update_idletasks()
                Notebook.calib_canvas.config(scrollregion=Notebook.Calib_Result2.bbox())
                Notebook.Calib_Result2.bind('<Configure>',
                    lambda e: Notebook.calib_canvas.configure(
                        scrollregion=Notebook.calib_canvas.bbox('all'), width=e.width))

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
                        tk.Label(Notebook.Mat_Result2, text=f'Material Trial {TabTracker[i]} - {TabList[i][1].Mat.get()}').grid(row=(4+size)*i+i+j, columnspan=2)
                        tk.Label(Notebook.Mat_Result2, text='Peak Centroid').grid(row=(4+size)*i+i+j+1, column=0)
                        tk.Label(Notebook.Mat_Result2, text='Thickness').grid(row=(4+size)*i+i+j+1, column=1)

                    # Peak energy and thickness
                    tk.Label(Notebook.Mat_Result2, text=f'{Peaks[j][0]:.1f}').grid(row=(4+size)*i+i+j+2, column=0)
                    tk.Label(Notebook.Mat_Result2, text='%.*f %s' % (int(Results[-1]), Results[j], units_list[index])).grid(row=(4+size)*i+i+j+2, column=1)

                # Display average and uncertainty
                tk.Label(Notebook.Mat_Result2, text='\nAverage').grid(row=(4+size)*i+i+j+3, column=0)
                tk.Label(Notebook.Mat_Result2, text='Uncertainty').grid(row=(4+size)*i+i+j+4, column=0)
                tk.Label(Notebook.Mat_Result2, text='\n%.*f %s' % (int(Results[-1]), Results[j+1], units_list[index])).grid(row=(4+size)*i+i+j+3, column=1)
                tk.Label(Notebook.Mat_Result2, text='%.*f %s' % (int(Results[-1]), Results[j+2], units_list[index])).grid(row=(4+size)*i+i+j+4, column=1)
                tk.Label(Notebook.Mat_Result2, text='').grid(row=(4+size)*i+i+j+5, columnspan=2)

                # Update scrollable canvas
                Notebook.mat_canvas.update_idletasks()
                Notebook.mat_canvas.config(scrollregion=Notebook.Mat_Result2.bbox())
                Notebook.Mat_Result2.bind('<Configure>',
                    lambda e: Notebook.mat_canvas.configure(
                        scrollregion=Notebook.mat_canvas.bbox('all'), width=e.width))
    return

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
        - Relies on tkinter for GUI elements and wng for popup management.

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
        wng.popup('No Linear Regressions detected')
        tk.Label(wng.warning, text=(
            'No linear regressions were detected.\n\n'
            'Please perform a Calibration Trial before calculating the film\'s thickness.\n\n')).pack()
        tk.Button(wng.warning, text='Return', command=lambda: wng.warning.destroy()).pack()

    # Otherwise, show selection popup for user to choose regressions
    else:
        wng.popup('Linear Regression Selection Menu')
        tk.Label(wng.warning, text=(
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
                            wng.warning,
                            text=f'Linear Regression of Calibration Trial {-validCalib[i]}',
                            variable=TabList[num][1].Regression_List[i],  # Correct indexing
                            onvalue=on_value,
                            offvalue=-1)
            TabList[num][1].Regression_List[i].set(-1)  # Ensure unchecked by default
            button_Choice.pack()
        
        # Store the mapping for later use (e.g., as an attribute of the tab or globally)
        TabList[num][1].regression_value_map = regression_value_map

        # Button to simply close the popup (could be extended with a Confirm button)
        tk.Button(wng.warning, text='Return', command=lambda: ClearWidget('Popup', 0)).pack()

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
            wng.popup('Invalid Linear Regression Configuration')
            tk.Label(wng.warning, 
                     text="Number of Radiation Decay does not match the number of Peaks detected.\n").pack()
            tk.Label(wng.warning, 
                     text="Please adjust the Searching Algorithms or the number of Decay Energy.\n\n").pack()
            tk.Button(wng.warning, text='Return', command=lambda: wng.warning.destroy()).pack()
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
            wng.popup('Invalid Linear Regression Configuration')
            tk.Label(wng.warning, 
                     text="Number of Radiation Decay does not match the number of Peaks detected.\n").pack()
            tk.Label(wng.warning, 
                     text="Please adjust the Searching Algorithms or the number of Decay Energy.\n\n").pack()
            tk.Button(wng.warning, text='Return', command=lambda: wng.warning.destroy()).pack()
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
# ROI_Select_Alg: Detects peaks within user-defined ROIs and records their      #
# centroids, uncertainties, and sigma/sqrt(N) for the ROI Select algorithm.     #
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
        - Relies on the wng.Images helper for displaying images and links.

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
                    wng.Images('Decay Chain of ' + alphas, picture, domain)  # Display image and link
                except:
                    # If clearing fails, still try to display the image
                    wng.Images('Decay Chain of ' + alphas, picture, domain)

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
    if Notebook.notebook.select() == Notebook.notebook.tabs()[-1] and len(Notebook.notebook.tabs()) < 15:
        # Show popup menu for tab selection
        wng.popup('Tab Selector Menu')

        tk.Label(wng.warning, text='Please Select New Type of Tab to Open \n').pack()
        tk.Button(wng.warning, command=lambda: Tabs.tab_change(1),
                  text='Calibration Trial').pack()
        tk.Button(wng.warning, command=lambda: Tabs.tab_change(2), 
                  text='Material Trial').pack()
        tk.Label(wng.warning, text='\n').pack()
        tk.Button(wng.warning, text='Return',
                  command=lambda: Tabs.tab_change(3)).pack()
        # The numbers passed to tab_change() indicate which type of tab to add

    # If there are already 15 or more tabs, hide the '+' tab to prevent adding more
    elif len(Notebook.notebook.tabs()) >= 15:
        Notebook.notebook.hide(14)

    return

############################################################################
# Esta funcao gere as fontes de radiacao
#############################################################################
def SourceReader(*args):

    num = Current_Tab()
    value = TabTracker[num]
    
    ClearWidget('Source', 0) # Reset dos widgets e da geometria
    TabList[num][1].SourceOptionsFrame.grid(row = 2, columnspan = 2)

    for i in range(len(TabList[num][1].DecayList)):
        TabList[num][1].DecayList[i].set(-1) # Aqui garantimos que a mudanca de fontes de radiacao
        # alpha, nao interfere com a regressao linear a ser efetuada
    
    Alpha = TabList[num][1].Source.get()
    Alpha = 'Files\Sources\Values\\' +  Alpha + '.txt' # Esta adicao garante que o programa 
                                                    #encontra o ficheiro pretendido
    with open(Alpha, 'r') as file: # Aqui vai-se buscar todas as opcoes possiveis contidas no ficheiro
                                    # das fontes de radiacao
        Decay = [float(line) for line in file]

    for i in range(len(Decay)): # Aqui esta a criacao dos checkbuttons para selecionar os valores que o 
                                # utilizador pretende empregrar nos calculos
        checkbutton = tk.Checkbutton(TabList[num][1].SourceOptionsFrame, text = str(Decay[i]) + ' MeV',
                       variable = TabList[num][1].DecayList[i], onvalue = Decay[i], offvalue = -1)
        checkbutton.grid(row = i, columnspan = 2)
        checkbutton.select()

    tk.Button(TabList[num][1].SourceOptionsFrame, # Mostra uma imagem da cadeia 
              text = 'Show Decay Chain', command = showimage).grid(row = i + 1, column = 0)
    tk.Button(TabList[num][1].SourceOptionsFrame, text = 'Linear Regression',  # Efetua a regressao linear e poe os resultados na primeira tab
              command = lambda: Final_Results(value)).grid(row = i + 1, column = 1)
    
##############################################################################
# Esta funcao altera a interface dos dados inputs para cada algoritmo
##############################################################################
def Method(*args):

    num = Current_Tab()

    ClearWidget('Algorithm', 0)
    TabList[num][1].AlgFrame.grid(row = 2, columnspan = 2)

    decider = TabList[num][1].Algorithm_Method.get() #Devolve a StringVar que decide qual o algoritmo em uso

    if decider == 'Manual Selection':   #Definicao dos controlos para o algoritmo ManSelec_Alg
        tk.Label(TabList[num][1].AlgFrame, 
                 text = 'Right-Click on/near the Peaks in the Graphic ').grid(row = 2, columnspan = 2)
        tk.Label(TabList[num][1].AlgFrame, 
                 text = 'For an automatic point detection: ').grid(row = 3, columnspan = 2)
        tk.Button(TabList[num][1].AlgFrame, text = 'Remove Unchecked', 
                  command = Unchecked_Results).grid(row = 4, column = 0)
        tk.Button(TabList[num][1].AlgFrame, text = 'Remove All',
                  command = lambda: ClearWidget('Results', 1)).grid(row = 4, column = 1)

        TabList[num][1].Algorithm.set(0)
        TabList[num][5].destroyer()
     
    elif decider == 'Threshold Input':  #Definicao dos controlos para o algoritmo Threshold_Alg
        tk.Label(TabList[num][1].AlgFrame, 
                 text = 'Please input Threshold: ').grid(row = 2, columnspan = 3)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].Algorithm, relief = 'sunken',
                  borderwidth = 2).grid(row = 3, columnspan= 3)
        tk.Button(TabList[num][1].AlgFrame,text = 'Search', 
                  command =  Threshold_Alg).grid(row = 4, column = 0)
        tk.Button(TabList[num][1].AlgFrame,text = 'Remove Unchecked',
                  command = Unchecked_Results).grid(row = 4, column = 1)
        tk.Button(TabList[num][1].AlgFrame,text = 'Remove All',
                  command = lambda: ClearWidget('Results', 1)).grid(row = 4, column = 2)
    
    elif decider == 'ROI Select':  #Definicao dos controlos para o algoritmo ROI Select
        tk.Label(TabList[num][1].AlgFrame, 
                 text = 'ROI Down: ').grid(row = 2, column = 0)
        tk.Label(TabList[num][1].AlgFrame, 
                 text = 'ROI Up: ').grid(row = 2, column = 1)

        tk.Label(TabList[num][1].AlgFrame, 
            text = 'Peak 1').grid(row = 3, column = 2)
        tk.Label(TabList[num][1].AlgFrame, 
            text = 'Peak 2').grid(row = 4, column = 2) 
        tk.Label(TabList[num][1].AlgFrame, 
            text = 'Peak 3').grid(row = 5, column = 2) 
        tk.Label(TabList[num][1].AlgFrame, 
            text = 'Peak 4').grid(row = 6, column = 2) 
        tk.Label(TabList[num][1].AlgFrame, 
            text = 'Peak 5').grid(row = 7, column = 2) 
        tk.Label(TabList[num][1].AlgFrame, 
            text = 'Peak 6').grid(row = 8, column = 2)          

        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIdown1, relief = 'sunken',
                  borderwidth = 2).grid(row = 3, column = 0)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIup1, relief = 'sunken',
                  borderwidth = 2).grid(row = 3, column = 1)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIdown2, relief = 'sunken',
                  borderwidth = 2).grid(row = 4, column = 0)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIup2, relief = 'sunken',
                  borderwidth = 2).grid(row = 4, column = 1)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIdown3, relief = 'sunken',
                  borderwidth = 2).grid(row = 5, column = 0)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIup3, relief = 'sunken',
                  borderwidth = 2).grid(row = 5, column = 1)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIdown4, relief = 'sunken',
                  borderwidth = 2).grid(row = 6, column = 0)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIup4, relief = 'sunken',
                  borderwidth = 2).grid(row = 6, column = 1)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIdown5, relief = 'sunken',
                  borderwidth = 2).grid(row = 7, column = 0)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIup5, relief = 'sunken',
                  borderwidth = 2).grid(row = 7, column = 1)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIdown6, relief = 'sunken',
                  borderwidth = 2).grid(row = 8, column = 0)
        tk.Entry(TabList[num][1].AlgFrame, textvariable = TabList[num][1].ROIup6, relief = 'sunken',
                  borderwidth = 2).grid(row = 8, column = 1)
        
        tk.Button(TabList[num][1].AlgFrame,text = 'Search', 
                  command = ROI_Select_Alg).grid(row = 9, column = 0)
        tk.Button(TabList[num][1].AlgFrame,text = 'Remove Unchecked',
                  command = Unchecked_Results).grid(row = 9, column = 1)
        tk.Button(TabList[num][1].AlgFrame,text = 'Remove All',
                  command = lambda: ClearWidget('Results', 1)).grid(row = 9, column = 2)

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

    wng.warning.destroy()

#############################################################################
# Esta funcao deixa dar upload ou apagar ficheiros para a pasta de dados
# permanentes que utiliza
############################################################################# 
def File_Manager(Choice, Nature, Action):

    if Choice == 'Source':
        
        if Nature == 1:
            if Action == 1:
                filename = fd.askopenfilename(filetypes = (('Text Files', '*.txt'), ('All Files', '*.*')), 
                                              title = 'Add Alpha Source Energy File')
                
                dir = os.getcwd()
                dir = dir + '\Files\Sources\Values'
                if not filename:
                    pass
                else:
                    copy2(filename, dir, follow_symlinks=True)

            else:
                domain = 'Files\Sources\Values'
                dir = os.scandir(domain)
                wng.popup('Delete Alpha Source Files')
                name_list = []
                delete_list = []

                for entry in dir:
                    if entry.is_file():
                        temp = (os.path.splitext(entry.name))
                        name_list.append(temp[0])
                        delete_list.append(tk.IntVar())

                tk.Label(wng.warning, text = 'Files available for deletion\n').pack()

                for i in range(0, len(name_list)):
                    tk.Checkbutton(wng.warning, text = name_list[i], variable = delete_list[i],
                                   onvalue = 1, offvalue = 0).pack()
                                        
                tk.Button(wng.warning, command = lambda: Delete(delete_list, name_list, domain, '.txt'),
                           text = 'Delete Files').pack()
                
                tk.Button(wng.warning, command =  lambda: wng.warning.destroy(), 
                          text = 'Return').pack()
                
            Dir = os.scandir('Files\Sources\Values')
            source_list.clear()
            for entry in Dir:
                if entry.is_file():
                    temp = (os.path.splitext(entry.name))
                    source_list.append(temp[0])
                
        elif Nature == 0:
            if Action == 1:
                filename = fd.askopenfilename(filetypes = (('Image Files', '.jpg .jpeg .pgn'), 
                                                           ('All Files', '*.*')), 
                                              title = 'Add Alpha Source Energy Decay Image')
                dir = os.getcwd()
                dir = dir + '\Files\Sources\Images\\'
                if not filename:
                    pass
                else:
                    copy2(filename, dir, follow_symlinks=True)
            
            else:
                domain = 'Files\Sources\Images'
                dir = os.scandir(domain)
                wng.popup('Delete Alpha Source Chain Images')
                name_list = []
                delete_list = []

                for entry in dir:
                    if entry.is_file():
                        temp = (os.path.splitext(entry.name))
                        name_list.append(temp[0])
                        delete_list.append(tk.IntVar())

                tk.Label(wng.warning, text = 'Files available for deletion\n').pack()

                for i in range(0, len(name_list)):
                    tk.Checkbutton(wng.warning, text = name_list[i], variable = delete_list[i],
                                   onvalue = 1, offvalue = 0).pack()
                                        
                tk.Button(wng.warning, command = lambda: Delete(delete_list, name_list, domain, '.txt'),
                           text = 'Delete Files').pack()
                
                tk.Button(wng.warning, command =  lambda: wng.warning.destroy(), 
                          text = 'Return').pack()

    elif Choice == 'Material':
        if Action == 1:
            filename = fd.askopenfilename(filetypes = (('Text Files', '*.txt'), ('All Files', '*.*')), 
                                            title = 'Add Material File')
            dir = os.getcwd()
            dir = dir + '\Files\Materials\\'
            if not filename:
                pass
            else:
                copy2(filename, dir, follow_symlinks=True)
        
        else:
            domain = 'Files\Materials'
            dir = os.scandir(domain)
            wng.popup('Delete Material Files')
            name_list = []
            delete_list = []

            for entry in dir:
                if entry.is_file():
                    temp = (os.path.splitext(entry.name))
                    name_list.append(temp[0])
                    delete_list.append(tk.IntVar())

            tk.Label(wng.warning, text = 'Files available for deletion\n').pack()

            for i in range(0, len(name_list)):
                tk.Checkbutton(wng.warning, text = name_list[i], variable = delete_list[i],
                                onvalue = 1, offvalue = 0).pack()
                                    
            tk.Button(wng.warning, command = lambda: Delete(delete_list, name_list, domain, '.txt'),
                        text = 'Delete Files').pack()
            
            tk.Button(wng.warning, command =  lambda: wng.warning.destroy(), 
                        text = 'Return').pack()
            
        Dir = os.scandir('Files\Materials')
        materials_list.clear()
        for entry in Dir:
            if entry.is_file():
                temp = (os.path.splitext(entry.name))
                materials_list.append(temp[0])

    for i in range(0, len(TabTracker)):
        if TabTracker[i] < 0:
            TabList[i][1].Source_Menu.destroy()
            TabList[i][1].Source_Menu = tk.OptionMenu(TabList[i][1].SourceFrame, TabList[i][1].Source, 
                          *source_list, command = SourceReader)
            TabList[i][1].Source_Menu.grid(row = 1, columnspan = 2)

        elif TabTracker[i] > 0:
            TabList[i][1].Mat_Menu.destroy()
            TabList[i][1].Mat_Menu = tk.OptionMenu(TabList[i][1].SourceFrame, 
                                                   TabList[i][1].Mat, *materials_list)
            TabList[i][1].Mat_Menu.grid(row = 1, columnspan = 2)

#############################################################################
# Permite guardar os resultados todos obtidos no programa para um txt
#############################################################################
def Save_Results():
    
    domain = (('Text Files', '*.txt'), ('All Files', '*.*'))
    file = fd.asksaveasfile(title = 'Save Results', initialdir = ",", filetypes = domain, defaultextension = ".txt")

    if file:
        file.write('The Results calculated by NUC-RIA\'s ARC-TF were the following:\n\n')
        for i in range(0, len(TabTracker)):
            if TabTracker[i] < 0:
                file.write('Calibration Trial ' + str(-TabTracker[i]) + '\n\n')
                file.write('The detected Peaks and Counts were: \n \n')
                Peaks = File_Reader(TabList[i][3], ',', 'Yes', 'No')

                for j in range(0, len(Peaks)):
                    file.write('Channel: ' + str("{:.1f}".format(Peaks[j][0])) + '\tCounts: ' + str("{:.1f}".format(Peaks[j][1])) + '\n')

                file.write('\nThe Radiation source used for this trial was: ' + 
                           TabList[i][1].Source.get() + '\n\n')
                Regression = File_Reader(TabList[i][4], '0', 'String', 'No', )
                file.write('The Linear Regression calculated was the following, using ' + 
                           Regression[0] + ' units.\n\n')
                
                file.write('The Slope: ' + '%.*f' % (int(Regression[5]), float(Regression[1])) + ' ' +
                           u"\u00B1" + ' ' + '%.*f' % (int(Regression[5]), float(Regression[2])) + '\n')
                file.write('The Intersect: ' + '%.*f' % (int(Regression[6]), float(Regression[3])) + ' ' +
                           u"\u00B1" + ' ' + '%.*f' % (int(Regression[6]), float(Regression[4])) + '\n')
                file.write('_______________________________________________________________\n\n')

        for i in range(0, len(TabTracker)):
            if TabTracker[i] > 0:
                file.write('Material Trial ' + str(TabTracker[i]) + '\n\n')
                file.write('The detected Channels, Counts and respectful Thickness approximation were:'
                            + '\n\n')
                Peaks = File_Reader(TabList[i][3], ',', 'Yes', 'No')
                Peaks.sort()
                thickness = File_Reader(TabList[i][4], '0', 'String', 'No')

                units_list = ['nm', '\u03bcm',
                  '\u03bcg' + ' cm' + '{}'.format('\u207B' + '\u00b2'),
                  '10' + '{}'.format('\u00b9' + '\u2075') + ' Atoms' 
                   + ' cm' + '{}'.format('\u207B' + '\u00b3')]

                units_values = [10.0**9, 10.0**6, 0.0, -1.0]
                index = units_values.index(TabList[i][1].units.get())

                for j in range(0, len(Peaks)):
                    file.write('Channel: ' + str("{:.1f}".format(Peaks[j][0])) + '\tCounts: ' + str("{:.1f}".format(Peaks[j][1])) + 
                               '\tThickness: ' + '%.*f' % (int(thickness[-1]), float(thickness[j])) + 
                               ' ' + units_list[index] + '\n')

                file.write('\nThe average thickness, of material ' + 
                           TabList[i][1].Mat.get() + ' was calculated to be: ' + '('
                           '%.*f' % (int(thickness[-1]), float(thickness[j + 1])) + 
                           ' ' + u"\u00B1 " + 
                           '%.*f' % (int(thickness[-1]), float(thickness[j + 2])) + ') ' + 
                           units_list[index] + '\n\n')
                
                '%.*f' % (int(Regression[6]), float(Regression[4]))
                
                file.write('_______________________________________________________________\n\n')
        
#############################################################################
# A classe do esqueleto, onde esta a barra de ferramentas e a janela principal
# do programa
#############################################################################
class Skeleton:

    def __init__(self):

        ############# A Janela Mae ##############
        self.main = tk.Tk()
        self.main.title('ARC_TF')
        self.main.state('zoomed')
        self.main.configure(background = 'dark grey')

        ############## A barra de Ferramentas e opcoes #################
        self.menu = tk.Menu(self.main)
        self.main.config(menu = self.menu)

            # O tipico File Menu. Hao de haver mais opcoes no futuro
        __file_menu = tk.Menu(self.menu, tearoff = False) 
        self.menu.add_cascade(label = 'File', menu = __file_menu)
        __file_menu.add_command(label = 'Plot Data', command = DataUploader)
        __file_menu.add_command(label = "Save Results", command = Save_Results)
        __file_menu.add_separator()
        __file_menu.add_command(label = 'Remove Current Plot',
                                command = lambda: ClearWidget('Graphic', 0))
        __file_menu.add_command(label = 'Remove Algorithm Results', 
                                command = lambda: ClearWidget('Results', 0))
        __file_menu.add_command(label = "Reset All Data from Current Tab",
                                command = lambda: ClearWidget('Everything', 1) )
        __file_menu.add_separator()
        __file_menu.add_command(label = 'Exit', command = self.main.quit)

            # O bloco de definicoes do programa
        self.menu.add_command(label = 'Settings', command = lambda: wng.Settings())

            # Gere as tabs do programa
        __tabs_menu = tk.Menu(self.menu, tearoff = False)
        self.menu.add_cascade(label = 'Manage Tabs', menu = __tabs_menu)
        __tabs_menu.add_command(label = 'Add Calibration Tab', command = lambda: Tabs.tab_change(1))
        __tabs_menu.add_command(label = 'Add Material Tab', command = lambda: Tabs.tab_change(2))
        __tabs_menu.add_separator()
        __tabs_menu.add_command(label = 'Remove Current Tab', command = lambda: Tabs.tab_change(4))
        
            # Gere os documentos da base de dados do programa
        __files_data = tk.Menu(self.menu, tearoff = False)
        self.menu.add_cascade(label = 'Manage Data Files', menu = __files_data)
        __files_data.add_command(label = 'Add Alpha Source File', 
                                 command = lambda: File_Manager('Source', 1, 1))
        __files_data.add_command(label = 'Remove Alpha Source File',
                                 command = lambda: File_Manager('Source', 1, 0))
        __files_data.add_separator()
        __files_data.add_command(label = 'Add Alpha Source Image',
                                 command = lambda: File_Manager('Source', 0, 1))
        __files_data.add_command(label = 'Remove Alpha Source Image', 
                                 command = lambda: File_Manager('Source', 0, 0))
        __files_data.add_separator()
        __files_data.add_command(label = 'Add Material File',
                                 command = lambda: File_Manager('Material', 0, 1))
        __files_data.add_command(label = 'Remove Material File',
                                 command = lambda: File_Manager('Material', 0, 0))
        
            # Abre o html de Help
        self.menu.add_command(label = 'Help', command = lambda: wng.Help())

            # Abre o link para o Readme
        self.menu.add_command(label = "About", command = lambda: webbrowser.open(
            'https://github.com/AlexVnGit/GUI_thin_films/blob/master/README.md', new = 1))
        
                
    def run(self):
        #### O metodo que permite o programa manter-se aberto
        self.main.mainloop()

#############################################################################
# A class dos avisos e popups. E facilmente reciclavel e versatil
#############################################################################
class Warnings:
    def popup(self, name):

        # A class cria poopups que tem que ser fechados antes de voltarem para o skeleton
        # funciona bem para avisos, mas tambem para menus com opcoes obrigatorias de submeter
        self.warning = tk.Toplevel(window.main)
        self.warning.title(name)
        self.warning.geometry('700x300')
        self.warning.grab_set()

    def Images(self, name, picture, site): # Para as imagens dos decaimentos

        self.decay = tk.Toplevel(window.main)
        self.decay.title(name)
        self.decay.geometry('1000x600')

        load = Image.open(picture)  
        load_resize = load.resize((700,500))
        img = ImageTk.PhotoImage(load_resize)
        label = tk.Label(self.decay, image = img)
        label.image = img
        label.pack()
        site_label = tk.Label(self.decay, text = 'From:  ' + site)
        site_label.pack()

    def Settings(self):

        num = Current_Tab()

        self.configuration = tk.Toplevel(window.main)
        self.configuration.title('Settings')
        self.configuration.geometry('700x400')
        self.configuration.grab_set()
        
        self.parameter = ttk.Notebook(self.configuration)
        self.parameter.pack(expand = True, fill = 'both')
        self.parameter.enable_traversal()
        self.general_tab = tk.Frame(self.parameter)
        self.algorithm_tab = tk.Frame(self.parameter)
        self.parameter.add(self.general_tab, text = 'General Settings')
        self.parameter.add(self.algorithm_tab, text = 'Algorithm Settings')

        self.general_tab.columnconfigure(0, weight = 1)
        self.general_tab.columnconfigure(1, weight = 1)
        self.general_tab.columnconfigure(2, weight = 1)

        self.algorithm_tab.columnconfigure(0, weight = 1)
        self.algorithm_tab.columnconfigure(1, weight = 1)
        self.algorithm_tab.columnconfigure(2, weight = 1)
        self.algorithm_tab.columnconfigure(3, weight = 1)        

        self.energy_value = TabList[num][1].energy.get()
        self.unit_value = TabList[num][1].units.get()
        self.peak_interval = TabList[num][1].peaks_widths.get()
        self.cut_low_energy = TabList[num][1].channels_cut.get()

        self.buttons_frame = tk.Frame(self.configuration)
        self.buttons_frame.pack()

        def close(choice):

            num = Current_Tab()

            if choice == 0:
                TabList[num][1].units.set(wng.unit_value)
                TabList[num][1].energy.set(wng.energy_value)
                TabList[num][1].peaks_widths.set(wng.peak_interval)
                TabList[num][1].channels_cut.set(wng.cut_low_energy)
                wng.configuration.destroy()

            elif choice == 1:
                try:
                    TabList[num][1].peaks_widths.set(wng.entry2.get())
                    TabList[num][1].channels_cut.set(wng.entry1.get())

                except ValueError:
                    TabList[num][1].peaks_widths.set(wng.peak_interval)
                    TabList[num][1].channels_cut.set(wng.cut_low_energy)   

                wng.configuration.destroy()

            elif choice == 2:   
                try:
                    TabList[num][1].peaks_widths.set(wng.entry2.get())
                    TabList[num][1].channels_cut.set(wng.entry1.get())

                except ValueError:
                    TabList[num][1].peaks_widths.set(wng.peak_interval)
                    TabList[num][1].channels_cut.set(wng.cut_low_energy)    
                    
                Energy_settings.set(TabList[num][1].energy.get())
                Unit_settings.set(TabList[num][1].units.get())

                Channel_cut.set(TabList[num][1].channels_cut.get())
                Peak_Width.set(TabList[num][1].peaks_widths.get())                

                wng.configuration.destroy()

        ################################################ 

        tk.Label(self.general_tab, text = ' ').grid(row = 0, 
                                                                       column = 0, pady = 10, padx = 10)
        tk.Label(self.general_tab, text = 'Specify Energy Units').grid(row = 1, 
                                                                       column = 0, pady = 10, padx = 10)
        tk.Label(self.general_tab, 
                 text = 'Specify Thickness Units').grid(row = 1, column = 2, pady = 10, padx = 10)

        tk.Radiobutton(self.general_tab, text = 'kev', 
                       variable = TabList[num][1].energy, value = 1).grid(row = 2, column = 0)
        tk.Radiobutton(self.general_tab, text = 'Mev', 
                       variable = TabList[num][1].energy, value = 1000).grid(row = 3, column = 0)
        
        tk.Radiobutton(self.general_tab, text = 'nm', 
                       variable= TabList[num][1].units, value = 10.0**9).grid(row = 2, column = 2)
        tk.Radiobutton(self.general_tab, text = '\u03bcm', 
                       variable= TabList[num][1].units, value = 10.0**6).grid(row = 3, column = 2)
        tk.Radiobutton(self.general_tab, text = '\u03bcg' + ' cm' + '{}'.format('\u207B' + '\u00b2'), 
                       variable= TabList[num][1].units, value = 0.0).grid(row = 4, column = 2)
        tk.Radiobutton(self.general_tab, text = '10' + '{}'.format('\u00b9' + '\u2075') + ' Atoms' +
                       ' cm' + '{}'.format('\u207B' + '\u00b3'), 
                       variable= TabList[num][1].units, value = -1.0).grid(row = 5, column = 2)
        
        ########################################################################################

        tk.Label(self.algorithm_tab, text = 'Threshold Input Algorithm Settings').grid(
            row = 0, columnspan = 3, pady = 10)
        tk.Label(self.algorithm_tab, text = 'Low Energy Cut: ').grid(
            row = 1, column = 1, pady = 10)
        self.entry1 = tk.Entry(self.algorithm_tab)
        self.entry1.grid(row = 1, column = 2, pady = 10)
        self.entry1.insert(0, TabList[num][1].channels_cut.get())
        tk.Label(self.algorithm_tab, text = 'Approximate width of peaks: ').grid(
            row = 2, column = 1, pady = 10)
        self.entry2 = tk.Entry(self.algorithm_tab)
        self.entry2.grid(row = 2, column = 2, pady = 10)
        self.entry2.insert(0, TabList[num][1].peaks_widths.get())

        ########################################################################################

        tk.Button(self.buttons_frame, text = 'Cancel and Return', command = lambda: 
                  close(0)).grid(row = 1, column = 2, padx = 5, pady = 5)
        tk.Button(self.buttons_frame, text = 'Apply to current Tab', command = lambda: 
                  close(1)).grid(row = 1, column = 0, padx = 5, pady = 5)
        tk.Button(self.buttons_frame, text = 'Apply to all Tabs', command = lambda: 
                  close(2)).grid(row = 1, column = 1, padx = 5, pady = 5)

    def Help(self):

        try:
            wng.helping.destroy()

        except:
            pass

        with open('Files\Help.txt', 'r') as OpenFile: # Abre (e fecha) o documento
            lines = OpenFile.read() # Le os dados como string

        self.helping = tk.Toplevel(window.main)
        self.helping.title('Help')
        self.helping.resizable(0,0)
        self.frame = tk.Frame(self.helping)
        self.frame.pack(expand= True, fill= 'both')
        self.canvas = tk.Canvas(self.frame)
        self.frame2 = tk.Frame(self.canvas)
        self.scrollbar = tk.Scrollbar(self.frame)
        self.canvas.config(yscrollcommand = self.scrollbar.set, highlightthickness = 0 )
        self.scrollbar.config(orient = tk.VERTICAL, command = self.canvas.yview)
        self.canvas.grid(row = 1, column = 0)
        self.scrollbar.grid(row = 1, column = 2, sticky = 'ns')
        self.scrollbar.place
        self.canvas.create_window(0, 0, window = self.frame2, anchor = tk.NW)
        self.frame2.bind('<Configure>', 
                              lambda e: self.canvas.configure(
                                  scrollregion = self.canvas.bbox('all'), width = e.width))
        
        menu = tk.Label(self.frame2, text = lines).grid()
        
############################################################################
# A class das Tabs. Inclui a estrutura propria do Notebook - widget de 
# separadores; as tabs dos resultados e de adicionar tabs; e a estrutura que
# muda a forma de tabs de calibracao e materiais
############################################################################
class Tabs:

    Counter_Mat = 0
    Counter_Calib = 0

    def First_Tabs(self):

        ####### A variavel do Notebook ##########
        self.notebook = ttk.Notebook(window.main)   
        self.notebook.pack(expand = True, fill = 'both') # Expande a frame da tab ate ao final
                                                         # da janela
        self.notebook.enable_traversal() # Permite o uso de Ctrl+Tab para circular entre tabs
        
        ########### As frames principais - Os resultados e a frame de adicionar
        self.CRFrame = tk.Frame(self.notebook, bg = 'dark grey')
        self.PlusFrame = tk.Frame(self.notebook, bg = 'dark grey')
        self.CRFrame.columnconfigure(0, weight = 3)
        self.CRFrame.columnconfigure(1, weight = 3)
        self.CRFrame.rowconfigure(0, weight = 3)
        self.CRFrame.rowconfigure(1, weight = 3)

        ############# Aqui adicionam-se as frames iniciadas acima
        self.notebook.add(self.CRFrame, text = 'Final Results')
        self.notebook.add(self.PlusFrame, text = '+')
        
        ############# Frames onde irao ser inseridos os resultados finais e scrollbars
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
        
        ########### Variavel para contar o numero de separadores 
        self.value = 0

        ############ Este comando adiciona o evento a funcao, para adicionar frames
        self.notebook.bind("<<NotebookTabChanged>>", handleTabChange)      

    def AnalysisTab(self, choice):

        ### Configuracao de geometria e Frames comuns a Calib e Material Tabs  ####
        TabList[Notebook.value][0].columnconfigure(0, weight = 4)
        TabList[Notebook.value][0].columnconfigure(1, weight = 1)
        TabList[Notebook.value][0].rowconfigure(0, weight = 1)
        TabList[Notebook.value][0].rowconfigure(1, weight = 1)

        self.GraphicFrame = tk.Frame(TabList[Notebook.value][0], borderwidth = 5, relief = 'ridge')
        self.GraphicFrame.grid(column = 0, row = 0, sticky = "nw", pady = 5, columnspan = 2)

        self.DataFrame = tk.Frame(TabList[Notebook.value][0], borderwidth = 5, relief = 'ridge')
        self.DataFrame.grid(column = 2, row = 0, sticky = "ne", pady = 5)

        self.SourceFrame = tk.Frame(TabList[Notebook.value][0], borderwidth = 5, relief = 'ridge')
        self.SourceFrame.grid(column = 1, row = 0, sticky = "ne", pady = 5)

        self.Extra_Frame = tk.Frame(TabList[Notebook.value][0], borderwidth = 5, relief = 'ridge')

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
            self.variable11, self.variable12, self.variable13, self.variable14, self.variable15
        ]
    
        tk.Label(self.DataFrame, text = 'Analysis Method Selected: ').grid(row = 0, columnspan = 2)
        Algs = ["Manual Selection", "Threshold Input", "ROI Select"]
        tk.OptionMenu(self.DataFrame, self.Algorithm_Method, *Algs, command = Method).grid(row = 1, columnspan = 2)

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
            
        if choice == 1:
            CalibTab(self)
            
        elif choice == 2:
            MatTab(self)

    @staticmethod
    def tab_change(num):

        value = Current_Tab()
        index = len(Notebook.notebook.tabs()) - 1             

        if num == 1:
            Tabs.Counter_Calib -= 1
            Data = "Temp\Data" + str(Tabs.Counter_Calib) + ".txt"
            Analysis = "Temp\Analysis" + str(Tabs.Counter_Calib) + ".txt"
            Result =  "Temp\Result" + str(Tabs.Counter_Calib) + ".txt"
            ROIs = "Temp\ROIs" + str(Tabs.Counter_Calib) + ".txt"
            TabList.append([tk.Frame(Notebook.notebook, bg = 'dark grey'), Tabs(),  Data,
                        Analysis,  Result, Plot(), ROIs]) 

            TabTracker.append(Tabs.Counter_Calib)
            TabList[Notebook.value][1].AnalysisTab(1)
            Notebook.notebook.insert(index, TabList[Notebook.value][0], 
                                     text = "Calibration Trial " + str(-Tabs.Counter_Calib))
            Notebook.notebook.select(index)
            Notebook.value += 1
            try:
                wng.warning.destroy()
            except:
                ()

        elif num == 2:
            Tabs.Counter_Mat += 1
            Data = "Temp\Data" + str(Tabs.Counter_Mat) + ".txt"
            Analysis = "Temp\Analysis" + str(Tabs.Counter_Mat) + ".txt"
            Result =  "Temp\Result" + str(Tabs.Counter_Mat) + ".txt"
            ROIs = "Temp\ROIs" + str(Tabs.Counter_Mat) + ".txt"
            TabList.append([tk.Frame(Notebook.notebook, bg = 'dark grey'), Tabs(),  Data,
                        Analysis,  Result, Plot(), ROIs]) 

            TabTracker.append(Tabs.Counter_Mat)
            TabList[Notebook.value][1].AnalysisTab(2)
            Notebook.notebook.insert(index, TabList[Notebook.value][0], 
                                     text = "Material Trial " + str(Tabs.Counter_Mat))
            Notebook.notebook.select(index)
            Notebook.value += 1
            try:
                wng.warning.destroy()
            except:
                ()
        
        elif num == 3:
            Notebook.notebook.select(index - 1)
            wng.warning.destroy()

        elif num == 4:
            if Notebook.notebook.select() == '.!notebook.!frame' or Notebook.notebook.select() == '.!notebook.!frame2':
                wng.popup('Bad Tab Deletion')
                tk.Label(wng.warning, text = '\n This Tab cannot be deleted.\nPlease delete Analysis Tabs').pack()
                tk.Label(wng.warning, text = '\n\n').pack()
                tk.Button(wng.warning, text = 'Return', command = lambda: wng.warning.destroy()).pack()

            else:
                Notebook.notebook.forget("current")
                Notebook.notebook.select(index - 2)

                if os.path.isfile(TabList[value][2]) == True:
                    os.remove(TabList[value][2])
                if os.path.isfile(TabList[value][3]) == True:
                    os.remove(TabList[value][3])
                if os.path.isfile(TabList[value][4]) == True:
                    os.remove(TabList[value][4])

                TabList.pop(value)

                Notebook.value -= 1
                TabTracker.pop(value)
                ClearWidget('Final', 0)
                Final_Results(0)

    def RenameTab(name):

        tab_num = Current_Tab()
        Notebook.notebook.tab(tab_num+1, text = str(name))

        return


       
###########################################################################
# Esta classe recebe os dados dos ficheiros externos
# e insere os graficos na frame grande do GUI.
# Ao mesmo tempo cria um txt para outras funções
# acederem aos dados
# Se estiver selecionado o threshold input, mostra uma linha do valor
############################################################################
class Plot:

    def Structure(self, File, Name):

        self.Channel = []    #Lista vazia para guardar o Channel
        self.Counts = []     #Lista vazia para guardar os Counts
        self.line = []

        num = Current_Tab()
        total_sum = 0
        j = 0

        ClearWidget('Graphic', 0)
        TabList[num][1].GraphicFrame.grid(column = 0, row = 0, sticky = "nw", pady = 5, columnspan = 2)

        Data = open(TabList[num][2], "w")

        if Name[-4:] == ".mca":         #Por enquanto esta configurado para os ficheiros 
                                        # da maquina para AEL. Se for configurado RBS
                                        #ha-de-se incluir outro if.
            for i in range(12, len(File) - 1):         #### FOI ALTERADO POR CAUSA DA ROI NOS FICHEIROS !!!!!!
                self.Counts.append(int(File[i]))
                total_sum = total_sum + self.Counts[j]
                self.Channel.append(i-11)               #### FOI ALTERADO POR CAUSA DA ROI NOS FICHEIROS !!!!!!
                Data.write(str(self.Counts[i-12])+"\n") #### FOI ALTERADO POR CAUSA DA ROI NOS FICHEIROS !!!!!!
                j += 1

            # Ciclo for para adquirir os valores dos dados

        Data.close()

        TabList[num][1].Total_Counts.set(total_sum)
        TabList[num][1].Extra_Frame.grid(column = 0, row = 1, sticky = "nw")
        tk.Label(TabList[num][1].Extra_Frame, text = TabList[num][1].Real_Time.get()).grid(row = 0)
        tk.Label(TabList[num][1].Extra_Frame, text = 'Total Sum of Counts is: ' + 
                 str(TabList[num][1].Total_Counts.get())).grid(row = 1)

        if TabTracker[num] < 0:
            self.Title = 'Calibration Trial ' + str(-TabTracker[num])

        elif TabTracker[num] > 0: 
            self.Title = 'Material Trial ' + str(TabTracker[num])


        self.figure = Figure(figsize = (6,4), dpi = 100) #A figura contem o grafico
        self.figure_canvas = FigureCanvasTkAgg(self.figure, TabList[num][1].GraphicFrame) #A class FigureCanvasTkAgg
        #liga o matplotlib ao tkinter

        NavigationToolbar2Tk(self.figure_canvas, TabList[num][1].GraphicFrame) # Esta linha permite que as ferramentas
        #do matplotlip aparecam na interface do tkinter

    def subplots(self):
        #Aqui inicia-se o grafico com os dados e os eixos

        self.axes = self.figure.add_subplot() 
        self.axes.plot(self.Channel, self.Counts, '.', markersize = 7, label = 'Run')
        self.axes.set_title(self.Title)
        self.axes.set_xlabel('Channel')
        self.axes.set_ylabel('Counts')

        self.figure.canvas.mpl_connect('button_press_event', onclick)

        #Por fim, acrescenta-se a geometria do tkinter
        self.figure_canvas.get_tk_widget().pack()

    def destroyer(self):

        if self.line:
            self.line.pop().remove()
            self.figure_canvas.draw()

    def threshold(self, height):

        if self.line:
            self.line.pop().remove()

        self.line.append(self.axes.axhline(y = height, color = 'r', linestyle = '-'))
        self.figure_canvas.draw()

############################################################################
Dir = os.scandir('Files\Sources\Values')
source_list = []
for entry in Dir:
    if entry.is_file():
        temp = (os.path.splitext(entry.name))
        source_list.append(temp[0])

Dir = os.scandir('Files\Materials')
materials_list = []
for entry in Dir:
    if entry.is_file():
        temp = (os.path.splitext(entry.name))
        materials_list.append(temp[0])

############ Variaveis Estruturais #############################
wng = Warnings()
window = Skeleton()
Notebook = Tabs()
Notebook.First_Tabs()

############## Tabs variaveis para serem criadas ###############
TabList = []
TabTracker = []

Energy_settings = tk.IntVar()
Energy_settings.set(1000)
Unit_settings = tk.DoubleVar()
Unit_settings.set((10**9))

Channel_cut = tk.IntVar()
Channel_cut.set(100)

Peak_Width = tk.IntVar()
Peak_Width.set(35)

#############################################################################################
Tabs.tab_change(1)
Notebook.notebook.select(1)

os.mkdir('Temp') # Pasta onde serao guardados os ficheiros temporarios
window.run()
##############################################################################################

for i in range(Notebook.value):
    if os.path.isfile(TabList[i][2]) == True:
        os.remove(TabList[i][2]) # Apaga os dados adquiridos quando se faz plot
    if os.path.isfile(TabList[i][3]) == True:
        os.remove(TabList[i][3]) # Apaga os dados dos resultados dos algoritmos
    if os.path.isfile(TabList[i][4]) == True:
        os.remove(TabList[i][4]) # Apaga os resultados das regressoes lineares

os.rmdir('Temp') # Apaga a pasta Temp