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

############################################################################################
# Le os resultados finais e exibe os na primeira tab
############################################################################################
def Final_Results(tracker):

    num = Current_Tab()
    method = TabList[num][1].Algorithm_Method.get()

    if tracker < 0 and method == 'ROI Select':
        LinearizeWithErrors()
    
    elif tracker < 0 and method != 'ROI Select':
        Linearize()

    elif tracker > 0 and method != 'ROI Select':
        Final_Calculation()

    elif tracker > 0 and method == 'ROI Select':
        ROI_Thick_Calculation()

    elif tracker == 0:
        pass        
        
    for i in range(0, len(TabTracker)):
        if (tracker < 0 or tracker == 0) and TabTracker[i] < 0:
            if os.path.isfile(TabList[i][4]) == True:

                if TabList[i][1].energy.get() == 1000:
                    unit_energy = 'MeV'

                elif TabList[i][1].energy.get() == 1:
                    unit_energy = 'keV'

                Results = File_Reader(TabList[i][4], '0', 'String', 'No')
                tk.Label(Notebook.Calib_Result2, text = 'Calibration Trial ' + 
                        str(-TabTracker[i]) + ' - ' +
                        TabList[i][1].Source.get()).grid(row = 4 * i + i, columnspan = 3)
                tk.Label(Notebook.Calib_Result2, 
                        text = '(' + unit_energy + ')').grid(row = 4 * i + i + 1, column = 0)
                tk.Label(Notebook.Calib_Result2, 
                        text = 'Values').grid(row = 4 * i + i + 1, column = 1)
                tk.Label(Notebook.Calib_Result2, 
                        text = 'Uncertainty').grid(row = 4 * i + i + 1, column = 2)
                tk.Label(Notebook.Calib_Result2, 
                        text = 'Slope').grid(row = 4 * i + i + 2, column = 0)
                tk.Label(Notebook.Calib_Result2, 
                        text = 'Intersect').grid(row = 4 * i + i + 3, column = 0)
                tk.Label(Notebook.Calib_Result2, 
                        text = '').grid(row = 4 * i + i + 4, column = 0)

                tk.Label(Notebook.Calib_Result2, 
                        text = '%.*f' %(int(Results[5]), float(Results[1]))).grid(
                            row = 4 * i + i + 2, column = 1)
                tk.Label(Notebook.Calib_Result2, 
                        text = '%.*f' %(int(Results[5]), float(Results[2]))).grid(
                            row = 4 * i + i + 2, column = 2)
                tk.Label(Notebook.Calib_Result2, 
                        text = '%.*f' %(int(Results[6]), float(Results[3]))).grid(
                            row = 4 * i + i + 3, column = 1)
                tk.Label(Notebook.Calib_Result2, 
                        text = '%.*f' %(int(Results[6]), float(Results[4]))).grid(
                            row = 4 * i + i + 3, column = 2)
            else:                    
                Notebook.calib_canvas.update_idletasks()    
                Notebook.calib_canvas.config(scrollregion = Notebook.Calib_Result2.bbox())
                Notebook.Calib_Result2.bind('<Configure>', 
                              lambda e: Notebook.calib_canvas.configure(
                                  scrollregion = Notebook.calib_canvas.bbox('all'), width = e.width)) 

        if (tracker > 0 or tracker == 0) and TabTracker[i] > 0:
            if os.path.isfile(TabList[i][4]) == True and os.path.isfile(TabList[i][3]) == True:

                Results = File_Reader(TabList[i][4], '0', 'Yes', 'No')
                Peaks = File_Reader(TabList[i][3], ',', 'Yes', 'No')
                Peaks.sort()

                units_list = ['nm', '\u03bcm',
                  '\u03bcg' + ' cm' + '{}'.format('\u207B' + '\u00b2'),
                  '10' + '{}'.format('\u00b9' + '\u2075') + ' Atoms' 
                   + ' cm' + '{}'.format('\u207B' + '\u00b3')]

                units_values = [10.0**9, 10.0**6, 0.0, -1.0]
                index = units_values.index(TabList[i][1].units.get())

                size = len(Peaks)

                for j in range(0, size):
                    if j == 0:

                        tk.Label(Notebook.Mat_Result2, text = 'Material Trial ' + 
                            str(TabTracker[i]) + ' - ' +
                            TabList[i][1].Mat.get()).grid(row = (4 + size) * i + i + j , columnspan = 2)
                        tk.Label(Notebook.Mat_Result2, 
                                 text = 'Peak Centroid').grid(row = (4 + size) * i + i + j + 1, column = 0)
                        tk.Label(Notebook.Mat_Result2, 
                                 text = 'Thickness').grid(row = (4 + size) * i + i + j + 1, column = 1)
                    
                    tk.Label(Notebook.Mat_Result2, 
                             text = str("{:.1f}".format(Peaks[j][0]))).grid(row = (4 + size) * i + i + j + 2, 
                                                          column = 0)
                    tk.Label(Notebook.Mat_Result2, 
                             text = '%.*f' % (int(Results[-1]), Results[j]) + 
                             ' ' + units_list[index] ).grid(row = (4 + size) * i + i + j + 2, 
                                                                  column = 1)
                    
                tk.Label(Notebook.Mat_Result2, 
                                 text = '\nAverage').grid(row = (4 + size) * i + i + j + 3, column = 0)
                tk.Label(Notebook.Mat_Result2, 
                                 text = 'Uncertainty').grid(row = (4 + size) * i + i + j + 4, column = 0)
                tk.Label(Notebook.Mat_Result2, 
                                 text = '\n' + '%.*f' % (int(Results[-1]), Results[j + 1]) 
                                 + ' ' + units_list[index]).grid(
                                     row = (4 + size) * i + i + j + 3, column = 1)
                tk.Label(Notebook.Mat_Result2, 
                                 text = '%.*f' % (int(Results[-1]), Results[j + 2]) 
                                 + ' ' + units_list[index]).grid(
                                     row = (4 + size) * i + i + j + 4, column = 1)
                tk.Label(Notebook.Mat_Result2, 
                                 text = '').grid(row = (4 + size) * i + i + j + 5, 
                                                                  columnspan = 2)
                
                Notebook.mat_canvas.update_idletasks()
                Notebook.mat_canvas.config(scrollregion = Notebook.Mat_Result2.bbox())
                Notebook.Mat_Result2.bind('<Configure>', 
                              lambda e: Notebook.mat_canvas.configure(
                                  scrollregion = Notebook.mat_canvas.bbox('all'), width = e.width)) 
            else:
                print('WDK')
        else:
            print('WTF')
###########################################################################################
# Permite a escolha de regressoes lineares por parte do utilizador
###########################################################################################
def Calib_Choice():

    num = Current_Tab()
    Measure = []

    for i in range(0, len(TabList[num][1].Regression_List)):
            TabList[num][1].Regression_List[i].set(-1)
    # Este for previne as escolhas antigas de interferir com a nova selecao de regressoes
    
    for i in range(0, len(TabTracker)):

        if TabTracker[i] < 0: # Certifica que so le ficheiros com regressoes lineares

            if os.path.isfile(TabList[i][4]) == True:
                Measure.append(TabTracker[i])
    # Neste ciclo, o measure regista quantas tabs de calibracao tem uma regressao linear,
    # ja que o tabtracker identifica as tabs de calibracao como sendo negativas,
    # o measure tera sempre valores entre -1 e -10

    if not Measure:
        wng.popup('No Linear Regressions detected')
        tk.Label(wng.warning, text = 'No linear Regressions were detected.\n\n' + 
                 'Please Perform a Calibration Trial before calculating the Film\'s Thickness.\n\n').pack()
        tk.Button(wng.warning, text = 'Return', command = lambda: wng.warning.destroy()).pack()
    # No caso do Measure estar vazio, aparece um aviso que nenhuma regressao linear foi efetuada

    else: # Caso contrario, entra o popup para selecionar quais as regressoes a utilizar
        wng.popup('Linear Regression Selection Menu')
        tk.Label(wng.warning, text = 'Please Select a Calibration Trial \n' +
                 'Choosing more than one calibration will average the slopes and intersects.\n\n').pack()

        for i in range(0, len(Measure)):
            button_Choice = tk.Checkbutton(wng.warning, 
                                           text = 'Linear Regression of Calibration Trial ' +
                                             str(-Measure[i]),
                                            variable = TabList[num][1].Regression_List[i], 
                                            onvalue = -Measure[i], offvalue = -1)
            button_Choice.pack()
        # Neste ciclo, por cada regressao linear efetuada, o utilizador pode esolher utilizar uma ou
        # mais, para a calibracao. A Regression_List guarda valores 1 ou -1 sendo que o indice desta
        # lista, e utilizado depois nos calculos finais

        tk.Button(wng.warning, text = 'Return', command = lambda: ClearWidget('Popup', 0)).pack()

################################################################################################
# Calcula a espessura por cada pico e a media das espessuras
################################################################################################
def Final_Calculation():

    num = Current_Tab()
    ClearWidget('Thickness', 0)
    TabList[num][1].ThicknessFrame.grid(row = 5, columnspan = 3, pady = 5)
    units_list = ['nm', '\u03bcm',
                  '\u03bcg' + ' cm' + '{}'.format('\u207B' + '\u00b2'),
                  '10' + '{}'.format('\u00b9' + '\u2075') + ' Atoms' 
                   + ' cm' + '{}'.format('\u207B' + '\u00b3')]

    units_values = [10.0**9, 10.0**6, 0.0, -1.0]
    index = units_values.index(TabList[num][1].units.get())

    Material_choice = TabList[num][1].Mat.get() # Determina qual o ficheiro do material a ler
    Material_choice = 'Files\Materials\\' +  Material_choice + '.txt'

    slope = 0
    intersect = 0
    points = []
    regressions_index = []
    material_data = File_Reader(Material_choice, '|', 'Yes', 'No')  #Daqui obtemos a lista que ira guardar
                                                                    # as energias e o stopping power do material em uso

# Este ciclo determina a tab de calibracao cuja regressao foi selecionada para uso
# o regressions_index guarda o indice do Tab Tracker
# Com este valor, podemos aceder ao indice do TabList para obter todos os dados que queiramos
    for i in range(0, len(TabList[num][1].Regression_List)):

        if TabList[num][1].Regression_List[i].get() != -1:
            regressions_index.append(TabTracker.index(-TabList[num][1].Regression_List[i].get()))

    i = 0
    j = 0
    Aux_Channel = []
    Temp = []

    for j in range(0, len(TabList[regressions_index[i]][1].DecayList)):
        if TabList[regressions_index[0]][1].DecayList[j].get() != -1: # Aqui, vamos buscar os valores de
                                                                        # decaimento que a utilizar para o intervalo
            Aux_Channel.append(TabList[regressions_index[0]][1].DecayList[j].get()) # de stopping powers
    # Como a fonte e a mesma para todas as calibracoes, nao importa qual delas e selecionada

    Aux_Channel.sort()
    peaks = len(Aux_Channel) # Para referencia do tamanho

    for i in range(0, len(regressions_index)):
        
        Aux = File_Reader(TabList[regressions_index[i]][4], '0', 'String', 'No') # Aqui 
        #vamos buscar as regressoes lineares para fazer uma media

        if Aux[0] == 'keV':
            placeholder1 = float(Aux[1]) * 1000
            Aux[1] = str(placeholder1)
            placeholder2 = float(Aux[3]) * 1000
            Aux[3] = str(placeholder2)

        slope = float(Aux[1]) + slope # Acumular o declive
        intersect = float(Aux[3]) + intersect # Acumular a ordenada na origem

    points = File_Reader(TabList[num][3], ',', 'No', 'No') # Analise dos picos de materiais em estudo
    points.sort() # Ordenar os picos por ordem crescente
    slope = slope / len(regressions_index)  # Buscar a media do declive
    intersect = intersect / len(regressions_index) # Buscar a ordenada na origem

    i = 0
    j = 0
    Aux.clear()     
    thickness = 0

    tk.Label(TabList[num][1].ThicknessFrame, 
             text = 'Thickness (' + units_list[index] + ')').grid(row = 0, column = 0)
    tk.Label(TabList[num][1].ThicknessFrame, text = 'Channel').grid(row = 0, column = 1)

    for j in range(0, peaks):

        Temp.append((slope * points[j][0]) + intersect )  # Calibracao dos picos do material
        uncertain = 0 # Ira devolver a incerteza da media da espessura
        summed_values = 0 # Faz o somatorio dos valores do stopping power
        i = 1

        for i in range(1, len(material_data)): # O ciclo comeca em 1 porque a linha 0 tem a densidade e 
                                                # o numero atomico
            if material_data[i][0] >= Temp[j] and material_data[i][0] <= Aux_Channel[j]:

                stopping_power = ( 1 / material_data[i][1]) # Stopping power a dividir pelo step
                summed_values = summed_values + stopping_power # O somatorio que 
                #resulta na aproximacao da espessura
                
        if index == 0 :
            summed_values = summed_values / material_data[0][1]
            summed_values = summed_values * 10000

        if index == 1:
            summed_values = summed_values / material_data[0][1]
            summed_values = summed_values * 10

        if index == 2:
            summed_values = summed_values *  0.001

        if index == 3:
            summed_values = summed_values * 1000 * ((6.02214076**(-23)) / material_data[0][0])

        Aux.append(summed_values) # Lista que guarda a espessura por perda de energia
        thickness = thickness + (Aux[j]) # Espessua media

        tk.Label(TabList[num][1].ThicknessFrame, text = '%.2f' % (Aux[j]) ).grid(
                row = j + 1, column = 0)
        tk.Label(TabList[num][1].ThicknessFrame, text = str(points[j][0])).grid(
                row = j + 1, column = 1)

        if j == peaks - 1: # No ultimo run do ciclo for
            i = 0
            thickness = thickness / len(Aux) # A fazer a media da espessura
 
            my_file = open(TabList[num][4], 'w')
            for i in range(0, len(Aux)): # Por fim faz se a incerteza da espessura media
                uncertain = (Aux[i] - thickness)**2 + uncertain
                my_file.write('%.3f' %(Aux[i]))
                my_file.write('\n')

            uncertain = uncertain / (peaks - 1)
            uncertain = math.sqrt(uncertain) # Ultimo passo que resulta na espessura certa

            sig_fig = Precision(('%.2g' % (uncertain)))

            my_file.write(str(thickness) + '\n')
            my_file.write(str(uncertain) + '\n')
            my_file.write(str(sig_fig))
            my_file.close()

    tk.Label(TabList[num][1].ThicknessFrame, 
             text = 'Average Thickness (' + units_list[index] + ')').grid(row = j + 2, column = 0)
    tk.Label(TabList[num][1].ThicknessFrame, 
             text = 'Uncertainty (' + units_list[index] + ')').grid(row = j + 2, column = 1)
    tk.Label(TabList[num][1].ThicknessFrame, 
             text = '%.*f' %(sig_fig, thickness)).grid(row = j + 3, column = 0)
    tk.Label(TabList[num][1].ThicknessFrame, 
             text = '%.*f' %(sig_fig, uncertain)).grid(row = j + 3, column = 1)
    tk.Button(TabList[num][1].ThicknessFrame, 
              command = lambda: ClearWidget('Thickness', 1),
              text = 'Reset Results').grid(row = j + 4, columnspan = 2)

#########################################################################################
# Calculo da espessura no método ROI Select
########################################################################################
def ROI_Thick_Calculation():

    num = Current_Tab() ## Get current tab's index nr
    ClearWidget('Thickness', 0)
    TabList[num][1].ThicknessFrame.grid(row = 5, columnspan = 3, pady = 5)
    units_list = ['nm', '\u03bcm',
                  '\u03bcg' + ' cm' + '{}'.format('\u207B' + '\u00b2'),
                  '10' + '{}'.format('\u00b9' + '\u2075') + ' Atoms' 
                   + ' cm' + '{}'.format('\u207B' + '\u00b3')]
    units_values = [10.0**9, 10.0**6, 0.0, -1.0]
    index = units_values.index(TabList[num][1].units.get())
    calibration = TabList[num][1].Regression_List[0].get() ## vamos selecionar apenas uma calibração
    
    ## Get the selected material and stop. pow.
    Material_choice = TabList[num][1].Mat.get() # Determina qual o ficheiro do material a ler
    Material_choice = 'Files\Materials\\' +  Material_choice + '.txt'
    material_data = File_Reader(Material_choice, '|', 'Yes', 'No') #Daqui obtemos a lista que ira guardar

    ## Get energies of selected source
    energies = [TabList[0][1].DecayList[k].get() for k in range(len(TabList[0][1].DecayList))]
    energies.remove(-1.0) ## Remove unselected energies

    ## Get slope and intercept of the selected calib
    calib_params = File_Reader(TabList[0][4], '0', 'String', 'No')
    ## Check energy units to match the stop. pow. ones
    if calib_params[0] == 'keV':
        m = str(float(calib_params[1])*1000) ## slope
        dm = str(float(calib_params[2])*1000) ## slope uncertainty
    else:  ###   !!!!!!!!!!    VERIFICAR ESTA PARTE   !!!!!!!!!!!   ###
        m = str(float(calib_params[1])*1000) 
        dm = str(float(calib_params[2])*1000) 

    ## Get calibration centroids
    calibCents = [File_Reader(TabList[0][3], ',', 'Yes', 'No')[k][0] for k in range(len(File_Reader(TabList[0][3], ',', 'Yes', 'No')))]
    ## Get calibration peak error
    calibErr = [File_Reader(TabList[0][3], ',', 'Yes', 'No')[k][1] for k in range(len(File_Reader(TabList[0][3], ',', 'Yes', 'No')))]
    ## Get film centroids
    filmCents = [File_Reader(TabList[num][3], ',', 'Yes', 'No')[k][0] for k in range(len(File_Reader(TabList[num][3], ',', 'Yes', 'No')))]
    ## Get film peak error
    filmErr = [File_Reader(TabList[num][3], ',', 'Yes', 'No')[k][1] for k in range(len(File_Reader(TabList[num][3], ',', 'Yes', 'No')))]
    
    ## Calculate energy loss and uncertainty, returning min and max energy of alphas after crossing the film
    Emin, Emax, eloss = Eloss(energies, calibCents, filmCents, calibErr, filmErr, m, dm)

    ## Get selected material stop. pow.
    Material_choice = TabList[num][1].Mat.get() 
    Material_choice = 'Files\Materials\\' +  Material_choice + '.txt'
    material_data = File_Reader(Material_choice, '|', 'Yes', 'No')

    ## Calculate thickness from energy loss
    ## for each peak
    thickPeak = Thickness(energies, Emin, Emax, material_data)
    ## the mean of all peaks
    meanThick = np.mean([thick for thick in thickPeak])
    ## Calculate mean thickness std deviation
    stdDevThick = np.std(thickPeak)

    ## Write results to file
    result_file = open(TabList[num][4], 'w')
    for i in range(0, len(thickPeak)): 
        result_file.write(str("{:.0f}".format(thickPeak[i]))+'\n')
    result_file.write(str(meanThick) + '\n')
    result_file.write(str(stdDevThick) + '\n')
    result_file.write(str(0))
    result_file.close()

    ## Display results in the frame
    ## header
    tk.Label(TabList[num][1].ThicknessFrame, text = 'Peak Energy\n(MeV)').grid(row = 0, column = 0)
    tk.Label(TabList[num][1].ThicknessFrame, text = ' ').grid(row = 0, column = 1) #spacer
    tk.Label(TabList[num][1].ThicknessFrame, text = 'Channel\n').grid(row = 0, column = 2)
    tk.Label(TabList[num][1].ThicknessFrame, text = ' ').grid(row = 0, column = 3) #spacer
    tk.Label(TabList[num][1].ThicknessFrame, text = 'Eloss\n(keV)').grid(row = 0, column = 4)
    tk.Label(TabList[num][1].ThicknessFrame, text = ' ').grid(row = 0, column = 5) #spacer
    tk.Label(TabList[num][1].ThicknessFrame, text = 'Thickness\n(' + units_list[index] + ')').grid(row = 0, column = 6)
    ## peak info
    for k in range(len(eloss)):
        tk.Label(TabList[num][1].ThicknessFrame, text = '%.3f' % (energies[k]) ).grid(row = k + 1, column = 0)
        tk.Label(TabList[num][1].ThicknessFrame, text = ' ').grid(row = 0, column = 1) #spacer
        tk.Label(TabList[num][1].ThicknessFrame, text = '%.1f' % (filmCents[k]) ).grid(row = k + 1, column = 2)
        tk.Label(TabList[num][1].ThicknessFrame, text = ' ').grid(row = 0, column = 3) #spacer
        tk.Label(TabList[num][1].ThicknessFrame, text = '%.0f' % (eloss[k]) ).grid(row = k + 1, column = 4)
        tk.Label(TabList[num][1].ThicknessFrame, text = ' ').grid(row = 0, column = 5) #spacer
        tk.Label(TabList[num][1].ThicknessFrame, text = '%.0f' % (thickPeak[k]) ).grid(row = k + 1, column = 6)
    ## thickness final result
    tk.Label(TabList[num][1].ThicknessFrame, 
             text = 'Average Thickness (' + units_list[index] + ')').grid(row = k + 2, column = 0)
    tk.Label(TabList[num][1].ThicknessFrame, 
             text = 'Uncertainty (' + units_list[index] + ')').grid(row = k + 2, column = 4)
    tk.Label(TabList[num][1].ThicknessFrame, 
             text = "{:.0f}".format(meanThick)).grid(row = k + 3, column = 0)
    tk.Label(TabList[num][1].ThicknessFrame, 
             text = "{:.0f}".format(stdDevThick)).grid(row = k + 3, column = 4)
    tk.Button(TabList[num][1].ThicknessFrame, 
              command = lambda: ClearWidget('Thickness', 1),
              text = 'Reset Results').grid(row = k + 4, columnspan = 2)
    
    return

#########################################################################################
# Recebe os resultados dos algoritmos e mostra no GUI
########################################################################################
def ResultManager():

    num = Current_Tab()

    ClearWidget('Results', 0) # A funcao e evocada para dar reset aos valores anteriores
    TabList[num][1].ResultFrame.grid(row = 3, columnspan = 2, pady = 5)

    values = File_Reader(TabList[num][3], ',', 'No', 'No') #Esta leitura devolve os valores de channel e counts

    for j in range(0, len(values)): # O ciclo for apenas cria os checkbuttons e as labels para depois guardar
                                    # os valores a serem apagados/usados nas funcoes seguintes
        Result_Button = tk.Checkbutton(TabList[num][1].ResultFrame, variable = TabList[num][1].Var_Data[j],
                onvalue = 1, offvalue = -1,
                text = 'Channel: ' + str(values[j][0]))
        Result_Button.grid(row = j , column = 0)
        Result_Button.select()
        tk.Label(TabList[num][1].ResultFrame, 
                text = '\t Counts: ' + str(values[j][1])).grid(row = j , column = 1)

#########################################################################################
# Recebe os resultados do algoritmo ROI Select e mostra no GUI
########################################################################################
def ROIResultManager():

    num = Current_Tab()

    ClearWidget('Results', 0) # A funcao e evocada para dar reset aos valores anteriores
    TabList[num][1].ResultFrame.grid(row = 3, columnspan = 2, pady = 5)

    values = File_Reader(TabList[num][3], ',', 'Yes', 'No') #Esta leitura devolve os valores de channel e counts

    for j in range(0, len(values)): # O ciclo for apenas cria os checkbuttons e as labels para depois guardar
                                    # os valores a serem apagados/usados nas funcoes seguintes
        Result_Button = tk.Checkbutton(TabList[num][1].ResultFrame, variable = TabList[num][1].Var_Data[j],
                onvalue = 1, offvalue = -1,
                text = 'Centroid: ' + str("{:.1f}".format(values[j][0])))
        Result_Button.grid(row = j , column = 0)
        Result_Button.select()
        tk.Label(TabList[num][1].ResultFrame, 
                text = '\t \u03C3 =  ' + str("{:.1f}".format(values[j][2]))).grid(row = j , column = 1)
        tk.Label(TabList[num][1].ResultFrame, 
                text = '\t \u03C3/\u221aN =  ' + str("{:.3f}".format(values[j][1]))).grid(row = j , column = 2)
        
##########################################################################################
# Retira os resultados que nao estao checked e atualiza o txt dos resultados 
##########################################################################################        
def Unchecked_Results():

    num = Current_Tab()
    i = 0
    j = 0
    k = 0
    Eraser = []
    Aux = []

    for widget in TabList[num][1].ResultFrame.winfo_children():

        Eraser.append(widget) # Aqui guardamos todos os widgets que exibem os resultados.
                                # E necessario guardar num vetor, para que depois sejam apagados os corretos

    values = File_Reader(TabList[num][3], '0', 'String', 'No') # Aqui vao se buscar os valores de channel e counts
                                                        # Como nao se fazem contas, apenas utilizamos a linha de
                                                        # string que contem ambos valores
    for i in range(len(TabList[num][1].Var_Data)):

        if TabList[num][1].Var_Data[i].get() == 1:
            Aux.append(values[k])   # Caso seja selecionada a opcao, guardamos a linha para depois reescrever
                                    # o documento que contem os resultados, de forma a guardar os pretendidos
    
        elif TabList[num][1].Var_Data[i].get() == -1: 
            # No caso do Var_data devolver um -1, significa que se removeu a selecao do valor
            # Nesse caso, tornamos esse valor num 0, destruimos o checkbutton no Eraser[j]
            # e destruimos a label dos counts no Eraser[j + 1]
            TabList[num][1].Var_Data[i].set(0)

            Eraser[j].destroy()
            Eraser[j + 1].destroy()

        elif TabList[num][1].Var_Data[i].get() == 0:
            # Caso haja um 0 no meio, na mudanca de dados, tiramos valores de iteracao do vetor Aux
            # e do vetor Eraser, para manter todas as contas certas
            k -= 1
            j -= 2

        j += 2
        k += 1

    for i in range(len(TabList[num][1].Var_Data)):
        # Este for esta separado para nao haver confusoes de indices no primeiro ciclo
        # Aqui, se for detetado um 0, pomos o IntVar no final do vetor Var_Data e
        # eliminamos os Var_Data[i] = 0 do meio dos valores selecionados/nao selecionados
        if TabList[num][1].Var_Data[i].get() == 0:
            TabList[num][1].Var_Data.append(TabList[num][1].Var_Data[i])
            TabList[num][1].Var_Data.pop(i)

    with open(TabList[num][3], "w") as file: # Por fim, reescrevemos o documento dos resultados
                                            # com os resultados unchecked removidos
        for i in range(len(Aux)):
            file.write(Aux[i] + '\n')

#########################################################################################
# Faz a regressao linear dos resultados
#########################################################################################
def Linearize():

    num = Current_Tab()
    
    xaxis = []
    yaxis = []
    values = File_Reader(TabList[num][3], ',', 'No', 'No') # Le os dados dos picos

    for i in range(0, len(values)):
        xaxis.append(values[i][0])  # Junto os canais apenas ao valor xaxis

    for i in range(0, len(TabList[num][1].DecayList)):

        if TabList[num][1].DecayList[i].get() != -1: 
# Na lista que guarda os valores dos alfas, juntam se aqueles que foram selecionados pelo utilizador
            yaxis.append(TabList[num][1].DecayList[i].get()) 
            
    
    xvalues = sorted(xaxis) #Organizam se ambos dados por ordem
    yvalues = sorted(yaxis)
    
    if len(xvalues) != len(yvalues): # Esta condicao verifica se para a regressao linear
                                     # existe uma relacao sobrejetiva
        wng.popup('Invalid Linear Regression Configuration')
        tk.Label(wng.warning, 
                 text = "Number of Radiation Decay does not " + 
                 "match the number of Peaks detected.\n").pack()
        tk.Label(wng.warning, text ="Please adjust the Searching Algorithms or the " +
                 "number of Decay Energy.\n\n").pack()
        tk.Button(wng.warning, text = 'Return',
                    command = lambda: wng.warning.destroy()).pack()
      
    else:
        ClearWidget('Linear', 0)
        TabList[num][1].LinearRegressionFrame.grid(row = 3, columnspan = 2, pady = 5)

        avgx = sum(xvalues) #Guardam se os valores das medias dos canais e da radiacao alfa
        avgy = sum(yvalues)
        avgx = avgx / len(xvalues)
        avgy = avgy / len(yvalues)

        Placeholder1 = 0
        Placeholder2 = 0

        for i in range(len(xvalues)):
            Placeholder1 = Placeholder1 + ((xvalues[i] - avgx) * (yvalues[i] - avgy))
            Placeholder2 = Placeholder2 + (xvalues[i] - avgx)**2

        m = Placeholder1 / Placeholder2 # Valor do declive
        b = avgy - m * avgx  # Valor da ordenada na origem

        sigma = 0
        Placeholder1 = 0
        Placeholder2 = 0

        # Estes proximos somatorios e contas servem para obter as incertezas dos valores do
        # declive e da ordenada de origem
        for i in range(0, len(xvalues)):
            sigma = (yvalues[i] - m * xvalues[i] - b)**2 + sigma
            Placeholder1 = xvalues[i]**2 + Placeholder1

        Placeholder2 = (sum(xvalues))**2
        sigma = sigma / (len(xvalues) - 2) 

        sigma_m = math.sqrt(sigma / ( Placeholder1 - (Placeholder2/len(xvalues))))
        sigma_b = math.sqrt((sigma * Placeholder1)/((len(xvalues) * Placeholder1) - Placeholder2))

        if TabList[num][1].energy.get() == 1000:
            unit_string = 'MeV'
            significant_digits_m = Precision('%.2g' % (sigma_m))
            significant_digits_b = Precision('%.2g' % (sigma_b))

        elif TabList[num][1].energy.get() == 1:
            unit_string = 'keV'
            m = m * 1000
            sigma_m = sigma_m * 1000
            b = b * 1000
            sigma_b = sigma_b * 1000
            if sigma_m > 1:
                significant_digits_m = 0
            else:
                significant_digits_m = Precision('%.2g' % (sigma_m))
            if sigma_b > 1:
                significant_digits_b = 0            
            else:
                significant_digits_b = Precision('%.2g' % (sigma_b))

        with open(TabList[num][4], 'w') as my_file:
            # Aqui, escrevem se os resultados num documento txt para outras funcoes
            # poderem aceder
            my_file.write(unit_string + '\n')
            my_file.write(str(m) + '\n')
            my_file.write(str(sigma_m) + '\n')
            my_file.write(str(b) + '\n')
            my_file.write(str(sigma_b) + '\n')
            my_file.write(str(significant_digits_m) + '\n')
            my_file.write(str(significant_digits_b))

        # Por fim, escreve se no GUI os resultados obtidos
        tk.Label(TabList[num][1].LinearRegressionFrame, text = '(' + unit_string + ')').grid(
            row = 0, column = 0)
        tk.Label(TabList[num][1].LinearRegressionFrame, text = 'Values').grid(row = 0, column = 1)
        tk.Label(TabList[num][1].LinearRegressionFrame, text = 'Uncertainty').grid(row = 0, column = 2)
        tk.Label(TabList[num][1].LinearRegressionFrame, text = 'Slope').grid(row = 1, column = 0)
        tk.Label(TabList[num][1].LinearRegressionFrame, text = 'Intersect').grid(row = 2, column = 0)

        tk.Label(TabList[num][1].LinearRegressionFrame, 
                 text = '%.*f' %(6, m)).grid(row = 1, column = 1)
        tk.Label(TabList[num][1].LinearRegressionFrame, 
                 text = '%.*f' % (6, sigma_m)).grid(row = 1, column = 2)
        tk.Label(TabList[num][1].LinearRegressionFrame, 
                 text = '%.*f' %(6, b)).grid(row = 2, column = 1)
        tk.Label(TabList[num][1].LinearRegressionFrame, 
                 text = '%.*f' % (6, sigma_b)).grid(row = 2, column = 2)

        tk.Button(TabList[num][1].LinearRegressionFrame, text = 'Clear Regression', 
                  command = lambda: ClearWidget('Linear', 1 )).grid(row = 3, column = 0, columnspan = 3)

#########################################################################################
# Faz a regressao linear dos resultados com input da incerteza nosa canais
#########################################################################################
def LinearizeWithErrors():

    num = Current_Tab()
    
    centroids = []
    errors = []
    energies = []
    values = File_Reader(TabList[num][3], ',', 'Yes', 'No') # Le os dados dos picos

    for i in range(0, len(values)):
        centroids.append(values[i][0])  # Junto os canais apenas ao valor xaxis
        errors.append(values[i][1])

    for i in range(0, len(TabList[num][1].DecayList)):
        if TabList[num][1].DecayList[i].get() != -1: 
            # Na lista que guarda os valores dos alfas, juntam se aqueles que foram selecionados pelo utilizador
            energies.append(TabList[num][1].DecayList[i].get()) 
            
    centroids = sorted(centroids) #Organizam se ambos dados por ordem
    errors = sorted(errors)
    energies = sorted(energies)
    
    if len(centroids) != len(energies): # Esta condicao verifica se para a regressao linear
                                        # existe uma relacao sobrejetiva
        wng.popup('Invalid Linear Regression Configuration')
        tk.Label(wng.warning, 
                 text = "Number of Radiation Decay does not " + 
                 "match the number of Peaks detected.\n").pack()
        tk.Label(wng.warning, text ="Please adjust the Searching Algorithms or the " +
                 "number of Decay Energy.\n\n").pack()
        tk.Button(wng.warning, text = 'Return',
                    command = lambda: wng.warning.destroy()).pack()
      
    else:
        ClearWidget('Linear', 0)
        TabList[num][1].LinearRegressionFrame.grid(row = 3, columnspan = 2, pady = 5)

        ## Faz a a regressao linear com a função do ficheiro Calibration
        m, b, sigma_m, sigma_b = Calib(energies, centroids, errors)

        if TabList[num][1].energy.get() == 1000:
            unit_string = 'MeV'

        elif TabList[num][1].energy.get() == 1:
            unit_string = 'keV'
            m = m * 1000
            sigma_m = sigma_m * 1000
            b = b * 1000
            sigma_b = sigma_b * 1000

        with open(TabList[num][4], 'w') as my_file:
            # Aqui, escrevem se os resultados num documento txt para outras funcoes
            # poderem aceder
            my_file.write(unit_string + '\n')
            my_file.write(str(m) + '\n')
            my_file.write(str(sigma_m) + '\n')
            my_file.write(str(b) + '\n')
            my_file.write(str(sigma_b) + '\n')
            my_file.write(str(6) + '\n')
            my_file.write(str(6))

        # Por fim, escreve se no GUI os resultados obtidos
        tk.Label(TabList[num][1].LinearRegressionFrame, text = '(' + unit_string + ')').grid(
            row = 0, column = 0)
        tk.Label(TabList[num][1].LinearRegressionFrame, text = 'Values').grid(row = 0, column = 1)
        tk.Label(TabList[num][1].LinearRegressionFrame, text = 'Uncertainty').grid(row = 0, column = 2)
        tk.Label(TabList[num][1].LinearRegressionFrame, text = 'Slope').grid(row = 1, column = 0)
        tk.Label(TabList[num][1].LinearRegressionFrame, text = 'Intersect').grid(row = 2, column = 0)

        tk.Label(TabList[num][1].LinearRegressionFrame, 
                 text = '%.*f' %(6, m)).grid(row = 1, column = 1)
        tk.Label(TabList[num][1].LinearRegressionFrame, 
                 text = '%.*f' % (6, sigma_m)).grid(row = 1, column = 2)
        tk.Label(TabList[num][1].LinearRegressionFrame, 
                 text = '%.*f' %(6, b)).grid(row = 2, column = 1)
        tk.Label(TabList[num][1].LinearRegressionFrame, 
                 text = '%.*f' % (6, sigma_b)).grid(row = 2, column = 2)

        tk.Button(TabList[num][1].LinearRegressionFrame, text = 'Clear Regression', 
                  command = lambda: ClearWidget('Linear', 1 )).grid(row = 3, column = 0, columnspan = 3)

###############################################################################
# Este e o algoritmo que determina a distancia quadrada minima entre pontos
# input, e pontos de dados
###############################################################################
def ManSelec_Alg(Valuex, Valuey):

    num = Current_Tab()
    Counts = File_Reader(TabList[num][2], '0', 'No', 'No')

    AuxList = []   #Lista para guardar minimos quadrados

    for i in range(0, len(Counts)):
        AuxList.append(math.sqrt(((Valuey - Counts[i])**2)+((Valuex - i + 1)**2)))

    # O ciclo for guarda toda a distancia quadrado dos valores input e dos valores dos dados

    Valuex = AuxList.index(min(AuxList)) + 1 #Corresponde ao channel com distancia minima
    Valuey = Counts[Valuex - 1] #Corresponde ao count do indice do xvalue

    if os.path.isfile(TabList[num][3]) == True:
        with open(TabList[num][3], 'a') as results:
            results.write(str(Valuex) + ',' + str(Valuey) + '\n')

    elif os.path.isfile(TabList[num][3]) == False:

        with open(TabList[num][3], 'w') as results:
            results.write(str(Valuex) + ',' + str(Valuey) + '\n')

    ResultManager()
    
###############################################################################
# Este e o algoritmo que regista os picos acima de um determinado threshold
###############################################################################
def Threshold_Alg():

    num = Current_Tab()
    Counts = File_Reader(TabList[num][2], '0', 'No', 'No')

    Threshold = TabList[num][1].Algorithm.get()
    yaxis = []
    current_peak_counts = 0
    current_peak = []
    counter = 0
    xaxis = []
    i = 0
    j = 0

    cut = TabList[num][1].channels_cut.get()
    width = TabList[num][1].peaks_widths.get()

#Just iterate of the values of the list itself. Each "count" is an element of Counts
    for count in (Counts):
       #Remember, if you hit a peak the condition 
       # will enter this part always. So the "current_peak_counts" will just 
       #add the total sum of the area under the peak 
       # (this may even be useful for FWHM calculations).
       
        if count > Threshold and i > cut:
                current_peak_counts += count
                current_peak.append(count)
        else:
            counter += 1
            #Just checking that our current_peak list is not empty, 
            #if it is then we're out of a peak and this will just skip this if statment
            #but if it is not, then we're at the end of a peak 
            #(count is now <= Threshold) and we should save the max value of
            #current_peak and put everything to zero (aka move to the next peak)
            if current_peak:
                    yaxis.append(max(current_peak))
                    xaxis.append(counter + current_peak.index(max(current_peak)))
                    if len(xaxis) > 1:
                        if xaxis[j] - xaxis[j-1] < width:  
                            decider = max(yaxis[j-1], yaxis[j])
                            if decider == yaxis[j]:
                                yaxis.pop(j-1)
                                xaxis.pop(j-1)
                            elif decider == yaxis[j-1]:
                                yaxis.pop(j)
                                xaxis.pop(j)
                            j -= 1
                    counter = counter + len(current_peak)
                    current_peak_counts = 0
                    current_peak = []
                    j += 1
        i += 1

    if os.path.isfile(TabList[num][3]) == True:
        with open(TabList[num][3], 'a') as results:
            for i in range(len(xaxis)):
                results.write(str(xaxis[i]) + ',' + str(yaxis[i]) + '\n')

    elif os.path.isfile(TabList[num][3]) == False:
        with open(TabList[num][3], 'w') as results:
            for i in range(len(xaxis)):
                results.write(str(xaxis[i]) + ',' + str(yaxis[i]) + '\n')

    #If the last peak is not added because the data ends with a peak (could be a problem, you should check it)
    ResultManager()

###############################################################################
# Funcao para o algoritmo ROI Select
################################################################################
def ROI_Select_Alg():
    
    num = Current_Tab()
    counts = File_Reader(TabList[num][2], '0', 'Yes', 'No')

    roi_down = [TabList[num][1].ROIdown1.get(),
                TabList[num][1].ROIdown2.get(),
                TabList[num][1].ROIdown3.get(),
                TabList[num][1].ROIdown4.get(),
                TabList[num][1].ROIdown5.get(),
                TabList[num][1].ROIdown6.get()]
    
    roi_up = [  TabList[num][1].ROIup1.get(),
                TabList[num][1].ROIup2.get(),
                TabList[num][1].ROIup3.get(),
                TabList[num][1].ROIup4.get(),
                TabList[num][1].ROIup5.get(),
                TabList[num][1].ROIup6.get()]

    cents, errs, sigmas  = Analyze(counts, roi_down, roi_up)

    if os.path.isfile(TabList[num][3]) == True:
        with open(TabList[num][3], 'a') as results:
            for i in range(len(cents)):
                results.write(str(cents[i]) + ',' + str(errs[i]) + ',' + str(sigmas[i]) + '\n')

    elif os.path.isfile(TabList[num][3]) == False:
        with open(TabList[num][3], 'w') as results:
            for i in range(len(cents)):
                results.write(str(cents[i]) + ',' + str(errs[i]) + ',' + str(sigmas[i]) + '\n')

    ROIResultManager()
    return
  
###############################################################################
# Esta e a funcao que abre as imagens da cadeia de decaimento
################################################################################
def showimage():
    num = Current_Tab()
    alphas = TabList[num][1].Source.get()
    Dir = os.scandir('Files\Sources\Images')

    for entry in Dir:
        if entry.is_file():
            temp = (os.path.splitext(entry.name))
            if alphas == temp[0]:
                picture = 'Files\Sources\Images\\' + alphas + temp[1]

                if alphas == '226Ra':
                    domain = 'https://www.nist.gov/image-23773'

                elif alphas == '232U':
                    domain = 'https://www.nuclear-power.com/nuclear-power-plant/nuclear-fuel/uranium/uranium-232/decay-half-life-uranium-232/'

                else:
                    domain = ''

                try:
                    ClearWidget('Image', 0)
                    wng.Images('Decay Chain of ' + alphas, picture, domain)

                except:
                    wng.Images('Decay Chain of ' + alphas, picture, domain)

################################################################################
# Gere o evento de obtencao de pontos diretamente do grafico
################################################################################
def onclick(event):

    num = Current_Tab()

    decider = TabList[num][1].Algorithm_Method.get()

    if decider == 'Manual Selection' and event.button == 3:
    # Estas duas opcoes verificam que so no modo Manual selection e so com o botao direito do rato
    # e que aparecem os dados
        xpoint = event.xdata
        ypoint = event.ydata
        ManSelec_Alg(xpoint, ypoint)

##############################################################################
#Esta funcao ira ler o ficheiro input.
##############################################################################
def DataUploader(): 

    num = Current_Tab()

    domain = (('Text Files', '*.mca'), ('All Files', '*.*')) # Aqui limitamos os ficheiros que podem ser abertos
    filename = fd.askopenfilename(title = 'Open a file', initialdir = ',', filetypes = domain) 
    # Faz a ligacao tkinter - janela do SO para abrir o ficheiro de dados

    if not filename:
        pass # Certifica que o programa nao queixa caso nao seja efetuado um upload

    else:
        file = File_Reader(filename, '0', 'string', 'Yes') # Aqui le se o ficheiro que foi feito o upload
        TabList[num][5].Structure(file, filename)    # Logo de seguida faz se o grafico
        TabList[num][5].subplots()
        value = TabList[num][1].Algorithm.get()
        if value != 0:
            TabList[num][5].threshold(value)

    ## Renames the tab with the filename
    try:
        Tabs.RenameTab(filename.split('/')[-1])
    except:
        ()

    return

##############################################################################
# Esta funcao gere o evento especifico de adicionar tabs, futuramente
# ira ter uma seccao especifica para esconder tabs e recuperar a tab
# de adicionar tabs, caso haja menos que 15 tabs
#############################################################################
def handleTabChange(event):

    if Notebook.notebook.select() == Notebook.notebook.tabs()[-1] and len(Notebook.notebook.tabs()) < 15:
        # Se a tab 'atual' for igual a ultima, que sera sempre a tab '+', e se o numero de tabs for menor que 15
        # Adiciona se outra tab
            
        wng.popup('Tab Selector Menu')

        tk.Label(wng.warning, text = 'Please Select New Type of Tab to Open \n').pack()
        tk.Button(wng.warning, command = lambda: Tabs.tab_change(1),
                        text = 'Calibration Trial').pack()
        tk.Button(wng.warning, command = lambda: Tabs.tab_change(2), 
                        text = 'Material Trial').pack()
        tk.Label(wng.warning, text ='\n').pack()
        tk.Button(wng.warning, text = 'Return',
                    command = lambda : Tabs.tab_change(3)).pack()
        # os numeros do tab_change indicam quais os tipos de tabs a adicionar

    elif len(Notebook.notebook.tabs()) >= 15:
        Notebook.notebook.hide(14)
        # Caso se chege ao numero 15 de tabs, o botao '+' e escondido

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