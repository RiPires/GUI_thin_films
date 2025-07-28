#RiP
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import csv

def PlotData(File):
    """
    Plots yield vs channel data from our .mca file

    INPUTS: "FileName.mca"
    OUTPUTS: yield vs channel plot
    """

    with open(File, 'r', encoding='latin1') as file:
        reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
        data = list(reader)
        #print(data)
    x = []
    y = []
    aux = []
    time=float(data[7][0][data[7][0].find('- ')+len('- '):])
    end=[i for i in range(len(data)) if data[i]== ['<<END>>']][-1]
    #print(end)
    for i in range(12,end): # len(data)-1
        aux.append(data[i][0].split())
    for i in range(len(aux)):
        #x.append(float(0.004482*i-0.105144)) ## energy calibraion needs to be corrected
        x.append(float(i)) ## axes in channel
        y.append(float(aux[i][0])/time)

    fig, ax = plt.subplots()
    ax.plot(x,y,'*-', color ='xkcd:black', label=(str(File)))
    #ax.semilogy(x,y,'^', color ='xkcd:purple', label=(str(File)))
    legend = ax.legend(loc="upper right",ncol=2, shadow=False,fancybox=True,framealpha = 0.0,fontsize=20)
    legend.get_frame().set_facecolor('#DAEBF2')
    tick_params(axis='both', which='major', labelsize=22)
    ##xlabel('Energy (MeV)',fontsize=22)
    xlabel('Channel',fontsize=22)
    ylabel('Count rate ($\\rm{s}^{\\rm{-1}}$)', fontsize=22)
    #yscale('log')
    #grid()
    show()

    return
##################################################


PlotData('Data/Combined_bkg_Internship.txt')
#PlotData('Data/NoFilm_for_C12_4h_2.mca')







