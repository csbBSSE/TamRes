import matplotlib.pyplot as plt
import numpy as np
import os
plt.rcParams.update({'font.size': 5})
from matplotlib import cm
import matplotlib as mpl
import seaborn as sb


path_to_files = "Enter/Path/to/output/files/"
path_to_savefig = "Enter/Path/to/savefigures/"

def calculate_mean_std(rep):
    rep = np.array(rep)
    mean_arr = []
    std_arr = []
    for i in range(len(rep[0])):
        mean_arr.append(np.mean(rep.transpose()[i]))
        std_arr.append(np.std(rep.transpose()[i]))
    
    return mean_arr, std_arr


def plot_parameter(initial_pop,heterogeneity,x_array,ylim,ylim_on_off,savefig_filename,parameter):
    rev_grid = [ele for ele in reversed(grid)]
    for i in initial_pop:
        os.chdir(path_to_files) #Put the path to the output files from the simulation
        figure, axes = plt.subplots(nrows=len(grid), ncols=len(grid),sharex='col', sharey='row')
        for j in heterogeneity:
            for p,k in enumerate(rev_grid):
                for q,l in enumerate(grid): 
                    filename = "out_tamres_"+str(i)
                    if j == 0:
                        filename += "_"+str(j)
                    else:
                        filename += "_"+("%.2f" % j).lstrip('0')
                    if k == 0:
                        filename += "_"+str(k)
                    else:
                        filename +="_"+("%.1f" % k).lstrip('0')
                    if l == 0:
                        filename += "_"+str(l)
                    else:
                        filename +="_"+("%.1f" % l).lstrip('0')

                    filename +="_summary.txt"
                    print(filename)
                    with open(filename) as file:
                        rep = []
                        for line in file:
                            if "Replicate" in line:
                                continue
                            if line == "\n":
                                continue
                            if "###" in line:
                                continue
                            if parameter in line:
                                temp = line.split("\t")[1:-1]
                                temp = list(map(int, temp)) 
                                rep.append(temp)

                    mean_arr, std_arr = calculate_mean_std(rep)
                    plot_y_array = [mean_arr[int(i)] for i in x_array]
                    plot_y_err_array = [std_arr[int(i)] for i in x_array]
                    axes[p,q].errorbar(x_array,plot_y_array,yerr = plot_y_err_array,linewidth=0.5,fmt='-o',markersize=0.75,capsize=1,capthick=0.5)
                    axes[p,q].set_xlim([0,100])
                    if ylim_on_off == "on":
                        axes[p,q].set_ylim([0,ylim])
            figure.set_facecolor("w")
    plt.savefig(path_to_savefig + "FigureS6.svg", dpi = 1200)
    plt.close()

initial_pop = [1]
heterogeneity= [0,0.25,0.5,0.75,1]
grid = [0,0.2,0.4,0.6,0.8,1]
time = [0,10,20,30,40,50,60,70,80,90,99]
x_array =time
ylim_1 = 2300
ylim_on_off = "off"
savefig_filename = "Population size Matrix"
parameter = "Population size"
plot_parameter(initial_pop,heterogeneity,x_array,ylim_1,ylim_on_off,savefig_filename,parameter)




