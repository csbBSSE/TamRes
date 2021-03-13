import matplotlib.pyplot as plt
import numpy as np
import os
plt.rcParams.update({'font.size': 5})
from matplotlib import cm
import matplotlib as mpl
import seaborn as sb
from matplotlib.colors import LogNorm

path_to_files = "Enter/Path/to/output/files"
path_to_savefig = "Enter/Path/to/savefigures"


def calculate_mean_std(rep):
    rep = np.array(rep)
    mean_arr = []
    std_arr = []
    for i in range(len(rep[0])):
        mean_arr.append(np.mean(rep.transpose()[i]))
        std_arr.append(np.std(rep.transpose()[i]))
    
    return mean_arr, std_arr



def plot_parameter(initial_pop,heterogeneity,x_array,ylim,ylim_on_off,savefig_filename,parameter):
    matrix = np.zeros((len(grid),len(grid)))
    rev_grid = [ele for ele in reversed(grid)]
    for i in initial_pop:
        os.chdir(path_to_files) #replace this with the path to the output files from the simulation
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
                    matrix[p,q] += mean_arr[-1]
    matrix = matrix/len(heterogeneity)
    matrix = np.matrix(matrix)
    matrix = matrix 
    svm = sb.heatmap(matrix,xticklabels=grid,yticklabels=rev_grid, cmap="coolwarm")
    figure = svm.get_figure()
    figure.set_facecolor("w")
    plt.xlabel("Backward Switching Probability")
    plt.ylabel("Forward Switching Probability")
    plt.title(parameter)
    figure.savefig(path_to_savefig+"Figure5D.svg", dpi=800)



initial_pop = [1]
heterogeneity= [0.0,0.25,0.5,0.75,1.0]
grid = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
time = [0,10,20,30,40,50,60,70,80,90,99]
#time = np.linspace(0,100,11)[:-1]
x_array =time
ylim_1 = 2300
ylim_on_off = "off"
savefig_filename = "FigureS6 Het = "+str(heterogeneity[0])+"  Pop size Heatmap"
parameter = "Population size"
plot_parameter(initial_pop,heterogeneity,x_array,ylim_1,ylim_on_off,savefig_filename,parameter)



