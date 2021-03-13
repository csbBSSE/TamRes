
import matplotlib.pyplot as plt
import numpy as np
import os
plt.rcParams.update({'font.size': 5})
from matplotlib import cm
import matplotlib as mpl
import seaborn as sb
import random

path_to_files = "Enter/Path/to/output/files/"
path_to_savefig = "Enter/Path/to/savefigures/"
filename = "EnterFileName.txt"

# %%
def calculate_mean_std(rep):
    rep = np.array(rep)
    mean_arr = []
    std_arr = []
    for i in range(len(rep[0])):
        mean_arr.append(np.mean(rep.transpose()[i]))
        std_arr.append(np.std(rep.transpose()[i]))
        #print(rep.transpose()[i])
    
    return mean_arr, std_arr


# %%
def sample_gaussian(mean,std):
    return np.random.normal(mean,std,100)


# %%
def plot_3d_hist(initial_pop,heterogeneity,x_array,y_array,ylim_1,ylim_on_off,savefig_filename,filename):
    os.chdir(path_to_files)
    initial_arr = np.histogram(sample_gaussian(-2,het),np.linspace(-6,6,121))[0]
    print(len(initial_arr))
    with open(filename) as file:
        rep = []
        for line in file:
            if "Replicate" in line:
                div_array = []
                continue
            elif line == "\n":
                continue
            elif "###" in line:
                rep.append(div_array)
                continue
            else:
                temp = line.split("\t")[:-1]
                temp = list(map(float, temp))
                sum_temp = np.sum(temp)
                if sum_temp == 0:
                    temp = np.array(temp)
                else:
                    temp = np.array(temp)#/sum_temp
                if len(temp) == 121:
                    temp = temp[:-1]
                div_array.append(temp)
            
        
        for p,i in enumerate(rep):
            for q,j in enumerate(i):
                if len(j) == 24:
                    rep[p][q] = [0 for x in range(120)]
        #for i in rep:
        #    for j in i:
        #        print(len(j))
        rep = np.array(rep)
        mean_arr = np.mean(rep,axis=0)
        std_arr = np.std(rep,axis=0)
        bottom = np.zeros(len(x_array))
        #print(len(mean_arr))
        matrix = []
        #matrix.append(initial_arr)
        for w in range(len(mean_arr[0])):
            plot_y_array = [mean_arr[int(k)][w] for k in x_array]
            matrix.append(plot_y_array)
        
    xp = np.array([x_array for i in range(len(y_array))]).flatten()
    yp = np.array([y_array for i in range(len(x_array))]).flatten("F")
    zp = np.array([0 for i in range(len(x_array)*len(y_array))])
    zq = np.array(matrix).flatten()
    xq = [0.1 for i in range(len(x_array)*len(y_array))]
    yq = [0.1 for i in range(len(x_array)*len(y_array))]
    cmap = cm.get_cmap('jet')
    max_height = np.max(zq)   # get range of colorbars so we can normalize
    min_height = np.min(zq)
    fig = plt.figure()
    fig.set_facecolor("w")
    ax = plt.axes(projection="3d")
    ax.bar3d(xp,yp,zp,xq,yq,zq, color=rgba)
    ax.view_init(elev=16., azim=44)
    ax.grid(False)
    plt.savefig(path_to_savefig+"Figure5F_3D_histogram.svg")
            


# %%
filename = path_to_files + filename
initial_pop = [1]
heterogeneity= [0,0.25,0.5,0.75,1]
het = 0.25
grid = [0,0.2,0.4,0.6,0.8,1]
#time = [0,20,40,60,80,99]
time = np.linspace(0,100,101)[:-1]
x_array =time
y_array = np.linspace(-6,6,120)
ylim_1 = 30000
ylim_2 = 40000
ylim_3 = 80000
ylim_on_off = "off"
savefig_filename = "Final Diversity Matrix"
plot_num = 1
plot_3d_hist(initial_pop,heterogeneity,x_array,y_array,ylim_1,ylim_on_off,savefig_filename,filename)


