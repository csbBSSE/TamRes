import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import math as math
from scipy import stats
from textwrap import wrap
import matplotlib.pyplot as plt

"""The functions in this file plot EM and resistance score on a scatter plot"""


#This function reads the CFG file to get information about the genes in the network.
def reading_the_cfg_file(path_to_folder): 
	for filename in os.listdir(path_to_folder):
		#print(filename[-4:])
		if filename[-4:] == ".cfg":
			name = filename[:-4]
	flag = -1
	d_genes = {}
	with open(path_to_folder+name+".cfg") as f:
		for line in f:
			a = line[:-1].split("\t")
			if a[0] == "NumberOfRACIPEModels":
				num_models = float(a[1])
			if a[0] == "NumberOfGenes":
				nodes = int(a[1])
				flag = 0
				continue
			if flag >= 0 and flag < nodes:
				flag += 1
				d_genes[a[0]] = a[1]

	return name,num_models,nodes,d_genes

#This function reads the Z normalized solution files generated from simulations upto a level of stability. 

def collating_all_runs_for_z_score_calculation(path_to_dat_files,network_name,genes,num_stability_to_consider):

    state_dataframe = pd.DataFrame(columns = np.sort(genes))

    # iterating over all the files with mono, bi, tri and tetra stable cases
    for i in range(1,num_stability_to_consider+1):
        
        # reading each file separately, getting the sub data structures and appending it to the state dataframe
        data = pd.read_csv(path_to_dat_files+network_name+"_solution_"+str(i)+"_z_norm.dat",delimiter="\t",header=None)
        
        for j in range(0,i):
            sub_dataframe = data[data.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = np.sort(genes)
            state_dataframe = state_dataframe.append(sub_dataframe,ignore_index = True)

    return state_dataframe

# This function reads the file containing the mean and standard deviation for each gaussian in the score distribution.
def read_scoring_file(path_to_score_file,phenotype):

	data = pd.read_csv(path_to_score_file+'_'+phenotype+'score_stats.txt',delimiter="\t",header=None)
	mean_of_the_columns = data.mean(axis = 0)
	sdev_of_the_columns = data.std(axis = 0)

	return mean_of_the_columns, sdev_of_the_columns

# This function uses mean and standard deviation for each gauassian in the score distribution determine a boundary condition
def return_inidividual_phenotypes_boundary_condition(mean_EMT,mean_TamRes,std_EMT,std_TamRes):
	boundary = []
	boundary.append(mean_EMT[6] -(abs(mean_EMT[6]-mean_EMT[3])/2))
	boundary.append(mean_EMT[6] +(abs(mean_EMT[6]-mean_EMT[0])/2))
	boundary.append(mean_TamRes[3] +(abs(mean_TamRes[0]-mean_TamRes[3])/2))
	return boundary

# This function generates score files from input z normalized dataframes
def collated_score_for_all_dataframes(state_dataframe):

	state_dataframe['EMT_score'] = state_dataframe["ZEB1"]- state_dataframe["MIR200"]
	state_dataframe['TamRes_score'] = state_dataframe["ERa36"]- state_dataframe["ERa66"]

	state_score_dataframe = pd.DataFrame(columns = ['EMT_score','TamRes_score'])
	state_score_dataframe['EMT_score'] = state_dataframe['EMT_score']
	state_score_dataframe['TamRes_score'] = state_dataframe['TamRes_score']

	return state_dataframe, state_score_dataframe

# This function generates scatter plots with boundary condition for identification of phenotype marked as red line on the plot
def plotting_scatter_plots_for_score(state_dataframe,path_to_plots,boundaries):

	sc = plt.scatter(state_dataframe['EMT_score'],state_dataframe['TamRes_score'],marker="o",s=0.1,c='black')
	plt.axvline(x=boundaries[0], color='r')
	plt.axvline(x=boundaries[1], color='r')
	plt.axhline(y=boundaries[2], color='r')
	plt.xlabel('EMT_score')
	plt.ylabel('TamRes_score')
	plt.ylim(-5, 4)
	plt.xlim(-4.1, 4.1)
	plt.savefig(path_to_plots+"score_scatter.png",dpi=500)
	plt.close()

num_stability_to_consider = 3
replicates = ['r1','r2','r3']

def run_for_all_analysis(replicate):
	core_path = "./../../"
	
	## running for self_activation perturbation
	network_name = "core"
	path_to_dat_files = core_path+"topo_file/"+replicate+"/"
	name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
	genes = list(d_genes.values())
	path_to_output_z_norm = path_to_dat_files+"Z_normed/"
	path_to_output_phenotype = path_to_dat_files+"phenotype/"
	path_to_plots = path_to_dat_files+"plots/scatters/"
	path_to_score_file = path_to_dat_files+"phenotype/"
	if not os.path.exists(path_to_plots):
		os.makedirs(path_to_plots)
	state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
	state_dataframe, state_score_dataframe = collated_score_for_all_dataframes(state_dataframe)
	mean_EMT, std_EMT = read_scoring_file(path_to_score_file,'EMT')
	mean_TamRes, std_TamRes = read_scoring_file(path_to_score_file,'TamRes')
	boundaries = return_inidividual_phenotypes_boundary_condition(mean_EMT,mean_TamRes,std_EMT,std_TamRes)
	plotting_scatter_plots_for_score(state_score_dataframe,path_to_plots,boundaries)
 