import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from textwrap import wrap
import matplotlib.pyplot as plt


"""This file contains the functions that generates the heatmap for the Resistance and EMT scores calculated from simulations data"""

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


#This function generates the score dataframes from the input z normalized data
def collated_score_for_all_dataframes(state_dataframe):

	state_dataframe['EMT_score'] = state_dataframe["ZEB1"]- state_dataframe["MIR200"]
	state_dataframe['TamRes_score'] = state_dataframe["ERa36"]- state_dataframe["ERa66"]

	state_score_dataframe = pd.DataFrame(columns = ['EMT_score','TamRes_score'])
	state_score_dataframe['EMT_score'] = state_dataframe['EMT_score']
	state_score_dataframe['TamRes_score'] = state_dataframe['TamRes_score']

	return state_dataframe, state_score_dataframe

#This function generates the heatmap 
def plotting_heatmaps_for_score(state_dataframe,path_to_heatmap_plots):
	# sampling 100000 solutions from the total solutions generated from simulations.
	
	state_dataframe_hm = state_dataframe.sample(n=100000)
	array_max_min = np.array([state_dataframe_hm.max().max(),-1*state_dataframe_hm.min().min()]).max()
	ax = sns.clustermap(state_dataframe_hm,cmap='seismic',vmin=-1*array_max_min,vmax=array_max_min)
	plt.savefig(path_to_heatmap_plots+"heatmap_score.png",dpi=500)
	plt.close()


def run_for_all_analysis(replicate):
	core_path = "./../../"
	
	## running for self_activation perturbation
	print("Running for self-activation perturbations: ")
	network_name = "core"
	path_to_dat_files = core_path+"topo_file/"+replicate+"/"
	name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
	genes = list(d_genes.values())
	path_to_output_z_norm = path_to_dat_files+"Z_normed/"
	path_to_plots = path_to_dat_files+"plots/heatmaps/"
	if not os.path.exists(path_to_plots):
		os.makedirs(path_to_plots)
	state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
	plotting_heatmaps_for_score(state_score_dataframe,path_to_heatmap_plots)
	
for replicate in replicates:
	run_for_all_analysis(replicate)