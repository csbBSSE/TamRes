import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from textwrap import wrap
import matplotlib.pyplot as plt


"""This file contains the functions that generates a clustermap for the simulated data"""


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

#This function generates a heatmap for the z normalized dataframe given to it as input.  

def plotting_overall_heatmaps(state_dataframe,network_name,path_to_plots):
	
	#sampling 100000 solutions from the total solutions generated from simulations.
	state_dataframe_hm = state_dataframe.sample(n=100000)

	#generating heatmap
	array_max_min = np.array([state_dataframe_hm.max().max(),-1*state_dataframe.min().min()]).max()
	sns.set(font_scale=1.4)
	ax = sns.clustermap(state_dataframe_hm,cmap='seismic',vmin=-1*array_max_min,vmax=array_max_min,annot_kws={"size": 16})

	#saving the plot
	plt.savefig(path_to_plots+network_name+"heatmap_all_considered.png",dpi=500)
	plt.close()

num_stability_to_consider = 4
replicates = ['r1','r2','r3']

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
	plotting_overall_heatmaps(state_dataframe,network_name,path_to_plots)
	
for replicate in replicates:
	run_for_all_analysis(replicate)