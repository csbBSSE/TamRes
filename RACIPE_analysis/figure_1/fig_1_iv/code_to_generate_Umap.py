import os
import umap
import itertools
import collections
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""This file contains the functions that generates the Umap"""

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
            sub_dataframe.columns = genes
            state_dataframe = state_dataframe.append(sub_dataframe,ignore_index = True)

    return state_dataframe

#This function performs the Umap analysis. 
def Umap_analysis(data,n_neighbors=15, min_dist=0.1, n_components=2, metric='euclidean'):
    fit = umap.UMAP(n_neighbors=n_neighbors,min_dist=min_dist,n_components=n_components,metric=metric)
    umap_data = fit.fit_transform(data)

    return umap_data

#This function plots the output from the Umap analysis. 
def Umap_scatter(state_dataframe,path_to_plots):
    sub_dataframe = state_dataframe.sample(n=10000)
    embedding = Umap_analysis(sub_dataframe,n_neighbors=100)
    Umap = pd.DataFrame()

    Umap['UMAP_1'] = embedding[:,0]
    Umap['UMAP_2'] = embedding[:,1]
    Umap.plot.scatter(x = 'UMAP_1',y='UMAP_2',s=1,c=sub_dataframe['ZEB1']-sub_dataframe['MIR200'],cmap=plt.cm.nipy_spectral)
    plt.title('EMT score UMAP',fontsize=14)
    plt.xlabel('UMAP_1',fontsize=14)
    plt.ylabel('UMAP_2',fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(path_to_plots+"ZEB1_MIR200_expression_diff_Umap_scatter.png",dpi=400)
    plt.close()

    Umap.plot.scatter(x='UMAP_1',y='UMAP_2',s=1,c=sub_dataframe['ERa36']-sub_dataframe['ERa66'],cmap=plt.cm.nipy_spectral)
    plt.title('Tamoxifen resistance score UMAP',fontsize=14)
    plt.xlabel('UMAP_1',fontsize=14)
    plt.ylabel('UMAP_2',fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(path_to_plots+"ERa36_ERa66_expression_diff_Umap_scatter.png",dpi=400)
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
    path_to_plots = path_to_dat_files+"plots/Umap_plots/"
    if not os.path.exists(path_to_plots):
        os.makedirs(path_to_plots)
    state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
    Umap_scatter(state_dataframe,path_to_plots)

    
for replicate in replicte_array:
	run_for_all_analysis(replicate)