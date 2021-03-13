import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import math as math
from scipy import stats
from textwrap import wrap
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from sklearn import metrics
from scipy.spatial.distance import cdist

"""The functions in this file generates a bar plot for expression of nodes in the network for each phenotype"""

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

# This function uses mean and standard deviation for each gauassian in the score distribution determine a boundary condition.
def return_inidividual_phenotypes_boundary_condition(mean_EMT,mean_TamRes,std_EMT,std_TamRes):
	boundary = []
	boundary.append(mean_EMT[6] -(abs(mean_EMT[6]-mean_EMT[3])/2))
	boundary.append(mean_EMT[6] +(abs(mean_EMT[6]-mean_EMT[0])/2))
	boundary.append(mean_TamRes[3] +(abs(mean_TamRes[0]-mean_TamRes[3])/2))
	return boundary

# This function generates score files from input z normalized dataframes.
def collated_score_for_all_dataframes(state_dataframe):

	state_dataframe['EMT_score'] = state_dataframe["ZEB1"]- state_dataframe["MIR200"]
	state_dataframe['TamRes_score'] = state_dataframe["ERa36"]- state_dataframe["ERa66"]

	state_score_dataframe = pd.DataFrame(columns = ['EMT_score','TamRes_score'])
	state_score_dataframe['EMT_score'] = state_dataframe['EMT_score']
	state_score_dataframe['TamRes_score'] = state_dataframe['TamRes_score']

	return state_dataframe, state_score_dataframe

# This function determines the phenotype based upon given EM and resistance score. 
def whats_my_phenotype(EMT_score,TamRes_score,boundaries):
	if EMT_score < boundaries[0]:
		s = "E"
	elif EMT_score >= boundaries[0] and EMT_score < boundaries[1]:
		s = "H"
	else:
		s = "M"
	if TamRes_score < boundaries[2]:
		s += "S"
	else:
		s += "R"
	return s

# This function estimates expression of each node for a given phenotype cluster.
def cluster_chacteristics(cluster_number,df,genes):
	print("\nCluster ",cluster_number)
	print("num\t",genes)
	# extracting a given phenotype cluster
	sub_cluster = df.loc[df['Phenotype'].isin([cluster_number])] 

	# estimating mean, standard deviation and standard error of mean for each node's expression
	mean_array = sub_cluster.mean()
	std_array = sub_cluster.std()
	sem_array = mean_array/(math.sqrt(len(sub_cluster)))

	s_mean = str(len(sub_cluster))
	s_std = str(len(sub_cluster))
    for i in genes:
        s_mean += "\t"+str(mean_array[i])
        s_std += "\t"+str(std_array[i])
	print(s_mean)
	print(s_std)

	mean_exp = []
	std_exp = []
	sem_exp = []
	for i in genes:
		mean_exp.append(mean_array[i])
		std_exp.append(std_array[i])
		sem_exp.append(sem_array[i])

	return mean_exp,std_exp,sem_exp


"""This function quantifies the number of ES, ER, HS, HR, MS, MR phenotypes. 
It also calls the cluster_characteristic function and plot the expression of each node in a phenotype cluster."""
def quantify_all_phenotypes(state_dataframe,boundaries,genes,path_to_plots):

	# counting the number of data points for each phenotype
	phe = []
	phenotype_list = ["ES","ER","HS","HR","MS","MR"]
	mean_exp_df = pd.DataFrame(index = genes)
	std_exp_df = pd.DataFrame(index = genes)
	sem_exp_df = pd.DataFrame(index = genes)
	phe_count = pd.DataFrame()
	for index, row in state_dataframe.iterrows():
		phe.append(whats_my_phenotype(row['EMT_score'],row['TamRes_score'],boundaries))
	state_dataframe["Phenotype"] = phe
	phe_count = state_dataframe["Phenotype"].value_counts().sort_index()
	total = sum(phe_count)
	phe_count_percent = (phe_count/total)*100


	# estimating the expression of each node 
	for cluster_number in phenotype_list:
		mean,std,sem = cluster_chacteristics(cluster_number,state_dataframe,genes)
		mean_exp_df[cluster_number] = mean
		std_exp_df[cluster_number] = std
		sem_exp_df[cluster_number] = sem
	# print(std_exp_df)
	mean_exp_df_t = mean_exp_df.T
	std_exp_df_t = std_exp_df.T
	sem_exp_df_t = sem_exp_df.T


	#Plotting the expression value of each node
	n_groups = 6
	fig, ax = plt.subplots()
	index = np.arange(n_groups)
	bar_width = 0.15
	for gene in genes:
		count = genes.index(gene)
		rects_count = plt.bar(index+(bar_width*count), mean_exp_df_t[gene], bar_width,label=gene,yerr=sem_exp_df_t[gene],error_kw=dict(lw=2,capsize=2,capthick=1))
	plt.xticks(index+bar_width*2, phenotype_list)
	plt.ylabel('expression')
	plt.legend(loc = "lower right")
	plt.savefig(path_to_plots+"phenotype_expression_grouped_by_genes_wo_label.png")
	plt.close()

	return phe_count_percent

num_stability_to_consider = 3
replicates = ['r1','r2','r3']

def run_for_all_analysis(replicate):
	core_path = "./../../"
	
	## running for self_activation perturbation
	d = {}
	
	network_name = "core"
	path_to_dat_files = core_path+"topo_file/"+replicate+"/"
	name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
	genes = list(d_genes.values())
	path_to_output_z_norm = path_to_dat_files+"Z_normed/"
	path_to_output_phenotype = path_to_dat_files+"phenotype/"
	path_to_plots = path_to_dat_files+"plots/scatters/"
	path_to_heatmap_plots = path_to_dat_files+"plots/heatmaps/"
	path_to_score_file = path_to_dat_files+"phenotype/"
	if not os.path.exists(path_to_plots):
		os.makedirs(path_to_plots)
	state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
	state_dataframe, state_score_dataframe = collated_score_for_all_dataframes(state_dataframe)
	mean_EMT, std_EMT = read_scoring_file(path_to_score_file,'EMT')
	mean_TamRes, std_TamRes = read_scoring_file(path_to_score_file,'TamRes')
	boundaries = return_inidividual_phenotypes_boundary_condition(mean_EMT,mean_TamRes,std_EMT,std_TamRes)
	quantify_all_phenotypes(state_dataframe,boundaries,genes,path_to_plots)

for i in replicates:
	run_for_all_analysis(i)