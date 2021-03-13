import os
import itertools
import numpy as np
import scipy.stats
import pandas as pd
import seaborn as sns
from textwrap import wrap
from random import choices
import matplotlib.pyplot as plt


"""The functions in this file generate bar plot for each phenotype count percentage in wild type and perturbed network."""
#This function reads the CFG file to get information about the genes in the network.
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

""" This function reads file for each steady state solution and generates scores for cell state/model (rows of each dataframe). 
The phenotypes are assigned to each cell state based upon this score and the count of that phenotype is then determined"""

def phenotype_counter(mean_EMT,mean_TamRes,path_to_output_phenotype,network_name,boundaries,solution_num):
	state_dataframe = pd.read_csv(path_to_output_phenotype+network_name+"_solution_"+str(solution_num)+"_score.dat",delimiter="\t",header = 0)
	if solution_num == 1:
		phe = []
		for index, row in state_dataframe.iterrows():
			phe.append(whats_my_phenotype(row['EMT_score_1'],row['TamRes_score_1'],boundaries))
		state_dataframe["Phenotype"] = pd.Series(phe)
		phe_count = state_dataframe["Phenotype"].value_counts()
		total = sum(phe_count)
		phe_count_percent = (phe_count/total)*100
		phe_count_percent = phe_count_percent.rename_axis('phenotypes').to_frame('counts')

	
	elif solution_num == 2:
		phe = []
		for index, row in state_dataframe.iterrows():
			state1 = whats_my_phenotype(row['EMT_score_1'],row['TamRes_score_1'],boundaries)
			state2 = whats_my_phenotype(row['EMT_score_2'],row['TamRes_score_2'],boundaries)
			x = [state1,state2]
			x.sort()
			state = x[0]+x[1]
			phe.append(state)
		state_dataframe["Phenotype"] = pd.Series(phe)
		phe_count = state_dataframe["Phenotype"].value_counts()
		total = sum(phe_count)
		phe_count_percent = (phe_count/total)*100
		phe_count_percent = phe_count_percent.rename_axis('phenotypes').to_frame('counts')

	elif solution_num == 3:
		phe = []
		for index, row in state_dataframe.iterrows():
			state1 = whats_my_phenotype(row['EMT_score_1'],row['TamRes_score_1'],boundaries)
			state2 = whats_my_phenotype(row['EMT_score_2'],row['TamRes_score_2'],boundaries)
			state3 = whats_my_phenotype(row['EMT_score_3'],row['TamRes_score_3'],boundaries)
			x = [state1,state2,state3]
			x.sort()
			state = x[0]+x[1]+x[2]
			phe.append(state)
		state_dataframe["Phenotype"] = pd.Series(phe)
		phe_count = state_dataframe["Phenotype"].value_counts()
		total = sum(phe_count)
		phe_count_percent = (phe_count/total)*100
		phe_count_percent = phe_count_percent.rename_axis('phenotypes').to_frame('counts')

	# print(phe_count_percent)
	state_dataframe.to_csv(path_to_output_phenotype+network_name+"_solution_"+str(solution_num)+"_phenotype.dat",sep="\t", header=False)
	return state_dataframe,phe_count_percent

# This function determines the percentage of each steady state in the total number of simulated solutions (100000)
def counting_solutions(path_to_dat_files,network_name,num_stability):
	d = {}
	total = 0
	for i in range(1,num_stability+1):
		state_dataframe = pd.read_csv(path_to_dat_files+network_name+"_solution_"+str(i)+".dat",delimiter="\t",header = 0)
		count = len(state_dataframe.index)
		d[i] = abs(count/1000)
		total = total + count
	d['4+'] = abs((100000 - total)/1000)
	solution_count = pd.DataFrame.from_dict(d,orient='index',columns = ['solutions'])

	return solution_count

"""This function plots horizontal stacked bar plots for phenotype proportions for each steady state
It clubs the phenotypes with low counts into the others. The number of phenotypes getting clubbed into others is 
decided by the top_phenotypes input. """

def plot_stacked_bar_plot(data,path_to_output_phenotype,network_name,stability,color_code,top_phenotypes):
	data_to_plot = data.nlargest(top_phenotypes,'mean')
	legend_labels = data_to_plot.index.to_list()
	legend_labels.append('others')
	np_array_data = data_to_plot['mean'].to_numpy()
	others = 100 - np.sum(np_array_data,axis = 0)
	np_array_data = np.append(np_array_data,[others],axis = 0)
	np_array_data_cum = np_array_data.cumsum(axis=0)
	std_dev_data = data['std']
	np_array_std = std_dev_data.to_numpy()

	category_colors = plt.get_cmap(color_code)(np.linspace(0.15, 0.85, len(np_array_data)))

	for (i,colors) in zip(range(len(np_array_data)),category_colors):
		widths = np_array_data[i]
		starts = np_array_data_cum[i] - widths
		plt.barh('y', widths, left=starts, height=0.1,xerr = np_array_std[i],color = colors,edgecolor='black',error_kw=dict(lw=2,capsize=2,capthick=1))

	# plt.legend(legend_labels)
	plt.savefig(path_to_output_phenotype+network_name+"_"+stability+"_stacked_bar_plot_wo_legend.png",dpi=500)
	plt.close()

# This function plots horizontal stacked bar plots for proportion of each steady state 
def plot_stacked_bar_plot_solutions(data,path_to_output_phenotype,network_name,stability,color_code):
	data_to_plot = pd.DataFrame()
	data_to_plot = data['mean']
	legend_labels = data_to_plot.index.to_list()
	np_array_data = data_to_plot.to_numpy()
	np_array_data_cum = np_array_data.cumsum(axis=0)
	std_dev_data = data['std']
	np_array_std = std_dev_data.to_numpy()

	category_colors = plt.get_cmap(color_code)(np.linspace(0.15, 0.85, len(np_array_data)))

	for (i,colors) in zip(range(len(np_array_data)),category_colors):
		widths = np_array_data[i]
		starts = np_array_data_cum[i] - widths
		plt.barh('y', widths, left=starts, height=0.1,xerr = np_array_std[i],edgecolor='black',error_kw=dict(lw=3,capsize=4,capthick=1))

	# plt.legend(legend_labels)
	plt.savefig(path_to_output_phenotype+network_name+"_"+stability+"_stacked_bar_plot_wo_legend.png",dpi=500)
	plt.close()


replicates_array = ['r1','r2','r3']
path_to_bar_plots = "./../../topo_file/bar_plots"
if not os.path.exists(path_to_bar_plots):
	os.makedirs(path_to_bar_plots)

def run_for_all_analysis(replicate):
	core_path = "./../../"

	## running for self_activation perturbation
	print("Running for self-activation perturbations: ")
	network_name = "core"
	path_to_dat_files = core_path+"topo_file/"+replicate+"/"
	name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
	genes = list(d_genes.values())
	path_to_output_phenotype = path_to_dat_files+"phenotype/"
	path_to_score_file = path_to_dat_files+"phenotype/"
	mean_EMT,std_EMT = read_scoring_file(path_to_score_file,'EMT')
	mean_TamRes, std_TamRes = read_scoring_file(path_to_score_file,'TamRes')
	boundaries = return_inidividual_phenotypes_boundary_condition(mean_EMT,mean_TamRes,std_EMT,std_TamRes)
	mono_state_dataframe,mono_phe_count_percent = phenotype_counter(mean_EMT,mean_TamRes,path_to_output_phenotype,network_name,boundaries,1)
	bi_state_dataframe,bi_phe_count_percent = phenotype_counter(mean_EMT,mean_TamRes,path_to_output_phenotype,network_name,boundaries,2)
	tri_state_dataframe,tri_phe_count_percent = phenotype_counter(mean_EMT,mean_TamRes,path_to_output_phenotype,network_name,boundaries,3)
	solution_count = counting_solutions(path_to_dat_files,network_name,3)

	return mono_phe_count_percent,bi_phe_count_percent,tri_phe_count_percent,solution_count

mono_phe_percent = pd.DataFrame()
bi_phe_percent = pd.DataFrame()
tri_phe_percent = pd.DataFrame()
solution_count_df = pd.DataFrame()

for replicate in replicates_array:
	mono_phe_count_percent,bi_phe_count_percent,tri_phe_count_percent,solution_count = run_for_all_analysis(replicate)
	mono_phe_percent[replicate] = mono_phe_count_percent['counts']
	bi_phe_percent[replicate] = bi_phe_count_percent['counts']
	tri_phe_percent[replicate] = tri_phe_count_percent['counts']
	solution_count_df[i] = solution_count['solutions']

count = 1
for x in [mono_phe_percent,bi_phe_percent,tri_phe_percent]:
	x['mean'] = x.mean(axis=1)
	x['std'] = x.std(axis=1)
	count = count + 1

plot_stacked_bar_plot(mono_phe_percent,path_to_horizontal_bar_plots,'core','1','viridis',2)
plot_stacked_bar_plot(bi_phe_percent,path_to_horizontal_bar_plots,'core','2','cividis',3)
plot_stacked_bar_plot(tri_phe_percent,path_to_horizontal_bar_plots,'core','3','magma',5)

mean_data = mean_sdev_for_df(solution_count_df)
solution_count_df['mean'] = mean_data['mean']
solution_count_df['std'] = mean_data['std']

plot_stacked_bar_plot_solutions(solution_count_df,path_to_horizontal_bar_plots,'core','sol','PiYG')
