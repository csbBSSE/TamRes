import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from textwrap import wrap
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from sklearn.cluster import AgglomerativeClustering
from sklearn.mixture import GaussianMixture


"""EM score and Resistance score show trimodality and bimodality in their distribution. The functions in this file performs
the gaussian fitting in EM and Resistance score distribution and also validates their tri/bimodality
through information criterion estimation."""


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

"""This function performs multiple gaussian fitting on score histograms and also estimates AIC/BIC to find the optimal number of gaussians.
It also generates a statistics file that contains information on mean and standard deviation for each gaussian in the distribution. These values
are used to discern phenotype based upon the score for a cell state/solution. """

def multiple_gaussian_fitting(array_to_be_plotted,gene_name,path_to_plots,path_to_score_files,num_gaussians,sample_size,max_clusts,phenotype):

	AIC_all = {}
	BIC_all = {}

	file = open(path_to_score_files+'_'+phenotype+'score_stats.txt','w')
	for i in range(1,max_clusts+1):
		AIC_all[i] = []
		BIC_all[i] = []
	for sample in range(1,sample_size):

		random_state = np.random.RandomState(seed=1)

		X = np.array(array_to_be_plotted).reshape(-1, 1)

		# fit models with 1-10 components
		N = np.arange(1, max_clusts+1)
		models = [None for i in range(len(N))]

		for i in range(len(N)):
			models[i] = GaussianMixture(N[i]).fit(X)

		# compute the AIC and the BIC
		AIC = [m.aic(X) for m in models]
		BIC = [m.bic(X) for m in models]

		for idx,val in enumerate(AIC):
			AIC_all[idx+1].append(AIC[idx])
			BIC_all[idx+1].append(BIC[idx]) 

		M_best = models[num_gaussians-1]              # change number of gaussians here
		s = ""
		mean_array = []
		if num_gaussians>2:
			score_array = np.empty((15,9))
		else:
			score_array = np.empty((15,9))
			
		# for x in range(num_gaussians):
		# 	s += str(M_best.means_[x][0])+"\t"+str(M_best.covariances_[x][0][0])+"\t"+str(M_best.weights_[x])+"\t"
		# 	#print(M_best.means_[0][0],"\t",M_best.means_[1][0],"\t",M_best.means_[2][0],"\t",M_best.covariances_[0][0][0],"\t",M_best.covariances_[1][0][0],"\t",M_best.covariances_[2][0][0],"\t",M_best.weights_[0],"\t",M_best.weights_[1],"\t",M_best.weights_[2])
		# print(s)
		for x in range(num_gaussians):
			mean_array.append(M_best.means_[x][0])

		max_mean = max(mean_array)
		min_mean = min(mean_array)
		y = mean_array.index(max_mean)
		z = mean_array.index(min_mean)

		s += str(M_best.means_[y][0])+"\t"+str(M_best.covariances_[y][0][0])+"\t"+str(M_best.weights_[y])+"\t"
		s += str(M_best.means_[z][0])+"\t"+str(M_best.covariances_[z][0][0])+"\t"+str(M_best.weights_[z])+"\t"
		
		if(num_gaussians>2):
			for k in range(num_gaussians):
				if (k != y and k != z):
					s += str(M_best.means_[k][0])+"\t"+str(M_best.covariances_[k][0][0])+"\t"+str(M_best.weights_[k])+"\t"
					hybrid = M_best.means_[k][0]
		file.write(s+"\n")
		#print(s)

	file.close()

	# print(AIC_all)
	# print(BIC_all)


	fig = plt.figure()
	#fig.subplots_adjust(left=0.12, right=0.97,bottom=0.21, top=0.9, wspace=0.5)


	# plot 1: data + best-fit Mixture
	fig, ax = plt.subplots(figsize=(10,5))

	x = np.linspace(-6, 6, 1000)
	logprob = M_best.score_samples(x.reshape(-1, 1))
	responsibilities = M_best.predict_proba(x.reshape(-1, 1))
	pdf = np.exp(logprob)
	pdf_individual = responsibilities * pdf[:, np.newaxis]

	ax.hist(X, 30, density=True, histtype='stepfilled', alpha=0.4)
	ax.plot(x, pdf, '-k')
	ax.plot(x, pdf_individual, '--k')
	#ax.text(0.04, 0.96, "Best-fit Mixture",ha='left', va='top', transform=ax.transAxes,fontsize=18)
	ax.set_xlim(-5.1,5.1)
	ax.set_xlabel('score', fontsize=16)
	ax.set_ylabel('Frequency', fontsize=16)
	if phenotype == 'EMT':
			plt.axvline(x = hybrid - (abs(hybrid - min_mean)/2), color ='r')
			plt.axvline(x = hybrid + (abs(hybrid - max_mean)/2), color ='r')
	else:
		plt.axvline(x = min_mean +(abs(min_mean - max_mean)/2), color ='r')
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.savefig(path_to_plots+gene_name+".png",dpi=500)
	plt.close()

	# removing the first occurance 
	# AIC_all.pop(1)
	# BIC_all.pop(1)

	AIC_dataframe = pd.DataFrame.from_dict(AIC_all,orient='columns')
	BIC_dataframe = pd.DataFrame.from_dict(BIC_all,orient='columns')

	AIC_dataframe.to_csv(path_to_plots+gene_name+"_AIC.dat",sep="\t",header=None)
	BIC_dataframe.to_csv(path_to_plots+gene_name+"_BIC.dat",sep="\t",header=None)


	# # plot 2: AIC plots
	fig, ax = plt.subplots(figsize=(10,5))
	ax.boxplot(AIC_all.values())
	ax.set_xticklabels(AIC_all.keys())
	plt.xlabel('Cluster Count')
	plt.ylabel('AIC metric')
	plt.title('AIC Values by Cluster Count')

	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.savefig(path_to_plots+gene_name+"_AIC.png",dpi=500)
	plt.close()
	

	# # plot 3: BIC plots
	fig, ax = plt.subplots(figsize=(10,5))
	ax.boxplot(BIC_all.values())
	ax.set_xticklabels(BIC_all.keys())
	plt.xlabel('Cluster Count')
	plt.ylabel('BIC metric')
	plt.title('BIC Values by Cluster Count')

	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.savefig(path_to_plots+gene_name+"_BIC.png",dpi=500)
	plt.close()

# This function generates the scores for EM and Resistance axis and calls the multiple_gaussian_fitting function to perform gaussian fitting. 

def scoring_relevant_features(state_dataframe,network_name,path_to_plots,path_to_score_files):
	state_dataframe["EMT_score"] = state_dataframe["ZEB1"] - state_dataframe["MIR200"]
	state_dataframe["TamRes_score"] = state_dataframe["ERa36"] - state_dataframe["ERa66"]
	# print("EMT_scores")
	sns.distplot(state_dataframe["EMT_score"])
	plt.title("EMT_scores   "+network_name)
	plt.savefig(path_to_plots+network_name+"_EMT_scoring_histogram.png")
	#plt.show()
	plt.close()
	multiple_gaussian_fitting(state_dataframe["EMT_score"],"EMT_score_dist",path_to_plots,path_to_score_files,3,15,6,'EMT')
	# print("TamRes_scores")
	sns.distplot(state_dataframe["TamRes_score"])
	plt.title("TamRes_score   "+network_name)
	plt.savefig(path_to_plots+network_name+"_TamRes_scoring_histogram.png")
	#plt.show()
	plt.close()
	multiple_gaussian_fitting(state_dataframe["TamRes_score"],"TamRes_score_dist",path_to_plots,path_to_score_files,2,15,6,'TamRes')


num_stability_to_consider = 4
replicate_array = ['r1','r2','r3']


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
	path_to_score_files = path_to_dat_files+"phenotype/"
	if not os.path.exists(path_to_plots):
		os.makedirs(path_to_plots)
	if not os.path.exists(path_to_score_files):
		os.makedirs(path_to_score_files)
	state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
	scoring_relevant_features(state_dataframe,network_name,path_to_plots,path_to_score_files)

for replicate in replicates:
	run_for_all_analysis(replicate)