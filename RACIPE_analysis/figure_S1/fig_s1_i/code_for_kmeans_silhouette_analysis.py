import os
import itertools
import numpy as np
import collections
import pandas as pd
import seaborn as sns
from gap_statistic import OptimalK
from textwrap import wrap
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.spatial.distance import cdist
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.mixture import GaussianMixture
from matplotlib.colors import ListedColormap

"""This function estimates the silhouette coefficient for the given data and plots it as a boxplot"""

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

#This function performs Kmeans clustering on data and then performs silhoutte analysis to confirm the number of clusters.  

def kmeans_silhouette(state_dataframe_PCA,num_clusters,path_to_plots,network_name):
    score_array = [num_clusters]
    for N in range(15):   #Kmeans is performed 15 times 
        range_n_clusters = num_clusters
        x = []
        for n_clusters in range_n_clusters:    #Kmeans is performed for the number of clusters defined in num_clusters array
            clusterer = KMeans(n_clusters=n_clusters)
            preds = clusterer.fit_predict(state_dataframe_PCA)
            centers = clusterer.cluster_centers_
            score = silhouette_score(state_dataframe_PCA, preds)
            x.append(score)
        score_array = np.vstack((score_array,x))
    #print(score_array)
    silhouette_score_data_frame = pd.DataFrame(score_array)
    silhouette_score_data_frame = silhouette_score_data_frame.drop([0])
    
    column_name=[]
    for i, j in zip(num_clusters, range(len(num_clusters))):
        name = str(i) + ' clusters'
        silhouette_score_data_frame = silhouette_score_data_frame.rename(columns={j:name})
        column_name.append(name)

    # silhouette_score_data_frame.to_csv(path_to_plots+network_name+"_silhoutte.dat",sep="\t", header=True, index=True)
    
    #plotting the boxplot to find the ideal number of clusters for k means 
    fig, ax = plt.subplots()
    plot = silhouette_score_data_frame.boxplot(grid=False)
    ax.set(xlabel='Clusters', ylabel='mean silhouette coeff', title='Boxplot grouped by cluster pca_clustering')
    plt.savefig(path_to_plots+"boxplot_for_silhouette_analysis_png")
    plt.close()

    return silhouette_score_data_frame_mean


#This function is used to generate box plots 
def plot_a_box_plot(data,path_to_plots,data_type):
    data = data.transpose()
    fig, ax = plt.subplots()
    plot = data.boxplot(grid=False)
    plt.savefig(path_to_plots+"boxplot_"+data_type+".png")
    plt.close()



num_stability_to_consider = 4
replicate_array = ['r1','r2','r3']
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
    path_to_output_z_norm = path_to_dat_files+"Z_normed/"
    path_to_plots = path_to_dat_files+"plots/scatter_cluster/"
    if not os.path.exists(path_to_plots):
        os.makedirs(path_to_plots) 
    state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
    sub_dataframe = state_dataframe.sample(n=10000)
    num_clusters = [2,3,4,5,6,7,8]
    silhouette_score_data_frame = kmeans_silhouette(sub_dataframe,num_clusters,path_to_plots,network_name)

    return silhouette_score_data_frame


silhouette_score_data_frame_mean = pd.DataFrame()

for replicate in replicates_array:
	silhouette_score_data_frame = run_for_all_analysis(replicate)
    silhouette_score_data_frame_mean[replicate] = silhouette_score_data_frame['mean']

plot_a_box_plot(silhouette_score_data_frame_mean,path_to_replicate_plots,'silhoutte')