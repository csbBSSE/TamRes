Before running the pipelines for the figure generation, run [`RACIPE`](https://github.com/simonhb1990/RACIPE-1.0) for the topo file given in topo file folder. It will generate the RACIPE data which can be z-normalized using the code in `Misc_code` folder and used for further analysis.
This directory contains programs written in python script that can generate the plots represented in respective figures. Each folder contains subfolders with the codes written for the generation of each figure. Each program is an isolated code with all necessary functions included in that file. Before running each script, Z score normalization of the output from RACIPE is necessary which can be done by the `code_to_generate_z_score_norm_files`.py present in `Misc_Code` folder. Besides that, `code_to_generate_score_distribuitions.py` is also necessary to generate statistics for the estimation of boundary condition required for the identification of phenotypes used in figure 1 F, 3 A , 3 D and 5 A.



