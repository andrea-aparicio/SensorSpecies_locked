This is a repository for the paper 
"Identifying sensor species that anticipate critical transitions in multispecies ecological systems" 
by Andrea Aparicio, Jorge X. Velasco, Claude H. Moog, Yang-Yu Liu, and Marco Tulio Angulo.

The repository contains the files:

	For the results of Fig 1: 
		-Toy_Mutualistic.py

	For the simulations of the general systems: 
		-Define_Parameters.py - Code used to define the network parameters according to Section S1.2 
		-Solve_EmpiricalNet_Mutualistic.py - Code used to simulate a critical transition using a
		 mutualistic model according to Section S1.2 
		-Early_Warnings_Score.py - Code used to calculate the variance change and Early-Warnings score 
		 according to Section S1.3 
		-Sensor_Score.nb - Code used to calculate the Sensor Score according to Section S2 
		-build_community_matrix.py - File containing functions implemented by Define_Parameters.py 
		-Network_ID.py - File where the working WoL Network ID is specified
		-General_Plots.nb - Builds the plots for the collection of networks

The repository contains the folders:

	-WoL_matrices - To store the network structure .csv file obtained from WoL website 
	-Community_parameters - To store the network parameters defined in Define_Parameters.py 
	-data - To store the simulation data obtained in Solve_EmpiricalNet_Mutualistic.py, 
	 and the Early-Warnings score data obtained in Early_Warnings_Score.py
	-General_data - contains the summarized sensor score and early-warning score data of the 
	 collection of networks. 

Requirements:
	-python3 and packages csv, numpy, sdeint, matplotlib, pandas, scipy.integrate, math, and scipy.stats 
	-Mathematica 12

Instructions to simulate the toy ecosystems:

	Run Toy_Mutualistic.py. This will simulate the perturbed and critical transition cases, and produce 
	 the plots in Figure 1.

Instructions to simulate the general systems:

	1. Download to the folder Community_parameters the network matrix in .csv format of a mutualistic, 
	 single-component ecological network from http://www.web-of-life.es/ (e.g. M_SD_020.csv)
	2. Change the network ID in the file Network_ID.py (e.g. fn = "M_SD_020")
	3. Run Define_Parameters.py. This will produce five files containing 
		-initial conditions (e.g., M_SD_020_x0.csv) 
		-growth rate (e.g., M_SD_020_r.csv) 
		-competition matrix (e.g., M_SD_020_Ac.csv) 
		-mutualistic interaction matrix (e.g., M_SD_020_A.csv)
		-structural mutualistic interaction matrix (0-1 structure) (e.g., M_SD_020_Abin.csv)
	 They will be stored in the folder Community_parameters
	4.Run Solve_EmpiricalNet_Mutualistic.py. This will produce six files containing, for every species and 
	 every value of \mu: 
		-the variance of the abundance with and without interspecific competition 
		 (e.g., M_SD_020_var1.csv, M_SD_020_var0.csv, respectively) 
		-the mean of the abundance with and without interspecific competition 
		 (e.g., M_SD_020_mean1.csv, M_SD_020_mean0.csv, respectively) 
		-a sample of the abundance with and without interspecific competition 
		 (e.g., M_SD_020_samp1.csv, M_SD_020_samp0.csv, respectively) 
	 They will be stored in the folder data
	5. Run Early_Warnings_Score.py. This will produce two files containing 
		-Early-Warnings score of every species with interspecific competition 
		 (e.g., M_SD_020_detection1.csv) 
		-Early-Warnings score of every species without interspecific competition 
		 (e.g., M_SD_020_detection0.csv) 
	 They will be stored in the folder data
	6. Run Sensor_Score.nb that will calculate the sensor score of every species and produce the plots of 
	 Figure 2
	7. Run General_Plots.nb that will produce the plots of Figure 3


Developed by: Andrea Aparicio