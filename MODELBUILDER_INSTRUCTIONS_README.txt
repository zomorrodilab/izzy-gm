ComPy: A MICROBIAL COMMUNITY MODELING PIPELINE 

created by Izzy Goodchild-Michelman in June 2022. contact igoodchildmichelman@college.harvard.edu with any questions

OVERVIEW: Compy is a pipeline written in python that is inspired by the microbiome modeling toolbox (https://doi.org/10.1093/bioinformatics/bty941). It combines single species genome scale metabolic models of metabolism (GEMs) into a community level GEM. The community level objective function maximizes the growth of the sum of the individual species's biomass reactions, constrained by the provided relative abundances of each species. Currently the pipeline only works with AGORA models, I am working on expanding it to accept other model types (Bigg, CarveMe, Kbase etc.) 



REQUIREMENTS:

1. A saved directory of AGORA models. I have copied all 7,200 AGORA2 models into the microbiome_modeling_pipeline folder (Dropbox/Zomorrodi_lab_shared/microbiome_modeling_pipeline/AGORA2_mat) that you can use. 

2. Python installed (I created it to work with Python 3.9.12)

3. Several python packages: cobra, pandas, scipy (these can be installed using pip on your command line)




INPUTS:

abun_path: path to the species abundance .csv file
* see file (/Dropbox/Zomorrodi_lab_shared/microbiome_modeling_pipeline/normCoverageReduced.csv) in shared dropbox for how to format
* make sure the column names match the names of the AGORA model exactly (do not include .mat)
* include an X as the header for the individual species model name column

modpath: path to folder with all AGORA models 
EX: '/Users/isabellagm/Dropbox/Zomorrodi_lab_shared/microbiome_modeling_pipeline/AGORA2_mat'

respath: path where the community models will be output 
EX: "/Users/isabellagm/Desktop/working_compy_models/"

dietpath: path to the AGORA compatible diet (for community model) .csv file
* see file (/Dropbox/Zomorrodi_lab_shared/microbiome_modeling_pipeline/WesternDietAGORA.txt) for how to format
* take the VMH metabolite name add an "EX_" to the beginning and an "[d]" to the end for proper formatting, for example: if you want to add xyl_D to your diet you would add EX_xyl_D[d]



OUTPUTS:
All community models for the number of samples you inputed in the species abundance file (the columns of the file) in .json format. 



HOW TO RUN THE MODEL BUILDER: 

1. Open terminal (or non Mac equivalent command line tool) 
cd to the ".../Dropbox/Zomorrodi_lab_shared/microbiome_modeling_pipeline" folder 
2. Start python in terminal
3. type "from compy import compy"
4. Run the pipeline, ex:

compy(abunpath="/Users/theogoodchild-michelman/Dropbox/BLOOM_metabolic_modeling_shared/bloom_microbiome_models/lateabun_w78.csv",
     modpath='/Users/theogoodchild-michelman/Desktop/research/MGH/AGORA-1.03/AGORA2_mat/',
     respath="/Users/theogoodchild-michelman/Desktop/")


compy(abunpath="/Users/theogoodchild-michelman/Dropbox/BLOOM_metabolic_modeling_shared/bloom_microbiome_models/lateabun_w78_short.csv", 
	modpath="/Users/theogoodchild-michelman/Desktop/research/summer_research/AGORA_2_new/", 
	respath="/Users/theogoodchild-michelman/Downloads/")


