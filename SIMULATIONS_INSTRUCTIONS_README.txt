ComPy: A MICROBIAL COMMUNITY MODELING PIPELINE 

created by Izzy Goodchild-Michelman in December 2022. contact igoodchildmichelman@college.harvard.edu with any questions

Requirements:

1. A .json community metabolic model, the output from the model builder pipeline (see MODELBUILDER_INSTRUCTIONS_README.txt file to run this first)

2. A properly formatted diet file ** (see AverageEU_diet_fluxes.txt file for example)

** take the VMH metabolite name add an "Diet_EX_" to the beginning and 	an "[d]" to the end for proper formatting, for example: if you want to add xyl_D to your diet you would add "Diet_EX_xyl_D[d]"




HOW TO SIMULATE COMMUNITY GROWTH:

1. Open your python IDE of choice and load in your community model: 

	import cobra
	import pandas as pd
	from cobra import Model, Reaction, Metabolite
	from cobra.flux_analysis import flux_variability_analysis
	import cplex
	import csv
	

	model = cobra.io.load_json_model(PATH TO YOUR COMMUNITY MODEL FILE)

2. Load in your diet file:

	AvEuroDiet_dic = {}
	import csv
	with open(PATH TO YOUR DIET FILE) as f:
    		next(f)
   		for key, *values in csv.reader(f,  delimiter="\t"):
        		AvEuroDiet_dic[key.strip()] = - float(values[0])

3. bound the community biomass based on reasoning from:
Heinken, A., Ravcheev, D.A., Baldini, F. et al. Systematic assessment of secondary bile acid metabolism in gut microbes reveals distinct metabolic capabilities in inflammatory bowel disease. Microbiome 7, 75 (2019).
https://doi.org/10.1186/s40168-019-0689-3

	model.reactions.get_by_id("communityBiomass").bounds = (0.4,1)

4. First maximize the flux through all fecal exchange reactions:
	
	loop_reactions = []
	for reac in [i.id for i in model.reactions if "UFEt" in i.id]:
   		loop_reactions.append(reac)

	final_fecalexchange = []
	for i in loop_reactions:
   		model.objective = i
  		solution = model.optimize()
    		print(i)
   		print(solution.objective_value)

5. Then iterating over each fecal exchange reaction, constrain the flux of that reaction to its calculated maximum from step 4. With this constraint in place, then minimize the flux through each of the associated species specific metabolite exchange reactions:
	
	final_EX_fluxes = {}
	for i in range(len(loop_reactions)):
    		if final_fecalexchange[i] != 0:
      			old_bounds = model.reactions.get_by_id(loop_reactions[i]).bounds
       			model.reactions.get_by_id(loop_reactions[i]).bounds = (final_fecalexchange[i],final_fecalexchange[i])
        		metabolite = loop_reactions[i].replace("UFEt_","") + str("[u]")
        		for reac in model.metabolites.get_by_id(metabolite).reactions:
            			if "IEX" in reac.id:
               				model.objective = model.reactions.get_by_id(reac.id)
                			solution = model.optimize(objective_sense = "minimize")
                			final_EX_fluxes[reac.id] = solution.objective_value
    
    		model.reactions.get_by_id(loop_reactions[i]).bounds = old_bounds


The dictionary final_EX_fluxes will now hold the flux through each species exchange metabolite

		