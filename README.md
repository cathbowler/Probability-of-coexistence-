# Probability-of-coexistence-
Determining the probability of coexistence: Adding demographic variation to modern coexistence modelling

Instructions for running code:
- First run Make.comp.datasets.pairwise.R in the ‘raw data’ folder 
- Run Stan_models_all_species.R which calls on the prepared datasets and the Stan models (see BH stan models)
- In 'SurvivalAndGermination' folder, run Germination.R and Survival.R to get posterior probabilites for germination and survival rates for each species 

Can then simulation coexistence using either the Invasion or Ratio approach as follows:
    1. Invasion approach 
        1. First run Residents.to.equilibrium.R 
        2. And run Residents.to.equilibrium.means.R
        3. Then run Invasion.approach.revisions.R which simulates coexistence and propagates uncertainty through posterior probabilities 
        4. Also run Invasion.approach.ldgr.means.R which simulates coexistence the traditional way using the mean of each posterior parameter in the calculations 
        
    2. Ratio approach
      1. BH.equalizing.stabilizing.R calculates the equalizing and stabilizing ratios using mean values for each posterior parameter 
      2. Run CB.equalizing.stabilizing.distributions.R to calculate the ratio using the full posterior probilities for each parameter 