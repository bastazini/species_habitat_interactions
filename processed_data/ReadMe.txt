Summary of the site occupancy model output
The file contains the following objects:

- cob_corals: coral cover, per species and site
- coordenadas: site coordinates
- coral_species: list of coral species
- fish_species: list of fish species
- extracted_data: list of model outputs. The lenght of this list is equal to the number of coral species. Each object of the list consist of a table
depicting the coral species which the model used, the fish (peixe) species, the age (1="adult",2="juvenile"), estimated regression coefficient (estimate.coral),
lower and upper credible intervals (low.coral and high.coral), estimated effect of turf cover (estimate.turf), psi_i estimate (psi), posterior exceedance probability of the coefficient being different from zero (p_exceedance, not used), and predicted occupancy (pred) 
- extracted_data_occ: Pearson's correlation between site occupancy probability of coral-associated fish and other fish.  The lenght of this list is equal to the number of coral species. 
- fish_size: fish traits
- network_data: processed data to build the network


Basic data contains the formatted data for network analysis.