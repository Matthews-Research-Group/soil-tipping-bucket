library(BioCro)
library(BioCroWater)

#Yufeng: I made up some of these values.
input_list = list(canopy_temperature = 298,      #Kelvin
                  gas_constant       = 8.314,    #J/mol/K
                  b                  = soybean$parameters$soil_b_coefficient,
                  soil_depth         = 1,        #m
                  soil_relative_water_content  = 0.5,   #m3/m3
                  plant_relative_water_content = 0.5,   #m3/m3
                  water_density      = 1000,            #kg/m3
                  structural_dry_matter_fraction = 0.05 #kg/(m^2 ground)
                  )
output_list = evaluate_module("BioCroWater:canopy_potential",input_list)
