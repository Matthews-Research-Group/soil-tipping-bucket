library(BioCro)
library(BioCroWater)
library(DEoptim)

# Cost function
source('obj_function.R')
source('set_modules_for_soil_water.R')

use_pre_opt_vars = FALSE 
# Names of parameters being fit
if(use_pre_opt_vars){
  pre_arg_names <- c('alphaLeaf','alphaRoot','alphaStem','betaLeaf','betaRoot','betaStem',
                 'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
                 'alphaSeneStem','betaSeneStem','alphaShell','betaShell') 
  
  arg_names <- c('phi2') 
  
  pre_par = readRDS(paste0('opt_results/opt_result_old_water_opt_v1.rds'))
  pre_par = pre_par$optim$bestmem
}else{
  arg_names_all <- c('alphaLeaf','alphaRoot','alphaStem','betaLeaf','betaRoot','betaStem',
                 'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
                 'alphaSeneStem','betaSeneStem','alphaShell','betaShell','mrc1','mrc2')
  arg_names <- c('alphaLeaf','alphaStem','betaLeaf','betaStem',
                 'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
                 'alphaSeneStem','betaSeneStem','alphaShell','betaShell','mrc2')
}

#these weights are picked without specific reasons
#the current values work fine for my case
wts_on_errors <- data.frame("Stem" = 1, "Leaf" = 2, "Shell" = 0.5, "Seed" = 1,"TotalLitter" = 0.1,"Root" = 0.1)

use_new_water   = FALSE 
output_filename = "opt_result_old_water_v11.rds"

print(arg_names)
print(wts_on_errors)
print(output_filename)

# years, sowing dates, and harvesting dates of growing seasons being fit to
# optimization uses these two years as Matthews et al (2022). 
# The other two years (2004&2006) can be considered as validation sets
year      <- c('2002', '2005')
sow.date  <- c(152, 148)
harv.date <- c(288, 270)

## Initialize variables for the cost function
soybean_optsolver <- list()
ExpBiomass <- list()
ExpBiomass.std <- list()
RootVals <- list()
weights <- list()
numrows <- vector()

ctl <- list()

if(use_new_water){
  soybean_steadystate_modules0 = set_direct_modules() 
  soybean_derivative_modules0  = set_differential_modules() 
  soybean_initial_state0       = set_init_values() 
  soybean_parameters0          = set_parameters() 
  soybean_parameters0$kd = soybean_parameters0$k_diffuse
}else{
  soybean_steadystate_modules0 = soybean$direct_modules 
  soybean_derivative_modules0  = soybean$differential_modules 
  soybean_initial_state0       = soybean$initial_values 
  soybean_parameters0          = soybean$parameters 
}

#pre-define the growth_respiration_fraction to reduce Assimilation
#soybean_parameters0$growth_respiration_fraction = 0.25 

if(use_pre_opt_vars){
  soybean_parameters0[pre_arg_names] = pre_par 
}


soybean_solver_params=soybean$ode_solver
#use euler to speed up the optimization process.
#you may use the non-euler in the evaluation process. 
#I have not seen any sig difference between solvers
soybean_solver_params$type="homemade_euler"

for (i in 1:length(year)) {
  yr <- year[i]
  weather <- read.csv(file = paste0('Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))

  sd.ind <- which(weather$doy == sow.date[i])[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather_growing_season <- weather[sd.ind:hd.ind,]
  weather_growing_season$time_zone_offset = -6
  
#calculate CTL for diagnostic purposes
  ctl[[i]] <- run_biocro(
    soybean$initial_values,
    soybean$parameters,
    weather_growing_season,
    soybean$direct_modules,
    soybean$differential_modules,
    soybean$ode_solver
  )

  soybean_optsolver[[i]] <- partial_run_biocro(soybean_initial_state0, 
 					       soybean_parameters0, 
					       weather_growing_season,
                                               soybean_steadystate_modules0,
					       soybean_derivative_modules0, 
                                               soybean_solver_params,arg_names)
  
  ExpBiomass[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_ambient_biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Shell0","Seed","Litter","CumLitter")
  Shell = ExpBiomass[[i]]$Shell0 - ExpBiomass[[i]]$Seed
#make Shell not decline. When shell reaches peak, we assume it does not decline after 
#because in BioCro, we do not have shell senescence 
  Shell[which.max(Shell):length(Shell)] = max(Shell) 
  ExpBiomass[[i]]$Shell = Shell 

  ExpBiomass[[i]]$Shell0 = NULL
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_ambient_biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Shell0","Seed","Litter","CumLitter")
  ExpBiomass.std[[i]]$Shell = sqrt((ExpBiomass.std[[i]]$Shell0)^2 + (ExpBiomass.std[[i]]$Seed)^2)
  ExpBiomass.std[[i]]$Shell0 = NULL
  
#why use the 5th DOY?
  row_to_use = 5
  Root_estimate = 0.17*sum(ExpBiomass[[i]][row_to_use,c('Leaf','Stem','Seed','Shell')])
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[row_to_use], "Root"=Root_estimate) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather_growing_season)
  invwts <- ExpBiomass.std[[i]]
  #this weight is to put importance based on std
  #the larger the error bar is, the less the weight it has
  weights[[i]] <- log(1/(invwts[,2:ncol(invwts)]+1e-5))

}


## Optimization settings
ul = 50 
ll = -50

# parameter upper limit
upperlim_all<-c(ul,ul,ul,
            0,0,0,
            0.0125,.005,
            ul,0,ul,0,
            ul,0,0.08,0.075)

# parameter lower limit
lowerlim_all<-c(0,0,0,
            ll,ll,ll,
            0,0,
            0,ll,0,ll,
            0,ll,8e-4,0.0025)

upperlim = upperlim_all[match(arg_names,arg_names_all)]
lowerlim = lowerlim_all[match(arg_names,arg_names_all)]

if(length(arg_names) != length(upperlim)) stop("number of vars does not match the number of bounds")

# cost function
cost_func <- function(x){
  multiyear_BioCro_optim_obj(x, soybean_optsolver, ExpBiomass, numrows, weights, wts_on_errors, RootVals)
}


#fix a random seed for repeatability 
rng.seed <- 1234 
set.seed(rng.seed)

# initialize parameters being fitted as randome values from a uniform distribution
opt_pars <- runif(length(arg_names), min = lowerlim, max = upperlim)
# maximum number of iterations
max.iter <- 2000
# number of cores for parallelization
# you probably should pick a number less than the total cores on your machine
nc = 8
cl <- makeCluster(nc)

# Call DEoptim function to run optimization
parVars <- c('multiyear_BioCro_optim_obj','soybean_optsolver','ExpBiomass','numrows','weights','wts_on_errors','RootVals')
clusterExport(cl, parVars,envir=environment())
optim_result<-DEoptim(fn=cost_func, lower=lowerlim, upper = upperlim, 
		     control=list(itermax=max.iter,parallelType=1,
                     packages=c('BioCro'),parVar=parVars,cluster=cl))
#optim_result<-DEoptim(fn=cost_func, lower=lowerlim, upper = upperlim, 
#		     control=list(itermax=max.iter,parallelType=0,
#                     packages=c('BioCro')))

#save the optimized parameters to a file, which will be used for plotting later
saveRDS(optim_result,paste0('opt_results/',output_filename))
