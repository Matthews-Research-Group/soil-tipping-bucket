library(BioCro)
library(DEoptim)

# Cost function
source('obj_function.R')
source('set_modules_for_soil_water.R')

# Names of parameters being fit
#arg_names <- c('alphaLeaf','alphaRoot','alphaStem','betaLeaf','betaRoot','betaStem',
#               'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
#               'alphaSeneStem','betaSeneStem','alphaShell','betaShell') 

# years, sowing dates, and harvesting dates of growing seasons being fit to
# optimization uses these two years as Matthews et al (2022). 
# The other two years (2004&2006) can be considered as validation sets
year <- c('2002', '2005')
sow.date <- c(152, 148)
harv.date <- c(288, 270)

## Initialize variables for the cost function
soybean_optsolver <- list()
ExpBiomass <- list()
ExpBiomass.std <- list()
RootVals <- list()
weights <- list()
numrows <- vector()

ctl <- list()

soybean_steadystate_modules0 = set_direct_modules() 
soybean_derivative_modules0  = set_differential_modules() 
soybean_initial_state0       = set_init_values() 
soybean_parameters0          = set_parameters() 


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
  
  weather.growingseason <- weather[sd.ind:hd.ind,]
  weather.growingseason$time_zone_offset = -6
  
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
					       weather.growingseason,
                                               soybean_steadystate_modules0,
					       soybean_derivative_modules0, 
                                               soybean_solver_params,arg_names)
  
  ExpBiomass[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_ambient_biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Shell0","Seed")
  Shell = ExpBiomass[[i]]$Shell0 - ExpBiomass[[i]]$Seed
#  Shell[which.max(Shell):length(Shell)] = max(Shell) #make Shell not decline
  ExpBiomass[[i]]$Shell = Shell 

  ExpBiomass[[i]]$Shell0 = NULL
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_ambient_biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Shell0","Seed")
  ExpBiomass.std[[i]]$Shell = sqrt((ExpBiomass.std[[i]]$Shell0)^2 + (ExpBiomass.std[[i]]$Seed)^2)
  ExpBiomass.std[[i]]$Shell0 = NULL
  
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[5], "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather.growingseason)
  invwts <- ExpBiomass.std[[i]]
  #this weight is to put importance based on std
  #the larger the error bar is, the less the weight it has
  weights[[i]] <- log(1/(invwts[,2:ncol(invwts)]+1e-5))
}

#these weights are picked without specific reasons
#the current values work fine for my case
wts2 <- data.frame("Stem" = 1, "Leaf" = 1, "Shell" = 1, "Seed" = 1,"Root" = 0.01)

## Optimization settings
ul = 50 
ll = -50

# parameter upper limit
upperlim<-c(ul,ul,ul,
            0,0,0,
            0.0125,.005,
            ul,0,ul,0,
            ul,0)

# parameter lower limit
lowerlim<-c(0,0,0,
            ll,ll,ll,
            0,0,
            0,ll,0,ll,
            0,ll)

# cost function
cost_func <- function(x){
  multiyear_BioCro_optim_obj(x, soybean_optsolver, ExpBiomass, numrows, weights, wts2, RootVals)
}


#fix a random seed for repeatability 
rng.seed <- 1234 
set.seed(rng.seed)

# initialize parameters being fitted as randome values from a uniform distribution
opt_pars <- runif(length(arg_names), min = lowerlim, max = upperlim)
# maximum number of iterations
max.iter <- 1000
# number of cores for parallelization
# you probably should pick a number less than the total cores on your machine
nc = 8
cl <- makeCluster(nc)

# Call DEoptim function to run optimization
parVars <- c('multiyear_BioCro_optim_obj','soybean_optsolver','ExpBiomass','numrows','weights','wts2','RootVals')
clusterExport(cl, parVars,envir=environment())
optim_result<-DEoptim(fn=cost_func, lower=lowerlim, upper = upperlim, 
		     control=list(itermax=max.iter,parallelType=1,
                     packages=c('BioCro'),parVar=parVars,cluster=cl))

#save the optimized parameters to a file, which will be used for plotting later
saveRDS(optim_result,'opt_results/opt_result_soil_water_r1.rds')
