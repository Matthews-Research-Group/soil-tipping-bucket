library(BioCro)
library(BioCroWater)
years = 2004:2011
obs_data = readRDS("data/soil_moisture_data_rearranged.rds")
obs_data$sm = obs_data$sm/100
results = list()
for (i in 1:length(years)){
  year = years[i]
  weatherData <- weather[[as.character(year)]]
  weather_growing_season <- get_growing_season_climate(weatherData, threshold_temperature = 0)
  soiltype = 3
  direct_modules_new = list("BioCroWater:soil_type_selector", "BioCroWater:soil_surface_runoff",
                            "BioCroWater:soil_water_downflow", "BioCroWater:soil_water_tiledrain", 
                            "BioCroWater:soil_water_upflow","BioCroWater:multilayer_soil_profile_avg")
  
  differential_modules_new = list("BioCroWater:soil_evaporation","BioCroWater:multi_layer_soil_profile")
  
  obs_year_i = obs_data[obs_data$year==year,]
  init_values =   within(miscanthus_x_giganteus$initial_values,{
    soil_water_content_1 = obs_year_i$sm[obs_year_i$depth==1][1]  #10cm
    soil_water_content_2 = obs_year_i$sm[obs_year_i$depth==2][1]  #20cm
    soil_water_content_3 = obs_year_i$sm[obs_year_i$depth==3][3]  #30cm
    soil_water_content_4 = obs_year_i$sm[obs_year_i$depth==5][3]  #50cm
    soil_water_content_5 = obs_year_i$sm[obs_year_i$depth==5][3]  #50cm
    soil_water_content_6 = obs_year_i$sm[obs_year_i$depth==5][3]  #50cm
    sumes1               = 0.125
    sumes2               = 0
    time_factor = 0
    soil_evaporation_rate = 0
  })
  init_values$soil_water_content=NULL
  parameters =   within(miscanthus_x_giganteus$parameters, {
    irrigation = 0.0
    swcon = 0.05/24
    curve_number = 61
    soil_type_indicator_1 = soiltype
    soil_type_indicator_2 = soiltype
    soil_type_indicator_3 = soiltype
    soil_type_indicator_4 = soiltype
    soil_type_indicator_5 = soiltype
    soil_type_indicator_6 = soiltype
    soil_depth_1 = 5
    soil_depth_2 = 15
    soil_depth_3 = 30
    soil_depth_4 = 50
    soil_depth_5 = 50
    soil_depth_6 = 50
    tile_drainage_rate          = 0.2 # dimensionless, as in DSSAT (1/d ?)
    tile_drain_depth            = 90  # cm (typical value 3-4 ft)
    skc = 0.3
    kcbmax = 0.25
    elevation                   = 219
    lai = 5
    bare_soil_albedo            = 0.15
  })
  parameters[c('soil_field_capacity','soil_saturated_conductivity','soil_saturation_capacity','soil_wilting_point')]=NULL
  
  results[[i]] <- run_biocro(
    init_values,
    parameters,
    weather_growing_season,
    direct_modules_new,
    differential_modules_new
  )
}

#plotting
par(mfrow=c(2,2),mar=c(4,4,3,3),mai=c(1,1,0.5,0.5),cex=1.5) #c(bottom, left, top, right)
for (i in 1:4){
  year   = years[i]
  result = results[[i]]
  obs_year_i = obs_data[obs_data$year==year,]
  weatherData <- weather[[as.character(year)]]
  weather_growing_season <- get_growing_season_climate(weatherData, threshold_temperature = 0)
  linewidth=2
  axislabelspace = 2
  plot(result$time, result$soil_water_content,type='l',
       lwd = linewidth,
       xlab='',ylab='',ylim=c(0,0.5),
       main = paste0("year ",year))
  obs = obs_year_i[obs_year_i$depth<=5,] #top 50cm average for rings and depths
  obs_avg = aggregate(sm ~ doy, data = obs, FUN = mean)
  lines(obs_avg$doy,obs_avg$sm,col="red",lwd = linewidth,)
  lines(weather_growing_season$doy,weather_growing_season$precip/max(weather_growing_season$precip)/2,col="blue")
  if(i==1){
  legend(x=130,y=0.6, legend=c("BioCro", "OBS","Precip"),
         col=c("black","red","blue"), lty=1, cex=0.8,bty="n",lwd = linewidth,
         seg.len=0.5, x.intersp = 0.2, y.intersp = 0.2)
  }
  title(xlab = "time", line = axislabelspace)            # Add x-axis text
  title(ylab = "soil water content", line = axislabelspace)            # Add y-axis text
  # lines(result$time,result$soil_evaporation_rate/2,col="blue")
}
