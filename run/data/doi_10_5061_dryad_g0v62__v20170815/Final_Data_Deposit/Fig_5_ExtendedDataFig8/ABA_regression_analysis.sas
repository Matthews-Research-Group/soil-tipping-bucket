Data ABA_regression_data;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Fig_5_Extended_Data_Fig_8\ABA_soilVWC_middaygasex09_11.csv' 
dlm=',' firstobs=2 lrecl=10000; input Obs block ring CO2 h2o$ DOY year xylem_ABA CO2_num Ci CO2xyl pH logABA_corr VWCthru75 CO2_VWC75 logxyl ;;

**********Figure 5A regression full model******************;
proc reg data=ABA_regression_data;
model Ci=xylem_ABA CO2_num CO2xyl;
test CO2_num=0;
test ABA_corr=0;
run;
*********Figure 5A treatment-specific regression to yield R^2 values for each treatment;
***regression of Ci to xylem ABA for elevated CO2 treatment;
data elevated_CO2;
set ABA_regression_data;
where CO2=585;
proc reg data=elevated_CO2;
model Ci=xylem_ABA;
run;
***regression of Ci to xylem ABA for ambient CO2 treatment;
data ambient_CO2;
set ABA_regression_data;
where CO2=385;
proc reg data=ambient_CO2;
model Ci=xylem_ABA;
run;

**********Extended Data Figure 8A regression full model*****************;
proc reg data=ABA_regression_data;
model logABA_corr=VWCthru75 CO2_num CO2_VWC75;
test CO2_num=0;
test VWCthru75=0;
run;
***Extended Data Figure 8A treatment-specific regressions;
***ambient CO2;
proc reg data=ambient_CO2;
model logABA_corr=VWCthru75;
run;
***elevated CO2;
proc reg data=elevated_CO2;
model logABA_corr=VWCthru75;
run;


**********Extended Data Figure 8B regression full model*****************;
proc reg data=ABA_regression_data;
model logxyl=VWCthru75 CO2_num CO2_VWC75;
test CO2_num=0;
test VWCthru75=0;
run;
***Extended Data Figure 8B treatment-specific regressions;
***ambient CO2;
proc reg data=ambient_CO2;
model logxyl=VWCthru75;
run;
***elevated CO2;
proc reg data=elevated_CO2;
model logxyl=VWCthru75;
run;

**********Extended Data Figure 8C regression full model*****************;
proc reg data=ABA_regression_data;
model pH=VWCthru75 CO2_num CO2_VWC75;
test CO2_num=0;
test VWCthru75=0;
run;
***Extended Data Figure 8C treatment-specific regressions;
***ambient CO2;
proc reg data=ambient_CO2;
model pH=VWCthru75;
run;
***elevated CO2;
proc reg data=elevated_CO2;
model pH=VWCthru75;
run;
