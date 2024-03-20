Data alldepths_2011;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Fig_S7_Table_S7_Rhizotron\2011_subplot_means_RLD.csv'
dlm=',' firstobs=2 lrecl=10000; input block ring DOY CO2$ H2O$ category cm_cm3 TotLengthcm @@;


************************DEPTH CATEGORY 1 (CORRESPONDING TO 5-15 CM SOIL DEPTH********************;
data category1;
set alldepths_2011;
where category = 1;
run;
proc sort data =category1;
by block DOY co2 h2o;
run ;
proc mixed data =category1 order =data covtest;
class block DOY co2 h2o ;
model cm_cm3 = co2|h2o|DOY/ ddfm =kr;
random block block*co2 ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = rhizopdiffs_cat1;
ods output lsmeans = rhizolsmeans_cat1;
run ;
proc print data =rhizopdiffs_cat1; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizopdiffs_cat1_revised2014_revised2014.csv';
proc print data =rhizopdiffs_cat1; run; ODS CSV close;
proc print data =rhizolsmeans_cat1; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizolsmeans_cat1_revised2014_revised2014.csv';
proc print data =rhizolsmeans_cat1; run ; ODS CSV close;
************************DEPTH CATEGORY 2 (CORRESPONDING TO 15-25 CM SOIL DEPTH********************;
data category2;
set alldepths_2011;
where category = 2;
run;
proc sort data =category2;
by block DOY co2 h2o;
run ;
proc mixed data =category2 order =data covtest;
class block DOY co2 h2o;
model cm_cm3 = co2|h2o|DOY/ ddfm =kr;
random block block*co2 ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = rhizopdiffs_cat2;
ods output lsmeans = rhizolsmeans_cat2;
run ;
proc print data =rhizopdiffs_cat2; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizopdiffs_cat2_revised2014.csv';
proc print data =rhizopdiffs_cat2; run; ODS CSV close;
proc print data =rhizolsmeans_cat2; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizolsmeans_cat2_revised2014.csv';
proc print data =rhizolsmeans_cat2; run ; ODS CSV close;
************************DEPTH CATEGORY 3 (CORRESPONDING TO 25-35 CM SOIL DEPTH********************;
data category3;
set alldepths_2011;
where category = 3;
run;
proc sort data =category3;
by block DOY co2 h2o ;
run ;
proc mixed data =category3 order =data covtest;
class block DOY co2 h2o ;
model cm_cm3 = co2|h2o|DOY/ ddfm =kr;
random block block*co2 ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = rhizopdiffs_cat3;
ods output lsmeans = rhizolsmeans_cat3;

run ;
proc print data =rhizopdiffs_cat3; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizopdiffs_cat3_revised2014.csv';
proc print data =rhizopdiffs_cat3; run; ODS CSV close;
proc print data =rhizolsmeans_cat3; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizolsmeans_cat3_revised2014.csv';
proc print data =rhizolsmeans_cat3; run ; ODS CSV close;
************************DEPTH CATEGORY 4 (CORRESPONDING TO 35-45 CM SOIL DEPTH********************;
data category4;
set alldepths_2011;
where category = 4;
run;
proc sort data =category4;
by block DOY co2 h2o ;
run ;
proc mixed data =category4 order =data covtest;
class block DOY co2 h2o ;
model cm_cm3 = co2|h2o|DOY/ ddfm =kr;
random block block*co2 ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = rhizopdiffs_cat4;
ods output lsmeans = rhizolsmeans_cat4;
run ;
proc print data =rhizopdiffs_cat4; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizopdiffs_cat4_revised2014.csv';
proc print data =rhizopdiffs_cat4; run; ODS CSV close;
proc print data =rhizolsmeans_cat4; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizolsmeans_cat4_revised2014.csv';
proc print data =rhizolsmeans_cat4; run ; ODS CSV close;
************************DEPTH CATEGORY 5 (CORRESPONDING TO 45-55 CM SOIL DEPTH********************;
data category5;
set alldepths_2011;
where category = 5;
run;
proc sort data =category5;
by block DOY co2 h2o ;
run ;
proc mixed data =category5 order =data covtest;
class block DOY co2 h2o ;
model cm_cm3 = co2|h2o|DOY/ ddfm =kr;
random block block*co2 ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = rhizopdiffs_cat5;
ods output lsmeans = rhizolsmeans_cat5;
run ;
proc print data =rhizopdiffs_cat5; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizopdiffs_cat5_revised2014.csv';
proc print data =rhizopdiffs_cat5; run; ODS CSV close;
proc print data =rhizolsmeans_cat5; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizolsmeans_cat5_revised2014.csv';
proc print data =rhizolsmeans_cat5; run ; ODS CSV close;
************************DEPTH CATEGORY 6 (CORRESPONDING TO 55-65 CM SOIL DEPTH********************;
data category6;
set alldepths_2011;
where category = 6;
run;
proc sort data =category6;
by block DOY co2 h2o ;
run ;
proc mixed data =category6 order =data covtest;
class block DOY co2 h2o ;
model cm_cm3 = co2|h2o|DOY/ ddfm =kr;
random block block*co2 ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = rhizopdiffs_cat6;
ods output lsmeans = rhizolsmeans_cat6;
run ;
proc print data =rhizopdiffs_cat6; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizopdiffs_cat6_revised2014.csv';
proc print data =rhizopdiffs_cat6; run; ODS CSV close;
proc print data =rhizolsmeans_cat6; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizolsmeans_cat6_revised2014.csv';
proc print data =rhizolsmeans_cat6; run ; ODS CSV close;
************************DEPTH CATEGORY 7 (CORRESPONDING TO 65-75 CM SOIL DEPTH********************;
data category7;
set alldepths_2011;
where category = 7;
run;
proc sort data =category7;
by block DOY co2 h2o ;
run ;
proc mixed data =category7 order =data covtest;
class block DOY co2 h2o ;
model cm_cm3 = co2|h2o|DOY/ ddfm =kr;
random block block*co2 ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = rhizopdiffs_cat7;
ods output lsmeans = rhizolsmeans_cat7;
run ;
proc print data =rhizopdiffs_cat7; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizopdiffs_cat7_revised2014.csv';
proc print data =rhizopdiffs_cat7; run; ODS CSV close;
proc print data =rhizolsmeans_cat7; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizolsmeans_cat7_revised2014.csv';
proc print data =rhizolsmeans_cat7; run ; ODS CSV close;
************************DEPTH CATEGORY 8 (CORRESPONDING TO 75-85 CM SOIL DEPTH********************;
data category8;
set alldepths_2011;
where category = 8;
run;
proc sort data =category8;
by block DOY co2 h2o ;
run ;
proc mixed data =category8 order =data covtest;
class block DOY co2 h2o ;
model cm_cm3 = co2|h2o|DOY/ ddfm =kr;
random block block*co2 ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = rhizopdiffs_cat8;
ods output lsmeans = rhizolsmeans_cat8;
run ;
proc print data =rhizopdiffs_cat8; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizopdiffs_cat8_revised2014.csv';
proc print data =rhizopdiffs_cat8; run; ODS CSV close;
proc print data =rhizolsmeans_cat8; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizolsmeans_cat8_revised2014.csv';
proc print data =rhizolsmeans_cat8; run ; ODS CSV close;
************************DEPTH CATEGORY 9 (CORRESPONDING TO 85-95 CM SOIL DEPTH********************;
data category9;
set alldepths_2011;
where category = 9;
run;
proc sort data =category9;
by block DOY co2 h2o ;
run ;
proc mixed data =category9 order =data covtest;
class block DOY co2 h2o ;
model cm_cm3 = co2|h2o|DOY/ ddfm =kr;
random block block*co2 ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = rhizopdiffs_cat9;
ods output lsmeans = rhizolsmeans_cat9;
run ;
proc print data =rhizopdiffs_cat9; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizopdiffs_cat9_revised2014.csv';
proc print data =rhizopdiffs_cat9; run; ODS CSV close;
proc print data =rhizolsmeans_cat9; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Rhizotron\Rhizotron Statistics\Root Length Density\rhizolsmeans_cat9_revised2014.csv';
proc print data =rhizolsmeans_cat9; run ; ODS CSV close;
