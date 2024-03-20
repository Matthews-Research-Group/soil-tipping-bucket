options ls=100 ps=8000 pageno=1;
title 'midday_gas_ex_2004thru2011';
Data ringmeans_2004thru2011;	
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\ExtendedDataFig4_Table_S4_gs_ring_means.csv'
dlm=',' firstobs=2 lrecl=10000;
input year DOY block CO2 Ring Cond @@;
run;

data ringmeans_midday2004;
set ringmeans_2004thru2011;
where year=2004;
run;
data ringmeans_midday2005;
set ringmeans_2004thru2011;
where year=2005;
run;
data ringmeans_midday2006;
set ringmeans_2004thru2011;
where year=2006;
run;
data ringmeans_midday2007;
set ringmeans_2004thru2011;
where year=2007;
run;
data ringmeans_midday2008;
set ringmeans_2004thru2011;
where year=2008;
run;
data ringmeans_midday2009;
set ringmeans_2004thru2011;
where year=2009;
run;
data ringmeans_midday2010;
set ringmeans_2004thru2011;
where year=2010;
run;
data ringmeans_midday2011;
set ringmeans_2004thru2011;
where year=2011;
run;

***THE FOLLOWING CODE IS FOR 2004;;;;;;;;;;;;
data ringmeans_midday2004;
set ringmeans_midday2004;
run;
proc mixed data =ringmeans_midday2004 order =data covtest;
class DOY block co2 ;
model Cond = co2|DOY / ddfm =kr;
random block  ;
repeated DOY / type =ar(1) subject =block*co2;
lsmeans DOY|co2/ pdiff;  
ods output diffs = Cond2004pdiffs;
ods output lsmeans= Cond2004lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2004lsmeans.csv';
proc print data =Cond2004lsmeans;  run;ODS CSV close;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2004pdiffs.csv';
proc print data = Cond2004pdiffs;  where DOY= _DOY; run;ODS CSV close;

***THE FOLLOWING CODE IS FOR 2005;;;;;;;;;;;;
proc mixed data =ringmeans_midday2005 order =data covtest;
class DOY block co2 ;
model Cond = co2|DOY / ddfm =kr;
random block  ;
repeated DOY / type =ar(1) subject =block*co2;
lsmeans DOY|co2/ pdiff;  
ods output diffs = Cond2005pdiffs;
ods output lsmeans= Cond2005lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2005lsmeans.csv';
proc print data =Cond2005lsmeans;  run;ODS CSV close;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2005pdiffs.csv';
proc print data = Cond2005pdiffs;  where DOY= _DOY; run;ODS CSV close;

***THE FOLLOWING CODE IS FOR 2006;;;;;;;;;;;;
proc mixed data =ringmeans_midday2006 order =data covtest;
class DOY block co2 ;
model Cond = co2|DOY / ddfm =kr;
random block  ;
repeated DOY / type =ar(1) subject =block*co2;
lsmeans DOY|co2/ pdiff;  
ods output diffs = Cond2006pdiffs;
ods output lsmeans= Cond2006lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2006lsmeans.csv';
proc print data =Cond2006lsmeans;  run;ODS CSV close;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2006pdiffs.csv';
proc print data = Cond2006pdiffs;  where DOY= _DOY; run;ODS CSV close;

***THE FOLLOWING CODE IS FOR 2007;;;;;;;;;;;;
proc mixed data =ringmeans_midday2007 order =data covtest;
class DOY block co2 ;
model Cond = co2|DOY / ddfm =kr;
random block  ;
repeated DOY / type =ar(1) subject =block*co2;
lsmeans DOY|co2/ pdiff;  
ods output diffs = Cond2007pdiffs;
ods output lsmeans= Cond2007lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2007lsmeans.csv';
proc print data =Cond2007lsmeans;  run;ODS CSV close;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2007pdiffs.csv';
proc print data = Cond2007pdiffs;  where DOY= _DOY; run;ODS CSV close;

***THE FOLLOWING CODE IS FOR 2008;;;;;;;;;;;;
proc mixed data =ringmeans_midday2008 order =data covtest;
class DOY block co2 ;
model Cond = co2|DOY / ddfm =kr;
random block  ;
repeated DOY / type =ar(1) subject =block*co2;
lsmeans DOY|co2/ pdiff;  
ods output diffs = Cond2008pdiffs;
ods output lsmeans= Cond2008lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2008lsmeans.csv';
proc print data =Cond2008lsmeans;  run;ODS CSV close;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2008pdiffs.csv';
proc print data = Cond2008pdiffs;  where DOY= _DOY; run;ODS CSV close;

***THE FOLLOWING CODE IS FOR 2009;;;;;;;;;;;;
proc mixed data =ringmeans_midday2009 order =data covtest;
class DOY block co2 ;
model Cond = co2|DOY / ddfm =kr;
random block  ;
repeated DOY / type =ar(1) subject =block*co2;
lsmeans DOY|co2/ pdiff;  
ods output diffs = Cond2009pdiffs;
ods output lsmeans= Cond2009lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2009lsmeans.csv';
proc print data =Cond2009lsmeans;  run;ODS CSV close;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2009pdiffs.csv';
proc print data = Cond2009pdiffs;  where DOY= _DOY; run;ODS CSV close;

***THE FOLLOWING CODE IS FOR 2010;;;;;;;;;;;;
proc mixed data =ringmeans_midday2010 order =data covtest;
class DOY block co2 ;
model Cond = co2|DOY / ddfm =kr;
random block  ;
repeated DOY / type =ar(1) subject =block*co2;
lsmeans DOY|co2/ pdiff;  
ods output diffs = Cond2010pdiffs;
ods output lsmeans= Cond2010lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2010lsmeans.csv';
proc print data =Cond2010lsmeans;  run;ODS CSV close;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2010pdiffs.csv';
proc print data = Cond2010pdiffs;  where DOY= _DOY; run;ODS CSV close;

***THE FOLLOWING CODE IS FOR 2011;;;;;;;;;;;;
proc mixed data =ringmeans_midday2011 order =data covtest;
class DOY block co2 ;
model Cond = co2|DOY / ddfm =kr;
random block  ;
repeated DOY / type =ar(1) subject =block*co2;
lsmeans DOY|co2/ pdiff;  
ods output diffs = Cond2011pdiffs;
ods output lsmeans= Cond2011lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2011lsmeans.csv';
proc print data =Cond2011lsmeans;  run;ODS CSV close;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Gas Exchange Data\SAS Library\Cond2011pdiffs.csv';
proc print data = Cond2011pdiffs;  where DOY= _DOY; run;ODS CSV close;


