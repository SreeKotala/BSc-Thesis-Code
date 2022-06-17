# BSc-Thesis-
BSc Data Science and Artifical Intelligence Thesis supporting code.

The supporting code is written in MATLAB R2021a

Required Packages: 

map_toolbox
matlab
signal_toolbox
statistics_toolbox



Code Guide: 

Download all the files and the supporting GW event files in the 'gravitational waves' folder. 

Open MATLAB and run gwSSDPipeline to produce the pipeline output as detailed in the report. 

dataCollectiongwSSDPipeline will generate and record the performance metrics for all 48 GW events in the 'gravitational waves' folder and output them into the GWEventOutput.txt file. 

statisticalAnalysis conducts the t-tests used to determine the (lack there of) significance in the results from the dataCollectiongwSSDPipeline file. 



Supporting Documents:

The gweventinformation.txt file contains information downloaded from the GWOSC query tool (https://www.gw-openscience.org/eventapi/html/query/) and stores the GPSEventTime for each corresponding GW event. 

The file needs to be updated with the information in the same structure for newer events if the pipeline is to work with events not included within the 48 examples detailed in the report. 
