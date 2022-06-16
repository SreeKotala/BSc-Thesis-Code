%% Statistical Signficance testing of the outputs of the gwSSDPipeline data collection process 

%% Load in the data from the GWEventoutput.txt file
GWEventOutput = readtable('GWEventoutput.txt');

%% Seperate the data into the 7 columns for the name and the various
%concentrations and absolute energies. 

GWNames = GWEventOutput.GWEvent;

SpectralConcSpectrogram = GWEventOutput.SCSpec;
SpectralConcHS_SSD = GWEventOutput.SCHS_SSD;
SpectralConcHS_gwSSD = GWEventOutput.SCHS_gwSSD;

AbsEnergySpectrogram = GWEventOutput.AbsESpec;
AbsEnergyHS_SSD = GWEventOutput.AbsEHS_SSD;
AbsEnergyHS_gwSSD = GWEventOutput.AbsEHS_gwSSD;

%% Setup the first run of experiments with the spectral concentration
%comparing the three methods against each other. 
% paired t-test and the spectral concentration of the three tuples to 
% compare the methods as such

%Spec vs HS_SSD

[h1,p1,ci1,stats1] = ttest(SpectralConcSpectrogram,SpectralConcHS_SSD, 'Alpha', 0.05);%,'Tail','right');

%Spec vs HS_gwSSD

[h2,p2,ci2,stats2] = ttest(SpectralConcSpectrogram,SpectralConcHS_gwSSD, 'Alpha', 0.05);%,'Tail','right');

%Hs_SSD vs HS_gwSSD

[h3,p3,ci3,stats3] = ttest(SpectralConcHS_SSD,SpectralConcHS_gwSSD, 'Alpha', 0.05);%,'Tail','right');


%% Setup the second run of experiments. 
% Look at whether there is a statistical significance between the absolute
% energies. 
% Test significance between the absolute energies of the Spec, HS_SSD and 
% HS_gwSSD outputs. 

%Spec vs HS_SSD
[h4,p4,ci4,stats4] = ttest(AbsEnergySpectrogram,AbsEnergyHS_SSD);%,'Tail','right');

[h5,p5,c5,stats5] = ttest(AbsEnergyHS_SSD,AbsEnergyHS_gwSSD);%,'Tail','right');

%Spec vs HS_gwSSD
[h6,p6,ci6,stats6] = ttest(AbsEnergySpectrogram,AbsEnergyHS_gwSSD);%,'Tail','right');

%HS_SSD vs HS_gwSSD










