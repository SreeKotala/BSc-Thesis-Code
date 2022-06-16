%% (gw)SSD pipeline 
% Generates the Visulisations including the three time - frequency plots,
% and their respective tfSRC plots. Also plots the SSD and gwSSD components
% of the signal being decomposed in one figure to simplify interpretation.


%% Basic setup: importing dataset and the meta data e.g. event time etc...
clear 

[fileName,gwEvent] = uigetfile('*.hdf5', 'Please Select the Gravitation Wave Event');
gwEventPath = strcat(gwEvent,fileName);

detectorname = string(h5read(gwEventPath,'/meta/Detector'));
strain = h5read(gwEventPath,'/strain/Strain');

meanstrain = mean(strain);
%centralizing the stream about the mean to normalize the values such that
%they are plotted about the mean
strain = (strain-meanstrain)*10^(19); 

%setting the sampling frequency and getting the time from total time of the
%sample from it 
%fs = 163484; %Hz
tstart = double(h5read(gwEventPath,'/meta/GPSstart'));
duration = double(h5read(gwEventPath,'/meta/Duration'));
tend = tstart + duration;

fs = double(length(strain)/duration);
time = (1:length(strain))/fs;
time = time';

close all

time_H1 = tstart:1/fs:tend;
time_H1(end) = [];
time_H1 = time_H1-tstart;

%Extracting the 'normal name' of the GW signal from the meta data
GWUTC = convertCharsToStrings(h5read(gwEventPath,'/meta/UTCstart'));
GWyear = extractBetween(GWUTC, 3, 4);
GWmonth = extractBetween(GWUTC, 6, 7);
GWday = extractBetween(GWUTC, 9, 10);

GWEventName = strcat('GW',GWyear,GWmonth,GWday); 

%Retrieving the appropriate event time from the gwevent information file
GWEventCatalogue = readtable('gweventinformation.txt');
GWEventList = GWEventCatalogue{:,2};
GWEventGPSList = GWEventCatalogue{:,5};
GWEventChirpStart = GWEventCatalogue{:,7};
GWEventChirpEnd = GWEventCatalogue{:,8};
[GWrow, ~] = find(GWEventList == GWEventName);
 

%Here need to check if GWrow is empty, that means that the signal also has
%the time in the common name 
if isempty(GWrow)
    %as GWUTC is the start of the recording and event names are recorded at
    %the time of occurence we need to add half the length of the recording
    %and then check whether the seconds, minutes, or hours counters tick
    %over. This is setup only for recordings that are 32 seconds long. as
    %the UTC start is then 15 seconds before the 'eventtime' as indicated
    %in the 'common name' of the GW event in the documentation.
    
    %Note this also only does the incrememnting for the second and hour and
    %hour, assuming that the events tested won't occur right next to a day 
    %, month or year change
    incrementMin=0;
    incrementHour=0;
    GWsecond = str2double(extractBetween(GWUTC, 18, 19));
   
    if (GWsecond+15 >= 60)
       GWsecond = rem(GWsecond+15 , 60);
       incrementMin = 1;
    else
       GWsecond = rem(GWsecond+15 , 60); 
    end
    
    GWminute = str2double(extractBetween(GWUTC, 15, 16));
    if (GWminute+incrementMin >= 60 && incrementMin == 1)
       GWminute = rem(GWminute+incrementMin , 60);
       incrementHour = 1;
    else 
       GWminute = rem(GWminute+incrementMin , 60);
    end
    
    GWhour = str2double(extractBetween(GWUTC, 12, 13));
    if (GWhour+incrementHour >= 24 && incrementHour == 1)
       GWhour = rem(GWhour+incrementHour , 24);
       incrementday = 1; % incrementing beyond the day may be excessive
    else 
       GWhour = rem(GWhour+incrementHour , 24); 
    end
    if GWsecond <= 9
        GWsecond = strcat('0',string(GWsecond));
    end
    if GWminute <= 9
        GWminute = strcat('0',string(GWminute));
    end
    if GWhour <= 9
        GWhour = strcat('0',string(GWhour));
    end
    GWEventName = strcat(GWEventName,'_',string(GWhour),string(GWminute),string(GWsecond));
    [GWrow, ~] = find(GWEventList == GWEventName);
end

%Instance where there is a discrepancy in the UTC and event times
if isempty(GWrow)
    GWsecond = double(GWsecond);
    GWEventName = strcat('GW',GWyear,GWmonth,GWday,'_',string(GWhour),string(GWminute),string(GWsecond+1));
    [GWrow, ~] = find(GWEventList == GWEventName);
end
%In the instance that the event is still not found we search a second
%before the time at UTC in case the discrepancy is the other direction
if isempty(GWrow)
    GWsecond = double(GWsecond);
    GWEventName = strcat('GW',GWyear,GWmonth,GWday,'_',string(GWhour),string(GWminute),string(GWsecond-1));
    [GWrow, ~] = find(GWEventList == GWEventName);
end

eventtime = GWEventGPSList(GWrow);
%eventtime = 1243533585;
tevent = eventtime-tstart;

%% Pre-processing of data

deltat = 5; % seconds around the event

% index into the strain time series for this time interval:
indxt = find((time_H1 >= tevent-deltat) & (time_H1 < tevent+deltat));
indxt3 = find((time_H1 >= tevent-0.1) & (time_H1 < tevent+0.05));

GWsignal = strain(indxt);

%% Whittening the data % (taking care to get normalization right)
% compute PSD by means of Pwelch
nfft = length(GWsignal);

%Pxx is the power of the signal at a given freq. 
[Pxx,freqs] = pwelch(GWsignal,[],[],nfft,fs,'twosided');%length(data),fs);

% We will use interpolations of the ASDs computed above for whitening:
psd_noise = interp1(freqs, Pxx, freqs);

%whitten: %transform to freq domain,
FTgw1 = fft(GWsignal,nfft); 
%divide by asd,
white_gwf = FTgw1./(sqrt(psd_noise)/(time_H1(2)-time_H1(1))/2);
% and then transform back
strain_whiten = real(ifft(white_gwf,nfft));

%% Band pass using butterworth filter
%creating a 8th order butterworth filter with a min and max band of 20Hz
%and 400Hz respectively. 
[bb, ab] = butter(4, [20/(fs/2), 400/(fs/2)]);

strain_whitenbp = filtfilt(bb, ab, strain_whiten);


%% Extract the information that is within 0.2 seconds around the event for further evaluation/processing
tt = time_H1(indxt)-tevent;
lim = 0.2;
indxt2 = find(tt>=-lim & tt<=lim);

%% Plot the Spectrogram of the processed signal 
plotXLim = [-0.20, 0.20];
plotYLim = [50, 475];
%X(signal input), window(length), NOVERLAP(overlapping points between
%windows), F(two sided specto. at these normalised frequencies, Fs:sampling
[S,F,T,P] = spectrogram(strain_whitenbp(indxt2),blackman(fs/32),fs/32*31/32,fs/8,fs,'reassigned','yaxis');
tvec = 2*lim/(T(end)-T(1))*(T-T(1))-lim;

% figure 
% imagesc(tvec,F,abs(S))
% %Setting the default tex interpreter to none so that the names don't have a
% %subscript 0 due to the naming convention of GW wave events. 
% set(0,'DefaultTextInterpreter','none')
% axis xy
% xlim(plotXLim)
% ylim(plotYLim);
% xlabel('time (s)')
% ylabel('frequency (Hz)')
% title(strcat('Spectrogram of ', {' '}, GWEventName),'Interpreter','none')
% movegui('southwest');
% 
% %Plotting the overlay of the ridge onto the spectrogram 
[fridge, iridge, lridge] = tfridge(S,F);
% 
% hold on
% plot3(tvec,fridge,abs(S(lridge)),'LineWidth',4)
% xlim(plotXLim)
% ylim(plotYLim);
% hold off

%% Initiating the tiled layout of the plots
f = figure;
tiledlayout(2,4)
f.WindowState = 'maximized';

%% Decompose with  regular SSD
[SSDcomponents]=SSD(strain_whitenbp(indxt2),fs,0.01);
nexttile
stackedplot([strain_whitenbp(indxt2)';SSDcomponents]','Title','SSD decomposition')
title(strcat('SSD Components of', {' '}, GWEventName))


% Select SSD components for time-frequency representation and generate Hilbert Spectrum
t = (0:length(SSDcomponents)-1)/fs;
f = (0:length(strain_whitenbp(indxt2))-1)*fs/length(strain_whitenbp(indxt2));

Xs2 = fft(SSDcomponents');
[~,max_index2] = max(abs(Xs2));
index_ssd_af = find(f(max_index2)>=20 & f(max_index2)<=400);

% Plotting the Hilbert Spectrum of the signal decomposed using regular SSD 
%index_ssd_af = [3,2,4];

%Summing presents an issue if there is only 1 components and thus we need
%to check if there is only 1 component. 
if size(SSDcomponents,1) == 1 
    [HSAA,HSA,HStt,HSf,HSff] = buildHS(SSDcomponents(index_ssd_af,:),fs,t,20,400);
else
    [HSAA,HSA,HStt,HSf,HSff] = buildHS(sum(SSDcomponents(index_ssd_af,:)),fs,t,20,400);
end
% figure
% %imagesc(t,HSff,HSAA)
tvec2 = 2*lim/(HStt(end)-HStt(1))*(HStt-HStt(1))-lim;
% imagesc(tvec2,HSff,HSAA)
% axis xy
% ylim(plotYLim);
% xlim(plotXLim)
% xlabel('time (s)')
% ylabel('frequency (Hz)')
% title(strcat('Hilbert Spectrum of', {' '}, GWEventName, 's SSD components between 20 and 400Hz'),'Interpreter','none')
% movegui('south');
% 
% %Plotting the overlay of the ridge onto the Hilbert Spectrum 
[fridge2, iridge2, lridge2] = tfridge(HSAA, HSff);
% hold on
% plot3(tvec2,fridge2,abs(HSAA(iridge2)),'LineWidth',0.5)
% xlim(plotXLim)
% ylim(plotYLim);
% hold off


%% Decompose with  gwSSD
[SSDcomponents2]=gwSSD(strain_whitenbp(indxt2),fs,0.01);
% nexttile
% stackedplot([strain_whitenbp(indxt2)';SSDcomponents2]','Title','gwSSD decomposition')
% title(strcat('gwSSD Components of', {' '}, GWEventName))


% Select SSD components for time-frequency representation and generate Hilbert Spectrum
t = (0:length(SSDcomponents2)-1)/fs;
f = (0:length(strain_whitenbp(indxt2))-1)*fs/length(strain_whitenbp(indxt2));

Xs3 = fft(SSDcomponents2');
[~,max_index3] = max(abs(Xs3));
index_ssd_af2 = find(f(max_index3)>=20 & f(max_index3)<=400);

% Plotting the Hilbert Spectrum of the signal decomposed using regular SSD 

%index_ssd_af = [3,2,4];

if size(SSDcomponents2,1) == 1 
    [HSAA2,HSA2,HStt2,HSf2,HSff2] = buildHS(SSDcomponents2(index_ssd_af2,:),fs,t,20,400);
else
    [HSAA2,HSA2,HStt2,HSf2,HSff2] = buildHS(sum(SSDcomponents2(index_ssd_af2,:)),fs,t,20,400);
end

%imagesc(t,HSff,HSAA)
tvec3 = 2*lim/(HStt2(end)-HStt2(1))*(HStt2-HStt2(1))-lim;
% figure
% imagesc(tvec3,HSff2,HSAA2)
% axis xy
% ylim(plotYLim);
% xlim(plotXLim)
% xlabel('time (s)')
% ylabel('frequency (Hz)')
% title(strcat('Hilbert Spectrum of', {' '}, GWEventName, 's gwSSD components between 20 and 400Hz'),'Interpreter','none')
% movegui('southeast');
% 
% %Plotting the overlay of the ridge onto the Hilbert Spectrum 
[fridge3, iridge3, lridge3] = tfridge(HSAA2, HSff2);
% hold on
% plot3(tvec3,fridge3,abs(HSAA2(iridge3)),'LineWidth',0.5)
% xlim(plotXLim)
% ylim(plotYLim);
% hold off


%% Determining the width of the peak of the HS+gwSSD output to process the 'chirp' window
%getting the tfSRC -0.20,0.20 seconds around the event time to get the
%location of the peak (chirp) with a frameband of 20Hz
%frameband = 2.503667481662592; %Hz as the width of the HS output is this, to make the comparison between Spectrogram and HS outputs fair. 
frameband = 1; %hz
[tfSRCChirpFinder,tvecPlotChirpFinder,~] = buildtfSRC(HSAA2, tvec3, fridge3, iridge3, frameband, fs, [-0.20,0.20]);

%Based on the intial chirp finding we determine the time of the event and
%thus the limits to consider when building the tfSRC for the Spectrogram,
%HS+SSD and HS+gwSS outputs. 
[~, locsChirp,wChirp,~] = findpeaks(tfSRCChirpFinder,tvecPlotChirpFinder, 'SortStr', 'descend');
relativeChirpTime = locsChirp(1);
chirpWidth = wChirp(1);
chirpLimits = [relativeChirpTime - chirpWidth, relativeChirpTime + chirpWidth];
%% Plotting the Spectrogram
plotXLim = [chirpLimits(1), chirpLimits(2)];

nexttile 
imagesc(tvec,F,abs(S))
%Setting the default tex interpreter to none so that the names don't have a
%subscript 0 due to the naming convention of GW wave events. 
set(0,'DefaultTextInterpreter','none')
axis xy
xlim(plotXLim)
ylim(plotYLim);
xlabel('time (s)')
ylabel('frequency (Hz)')
title(strcat('Spectrogram of ', {' '}, GWEventName),'Interpreter','none')


%Plotting the overlay of the ridge onto the spectrogram 
%[fridge, iridge, lridge] = tfridge(S,F);

hold on
plot3(tvec,fridge,abs(S(lridge)),'LineWidth',4)
xlim(plotXLim)
ylim(plotYLim);
hold off

%% Plotting the HS+SSD output

nexttile
%imagesc(t,HSff,HSAA)
%tvec2 = 2*lim/(HStt(end)-HStt(1))*(HStt-HStt(1))-lim;
imagesc(tvec2,HSff,HSAA)
axis xy
ylim(plotYLim);
xlim(plotXLim)
xlabel('time (s)')
ylabel('frequency (Hz)')
title(strcat('HS of', {' '}, GWEventName, 's selected SSD components'),'Interpreter','none')


%Plotting the overlay of the ridge onto the Hilbert Spectrum 
%[fridge2, iridge2, lridge2] = tfridge(HSAA, HSff);
hold on
plot3(tvec2,fridge2,abs(HSAA(iridge2)),'LineWidth',0.5)
xlim(plotXLim)
ylim(plotYLim);
hold off

%% Plotting the HS + gwSSD
nexttile
imagesc(tvec3,HSff2,HSAA2)
axis xy
ylim(plotYLim);
xlim(plotXLim)
xlabel('time (s)')
ylabel('frequency (Hz)')
title(strcat('HS of', {' '}, GWEventName, 's selected gwSSD components'),'Interpreter','none')


%Plotting the overlay of the ridge onto the Hilbert Spectrum 
%[fridge3, iridge3, lridge3] = tfridge(HSAA2, HSff2);
hold on
plot3(tvec3,fridge3,abs(HSAA2(iridge3)),'LineWidth',0.5)
xlim(plotXLim)
ylim(plotYLim);
hold off

%% Plotting the gwSSD Components of the signal due to the tiled layout
nexttile
stackedplot([strain_whitenbp(indxt2)';SSDcomponents2]','Title','gwSSD decomposition')
title(strcat('gwSSD Components of', {' '}, GWEventName))


%% Determining the tf-Ridge Spectral Concentration of the HS+gwSSD plot. 

frameband2 = 0; %Setting to zero ignores frame bands and only considers the
%energy within the window containing the 'loudest' frequency at each
%instance
[tfSRC3,tvecPlot3,specConcHSgw, absEHSgw] = buildtfSRC(HSAA2, tvec3, fridge3, iridge3, frameband, fs, chirpLimits);
[~, locs3,w3,p3] = findpeaks(tfSRC3,tvecPlot3,'Annotate','extents', 'SortStr', 'descend');
plotPeakIndex = find(tvecPlot3 == locs3(1), 1);
tfSRCYLim = tfSRC3(plotPeakIndex);

%% Determining the tf-Ridge Spectral Concentration of the HS+SSD plot. 

frameband2 = 0; %Setting to zero ignores frame bands and only considers the
%energy within the window containing the 'loudest' frequency at each
%instance
[tfSRC2,tvecPlot2,specConcHS, absEHS] = buildtfSRC(HSAA, tvec2, fridge2, iridge2, frameband, fs, chirpLimits);
[~, locs2,w2,~] = findpeaks(tfSRC2,tvecPlot2,'Annotate','extents', 'SortStr', 'descend');

%% Determining the tf-SRC of the Spectrogram Output

[tfSRC,tvecPlot,specConcSpectro,absESpectro] = buildtfSRC(S, tvec, fridge, iridge, frameband, fs, chirpLimits);
%chirpLimits = [-0.25, 0.25];

%% Plotting the tf-Ridge Spectral Concentration of the Spectrogram about the chirp 
nexttile
%[~, locs,w,~] = findpeaks(tfSRC,tvecPlot,'Annotate','extents', 'SortStr', 'descend');
findpeaks(tfSRC,tvecPlot,'Annotate','extents', 'SortStr', 'descend')
ylim([-0.0001,tfSRCYLim*1.1])
title(strcat('tfSRC progression of', {' '}, GWEventName, 's spectrogram'),'Interpreter','none')


%% Plotting the tf-Ridge Spectral Concentration of the HS+SSD output about the chirp 
nexttile
findpeaks(tfSRC2,tvecPlot2,'Annotate','extents', 'SortStr', 'descend')
ylim([-0.0001, tfSRCYLim*1.1])
title(strcat('tfSRC progression of', {' '}, GWEventName, 's Hilbert Spectrum (SSD)'),'Interpreter','none')

%% Determining and plotting the tfRidge Concentration of the gwSSD HS plots. 
frameband2 = 0; %Setting to zero ignores frame bands and only considers the
%energy within the window containing the 'loudest' frequency at each
%instance
%[tfSRC3,tvecPlot3,specConcHSgw, absEHSgw] = buildtfSRC(HSAA2, tvec3, fridge3, iridge3, frameband, fs, chirpLimits);

nexttile
findpeaks(tfSRC3,tvecPlot3,'Annotate','extents', 'SortStr', 'descend')
ylim([-0.0001, tfSRCYLim*1.1])
title(strcat('tfSRC progression of', {' '}, GWEventName, 's Hilbert Spectrum (gwSSD)'),'Interpreter','none')


 %% Determining the dominant frequencies in each SSD Component to identify during selection
[numComponents,complength] = size(SSDcomponents);
DominantFrequency = zeros(numComponents,1);

for i = 1:numComponents
    frequency_vector = (0:floor(length(SSDcomponents(i,:))/2)-1)*(fs)/(length(SSDcomponents(i,:)));
    Spectrumcomponent = abs(fft(SSDcomponents(i,:)));
    Spectrumcomponent(floor(length(Spectrumcomponent)/2)+1:end) = [];
    [~, PositionMax] = max(Spectrumcomponent);
    DominantFrequency(i,1) = frequency_vector(PositionMax);
    
end 

%% Displaying the Spectral Concentration of the three outputs into the print file

disp("The Spectral Concentrations of the three outputs");
disp("Spectrogram: " + specConcSpectro);
disp("HS + SSD: " + specConcHS);
disp("HS + gwSSD: " + specConcHSgw);

disp("The Absolute Energies of the three outputs");
disp("Spectrogram: " + absESpectro);
disp("HS + SSD: " + absEHS);
disp("HS + gwSSD: " + absEHSgw);
