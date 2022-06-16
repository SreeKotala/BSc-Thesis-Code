%% Data Collection file for the gwSSD pipeline propsoed 
% Processes all the hdf5 files of GW events containing chirps from BBH
% mergers located in the 'gravitational waves' folder. Collects some
% information from the meta data within the file and partially from the
% gweventinformation.txt file. 

% Outputs the GWEventOutput.txt file which includes the six performance
% metrics described below. 

%%
clear

GWSignalFiles = dir('.\gravitational waves\*.hdf5');
%Defining the dimensions and headings for the table
varNames = ["GWEvent","SC Spec","SC HS+SSD", "SC HS+gwSSD", "AbsE Spec", "AbsE HS+SSD", "AbsE HS+gwSSD"];
varTypes = ["string","double","double","double","double","double","double"];
sz = [0, length(varNames)];

%Creating the table based upon the heading and datatypes above. 
GWEventResults = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

for gwEventFile = GWSignalFiles'
    
    gwEvent = string(extractfield(gwEventFile, 'folder'));
    gwEventPath = strcat(gwEvent, '//', string(extractfield(gwEventFile, 'name')));
     
    %Need to setup a loop to get files from the location

    %[fileName,gwEvent] = uigetfile('*.hdf5', 'Please Select the Gravitation Wave Event');
    %gwEventPath = strcat(gwEvent,fileName);

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
    %eventtime = 1126259462.4;
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
    %creating a 4th order butterworth filter with a min and max band of 20Hz
    %and 300Hz respectively. 
    [bb, ab] = butter(4, [20/(fs/2), 400/(fs/2)]);

    strain_whitenbp = filtfilt(bb, ab, strain_whiten);


    %% Extract the information that is within 0.2 seconds around the event for further evaluation/processing
    tt = time_H1(indxt)-tevent;
    lim = 0.2;
    indxt2 = find(tt>=-lim & tt<=lim);

    %% Generate the Spectrogram of the processed signal 
    plotXLim = [-0.05, 0.05];
    plotYLim = [50, 475];
    %X(signal input), window(length), NOVERLAP(overlapping points between
    %windows), F(two sided specto. at these normalised frequencies, Fs:sampling
    [S,F,T,P] = spectrogram(strain_whitenbp(indxt2),blackman(fs/32),fs/32*31/32,fs/8,fs,'reassigned','yaxis');
    tvec = 2*lim/(T(end)-T(1))*(T-T(1))-lim;
    %imagesc(tvec,F,abs(S))

    %Setting the default tex interpreter to none so that the names don't have a
    %subscript 0 due to the naming convention of GW wave events. 
    set(0,'DefaultTextInterpreter','none')

    %Plotting the overlay of the ridge onto the spectrogram 
    [fridge, iridge, lridge] = tfridge(S,F);

    %% Decompose with  regular SSD
    [SSDcomponents]=SSD(strain_whitenbp(indxt2),fs,0.01);
    
    % Select SSD components for time-frequency representation and generate Hilbert Spectrum
    t = (0:length(SSDcomponents)-1)/fs;
    f = (0:length(strain_whitenbp(indxt2))-1)*fs/length(strain_whitenbp(indxt2));

    Xs2 = fft(SSDcomponents');
    [~,max_index2] = max(abs(Xs2));
    index_ssd_af = find(f(max_index2)>=20 & f(max_index2)<=400);

    % Plotting the Hilbert Spectrum of the signal decomposed using regular SSD 

    %index_ssd_af = [3,2,4];
    % figure
    if size(SSDcomponents,1) == 1
        [HSAA,HSA,HStt,HSf,HSff] = buildHS(SSDcomponents(index_ssd_af,:),fs,t,20,400);
    else
        [HSAA,HSA,HStt,HSf,HSff] = buildHS(sum(SSDcomponents(index_ssd_af,:)),fs,t,20,400);
    end
    %imagesc(t,HSff,HSAA)
    tvec2 = 2*lim/(HStt(end)-HStt(1))*(HStt-HStt(1))-lim;

    %Generating the ridge of the Hilbert Spectrum 
    [fridge2, iridge2, lridge2] = tfridge(HSAA, HSff);


    %% Decompose with  gwSSD
    [SSDcomponents2]=gwSSD(strain_whitenbp(indxt2),fs,0.01);
    % figure
    % stackedplot([strain_whitenbp(indxt2)';SSDcomponents2]','Title','gwSSD decomposition')
    % title(strcat('gwSSD Components of', {' '}, GWEventName))
    % movegui('east');


    % Select SSD components for time-frequency representation and generate Hilbert Spectrum
    t = (0:length(SSDcomponents2)-1)/fs;
    f = (0:length(strain_whitenbp(indxt2))-1)*fs/length(strain_whitenbp(indxt2));

    Xs3 = fft(SSDcomponents2');
    [~,max_index3] = max(abs(Xs3));
    index_ssd_af2 = find(f(max_index3)>=20 & f(max_index3)<=400);

    %index_ssd_af = [3,2,4];
    %figure
    if size(SSDcomponents2,1) == 1
        [HSAA2,HSA2,HStt2,HSf2,HSff2] = buildHS(SSDcomponents2(index_ssd_af2,:),fs,t,20,400);
    else
        [HSAA2,HSA2,HStt2,HSf2,HSff2] = buildHS(sum(SSDcomponents2(index_ssd_af2,:)),fs,t,20,400);
    end

    tvec3 = 2*lim/(HStt2(end)-HStt2(1))*(HStt2-HStt2(1))-lim;


    %Generating the ridge onto the Hilbert Spectrum 
    [fridge3, iridge3, lridge3] = tfridge(HSAA2, HSff2);


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
    
    
    %% Determining the tf-Ridge Spectral Concentration of the HS+gwSSD plot.
    
    frameband2 = 0; %Setting to zero ignores frame bands and only considers the
    %energy within the window containing the 'loudest' frequency at each
    %instance
    [tfSRC3,tvecPlot3,specConcHSgw, absEHSgw] = buildtfSRC(HSAA2, tvec3, fridge3, iridge3, frameband, fs, chirpLimits);
    [~, locs3,w3,p3] = findpeaks(tfSRC3,tvecPlot3,'Annotate','extents', 'SortStr', 'descend');
    plotPeakIndex = find(tvecPlot3 == locs3(1), 1);
    tfSRCYLim = tfSRC3(plotPeakIndex);
    
    %% Determining the tf-Ridge Spectral Concentration of the HS+SSD plot
    [tfSRC2,tvecPlot2,specConcHS, absEHS] = buildtfSRC(HSAA, tvec2, fridge2, iridge2, frameband, fs, chirpLimits);
    [~, locs2,w2,~] = findpeaks(tfSRC2,tvecPlot2,'Annotate','extents', 'SortStr', 'descend');
    
    %% Determining the tf-SRC of the Spectrogram Output
    
    [tfSRC,tvecPlot,specConcSpectro,absESpectro] = buildtfSRC(S, tvec, fridge, iridge, frameband, fs, chirpLimits);
    %chirpLimits = [-0.25, 0.25];
    

    %% Determining and plotting the tfRidge Concentration of the gwSSD HS plots.
    frameband2 = 0; %Setting to zero ignores frame bands and only considers the
    %energy within the window containing the 'loudest' frequency at each
    %instance
    %[tfSRC3,tvecPlot3,specConcHSgw, absEHSgw] = buildtfSRC(HSAA2, tvec3, fridge3, iridge3, frameband, fs, chirpLimits);
    %% Storing the outputs of each signal into a table such that it can be written later. 
    eventinfo = {GWEventName,specConcSpectro, specConcHS, specConcHSgw, absESpectro, absEHS, absEHSgw};
    GWEventResults = [GWEventResults;eventinfo];
    disp("Finished Wave: " + GWEventName);
    
end

%% Outputting the Spectral Concentration of the three outputs into the print file

% fileID = fopen('output2.txt','w');
% 
% fprintf(fileID,'The Spectral Concentrations of the three outputs \r\n');
% fprintf(fileID,'Spectrogram: %f \r\n', specConcSpectro);
% fprintf(fileID,'HS + SSD: %f \r\n', specConcHS);
% fprintf(fileID,'HS + gwSSD: %f \r\n' , specConcHSgw);
% 
% fprintf(fileID,'The Absolute Energies of the three outputs \r\n');
% fprintf(fileID,'Spectrogram: %f \r\n', absESpectro);
% fprintf(fileID,'HS + SSD: %f \r\n', absEHS);
% fprintf(fileID,'HS + gwSSD: %f \r\n' , absEHSgw);
% fclose(fileID);

%Use the writetable function to write the table of all the results out to 
%the GWEventoutput.txt file. 

writetable(GWEventResults, 'GWEventoutput.txt');
