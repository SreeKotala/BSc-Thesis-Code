%% Signal processing of Gravitational Wave recordings from LIGO

% recreate partially tutorial on GVs signal processing available at:

% https://losc.ligo.org/s/events/GW150914/GW150914_tutorial.html


clear

%

% load('H-H1_LOSC_4_V2-1126257414-4096.mat')

% plot(VarName1)


% analyse 32 sec. recording, sampled at 4096Hz, version 1

data1 = h5read('H-H1_LOSC_4_V1-1126259446-32.hdf5','/strain/Strain');

data2 = h5read('L-L1_LOSC_4_V1-1126259446-32.hdf5','/strain/Strain');

fs = 4096;

tstart = 1126259446; % GPS time start recording

tend = 1126259478; % GPS time end recording


%% V1 different from V2 - why? - ask Gideon

% load('H-H1_LOSC_4_V2-1126259446-32.mat')

% plot(data)

% % hold on

% % plot(VarName1,'r')

%% start with pre-processing

% select 10 sec. window centered around the event

close all

time_H1 = tstart:1/fs:tend;

time_H1(end) = [];

time_H1 = time_H1-tstart;

tevent = 1126259462.422-tstart; % Mon Sep 14 09:50:45 GMT 2015

% tevent = 1126259462.435-tstart;%1126259462.44

deltat = 5; % seconds around the event

% index into the strain time series for this time interval:

indxt = find((time_H1 >= tevent-deltat) & (time_H1 < tevent+deltat));

indxt3 = find((time_H1 >= tevent-0.1) & (time_H1 < tevent+0.05));


% plot(time_H1(indxt)-tevent,data(indxt))

% hold on

% plot([time_H1(indxt3(1)) time_H1(indxt3(1))]-tevent,[min(data) max(data)],'r')

% hold on

% plot([time_H1(indxt3(end)) time_H1(indxt3(end))]-tevent,[min(data) max(data)],'r')

% xlabel('time (s) +/-5 sec. around event')

% ylabel('strain')

% title('Advanced LIGO strain data near GW150914')

%

GVsignal1 = data1(indxt);

GVsignal2 = data2(indxt);


% figure,

% tt = time_H1(indxt)-tevent;

% indxt2 = find(tt>=-0.1 & tt<=0.05);

% plot(tt(indxt2),GVsignal(indxt2))


%% whitening

% compute PSD by means of Pwelch

nfft = length(GVsignal1);

[Pxx1,freqs1] = pwelch(GVsignal1,[],[],nfft,fs,'twosided');%length(data),fs);

[Pxx2,freqs2] = pwelch(GVsignal2,[],[],nfft,fs,'twosided');%length(data),fs);


% plot(freqs,10*log10(Pxx))

% semilogx(freqs(1:nfft/2),10*log10(Pxx(1:nfft/2)))

% xlabel('freqs. (Hz)')

% ylabel('dBs')

% indfmin = find(freqs>=10,1,'first');

% indfmax = find(freqs<=2000,1,'last');

% plot(freqs(indfmin:indfmax),Pxx(indfmin:indfmax))


% The sample rate is fs = 4096 Hz (2^12 Hz), so the data cannot capture

% frequency content above the Nyquist frequency = fs/2 = 2048 Hz.

% That's OK, because GW150914 only has detectable frequency content in the range 20 Hz - 300 Hz.


% You can see strong spectral lines in the data; they are all of instrumental

% origin. Some are engineered into the detectors (mirror suspension resonances

% at ~500 Hz and harmonics, calibration lines, control dither lines, etc)

% and some (60 Hz and harmonics) are unwanted. We'll return to these, later.


% You can't see the signal in this plot, since it is relatively weak and less

% than a second long, while this plot averages over 32 seconds of data.

% So this plot is entirely dominated by instrumental noise.


%% whiten the data in the Fourier domain

% We will use interpolations of the ASDs computed above for whitening:

psd_noise1 = interp1(freqs1, Pxx1, freqs1);

psd_noise2 = interp1(freqs2, Pxx2, freqs2);


% whitening: transform to freq domain, divide by asd, then transform back,

% (taking care to get normalization right)

FTgv1 = fft(GVsignal1,nfft);

FTgv2 = fft(GVsignal2,nfft);

% white_gvf = FTgv ./  (sqrt(psd_noise)/(freqs(2)-freqs(1))/2); %(np.sqrt(interp_psd(freqs) /dt/2.))

white_gvf1 = FTgv1 ./  (sqrt(psd_noise1)/(time_H1(2)-time_H1(1))/2);

white_gvf2 = FTgv2 ./  (sqrt(psd_noise2)/(time_H1(2)-time_H1(1))/2);


% a = abs(psd_noise)'*abs(FTgv)/(abs(psd_noise)'*abs(psd_noise));

% white_gvf = FTgv./(a*psd_noise);


strain_H1_whiten = real(ifft(white_gvf1,nfft));

strain_L1_whiten = real(ifft(white_gvf2,nfft));


%% band-passing + notching

% To get rid of remaining high frequency noise, we will also bandpass the

% data (see bandpassing, below).


 

[bb, ab] = butter(4, [20/(fs/2), 300/(fs/2)]); % butter(4, [43/(fs/2), 260/(fs/2)]);

strain_H1_whitenbp1 = filtfilt(bb, ab, strain_H1_whiten);

strain_L1_whitenbp1 = filtfilt(bb, ab, strain_L1_whiten);


% strain_H1_whitenbp = strain_H1_whitenbp1;

%

% notchesAbsolute = [14.0,34.70, 35.30, 35.90, 36.70, 37.30, 40.95, 60.00, ...

%          120.00, 179.99, 304.99, 331.49, 510.02, 1009.99];

%

% for indfreq = 1:length(notchesAbsolute)

%     clear B A

%     omega = 2*pi*notchesAbsolute(indfreq)/fs;

%     z1 = exp(1i*omega);

%     z2 = exp(-1i*omega);

%     p1 = 0.999*exp(1i*omega);

%     p2 = 0.999*exp(-1i*omega);

%     B = [1 -(z1+z2) z1*z2];

%     A = [1 -(p1+p2) p1*p2];

% %     freqz(B,A,[],fs)

% %     pause

% %     close

%     strain_H1_whitenbp = filtfilt(B, A, strain_H1_whitenbp);

%

% end


% plot the data after whitening:

% plot(time_H1(indxt)-tevent,data(indxt))

% hold on

subplot(221)

plot(time_H1(indxt)-tevent,strain_H1_whitenbp1)

xlim([-0.1,0.05])

xlabel('time (s) around event')

ylabel('whitented strain')

title('Hanford strain data near GW150914')


subplot(222)

plot(time_H1(indxt)-tevent,strain_L1_whitenbp1)

xlim([-0.1,0.05])

% ylim([-1.2,1.2])


% plt.xlim([-0.1,0.05])

% plt.ylim([-4,4])

xlabel('time (s) around event')

ylabel('whitented strain')

title('Livingston strain data near GW150914')


c = gray;

c = flipud(c);


subplot(223)

tt = time_H1(indxt)-tevent;

lim = 1;

indxt2 = find(tt>=-lim & tt<=lim);

[S,F,T] = spectrogram(strain_H1_whitenbp1(indxt2),blackman(fs/16),fs/16*15/16,fs/16,fs,'reassigned','yaxis');

tvec = 2*lim/(T(end)-T(1))*(T-T(1))-lim;

% S(abs(S)<=0.10) = 0;

imagesc(tvec,F,abs(S))


% colorbar

% hcb=colorbar;

% set(hcb,'YTick',[])

% cmap = jet(256);

% cmap(1,:) = 1;

% colormap(cmap);


axis xy

xlim([-0.1,0.05])

ylim([0 500]);

xlabel('time (s)')

ylabel('frequency (Hz)')

colormap(c);

% colormap bone


subplot(224)

tt = time_H1(indxt)-tevent;

lim = 0.5;

indxt2 = find(tt>=-lim & tt<=lim);

[S,F,T] = spectrogram(strain_L1_whitenbp1(indxt2),blackman(fs/16),fs/16*15/16,fs/16,fs,'reassigned','yaxis');

tvec = 2*lim/(T(end)-T(1))*(T-T(1))-lim;

imagesc(tvec,F,abs(S))


% colorbar

% hcb=colorbar;

% set(hcb,'YTick',[])

% cmap = jet(256);

% cmap(1,:) = 1;

% colormap(cmap);


axis xy

xlim([-0.1,0.05])

ylim([0 500]);

xlabel('time (s)')

ylabel('frequency (Hz)')

colormap(c);

% colormap bone


 

% subplot(212)

% f = (0:length(strain_H1_whitenbp)-1)*fs/length(strain_H1_whitenbp);

% plot(f,abs(fft(strain_H1_whitenbp1)),f,abs(fft(strain_H1_whitenbp)))


%% SSD

% % [RR1] = fastSSD(strain_H1_whitenbp1,fs,0.001,[],100);

% % [RR2] = fastSSD(strain_L1_whitenbp1,fs,0.001,[],100);

%

% load decGW

%

% % load SSD_L-L1_LOSC_4_V2-1126259446-32

% % tt = time_H1(indxt)-tevent;

% indxt2 = find(tt>=-0.1 & tt<=0.05);

% % indxt2 = find(tt>=-0.5 & tt<=0.5);

% subplot(211)

% plot(tt(indxt2),strain_H1_whitenbp1(indxt2),tt(indxt2),RR1(:,indxt2)')

% xlabel('time (s) around event')

% ylabel('whitented strain')

% title('SSD components of Hanford strain data near GW150914')

% subplot(212)

% plot(tt(indxt2),strain_H1_whitenbp1(indxt2),tt(indxt2),RR2(:,indxt2)')

% xlabel('time (s) around event')

% ylabel('whitented strain')

% title('SSD components of Livingston strain data near GW150914')


%% SSD on time interval around event

tindex = 0.15;

% intzero = find(tt>=-tindex & tt<=tindex);

intzero = find(tt>=-tindex & tt<=0.05);

% strain_H1_whitenbp1(intzero) = 0;

y1 = strain_H1_whitenbp1(intzero);

y2 = strain_L1_whitenbp1(intzero);


% SSD

% indxt2 = find(tt>=-tindex & tt<=tindex);

indxt2 = find(tt>=-tindex & tt<=0.05);

[RR1] = peakSSD(y1,fs,0.01,100);

[RR2] = peakSSD(y2,fs,0.001,100);


% subplot(211)

% plot(tt(indxt2),y1,tt(indxt2),RR1')

% % xlim([-0.1, 0.05])

% xlabel('time (s) around event')

% ylabel('whitented strain')

% title('SSD components of Hanford strain data near GW150914')

% subplot(212)

% plot(tt(indxt2),y2,tt(indxt2),RR2')

% % xlim([-0.1, 0.05])

% xlabel('time (s) around event')

% ylabel('whitented strain')

% title('SSD components of Livingston strain data near GW150914')


f = (0:length(RR1)-1)*fs/length(RR1);

subplot(321)

plot(tt(indxt2),y1)

ylim([-0.02 0.02]);

title('Hanford strain data near GW150914')

subplot(322)

plot(f,abs(fft(y1)))

xlim([0 500]);

title('Amplitude spectrum')

subplot(323)

plot(tt(indxt2),RR1(2,:))

ylim([-0.02 0.02]);

title('SSD component 1')

subplot(324)

plot(f,abs(fft(RR1(2,:))))

xlim([0 500]);

title('Amplitude spectrum')

subplot(325)

plot(tt(indxt2),RR1(3,:))

ylim([-0.02 0.02]);

title('SSD component 2')

xlabel('time (s) around event')

subplot(326)

plot(f,abs(fft(RR1(3,:))))

xlim([0 500]);

title('Amplitude spectrum')

xlabel('frequency (Hz)')

% subplot(414)

% plot(tt(indxt2),RR1(3,:))

% ylim([-0.02 0.02]);

% title('SSD component 3')