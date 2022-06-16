function [AA,A,tt,f,ff] = buildHS(imf,fs,t,fmin,fmax)

% function to generate the Hilbert spectrum of a signal starting from its
% intrinsic mode functions
% inputs:
% imf: mxN matrix, collecting the m imfs of the signals, each N samples
% long
% t: time vector (in seconds)
% fs: sampling frequency
% fmin: minimum frequency to consider for building the Hilbert spectrum
% fmax: maximum frequency to consider for building the Hilbert spectrum

[A,f,tt] = hhspectrum(imf);
ff = (0:floor(length(imf)/2)-1)*(fs)/(length(imf));
AA = zeros(length(ff),length(tt));

for i = 1:size(imf,1),
    clear ii ftemp
    ii = find(f(i,:)*fs>fmax);
    ftemp = f(i,:);
    ftemp(ii) = NaN;
    f(i,:) = fillmissing(ftemp,'spline');
end


for c = 1:size(imf,1),
    Yimf = abs(fft(imf(c,:)));
    Yimf(floor(length(Yimf)/2)+1:end) = [];
    [~, Im] = max(Yimf);
    
    if ff(Im)>=fmin && ff(Im)<=fmax,
        for p = 1:size(imf,2)-2
            [~,ii] = min(abs(ff-f(c,p)*fs));
            AA(ii,p) = AA(ii,p) + A(c,p);
        end
    end
end

%imagesc(t,ff,AA)
%axis xy
%ylim([0 500]);
%xlim([0.19,0.21])
%xlabel('time (s)')
%ylabel('frequency (Hz)')