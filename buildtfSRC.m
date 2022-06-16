function [tfSRC,tvecPlot, spectralConc, absEnergy] = buildtfSRC(S,tvec, fridge, iridge, ridgewidth, fs, limits, plotwindow)
%TFSRC returns the time-frequency Spectral Ridge Concentration of a given
%plot

%   The spectral ridge concentration is the concentration of energy within
%   a band of frequencies above and below the ridge frequency normalised
%   relative to the energy in that time window. It outputs a vector that
%   contains the tf-SRC progression through the time window of interest. 

%   INPUTS: 
%   S: a time-frequency representation of a monovariate signal 
%   that has undergone preprocessing and a chirp needs to be detected
%   within
%   tvec: vector containing time points of the signal
%   fridge: the frequency of greatest energy within a time window
%   iridge: the row index of S that the peak frequency lies within
%   ridgewidth: frequency band around the ridge that we want to include.
%   band is centered around the ridge. e.g. band = 10Hz, if let fridge be
%   an arbitrary frequency = 254Hz. The band of frequencies who's energies
%   are considered are [249Hz, 259Hz]. If set to zero then only the peak
%   frequency energy is considered. (Useful if the plot only has 1 energy
%   output per time window such as in a Hilbert-Spectrum output based on a
%   Hilbert Huang Transform. 
%   fs: sampling frequency, defaults to 4096
%   limits: vector containing the lower and upper limit of time to consider
%   before and after the event. defaults to [-0.02, 0.05]
%   plotwindow: the overall window of time around the event to consider.
%   defaults to -0.25, 0.25. (To compare various BBH events against each
%   other this needs to stay consistent.)

%   OUTPUTS:
%   tfSRCConc: vector containing the tfSRC for each time window. 
%   tvecConc: vector containing the time intervals for which the tfSRC was
%   computed to allow plotting with a centered scale.
    
    
    if ridgewidth == 0 || isempty(ridgewidth)
        considerFrame = 0;
    else
        considerFrame = 1;
    end
    
    if nargin < 6 || isempty(fs)
        fs = 4096;
    end
    
    if nargin < 7 || isempty(limits)
        lowerlim = -0.02; upperlim = 0.05;
    else
        lowerlim = limits(1); upperlim = limits(2);
    end
    
    if nargin < 8 || isempty(plotwindow)
        lowerplotwindow = -0.25; upperplotwindow = 0.25;
    else
        lowerplotwindow = window(1); upperplotwindow = window(2);
    end
    
    [~, eventcols] = find(tvec > lowerlim & tvec < upperlim);
    [~, windowcols] = find(tvec > lowerplotwindow & tvec < upperplotwindow);
    
    
    concWindowS = S(:,eventcols);
    concWindowI = iridge(eventcols,:);
    
    concPlotWindowS = S(:,windowcols);
    
    %inner sum get the amount of energy based on the amplitude per column and
    %the second sum adds it all together to get the sum of energy for the
    %window
    totalEPlotWindow = sum(sum(abs(concPlotWindowS).^2));
    tfSRC = zeros(length(eventcols),1);
    
    %In the case for plots such as output of buildHS that output a
    %single energy output per time window and thus no window needs to be
    %considered. 
    if(~considerFrame) %we don't consider the ridgewidth
        for i = 1:length(eventcols)
            %resHS = (fs/2)/ (length(concWindowS)-1); 
            %disp("The res of HS: " + resHS);
            ridgeE = abs(concWindowS(concWindowI(i),i)).^2;
            tfSRC(i,1) = ridgeE / totalEPlotWindow;
        end 
    elseif(considerFrame)
        res = (fs/2)/ (length(concWindowS)-1); 
        if(res >= ridgewidth)
            scaleFactor = ridgewidth/res;
            for i = 1:length(eventcols)
                ridgeE = scaleFactor*sum(abs(concWindowS(concWindowI(i),i)).^2);
                tfSRC(i,1) = ridgeE / totalEPlotWindow;
            end
        else %In this case we need to then begin checking if the reso fits an odd or even number or time within the ridgewidth
            floorcount = floor(ridgewidth/res);
            
            %Cases when there is an odd number of windows in the ridgewidth
            if(rem(floorcount, 2) == 1)
                sideindexs = (floorcount-1)/2;%Number of indexs to go above and below in the plot
                if(floorcount == ridgewidth/res) %Case where it neatly fits into the res. E.g ridgewidth of 15hz with a reso of 5Hz fits 3 'blocks'
                    for i = 1+sideindexs:length(eventcols)-sideindexs
                        ridgeE = sum(abs(concWindowS(concWindowI(i),(i-sideindexs:i+sideindexs)).^2));
                        tfSRC(i,1) = ridgeE / totalEPlotWindow;
                    end
                else
                    %Case where the number fits wholely an odd number of times, then we get the
                    %remainder of the signal then split it into two for each side
                    excessScalar = (ridgewidth/res - floorcount)/2; 
                    for i = 1+sideindexs+1:length(eventcols)-sideindexs-1
                        %ridgeEleft: the split relative partition multiplied by the leftmost index, 
                        ridgeEleft = excessScalar*sum(abs(concWindowS(concWindowI(i),(i-sideindexs-1)).^2));
                        %ridge center(center parts), and 
                        ridgeE = sum(abs(concWindowS(concWindowI(i),(i-sideindexs:i+sideindexs)).^2));
                        %ridgeEright: same as ridgeEleft but for the right end 
                        ridgeEright = excessScalar*sum(abs(concWindowS(concWindowI(i),(i+sideindexs+1)).^2)); 
                        tfSRC(i,1) = (ridgeEleft + ridgeE + ridgeEright) / totalEPlotWindow;
                    end

                end

            elseif(rem(floorcount, 2) == 0) %reso fits wholely an even number of times within the ridgewidth
                sideindexs = (floorcount-2)/2;
                if(floorcount == ridgewidth/res) %case where the even 'blocks' fit neatly in the window. e.g. 20Hz band with 5Hz reso fits 4 blocks
                    excessScalar = 0.5; %as we can split the block into a perfect half on each side we use 0.5 of the lowest and highest window
                    for i = 1+sideindexs+1:length(eventcols)-sideindexs-1
                        %ridgeEleft: the split part times the leftmost index, 
                        ridgeEleft = excessScalar*sum(abs(concWindowS(concWindowI(i),(i-sideindexs-1)).^2));
                        %ridge center(center parts), and 
                        ridgeE = sum(abs(concWindowS(concWindowI(i),(i-sideindexs:i+sideindexs)).^2));
                        %ridgeEright: same as ridgeEleft but for the right end 
                        ridgeEright = excessScalar*sum(abs(concWindowS(concWindowI(i),(i+sideindexs+1)).^2)); 
                        tfSRC(i,1) = (ridgeEleft + ridgeE + ridgeEright) / totalEPlotWindow;
                    end
                else
                %Do some stuff, here the reso fits at most an even number of times and thus
                %needs some scaling to centralise the even number of reso into the
                %ridgewidth 
                %getting the residual excess for each side
                    residualExcess = ((ridgewidth/res - floorcount)/2) ; 
                    excessScalar = 0.5+residualExcess;
                    for i = 1+sideindexs+1:length(eventcols)-sideindexs-1
                        %ridgeEleft: the split part times the leftmost index, 
                        ridgeEleft = excessScalar*sum(abs(concWindowS(concWindowI(i),(i-sideindexs-1)).^2));
                        %ridge center(center parts), and 
                        ridgeE = sum(abs(concWindowS(concWindowI(i),(i-sideindexs:i+sideindexs)).^2));
                        %ridgeEright: same as ridgeEleft but for the right end 
                        ridgeEright = excessScalar*sum(abs(concWindowS(concWindowI(i),(i+sideindexs+1)).^2)); 
                        tfSRC(i,1) = (ridgeEleft + ridgeE + ridgeEright) / totalEPlotWindow;
                    end
                end
           end 
        end
    end 
    %Returning the vector of times to help with plotting the tfSRC
    tvecPlot = tvec(eventcols);
    spectralConc = sum(tfSRC);
    absEnergy = sum(tfSRC)*totalEPlotWindow;
    %disp("This is the total E in the window:" + totalEPlotWindow);
    
end

