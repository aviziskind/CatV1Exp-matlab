function [x_cent, N_spikesPerSec] = calcPSTH( spkTsRelToWindow, windowSize, nPres, timeUnitsPerSec, nBins)

    if (~exist('nPres', 'var') || isempty(nPres)) && iscell(spkTsRelToWindow)
        nPres = length(spkTsRelToWindow); 
    end
    if ~exist('timeUnitsPerSec', 'var') || isempty(timeUnitsPerSec)
        timeUnitsPerSec = 1000;  % assume input is in ms unless otherwise specified
    end
    if ~exist('nBins', 'var') || isempty(nBins)
        nBins = 20;              % default value
    end
%     [binSize_ms, nBins] = getBestBinSize(spkTsRelToWindow, windowSize);

    if iscell(spkTsRelToWindow)
        spkTsRelToWindow = [spkTsRelToWindow{:}];
    end
    
    if length(windowSize) == 1
        [windowStart, windowEnd] = deal(0, windowSize);
    elseif length(windowSize) == 2
        [windowStart, windowEnd] = deal(windowSize(1), windowSize(2));        
    end
    windowSize_ms = windowEnd-windowStart;
        
    binSize_ms = windowSize_ms / nBins;
    x_edges = [windowStart : binSize_ms : windowEnd]';
    x_cent = x_edges(1:nBins) + binSize_ms/2;

    binTime_sec = (windowSize_ms / nBins) / timeUnitsPerSec;

    N_hist = histcnt( spkTsRelToWindow , x_edges );

    N_spikesPerSec = (N_hist / nPres) / binTime_sec;  % convert rate to spk/sec
end  


%     N_hist0 = histc( spkTsRelToWindow , x_edges );
%     if isempty(N_hist0),
%         N_hist0(end,1) = 0; % make all zeros, instead of empty.
%     end
%     N_hist = N_hist0(1:nBins);
%     N_hist(nBins) = N_hist(nBins) + N_hist0(nBins+1); % add the last bin (The last bin counts any values of x that match edges(end).)



% function [x_cent, N_hist] = calcPSTH( spk_t_ms, N_bins, X_size_ms, Y_scale)
% %   X_size_ms = 100;
% %   N_bins = 20;
% %   spk_t_ms = rel_t;
% %   Y_scale = 1/N_frm;
%   X_step = X_size_ms / N_bins;
%   x_edges = [0: X_step : X_size_ms];
%   x_cent = x_edges(1:N_bins) + X_step/2;
% 	N_hist0 = histc( spk_t_ms, x_edges );
%   
%   N_hist = N_hist0(1:N_bins);
%   %  The last bin counts any values of x that match edges(end). 
%   N_hist(N_bins) = N_hist(N_bins) + N_hist0(N_bins+1); % add the last bin
%   bin_t_sec = X_size_ms / N_bins / 1000;
% 	N_hist = N_hist * Y_scale / bin_t_sec;  % convert to rate spk/sec
% 	
% end  
