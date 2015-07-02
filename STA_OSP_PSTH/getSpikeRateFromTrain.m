function [r, T_binStarts] = getSpikeRateFromTrain(t_sec, binsize_ms)

    if nargin < 2
        binsize_ms = 50;
    end
    msPerSec = 1000;
    binsize_sec = binsize_ms / msPerSec;

    T_start = min(t_sec) - binsize_sec/2;
    T_end   = max(t_sec) + binsize_sec/2;
    
    T_end = T_end + mod( T_end - T_start , binsize_sec);
%     nBins = (T_end - T_start)/binsize_sec;

    T_bins = T_start : binsize_sec : T_end;
    T_binStarts = T_bins(1:end-1);    
        
    r = elementsInRange(t_sec, T_bins, 'count') / binsize_sec;


        
end