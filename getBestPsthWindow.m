function [L_bin, R_bin, nStim] = getBestPsthWindow(data, W)

    [pvals, slopes, areas] = deal(data.pvals, data.slopes, data.areas);
    nBins = size(pvals,1);    
    bins = (25/6)*[.5:1:nBins];

    singleBins = eye(nBins);
    tooEarly = zeros(nBins); tooEarly( :, bins < 12 ) = 1;  %
    tooLate  = zeros(nBins); tooLate ( bins > 130, :) = 1;
    [X,Y] = meshgrid(1:nBins, 1:nBins);            
    tooLong  = arrayfun(@(x,y) y > x + 17, X, Y);  % >70ms

    ok = ~(singleBins + tooEarly + tooLate + tooLong);
    
    pvals_n = slopes / max(slopes(:));
    psth_area_n = areas / max(areas(:));
    T = (W*pvals_n + (1-W)*psth_area_n) .* ok;
    
    [max_T, idx_maxT] = maxElement(T);
    if ndims(T) == 2
        [R_bin, L_bin] = dealV( idx_maxT );            
    elseif ndims(T) == 3
        [R_bin, L_bin, nStim] = dealV( idx_maxT );            
    end
    
end
