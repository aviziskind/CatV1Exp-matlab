function [l_bin, r_bin, T] = getBestLbinRbin(bins, pvals, areas, W)
    persistent ok_lr

    nBins = length(bins);
    binWidth = diff(bins(1:2));

    if isempty(ok_lr)
    
        singleBins = eye(nBins);
        tooLong_ms = 90;
        tooLong_nbins = floor(tooLong_ms / binWidth);
    %   doubleBins = diag(ones(nBins-1, 1), -1);
        tooEarly = zeros(nBins); tooEarly( :, bins < 15 ) = 1;
        tooLate  = zeros(nBins); tooLate ( bins > 150, :) = 1;
          [X,Y] = meshgrid(1:nBins, 1:nBins);            
        tooLong  = arrayfun(@(x,y) y > x + tooLong_nbins, X, Y);
        
        ok_lr = ~(singleBins + tooEarly + tooLate + tooLong);
    end
    
    p_max = max(pvals(:)); p_min = min(pvals(:));
    pvals_n = (pvals-p_min)/(p_max - p_min);
    
    a_max = max(areas(:)); a_min = min(areas(:));
    areas_n = (areas-a_min)/(a_max - a_min);

    T = (W*pvals_n + (1-W)*areas_n) .* ok_lr;
            
    [max_T, idx_maxT] = maxElement(T);
    [r_bin, l_bin] = dealV( idx_maxT );
end