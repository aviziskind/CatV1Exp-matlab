function ctrl_dist = getCtrlDist(X, binEdges)
    X = X(~isnan(X));
    if nargin < 2
        binEdges = [-inf ; unique(X(:)); inf];    
    end
    binCounts  =  histc (X(:), binEdges, 1);
    sumCounts  =  cumsum(binCounts)./sum(binCounts);
    
    CDF  =  sumCounts(1:end-1);    
    ctrl_dist = struct('binEdges', binEdges, 'CDF', CDF);
end
   
