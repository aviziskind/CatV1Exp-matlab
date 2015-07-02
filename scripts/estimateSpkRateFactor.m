function [f, X2] = estimateSpkRateFactor(X)
    vals = unique(X(:));
    v_max = max(vals);
    
    diff_min = 1e-3;
    allValDiffs = diff(vals);
    allValDiffs = allValDiffs(allValDiffs > diff_min);
        
    i_max = ceil(v_max / min(allValDiffs)*3);
    
    is = 1:i_max;
    maxRelDiffs = zeros(1, length(is));
    for i = is
        v_step = v_max/i;
        ns = vals / v_step;
        diffs = abs(ns - round(ns));
        
        maxRelDiffs(i) = max(diffs);
        
    end

    idx_best = find(maxRelDiffs < .01, 1);
    
    f = v_max / is(idx_best);
    
    X2 = round(X/f);
    3;


end