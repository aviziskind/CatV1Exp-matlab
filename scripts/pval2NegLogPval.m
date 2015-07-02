function p = pval2NegLogPval(p)
    global psthStatsSettings
    zerosHaveBeenNudged = psthStatsSettings.nudge0pval_amt > 0;
           % if true: any highly significant p-values (ie. p=0) have been set to 1e-50 with the nudge0pval function.
           % thus, any zeros in the array are 'unassigned' values, and should be ignored (and not converted to
           % highly significant values.)  (useful with sparse p-value matrics).
           % false: zeros are highly significant, so convert them to highly
           % sigificant numbers.

    non_zero_idx = (p ~= 0) & ~isnan(p);
    non_zero_ps = p(non_zero_idx);
    p(non_zero_idx)  = -log10(non_zero_ps); 

    zero_idx     = (p == 0);
    if ~zerosHaveBeenNudged % so any zeros are highly significant: assign the highest value .
        if ~isempty(non_zero_ps)
            m = min(non_zero_ps);
            p(zero_idx) = ceil(-log10(m))+1;  
        else  % ie. p consists only of zeros and/or nans.
            p(zero_idx) = 50;        
        end                
        
    end % otherwise, any zeros are unassigned: leave them.        
        
    
end
