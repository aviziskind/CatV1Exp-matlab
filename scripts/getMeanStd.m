function [mn, L, U] = getMeanStd(vals, dim, abs_flag)
    if ~exist('dim', 'var') || isempty(dim)
        dim = find(size(vals) > 1, 1);
    end        
    do_absPos = exist('abs_flag', 'var')  && isequal(abs_flag, 1);
    
    mn = nanmean(vals, dim);
    
    stdev = nanstd(vals, [], dim);
    
    
    if do_absPos    
        L = mn+stdev;
        U = mn-stdev;
    else
        L = stdev;
        U = stdev;
    end
end