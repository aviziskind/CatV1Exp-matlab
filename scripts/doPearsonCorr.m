function [cc, cc_p] = doPearsonCorr(x,y)
    if nargout == 1
        try
            cc = pearsonR(x(:), y(:) );
        catch
            cc = corr(x(:), y(:) );
        end  
    elseif nargout == 2
        try
            [cc, cc_p] = pearsonR(x(:), y(:) );
        catch
            [cc, cc_p] = corr(x(:), y(:) );
        end                  
    end
    if isnan(cc)
        3;
    end
end