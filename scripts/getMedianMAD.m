function [med, L, U] = getMedianMAD(vals, dim, abs_flag)
    if ~exist('dim', 'var') || isempty(dim)
        dim = find(size(vals) > 1, 1);
    end        
    do_absPos = exist('abs_flag', 'var')  && isequal(abs_flag, 1);
        
    med = nanmedian(vals, dim);

    abs_deviations = abs(bsxfun(@minus, vals, med));
    MAD = nanmedian( abs_deviations, dim);
    
    if do_absPos
        L = med+MAD;
        U = med-MAD;
    else
        L = MAD;
        U = MAD;
    end
end