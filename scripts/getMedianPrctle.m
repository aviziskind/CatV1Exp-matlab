function [med, L, U] = getMedianPrctle(vals, dim, abs_flag)
    if ~exist('dim', 'var') || isempty(dim)
        dim = find(size(vals) > 1, 1);
    end
    do_absPos = exist('abs_flag', 'var')  && isequal(abs_flag, 1);

    med = prctile(vals, 50, dim);

    std_dev_pct = 34.1;
    half_vals = 25;
%     pctle_use = std_dev_pct;
    pctle_use = half_vals;
    
    p_lo = prctile(vals, 50-pctle_use, dim);
    p_hi = prctile(vals, 50+pctle_use, dim);
    
    if do_absPos
        L = p_lo;
        U = p_hi;
    else
        L = med-p_lo;
        U = p_hi-med;    
    end
end