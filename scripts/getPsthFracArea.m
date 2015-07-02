function A = getPsthFracArea(cur_psth, l_bin, r_bin)
    nBins = length(cur_psth);
    psth_total = sum(cur_psth);
    % PSTH should be mean-subtracted and rectified *before* being passed to this function.
    frac_area = sum(cur_psth(l_bin:r_bin))/psth_total;
    frac_bins = (r_bin-l_bin+1)/nBins;
    
    A = frac_area/sqrt(frac_bins);
end