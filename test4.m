for i = 1:length(id)
    Si = allOSPs(id(i));
    R = Si.R;
    figure(89);
    R_ori_pref = squeeze( R(Si.oriSp_maxR_av(1), :,:) );
    imagesc(R_ori_pref)
    title(sprintf('ph_max = %d', Si.ph_max));
    3;    
end