
n = length(allOSPs);

nPh = arrayfun(@(s) length(s.ph), allOSPs);
idx_40ph = find(nPh == 40);

% osp_inds = 1:n;
osp_inds = idx_40ph;

ph_max_av = zeros(size(osp_inds));
ph_max_no_av = zeros(size(osp_inds));


for i = 1:length(osp_inds)
    idx = osp_inds(i);
    ori_sp_av = allOSPs(idx).oriSp_maxR_av;
    ori_sp_no_av = allOSPs(idx).oriSp_maxR_no_av;
    ph = allOSPs(idx).ph;
    R = allOSPs(idx).R;
    [nOri, nSp, nPh] = size(R);
    
    ori_sp_av = [randi(nOri), randi(nSp)];
    ori_sp_no_av = [randi(nOri), randi(nSp)];
    
    tc_av = squeeze(R(ori_sp_av(1), ori_sp_av(2), :));
    tc_no_av = squeeze(R(ori_sp_no_av(1), ori_sp_no_av(2), :));
    [mx] = max( tc_av );
    if any(tc_av > 0) 
        idx = find(tc_av==mx); 
        ph_max_av_i = ph(idx(randi(length(idx))));
    else
        ph_max_av_i = nan;
    end    
    
    [mx] = max( tc_no_av  );
    if any(tc_no_av > 0)
        idx = find(tc_no_av==mx); 
        ph_max_no_av_i = ph(idx(randi(length(idx))));
    else
        ph_max_no_av_i = nan;
    end    
    
    ph_max_av(i) = ph_max_av_i;
    ph_max_no_av(i) = ph_max_no_av_i;    
end
str = 'skipping 1st cycle';
nbin = 20;
figure(101); clf; hist(ph_max_av, nbin); title(['peak (av): ' str '. mean = ' num2str(nanmean(ph_max_av), '%.1f')]); xlim([0 360]);
figure(102); clf; hist(ph_max_no_av, nbin); title(['peak (no av): ' str '. mean = ' num2str(nanmean(ph_max_no_av), '%.1f')]); xlim([0 360]);