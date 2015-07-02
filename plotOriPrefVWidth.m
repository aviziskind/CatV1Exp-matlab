S = load('driftingGratingCells_GLFcuw8_degree_ori.mat');

idx_cell = [S.allCells.cellId] > 0;
oriStats = nestedFields(S.allCells(idx_cell), 'stats', 'tuningStats', 'oriStats_ss');
w_loc = [oriStats.w_ori_local];
w_glob = [oriStats.w_ori_global];
ori_pref = [oriStats.ori_pref_deg];

%%
ori_bins = [0:5:180];
[n, ori_bin_idx] = histcnt(ori_pref, ori_bins);

nbins = length(ori_bins)-1;
for i = 1:nbins
    meanLoc(i) = mean( w_loc(ori_bin_idx == i) );
    stdLoc(i) = std( w_loc(ori_bin_idx == i) );
    meanGlob(i) = mean( w_glob(ori_bin_idx == i) );
    stdGlob(i) = std( w_glob(ori_bin_idx == i) );
end
%%
figure(1); clf;
plot(ori_pref, w_loc, '.'); hold on;
stairs(ori_bins, meanLoc([1:end, end]), 'r', 'linewidth', 2)
 h = errorbar(binEdge2cent(ori_bins), meanLoc, stdLoc, 'r.', 'linestyle', 'none');
xlabel('Preferred Orientation'); ylabel('Local Orientation Width');

%%
figure(2); clf;
plot(ori_pref, w_glob, '.'); hold on;
stairs(ori_bins, meanGlob([1:end, end]), 'r', 'linewidth', 2)
 h = errorbar(binEdge2cent(ori_bins), meanGlob, stdGlob, 'r.', 'linestyle', 'none');
xlabel('Preferred Orientation'); ylabel('Global Orientation Width');

