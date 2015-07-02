function plotOriWidthVsOriPref

%     gratType = 'Flashed';
    gratType = 'Drifting';

    switch gratType
        case 'Flashed', S_ori = load('flashedGratingCells_GLFcuw8_degree_ori.mat');
        case 'Drifting', S_ori = load('driftingGratingCells_GLFcuw8_degree_ori');
    end
    oriStats = nestedFields(S_ori, 'allCells', 'stats', 'tuningStats', 'oriStats_si');
    %%
    w_glob = [oriStats.w_ori_global];
    w_loc = [oriStats.w_ori_local];
    ori_pref = [oriStats.ori_pref_deg];

    binE = [0:15:180]; binC = binEdge2cent(binE);    
    [n, binIds] = histcnt(ori_pref, binE);
    [uBinIds, binIdxs] = uniqueList(binIds);
    
    binned_w_glob_mean = cellfun(@(idx) nanmean(w_glob(idx)), binIdxs);
    binned_w_glob_std = cellfun(@(idx) nanstderr(w_glob(idx)), binIdxs);
    binned_w_loc_mean = cellfun(@(idx) nanmean(w_loc(idx)), binIdxs);
    binned_w_loc_std = cellfun(@(idx) nanstderr(w_loc(idx)), binIdxs);

    %%
    figure(10); clf;
    hist(ori_pref, length(binC)*1.5);    
    xlabel('Preferred Orientation'); 
    set(gca, 'xtick', binE);
    title(sprintf('%s Gratings : Preferred Orientation', gratType));
    3;
    %%
    figure(11); clf;
    plot(ori_pref, w_glob, '.')
    hold on
%     plot(binC, binned_w_glob_mean, 'ro-');
    errorbar(binC, binned_w_glob_mean, binned_w_glob_std, 'ro-');
    xlabel('Preferred Orientation'); ylabel('Global ORI tuning width');
    set(gca, 'xtick', binE);
    title(sprintf('%s Gratings : Global ORI width vs Preferred ori', gratType));
    3;

    figure(12); clf;
    plot(ori_pref, w_loc, '.');
    hold on
%     plot(binC, binned_w_loc_mean, 'ro-');
    errorbar(binC, binned_w_loc_mean, binned_w_loc_std, 'ro-');
    xlabel('Preferred Orientation'); ylabel('Local ORI tuning width');
    title(sprintf('%s Gratings : Local ORI width vs Preferred ori', gratType));
    set(gca, 'xtick', binE);
    3;
    
    

end


% function binned