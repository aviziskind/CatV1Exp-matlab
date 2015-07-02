%%
[allGids, allCellIds] = getAllGids('f');
%%
for i = 1:length(allCellIds)
    %%
    Gid = allGids(i);
    cellId = allCellIds(i);
    if flashedOrDrifting(Gid) == 1
        [PSTHdata, stats] = getPSTHforCell(Gid, cellId);
        LR_bins = PSTHdata.timeWindow_bins;
        [L_bin_best, R_bin_best, windowProfile_best] = deal(LR_bins(1), LR_bins(2), PSTHdata.windowProfile);
    else
        [L_bin_best, R_bin_best, windowProfile_best] = deal(1,1,[]);
    end
    
    fprintf(' (%d/%d), Gid = %d, cellId = %d [%d, %d]\n', i, length(allCellIds), Gid, cellId, L_bin_best, R_bin_best)
    getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_best, R_bin_best, windowProfile_best, {'osp_ph'}, 0);

end


