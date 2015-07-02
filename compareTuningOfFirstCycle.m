function compareTuningOfFirstCycle


    curGratingType('d');
    
%     [Gids_all, cellIds_all] = getAllGids('o'); 
    S = load('driftingGratingCells_GLFcuw8_degree_ori.mat');
    idx_cell = [S.allCells.cellId] > 0;    
    idx_rep = nestedFields(S.allCells, 'stats', 'isRep');
    idx_use = find(idx_cell & idx_rep);
    
    Gids = [S.allCells(idx_use).Gid];
    cellIds = [S.allCells(idx_use).cellId];
    nCells = length(cellIds);
    %%

    GC = [619,2;
          2296,2;
          430,3;
          1118,1;
          742,3;2264,4;2264,5;2294,3];
            
    Gids    = GC(:,1);
    cellIds = GC(:,2);
        
    
    
    for ci = 1:nCells
        
        [otc_c1, otc_rest] = getOriTuningCurves(Gids(ci), cellIds(ci));
        3;

    end
    
end

function [otc_cycle1, otc_rest] = getOriTuningCurves(Gid, cellId)

    opt.keepFirstDriftingGratingCycle = 1;

    [binCent, allHistVals, meanRate] = dbGetCellSpkStimHists(Gid, cellId, opt);
%%
    [nOri, nSp, nPh, nCycles, nRep] = dbGetUniqueOriSpPh('Gid', Gid, 'length');
    [oris, spf, phs] = dbGetUniqueOriSpPh('Gid', Gid);

    OSP_full = reshape(allHistVals, [nOri, nSp, nPh, nCycles, nRep]);
    OSP_full = OSP_full/mean(OSP_full(:))*meanRate;
    
    OSP_cycle1 = OSP_full(:,:,:,1,:);
    OSP_rest = OSP_full(:,:,:,2:end,:);
    
    ori_ph_cycle1 = reshape( apply2dims(@mean, OSP_cycle1, [2,4,5]), [nOri, nPh]);
    ori_ph_rest   = reshape( apply2dims(@mean, OSP_rest,   [2,4,5]), [nOri, nPh]);
    ori_ph_all    = reshape( apply2dims(@mean, OSP_full,   [2,4,5]), [nOri, nPh]);
        
    otc_cycle1 = mean(ori_ph_cycle1, 2);
    otc_rest = mean(ori_ph_rest, 2);
    otc_all = mean(ori_ph_all, 2);
    
    %%
    rescaleColors = 0;
    addColorBars = 1;
    
    figure(1); clf;
    h = subplotGap(1,2,1, 1);
    plot(oris, otc_cycle1, 'bo-', oris, otc_rest, 'r.-', oris, otc_all, 'k.-');
    title(sprintf('Phase averaged Ori-tuning curve [%d : %d]', Gid, cellId));
    xlabel('Orientation'); set(gca, 'xtick', [0:90:360]);
    xlim([0 360]);
    ylabel('Average firing rate');
    legend({'First cycle', 'other cycles', 'all cycles'}, 'location', 'best');
    
    h(1) = subplotGap(2,2,1,2); imagesc(oris, phs, ori_ph_cycle1);
    set(gca, 'xtick', [0:90:360], 'ytick', [0:90:360]);
    ylabel('Orientation'); 
    xlabel(' ');
    title('First cycle');
    if addColorBars
        colorbar;
    end
    
    h(2) = subplotGap(2,2,2,2); imagesc(oris, phs, ori_ph_rest);
    ylabel('Orientation'); 
    set(gca, 'xtick', [0:90:360], 'ytick', [0:90:360]);
    xlabel('Phase');
    title(sprintf('Cycles %d-%d', 2, nCycles));
    if addColorBars
        colorbar;
    end
    
    if rescaleColors
        CL = lims([get(h(1), 'clim'), get(h(2), 'clim')]);
        set(h, 'clim', CL);
    end
    
    3;

end