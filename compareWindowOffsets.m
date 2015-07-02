
    %%
    S = load('flashedGratingCells_GLFcuw8_degree_SS_ori.mat');
    idx_cells = [S.allCells.cellId] > 0;
    allCells = S.allCells(idx_cells);

    allGids = [allCells.Gid];
    allCellIds = [allCells.cellId];
    nCells = length(allCellIds);

    offsets = [-2, -1, 0, 1, 2, nan, 3];
    nOffsets = length(offsets);

    % allFrameLengths = arrayfun(@(gid) 1000/dbGetFramesPerSecond('Gid', gid), allGids);
    fps = arrayfun(@(gid) dbGetFramesPerSecond('Gid', gid), allGids);

    [ori_pref,d_ori_pref, w_ori_global, d_w_ori_global, d_w_ori_local, w_ori_local, d_w_spf, w_spf] = deal(nan(nCells, nOffsets));

    progressBar('init-', nCells*nOffsets);
    for i = 1:nCells

        for j = 1:nOffsets
            progressBar;
            if isnan(offsets(j))
                continue;
            end
            if offsets(j) == 3
                 LR_j = allCells(i).PSTH.timeWindow_stimw;
                 windowProfile_j = allCells(i).PSTH.windowProfile_stimw;
                 3;
            else
                LR = allCells(i).PSTH.timeWindow_bins;
                LR_j = LR + offsets(j)*[-1, 1];
                windowProfile_j = allCells(i).PSTH.vals(LR_j(1):LR_j(2));
                if diff(LR_j) < 2
                    continue;
                end
            end
            tuningStats = getOspDataForPsthWindow(allGids(i), allCellIds(i), [], [], LR_j(1), LR_j(2), windowProfile_j, {'tuningStats'}, 0);
            ori_pref(i,j) = tuningStats.oriStats_ss.ori_pref_deg;
            w_ori_global(i,j) = tuningStats.oriStats_ss.w_ori_global;
            w_ori_local(i,j) = tuningStats.oriStats_ss.w_ori_local;
            w_spf(i,j) = tuningStats.spfStats_ss.w_spf;
            if w_spf(i,j) == 0
                3;
            end
            3;

            
        end        
    end
    %%
    for i = 1:nCells
        d_ori_pref(i,:) = circDist( ori_pref(i,:), ori_pref(i,3), 180);
        d_w_ori_global(i,:) = (w_ori_global(i,:) - w_ori_global(i,3))/w_ori_global(i,3);
        d_w_ori_local(i,:) = (w_ori_local(i,:) - w_ori_local(i,3))/w_ori_local(i,3);
        d_w_spf(i,:) = (w_spf(i,:) - w_spf(i,3))/w_spf(i,3);



    end
    %%
    all_fps = [10, 24, 60];
    [uFPS, fps_idx] = uniqueList(fps);
    nfps = 3;
    figure(1); clf; hold on; box on;
    for j = 1:nfps
        X = d_w_ori_local(fps_idx{j},:)*100;
        errorbar(offsets, nanmedian(X, 1), prctile(X, 33), prctile(X, 66), [color_s(j) 'o-'], 'markersize', 7)

    end
    xticklabels = cellnum2cellstr( num2cell( nonnans(   offsets )) );
    xticklabels{end} = 'frameLength';
    xlabel('Change in size of analysis window on each side');
    ylabel('Percent change in tuning width');
    title('SpfFreq. width');
    title('SpfFreq. width');
%     title('ORI width (Local)');
    set(gca, 'xtick', nonnans(offsets), 'xticklabel', xticklabels)
    legend(legendarray('fps = ', uFPS, ' Hz' ), 'location', 'SE');
    drawHorizontalLine(0, 'color', 'k', 'linestyle', ':');

    % tfigure(3);


    3;

