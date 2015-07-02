function explorePsthWindowStats

    global psthVarHandles statVarHandles allGidsX allCellIdsX
    global globOmitFlag globW

    persistent allPsthWindowData allCurOsps 
    
    
    % load gid/cellids of all cells
    S = load('flashedGratingCells_all');
    allCells = S.allCells;
    allGids = [allCells.Gid];
    allCellIds = [allCells.cellId];    
    [tmp, idx] = sortrows( [[allGids]', [allCellIds]' ]);    

    
    % load my preferences
    allPrefs = getPrefs;
    
    
    nCells = length(allCells);
    psths = [allCells.PSTH];
    firingRates = [psths.meanRate];
    nspikes = [allCells.nspikes]; %#ok<NASGU>
    windowWidths = zeros(nCells, 1);
    frameLengthMs = round([psths.frameLength_ms]); %#ok<NASGU>
    nonZeroFracs = zeros(nCells, 1);

    [idx_nonZero_ovlp, idx_nonZero_sep, idx_autoZero, idx_savedZero, idx_bothZero, idx_undecided, Lbins_auto_tmp, Rbins_auto_tmp, ...
    idx_nonZero_ovlp_prev, idx_nonZero_sep_prev, idx_autoZero_prev, idx_savedZero_prev, idx_bothZero_prev, idx_undecided_prev] = deal(0);
    
%     numSpikes = ;
    
    binW = 25/6;
    all_bins = [-300 + binW/2 : binW : 200];

%     allPsthWindowData = [];
%     methodName = 'stimulus';
%     methodName = methodName(1:2);
%     psthWindowDataFile = [CatV1Path 'MatLabDB_avi' filesep 'psthWindowData_' methodName];            
%     name_fun = @(Gid, cellId, methName) sprintf('Gid_%04d_cell_%d_%s', Gid, cellId, lower(methName(1:2)));
    
%     allPsthWindowData = load(psthWindowDataFile);    
    redo = true;
    
    [Lbins_saved, Rbins_saved] = deal([allPrefs.L_bin], [allPrefs.R_bin]);        
        
    allStatNames = {'cc_p', 'rho_p', 'rho_p_nz', 'rho_p_nznz', 'r_entropy', 'rep_p'};
    if isempty(allPsthWindowData) || redo
        fprintf('loading ... '); tic;
        tmp_file = [CatV1Path 'MatlabDB_avi' filesep 'psthWindowStatsTmp.mat'];
        if exist(tmp_file, 'file') && ~redo
            S_tmp = load(tmp_file);
            allPsthWindowData = S_tmp.allPsthWindowData;
            allCurOsps = S_tmp.allCurOsps;
        else            
%             allPsthWindowData = arrayfun(@(gid,cellid) getPSTHwindowData(gid, cellid, allStatNames), allGids, allCellIds, 'un', 0);
%             allPsthWindowData = cellfun( @(s) structfun(@pval2NegLogPval, s, 'un', 0), allPsthWindowData, 'un', 0);
            allPsthWindowData = getPSTHwindowData(allGids, allCellIds, allStatNames);

            [Lbins_saved_tmp, Rbins_saved_tmp] = deal(Lbins_saved, Rbins_saved);
            idx = (Lbins_saved_tmp == 0) | isnan(Lbins_saved_tmp);    
            Lbins_saved_tmp(idx) = 78;
            Rbins_saved_tmp(idx) = 101;
%             allCurOsps = getOspsFromLRbins([], allGids, allCellIds, Lbins_saved_tmp, Rbins_saved_tmp);
%             allCurOsps = getOspsFromLRbins([], allGids, allCellIds, Lbins_saved_tmp, Rbins_saved_tmp);
%             [osp, osp_odd, osp_even] = getOspForCellWindow(Gid, cellId, allBins, psthVals, l_bin, r_bin)
            allCurOsps = arrayfun(@(gid,cellid,L,R) getOspDataForPsthWindow(gid,cellid, [],[], L,R, [], 'osp'), ...
                allGids, allCellIds, Lbins_saved_tmp, Rbins_saved_tmp, 'un', 0);
            
            save(tmp_file, 'allPsthWindowData', 'allCurOsps', '-v6');
            
        end
        fprintf(' done. '); toc;
    end
    
    dispFigId = 151;
    figure(dispFigId);  clf; 
%     N = nan(nCells,1);
    h_varxy = plot(0,0, 'bo', 0,0, 'ro', 0,0, 'go', 0,0, 'mo', 0,0, 'ko');  hold on;
    h_varxy_prev = plot(0,0, 'b.', 0,0, 'r.', 0,0, 'g.', 0,0, 'm.', 0,0, 'k.'); 
    h_varxy_xlab = xlabel(' ');
    h_varxy_ylab = ylabel(' abs(diff(bins))');
    h_varxy_tit = title(' ');
    glob_idx= [];
    
    figure(2);  clf; 
    h_fig2     = plot(0,0, 'bo', 0, 0, 'ko'); 
    h_fig2_x = xlabel(' ');
    h_fig2_y = ylabel(' ');
    h_fig2_t = title(' ');

    %     figure(156); clf; h_lrdist = plot(0,0, 'bo'); 
%     selectMethod = 'std';    
    alwaysRedoOsps = false;
    
    allvars = {'binAbsDiffs', 'meanStatVals', 'firingRates', 'numUnique', 'noisiness', 'nspikes', 'windowWidths', 'frameLengthMs', 'nonZeroFracs'};
    osp_vars = {'numUnique', 'noisiness', 'nspikes', 'nonZeroFracs'}; 
    binW = 25/6;
    allBins_ms = [-300+binW/2:binW:200];
    bins_ms_default = [20 120];
    Lbin_default = indmin(abs(allBins_ms-bins_ms_default(1)));
    Rbin_default = indmin(abs(allBins_ms-bins_ms_default(2)));
    
    [Lbins_auto, Rbins_auto, meanStatVals, binAbsDiffs, meanDiffFrom0s] = deal( zeros(1, nCells) );
    doOverlap = @(xL, xR, yL, yR) any(ibetween([xL xR], yL, yR)) || any(ibetween([yL yR], xL, xR));
    
    function updatePlot(doFigs, stat1, intervalSize, nStdTh, objTh, x_var, y_var, w, omitCenter, showPrev)
        allPsthWindowData = getPSTHwindowData(allGids, allCellIds, allStatNames);
        alwaysRedoOsps = false;
        allPrefs = getPrefs;
        [Lbins_saved, Rbins_saved] = deal([allPrefs.L_bin], [allPrefs.R_bin]);

        if doFigs(1)
            randn('state', 0)
            %%% for figure 1:        
            for ci = 1:nCells  
                GC = [allGids(ci), allCellIds(ci)];
                if all(GC == [1929, 1])
                    3;
                end
                stimType = getGratingStimType(allGids(ci));

                [l_bin_auto, r_bin_auto, stat_bin_x, stat_val_y, stat_m, stat_s, meanWindowStat] = ...
                    getLRbin(allGids(ci), allCellIds(ci), all_bins, allPsthWindowData{ci}, intervalSize, stat1, stimType, nStdTh, objTh, firingRates(ci) );
                if ~isnan(meanWindowStat)
                    meanWindowStat = min(meanWindowStat, 15+randn);                        
                end
                                                        
                [Lbins_auto(ci), Rbins_auto(ci), meanStatVals(ci)] = deal(l_bin_auto, r_bin_auto, meanWindowStat);
            end
            updateOSPs = any( strcmp(x_var, osp_vars)) || any( strcmp(y_var, osp_vars)) || alwaysRedoOsps;            
            if updateOSPs || doFigs(2)
                [Lbins_auto_tmp, Rbins_auto_tmp] = deal(Lbins_auto, Rbins_auto);
                idx = (Lbins_auto_tmp == 0) | isnan(Lbins_auto_tmp);
                Lbins_auto_tmp(idx) = Lbin_default;
                Rbins_auto_tmp(idx) = Rbin_default;                
                [allCurOsps] = arrayfun(@(gid,cellid,L,R) getOspDataForPsthWindow(gid,cellid,[], [], L,R, [], 'osp'), ...
                      allGids, allCellIds, Lbins_auto_tmp, Rbins_auto_tmp, 'un', 0);                                    
            end
%             if any(strcmp('numUnique', xy_vars))
%                 numUnique = cellfun(@(s) length(unique(s(:))), allCurOsps); 
%             end
%             numUniqueTh = 4;
    
            overlaps = arrayfun(@(xl, xr, yl, yr) doOverlap(xl, xr, yl, yr), Lbins_auto, Rbins_auto, Lbins_saved, Rbins_saved);
            %nan = undecided. 0 = not reproducible.
                        
            algorithmFound =     (Lbins_auto > 0);% & (numUnique > numUniqueTh);            
            algorithmNotFound = (Lbins_auto == 0);% | ~(numUnique > numUniqueTh);
            manualFound = (Lbins_saved > 0);
            manualNotFound = (Lbins_saved == 0);
            manualUndecided = isnan(Lbins_saved);
            
                nonZeros      = (algorithmFound) & (manualFound);
                idx_nonZero      = find( nonZeros );
            idx_nonZero_ovlp = find( nonZeros & overlaps );        
            idx_nonZero_sep  = find( nonZeros & ~overlaps );
            idx_autoZero     = find(  algorithmNotFound & manualFound );
            idx_savedZero    = find(  algorithmFound & manualNotFound );        
            idx_bothZero     = find(  algorithmNotFound & manualNotFound );
            idx_undecided    = find(  manualUndecided  );
            glob_idx = {idx_nonZero_ovlp, idx_nonZero_sep, idx_autoZero, idx_savedZero, idx_bothZero, idx_undecided};
            assert( all(sort([glob_idx{:}]) == 1:nCells )); % make sure all appear exactly once.
            

            xy_vars = {x_var, y_var};
            if any(strcmp('binAbsDiffs', xy_vars))
                binAbsDiffs(idx_nonZero) = abs(Lbins_auto(idx_nonZero)-Lbins_saved(idx_nonZero)) + abs(Rbins_auto(idx_nonZero)-Rbins_saved(idx_nonZero)) ...
                    + randn(1,nnz(idx_nonZero))/10; %#ok<SETNU>
        %         binAbsDiffs(idx_bothZero) = abs(Lbins_auto(idx_nonZero)-Lbins_saved(idx_nonZero)) + abs(Rbins_auto(idx_nonZero)-Rbins_saved(idx_nonZero)) ...
        %             + randn(1,length(idx_nonZero))/10;
                spc = 5;
                binAbsDiffs(idx_autoZero)  = -1*spc + randn(1,nnz(idx_autoZero))/2;
                binAbsDiffs(idx_savedZero) = -2*spc + randn(1,nnz(idx_savedZero))/2;
                binAbsDiffs(idx_bothZero)  = -3*spc + randn(1,nnz(idx_bothZero))/2;                
            end
            if any(strcmp('numUnique', xy_vars))
                3;
                numUnique = cellfun(@(s) length(unique(s(:))), allCurOsps); %#ok<NASGU>
            end
            if any(strcmp('noisiness', xy_vars))
                globW = w;
                globOmitFlag = omitCenter;
                noisiness = cellfun(@(s) ospNoisiness(s), allCurOsps); %#ok<NASGU>
%                 disp(mean(noisiness))
                noisiness = noisiness/mean(noisiness);                
            end
            if any(strcmp('windowWidths', xy_vars))                                
                spc = 5;
                windowWidths(idx_nonZero)   = Rbins_auto(idx_nonZero)-Lbins_auto(idx_nonZero)+1; %#ok<SETNU>
                windowWidths(idx_autoZero)  = -1*spc + randn(1,nnz(idx_autoZero))/2;
                windowWidths(idx_savedZero) = -2*spc + randn(1,nnz(idx_savedZero))/2;
                windowWidths(idx_bothZero)  = -3*spc + randn(1,nnz(idx_bothZero))/2;                                
            end
            if any(strcmp('nonZeroFracs', xy_vars))                
                nonZeroFracs = cellfun(@nonZeroFrac_fun, allCurOsps);  %#ok<SETNU>
            end
            
            
            Xdata = eval(x_var);
            Ydata = eval(y_var);
            
%             set(h_varxy_xlab, 'string', x_var);
            set(h_varxy_ylab, 'string', y_var);
            
            
            set(h_varxy(1),  'xdata', Xdata(idx_nonZero_ovlp), 'ydata', Ydata(idx_nonZero_ovlp), 'color', 'b' ); 
            set(h_varxy(2),  'xdata', Xdata(idx_nonZero_sep), 'ydata', Ydata(idx_nonZero_sep), 'color', 'r' ); 
            set(h_varxy(3),  'xdata', Xdata(idx_autoZero),  'ydata', Ydata(idx_autoZero), 'color', 'g'  ); 
            set(h_varxy(4),  'xdata', Xdata(idx_savedZero), 'ydata', Ydata(idx_savedZero), 'color', 'm' ); 
            set(h_varxy(5),  'xdata', Xdata(idx_bothZero),  'ydata', Ydata(idx_bothZero), 'color', 'k' ); 

            
            set(h_varxy_prev, 'visible', iff(showPrev, 'on', 'off'));
            if showPrev && (length(idx_nonZero_ovlp_prev) > 1)
                set(h_varxy_prev(1),  'xdata', Xdata(idx_nonZero_ovlp_prev), 'ydata', Ydata(idx_nonZero_ovlp_prev), 'color', 'b' ); 
                set(h_varxy_prev(2),  'xdata', Xdata(idx_nonZero_sep_prev), 'ydata', Ydata(idx_nonZero_sep_prev), 'color', 'r' ); 
                set(h_varxy_prev(3),  'xdata', Xdata(idx_autoZero_prev),  'ydata', Ydata(idx_autoZero_prev), 'color', 'g'  ); 
                set(h_varxy_prev(4),  'xdata', Xdata(idx_savedZero_prev), 'ydata', Ydata(idx_savedZero_prev), 'color', 'm' ); 
                set(h_varxy_prev(5),  'xdata', Xdata(idx_bothZero_prev),  'ydata', Ydata(idx_bothZero_prev), 'color', 'k' );                         
            end
            
            [idx_nonZero_ovlp_prev, idx_nonZero_sep_prev, idx_autoZero_prev, idx_savedZero_prev, idx_bothZero_prev, idx_undecided_prev] = ...
                deal(idx_nonZero_ovlp, idx_nonZero_sep, idx_autoZero, idx_savedZero, idx_bothZero, idx_undecided);

            
            nB = length(idx_nonZero_ovlp); nR = length(idx_nonZero_sep); nG = length(idx_autoZero); nM = length(idx_savedZero); nK = length(idx_bothZero);
            nT = length(meanStatVals);
            str1 = sprintf('R: %d (%.1f%%) B: %d (%.1f%%) G: %d (%.1f%%)', nR, nR/nT*100,   nB, nB/nT*100,    nG, nG/nT*100);
            str2 = sprintf('M: %d (%.1f%%) K: %d (%.1f%%) ', nM, nM/nT*100,    nK, nK/nT*100);

            set(h_varxy_xlab, 'string', {x_var, str1, str2})
            
            
            3;
        end
        
        if doFigs(2)
            % phase cc   vs   f1/dc
            idx = find(algorithmFound);
            
            allStimPhaseTcCCs = zeros(1, nCells);
            F1oDCs80 = zeros(1,nCells);
            progressBar('init-', nCells, 40);
            for i = 1:nCells
%                 fprintf('*');
                progressBar;                
                L = Lbins_auto_tmp(i); R = Rbins_auto_tmp(i);
                Gid = allGids(i);
                cellId = allCellIds(i);
                th = 0.80;
                
%                 [tgtWind_osp, phaseTC_CCs, stimF1oDCs] = getOspDataForPsthWindow(Gid, cellId, [], [], L, R, [], {'osp', 'phaseTC_CCs', 'stimF1oDCs'});
                [tgtWind_osp, phaseTC_CCs, phaseTC_CC_ps, phaseTC_CC_dots, stimF1oDCs] = getOspDataForPsthWindow(Gid, cellId, [], [], L, R, [], {'osp', 'phaseTC_CCs', 'phaseTC_CC_ps', 'phaseTC_Dots', 'stimF1oDCs'});

                % x_axis: mean f1odc of stim > 80%                            
                F1oDCs80(i) = getAllAboveTh(stimF1oDCs, tgtWind_osp, th);                
                
                % y axis: cc of phase tuning cc's & windows.                
                allStimPhaseTcCCs(i) = getAllAboveTh(phaseTC_CCs, tgtWind_osp, th);                
                
%                 tgtWind_osp_norm = tgtWind_osp / max(tgtWind_osp(:));
%                 ph_ccs_wgt = phaseTC_CCs .* tgtWind_osp_norm;                
%                 allStimPhaseTcCCs(i) = doPearsonCorr(tgtWind_osp_norm(:), ph_ccs_wgt(:));
            end
            
            tmpNans = nan(1, length([idx_autoZero_prev idx_savedZero_prev]));
            idx_good = [idx_nonZero_ovlp_prev idx_nonZero_sep_prev];
            idx_bad = idx_bothZero_prev;
%             
%                 set(h_varxy_prev(1),  'xdata', Xdata(), 'ydata', Ydata(idx_nonZero_ovlp_prev), 'color', 'b' );
%                 set(h_varxy_prev(2),  'xdata', Xdata(idx_nonZero_sep_prev), 'ydata', Ydata(idx_nonZero_sep_prev), 'color', 'r' ); 
%                 set(h_varxy_prev(3),  'xdata', Xdata(idx_autoZero_prev),  'ydata', Ydata(idx_autoZero_prev), 'color', 'g'  ); 
%                 set(h_varxy_prev(4),  'xdata', Xdata(idx_savedZero_prev), 'ydata', Ydata(idx_savedZero_prev), 'color', 'm' ); 
%                 set(h_varxy_prev(5),  'xdata', Xdata(idx_bothZero_prev),  'ydata', Ydata(idx_bothZero_prev), 'color', 'k' );                                     
                        
            set(h_fig2(1), 'xdata', [F1oDCs80(idx_good), tmpNans], 'ydata', [allStimPhaseTcCCs(idx_good), tmpNans]);            
            set(h_fig2(2), 'xdata', F1oDCs80(idx_bad), 'ydata', allStimPhaseTcCCs(idx_bad), 'color', 'r');
            
            set(h_fig2_x, 'string', sprintf('F1/DC (averaged over stimuli %d%%+ of max)', th*100));
            set(h_fig2_y, 'string', sprintf('phase-TC CC (av. over stimuli %d%%+ of max)', th*100));
            
        end
        
        
%         if (statNameChanged || intervalSizeChanged || nStdThChanged || objThChanged) && ~isempty(psthVarHandles)
%             manipulateSet(psthVarHandles, {'statName', 'intSize', 'nStdTh', 'objTh'}, {statName, intSize, nStd, objTh});
%         end
        
    end

    doFigs0 = [true, false];
    w0 = .6;
    omitCenter0 = 0;
    objTh0 = 2;
    intSize0 = 5;
    nStdTh0 = 2;
    updatePlot(doFigs0, allStatNames{1}, intSize0, nStdTh0, objTh0, allvars{1}, allvars{2}, w0, omitCenter0, false)
    set(dispFigId, 'windowButtonDownFcn', {@selectPointsInFigure, @updateStatsOfChosenCell_grp});
%     set(nzFigId,   'windowButtonDownFcn', {@selectPointsInFigure, @updateStatsOfChosenCell_grp});
    
%     set(dispFigId, 'buttonDownFcn', @selectPointsInFigure);
%     set(dispFigId, 'windowButtonDownFcn', @selectPointsInFigure);
    
%     set selectPointsInFigure(src, evnt, click1Function)
    
    
    3;

    
%     nStimMethods = {'top n', 'top x%'};
        
    args = { {'doFigs', {[false, true], [false, true]}, doFigs0}, ...
             {'statName', allStatNames}, ...
...              {'stat2', allStatNames}, ...
             {'intSize', [1:20], intSize0}, ...             
             {'nStdTh', [.1:.1:4], nStdTh0}, ...
             {'objTh', [.1:.1:4], objTh0}, ...
             {'x_var', allvars, allvars{1}}, ...
             {'y_var', allvars, allvars{2}}, ...
             {'w', [.6:.05:1], w0}, ...
             {'omitCenter0', [false, true], false}, ...
             {'showPrev', [false, true], false}, ...
           };

    statVarHandles = manipulate(@updatePlot, args, 'FigId', 100);

    
    function updateStatsOfChosenCell_grp(glob_id, grp_id, loc_id, x, y) %#ok<INUSL>

        glob_id = glob_idx{grp_id}(loc_id);
        Gid = allGids(glob_id);
        cellId = allCellIds(glob_id);
        set(h_varxy_tit, 'string', sprintf('(%d) Gid = %d, cellId = %d. [X = %.2f, Y = %.2f]', glob_id, Gid, cellId, x, y));
        cellIdx = find(allGidsX == Gid & allCellIdsX == cellId);
        if ~isempty(cellIdx)
            manipulateSet(psthVarHandles, 'Cell Index', cellIdx );
        else
            beep;
            warning('This cell not found');
        end
    end


    function updateStatsOfChosenCell_onlyFound(glob_id, grp_id, loc_id, x, y) %#ok<INUSL>

        glob_id = glob_idx{grp_id}(loc_id);
        Gid = allGids(glob_id);
        cellId = allCellIds(glob_id);
        set(h_varxy_tit, 'string', sprintf('(%d) Gid = %d, cellId = %d. [X = %.2f, Y = %.2f]', glob_id, Gid, cellId, x, y));
        cellIdx = find(allGidsX == Gid & allCellIdsX == cellId);
        if ~isempty(cellIdx)
            manipulateSet(psthVarHandles, 'Cell Index', cellIdx );
        else
            beep;
            warning('This cell not found');
        end
    end

%     function updateStatsOfChosenCell_ungrp(glob_id, grp_id, loc_id, x, y)
% 
%         Gid = allGids(glob_id);
%         cellId = allCellIds(glob_id);
%         set(h_varxy_tit, 'string', sprintf('(%d) Gid = %d, cellId = %d', p_idx, Gid, cellId));
%         
%         cellIdx = find(allGidsX == Gid & allCellIdsX == cellId);
%         manipulateSet(psthVarHandles, 'Cell Index', cellIdx );
%     end    


    function allPrefs = getPrefs
        psthWindowPrefsFile = [CatV1Path 'MatLabDB_avi' filesep 'psthWindowPrefsFile.mat'];
        S = load(psthWindowPrefsFile);
        allPrefs = S.allPrefs;
%         idx_haveSavedPref = [allPrefs.L_bin]>0;
%         allPrefs = allPrefs(idx_haveSavedPref);    
        [tmp2, idx] = sortrows( [[allPrefs.Gid]', [allPrefs.cellId]' ]);
        allPrefs = allPrefs(idx);
        
        idxWant = arrayfun(@(s) any(s.Gid == allGids & s.cellId == allCellIds), allPrefs);
        allPrefs = allPrefs(idxWant);
    end
    
    function vals = getDiag(statsM, intSize, testPeriod)
        allvals = diag(statsM, -(intSize-1));
        idx_offset = ceil((intSize-1)/2); % [0, 1, 1, 2, 2, ...]

        stat_bin_x = all_bins(1+idx_offset:end-idx_offset);
        fromBin = indmin( abs(stat_bin_x-testPeriod(1)));
        toBin   = indmin( abs(stat_bin_x-testPeriod(2)));

        vals = allvals(fromBin:toBin);        
    end
    
end


%{
function allCurOsps = getOspsFromLRbins(allCurOsps, Gids, cellIds, Lbins, Rbins, redoFlag)
    fprintf('updating osps... ');  tic;
    redo = exist('redoFlag', 'var') && ~isempty(redoFlag);
    if isempty(allCurOsps)
        allCurOsps = struct('osp', [], 'Gid', num2cell(Gids), 'cellId', num2cell(cellIds), 'L_bin', -1, 'R_bin', -1, 'oe_stats', 0 );
    end
    count = 0;
    for i = 1:length(Gids)                

        if ( allCurOsps(i).L_bin ~= Lbins(i) || allCurOsps(i).R_bin ~= Rbins(i) ) || redo
            assert(allCurOsps(i).Gid == Gids(i))
            assert(allCurOsps(i).cellId == cellIds(i))
%             if Gids(i) == 1666 && cellIds(i) == 3
%                 3;
%             end
            [psthBins, allPsthVals] = dbGetCellSpkStimHists(Gids(i), cellIds(i));
            shiftCent = false; curPsth = []; outputVars = {};            
            [tmp1, r, r_odd, r_even] = getOspStatVsPsthBinning(psthBins, allPsthVals, Lbins(i), Rbins(i), shiftCent, getGratingStimType(Gids(i)), curPsth, outputVars);
            oe_stats = getOEstats(r_odd, r_even);
            osp = reshape(r, [36 10]);
            allCurOsps(i).osp = osp;
            allCurOsps(i).L_bin = Lbins(i);
            allCurOsps(i).R_bin = Rbins(i);
            allCurOsps(i).oe_stats = oe_stats;
            count = count + 1;
        end

    end
    fprintf(' done (updated %d/%d). ', count, length(Gids));
    toc;
    

end
%}

function oe_stats = getOEstats(r_odd, r_even)
    z = nnz(     r_odd == 0 & r_even == 0);
    h = nnz( xor(r_odd == 0, r_even == 0) );
    f = nnz(     r_odd > 0 & r_even > 0);
    oe_stats = [z h f];
end


function y = nonZeroFrac_fun(osp)
    st = osp.oe_stats;
    z = st(1); h = st(2); f = st(3);
    y = f/(h);
    
end

function y = getAllAboveTh(X, osp, th)
    idx = (osp(:) > max(osp(:))*th);
    y = mean(X(idx));
end


% function getOsp(id, l_bin, r_bin)
% 
%     persistent allOSPs
% 
%         if isempty( allOSPs{id}{r_bin, l_bin} )
%             [psthBins, allPsthVals] = dbGetCellSpkStimHists(Gids(i), cellIds(i), [], 'stimulus');
%                          
%              [tmp, r, r_odd, r_even] = getOspStatVsPsthBinning(psthBins, allPsthVals, l_bin, r_bin_orig, shiftForEvenNBins, stimType, curPsth, outputVars)
%              v = getValFor1Spk(r(:));
%              r = uint8(round(r/v));
%             
%             
%     
%         else
% 
% 
% 
% end

%{
            %%% for fig 2
            testPeriod = [20, 120];        
            stat_with_0 = 'rho_p';
            stat_without_0 = 'rho_p_nz';
            for ci = 1:nCells  
    %             vals_with_0    = getDiag(allPsthWindowData{ci}.(stat_with_0), intervalSize, testPeriod);
                vals_without_0 = getDiag(allPsthWindowData{ci}.(stat_without_0), intervalSize, testPeriod);            
    %             meanDiffFrom0s(ci) = mean(vals_with_0-vals_without_0);
                meanDiffFrom0s(ci) = max(abs(vals_without_0));
            end
            Xdata = firingRates;
            Ydata = meanDiffFrom0s;            

            set(h_nz(1),  'xdata', Xdata(idx_nonZero_ovlp), 'ydata', Ydata(idx_nonZero_ovlp) ); 
            set(h_nz(2),  'xdata', Xdata(idx_nonZero_sep),  'ydata', Ydata(idx_nonZero_sep) ); 
            set(h_nz(3),  'xdata', Xdata(idx_autoZero),     'ydata', Ydata(idx_autoZero)  ); 
            set(h_nz(4),  'xdata', Xdata(idx_savedZero),    'ydata', Ydata(idx_savedZero) ); 
            set(h_nz(5),  'xdata', Xdata(idx_bothZero),     'ydata', Ydata(idx_bothZero)  );                         
%}