function explorePSTHeffectOnReproducibility(Gid_start, cellId_start)

    global psthVarHandles allGidsX allCellIdsX
%     global statVarHandles 
    global psthStatsSettings
    
    if isempty(psthStatsSettings)
        setPsthGlobalSettings;
    end
    
    ext_rep   = @(x) [x(:); x(end)];        
    wrp = @(x) [x(:); x(1)];

    S1 = load('flashedGratingCells_all');
    allCells = S1.allCells;
    
    S2 = load('cellsGroups_movie_fg.mat');
    mcg = S2.movieGroups_fg;
    mcg_Gids = [mcg.Gid];
    
    
    doPsth1CurAll = false;
    ignoreMUs = false;
    sortByRep = true;
        sortByFR = false;
    doCheckStimOrder = false;
    recordPreferences = false;
        doOnlyCellsWithPref = [];
    selectStimType = '';
    cph_separateBeforeAfterTrials = true;
        cphSplitTime_ms = 70;    
    showStatImages = true;
%         showOnlyCurStat = true;
    showPhRepAllStim = false;
    showPhRepAllStim_sampled = false;
    showF1oDCallStim = false;
    
    
    doOSP8 = false;
    doWhichStim = false;
    doPhaseTC = false;
    doRvsNStim = false;
    psthWindow = [-300, 200];                
    nStimMax = 50;

    doFullWindowProfile = true;
    doTgtWindowProfile = true;
    doTgtWgtWindowProfile = false;
    doAllProfilesIn1Fig = true;
    
    doProfiles = [doFullWindowProfile, doTgtWindowProfile, doTgtWgtWindowProfile];
    
    
    if recordPreferences 
       psthWindowPrefsFile = [CatV1Path 'MatLabDB_avi' filesep 'psthWindowPrefsFile.mat'];
       if exist(psthWindowPrefsFile, 'file')
           S_pref = load(psthWindowPrefsFile);
           allPrefs = S_pref.allPrefs;
       else
           allPrefs = struct('Gid', {allCells.Gid}, 'cellId', {allCells.cellId}, 'L_bin', 0, 'R_bin', 0, 'intSize', 0);
       end
       
       if ~isempty(doOnlyCellsWithPref)
            if isnan(doOnlyCellsWithPref)
                hasPref = isnan([allPrefs.L_bin]);
            elseif doOnlyCellsWithPref == 0 %#ok<BDSCI>
                hasPref = [allPrefs.L_bin] == 0;
            else
                hasPref = [allPrefs.L_bin] > 0;
            end
            assert(length(hasPref) == length(allPrefs)); %ie. no empty entries
           
            Gids_withPref    = [allPrefs(hasPref).Gid];
            cellIds_withPref = [allPrefs(hasPref).cellId];
            
            idx = arrayfun(@(c) any(c.Gid == Gids_withPref & c.cellId == cellIds_withPref), allCells);
            allCells = allCells(idx);                        
            
%             ss = struct('Gid', {allCells.Gid}, 'cellId', {allCells.cellId}, 'L_bin', 0, 'R_bin', 0, 'intSize', 0);
%             allPrefs = allPrefs;
%             allPrefs(end+1:end+length(find(idx_no_pref))) = ss;                
            
       end
       
    end
       
       
    
    if sortByRep
        rep = arrayfun(@(s) s.stats.allWindowStats.cc_p, allCells);
        allCells = allCells(ord(rep, 'descend'));
    elseif sortByFR
        psths = [allCells.PSTH];
        allFRs = [psths.meanRate];
        allCells = allCells(ord(allFRs, 'ascend'));
    end

    
    if ignoreMUs
        allCells = allCells( [allCells.cellId] > 0);
    end
    if ~isempty(selectStimType);
        indStimTypeMatch = find( arrayfun(@(s) ~isempty(strfind(s.stimType, selectStimType)), mcg) );
        ok_gids = [mcg(indStimTypeMatch).Gid]; %#ok<FNDSB>
        ok_ids = arrayfun(@(gid) any(gid == ok_gids), [allCells.Gid]);
        allCells = allCells(ok_ids);        
    end
    
    
    allGidsX = [allCells.Gid];
    allCellIdsX = [allCells.cellId];            
    nCellsTot = length(allGidsX);        
    
    if nargin < 2
        Gid_start = allGidsX(1);
        cellId_start = allCellIdsX(1);
    end        
    
    cell_idx_start = find(allGidsX == Gid_start & allCellIdsX == cellId_start, 1);    
    if isempty(cell_idx_start)
        error('cell not found');
    end    
    
        
    
    [bins, stimPSTH_vals_oe, stimPSTH_vals, stimPSTH_vals_ordered, stimPSTH_vals_cum, new_order, ...
        nBins, binWidth, binE, nOri, nSp, nPh, cur_psth, ori_sp_ind, ...
        allStats, stat_bin_x, stat_val_y, statVals, statVals_m, statVals_s, targetWindowMode, autoMethod, ...
        frameStimIds, uStimIds, stimIdsIdx, tgtWind_osp_ph ...
        cur_L_bin_wind, cur_R_bin_wind, psthMethod, psth_title_str, ...
        nBinsPerStim, tgtFracBins, nSpikesInPres, fullPsth_bins, phaseCCType, ...
        stimPSTH_vals_allTrials, binOffset] = deal([]);
    [cell_idx, cur_nStim, L_bin, R_bin, L_bin_auto, R_bin_auto, intervalSize, Gid, cellId, statName, nStdTh, objTh, meanRate] = deal(0);    
    [max_p_str, max_t_str, stimType] = deal({});
    
    
    figure(1); clf; 
    [h_cur_psth_ax, h_cur_psth, h_cur_stat] = plotyy(0,0, 0,0, @stairs, @plot); 
    set(h_cur_psth_ax, 'nextplot', 'add'); %"hold on" only does first axes.
    [h_cur_stat_sm] = plot(h_cur_psth_ax(2), 0,0, 'k', 0,0, 'k*', 0,0, 'gs');    
    set(h_cur_psth_ax, 'activeposition', 'outer');
    set(h_cur_psth, 'color', 'r');
    h_psth_tit = title('ABC ', 'fontsize', 13);  
    ax_pos = get(h_cur_psth_ax(1), 'position');
    h_psth_ann = annotation('textbox', [ax_pos(1)+.01, ax_pos(2)+ax_pos(4)*.8, .1, .1], 'string', 'X', 'edgecolor', 'none');
    xlabel('  ');  
    box off;
    
    if doPsth1CurAll
        figure(2); clf; hold on;        
        h_cum_psth_1   = stairs(0, 0, 'r'); hold on;
        h_cum_psth_end = stairs(0, 0, 'k:');
        h_cum_psth_n   = stairs(0, 0, 'b', 'linewidth', 2); %#ok<NASGU>
        h_cum_psth_ax = gca;    
    end    
%     h_bst_wgts = stairs(0,0, 'g');
    % placeholder figure
%     figure(3);
%     figure(4);

    L0 = [0;0];
    h_l_margin = line( L0 , L0, 'parent', h_cur_psth_ax(1), 'linewidth', 2);
    h_r_margin = line( L0 , L0, 'parent', h_cur_psth_ax(1), 'linewidth', 2);
    h_l_margin_full = line( L0 , L0, 'parent', h_cur_psth_ax(1), 'linewidth', 2);
    h_r_margin_full = line( L0 , L0, 'parent', h_cur_psth_ax(1), 'linewidth', 2);
    h_l_margin_auto = line( L0, L0, 'color', 'm', 'linestyle', ':', 'linewidth', 3, 'parent', h_cur_psth_ax(1));
    h_r_margin_auto = line( L0, L0, 'color', 'm', 'linestyle', ':', 'linewidth', 3, 'parent', h_cur_psth_ax(1));
    show_saved = iff(recordPreferences, 'on', 'off');
    h_l_margin_saved = line( L0, L0, 'color', 'g', 'linestyle', ':', 'parent', h_cur_psth_ax(1), 'visible', show_saved);
    h_r_margin_saved = line( L0, L0, 'color', 'g', 'linestyle', ':', 'parent', h_cur_psth_ax(1), 'visible', show_saved);
    h_std_lo = line( L0, L0, 'color', 'g', 'linestyle', ':', 'parent', h_cur_psth_ax(2));
    h_mean   = line( L0, L0, 'color', 'g', 'linestyle', ':', 'parent', h_cur_psth_ax(2));
    h_std_hi = line( L0, L0, 'color', 'g', 'linestyle', ':', 'parent', h_cur_psth_ax(2));

    tmp = zeros(36, 10, 8);    
    if doWhichStim
        figure(3);    
        ext_rep   = @(x) [x(:); x(end)];        
        [h, h_count_im] = imageOSP( tmp, 'mean:ph' );
        title('Stimuli in PSTH')    
    end    

%     allStatNames = {'rho', 'rho_t', 'rho_p',   'tau', 'tau_t', 'tau_p',  'rep_tstat', 'rep_p',   'r_var', 'r_entropy'};
%     allStatNames = {'rho', 'rho_t', 'rho_p',   'rep_tstat', 'rep_p',   'r_var', 'r_entropy'};
%     allStatNames = {'cc_p', 'rho_p', 'rho_p_nz', 'rho_p_nznz', 'r_entropy'};
    allStatNames = {'cc_p'};
    nstats = length(allStatNames);
    
%     allPsthMethods = {'stimulus', 'spikes'};
    allPsthMethods = {'stimulus'};

    
    if showStatImages
        [h_stat_im, h_stat_im_ax, h_stat_im_mx, h_stat_im_cur, h_stat_im_tit] = deal( zeros(1,nstats) );
        for stat_j = 1:nstats
            figure(400+stat_j); clf;
            h_stat_im(stat_j) = imagesc(zeros(2));  hold on;
            h_stat_im_ax(stat_j) = gca;
    %         h_stat_im_mx(stat_j) = plot(0,0, 'wo', 'markersize', 6, 'markerfacecolor', 'r');
            h_stat_im_cur(stat_j) = plot(0,0, 'wo', 'markersize', 6, 'markerfacecolor', 'b');        
            h_stat_im_tit(stat_j) = title('', 'interpreter', 'none');
    %         h_stat_int(stat_j) = plot(0,0, '-o', marker(stat_j), 'markersize', 3, 'markerfacecolor', 'b', 'parent', h_cur_psth_ax);
            xlabel('Left bin'); ylabel('right bin');    
            axis tight xy;
            set(h_stat_im_ax(stat_j), 'xlimMode', 'auto', 'ylimMode', 'auto');
            colorbar;
        end
    end
    
    getIdx = @(gid, cellid) find(gid == allGidsX & cellid == allCellIdsX,1); %#ok<NASGU>

    
%     figure(402); clf;
%     h_rep_t_im = imagesc(zeros(2));  hold on;
%     h_rep_t_im_ax = gca;
%     h_rep_t_im_mx = plot(0,0, 'wo', 'markersize', 6, 'markerfacecolor', 'r');
%     h_rep_t_im_cur = plot(0,0, 'wo', 'markersize', 6, 'markerfacecolor', 'b');
%     h_rep_t_im_tit = title('');
%     xlabel('Left bin'); ylabel('right bin');
%     axis tight xy;
%     colorbar;
    

%     figure(501); clf;
%     h_area = imagesc(zeros(2));  hold on;
%     h_area_ax = gca;
%     h_area_mx = plot(0,0, 'wo', 'markersize', 6, 'markerfacecolor', 'r');
%     h_area_cur = plot(0,0, 'wo', 'markersize', 6, 'markerfacecolor', 'b');
%     h_area_tit = title('Area');
%     xlabel('Left bin'); ylabel('right bin');
%     axis tight xy;
%     colorbar;

%     figure(502); clf;
%     h_T = imagesc(zeros(2));  hold on;
%     h_T_ax = gca;
%     h_T_mx = plot(0,0, 'wo', 'markersize', 6, 'markerfacecolor', 'r');
%     h_T_cur = plot(0,0, 'wo', 'markersize', 6, 'markerfacecolor', 'b');
%     h_T_tit = title('');
%     xlabel('Left bin'); ylabel('right bin');
%     axis tight xy;
%     colorbar;
    
    
    
%     loc = 'southOutside';

    doFullWindowProfile = true;
    doTgtWindowProfile = true;
    doTgtWgtWindowProfile = false;
    doAllProfilesIn1Fig = true;
    
    doProfiles = [doFullWindowProfile, doTgtWindowProfile, doTgtWgtWindowProfile];
    profileNum = cumsum(doProfiles);

    
    if doFullWindowProfile
        if doAllProfilesIn1Fig
            if (profileNum(1)==1), figure(5); clf; end;
            subplot(1, nnz(doProfiles), profileNum(1));
        else
            figure(5); clf;
        end
        h_fullWind_osp_ax = imageOSP( tmp, 'mean:ph', 'OSP', 'noLabels' ); colorbar();
        title('Full PSTH');  xlabel(' ');    
    end
    
    if doTgtWindowProfile
        if doAllProfilesIn1Fig
            if (profileNum(2)==1), figure(5); clf; end;
            subplot(1, nnz(doProfiles), profileNum(2));
        else
            figure(6); clf;
        end
        
        h_os_wind_ax = imageOSP(tmp, 'mean:ph', 'OSP', 'noLabels'); colorbar();
        title('Tgt PSTH'); xlabel(' ');
    end
    
    if doTgtWgtWindowProfile
        if doAllProfilesIn1Fig
            if (profileNum(3)==1), figure(5); clf; end;
            subplot(1, nnz(doProfiles), profileNum(3));
        else
            figure(7); clf;
        end
        h_os_wgt_wind_ax = imageOSP(tmp, 'mean:ph', 'OSP', 'noLabels'); colorbar();
        title('Wgt PSTH'); xlabel(' ');
    end
    
    if showPhRepAllStim
        figure(8); clf;
        [h_ph_rep_ax, h_ph_rep]= imageOSP(tmp, 'mean:ph', 'OSP', 'noLabels');
        caxis([-1 1]); colorbar;
        title('Phase cc');
    end
    if showPhRepAllStim_sampled
        figure(81); clf; 
        subplot(1,2,1);
        [h_ph_rep_smp_ax, h_ph_rep_smp]= imageOSP(tmp, 'mean:ph', 'OSP', 'noLabels'); 
        caxis([-1 1]); colorbar;        
        h_ph_rep_smp_t = title('Sampled dot');

%         figure(82); clf; 
        subplot(1,2,2);
        [h_ph_rep_wgt_ax, h_ph_rep_wgt]= imageOSP(tmp, 'mean:ph', 'OSP', 'noLabels'); 
        caxis([-1 1]); colorbar;        
        h_ph_rep_wgt_t = title({'Sampled, ', 'weighted dot'});
        h_ph_rep_wgt_x = xlabel(' ');
        
    end
    if showF1oDCallStim 
        figure(83); clf;
        subplot(1,2,3);
        [h_tmp, h_f1odc]= imageOSP(tmp, 'mean:ph', 'OSP', 'noLabels'); 
        caxis([0 2]); colorbar;        
        title('F1/DC');
        h_f1odc_x = xlabel(' ');
        
%         subplot(1,4,4);
%         [h_tmp, h_f1odc_wgt]= imageOSP(tmp, 'mean:ph', 'OSP', 'noLabels'); 
%         caxis([0 2]); colorbar;        
%         title('F1/DC (weighted)');
%         h_f1odc_wgt_x = xlabel(' ');
                
    end
    
    figure(9); clf;
    h_rep = plot(0, 0, 'o');  hold on
    h_oe_title = title('Even vs Odd');
    h_oe_xlab = xlabel('');
    h_oe_ylab = ylabel('');
    h_oe_ax = gca;
    
    
    if doPhaseTC
        figure(10); clf;
        h_phase_tc = plot(0, 0, 'bo-', 0, 0, 'gs:');
        xlim([0 360]);
        ylim([0 inf]);
        set(gca, 'xtick', [0:90:360]);
        h_phase_tc_tit = title(' ');
        h_phase_tc_xlab = xlabel(' ');
    end

    if doOSP8 
        figure(11); clf;
        [h_wgtWind_osp_ax, h_wgtWind_osp_im] = imageOSP(tmp);
    end
    
    if doRvsNStim
       figure(12); clf;
       h_rvsn = plot(0,0, 'bo-', 'markersize', 3);
       th = [1, .9, .8, .7];
       h_rvsn_ax = gca;
       h_rvsn_lines = line(nan(2,length(th)*2), nan(2,length(th)*2), 'linestyle', ':', 'color', 'k');
    end

    if doCheckStimOrder
       figure(15); clf;
       nChks = 3;
       for chk_j = 1:nChks
           subplot(1,nChks,chk_j);
           cso_L(chk_j) = plot(0,0, 'bo', 'markersize', 3); %#ok<AGROW>
           cso_tit(chk_j) = title(' '); %#ok<AGROW>
       end
    end
    
    if recordPreferences
        pref_fig = 16;
       figure(pref_fig); clf;
       hgt = 20;
       spc = 5;
       lab_w = 80;
       LB = [5 5];
       h = 8; s_gc_labl    = uicontrol('units', 'pixels', 'position', [LB+[0,         hgt*(h-1)+spc*(h-2)], 180 hgt], 'String', sprintf('Gid = __. cellId = __'), 'style', 'text' );
       h = 7;                uicontrol('units', 'pixels', 'position', [LB+[0,         hgt*(h-1)+spc*(h-2)], lab_w, hgt], 'String', 'L_bin:', 'style', 'text');
       h = 7; s_L_val      = uicontrol('units', 'pixels', 'position', [LB+[lab_w+spc, hgt*(h-1)+spc*(h-2)], lab_w, hgt], 'String', '3 [4]', 'style', 'text');
       h = 6;                uicontrol('units', 'pixels', 'position', [LB+[0,         hgt*(h-1)+spc*(h-2)], lab_w, hgt], 'String', 'R_bin:', 'style', 'text');
       h = 6; s_R_val      = uicontrol('units', 'pixels', 'position', [LB+[lab_w+spc, hgt*(h-1)+spc*(h-2)], lab_w, hgt], 'String', '3 [4]', 'style', 'text');
       h = 5; 	             uicontrol('units', 'pixels', 'position', [LB+[0,         hgt*(h-1)+spc*(h-2)], lab_w, hgt], 'String', 'intSize:', 'style', 'text');
       h = 5; s_nbin_val   = uicontrol('units', 'pixels', 'position', [LB+[lab_w+spc, hgt*(h-1)+spc*(h-2)], lab_w, hgt], 'String', '3 [4]', 'style', 'text');
       h = 4;                uicontrol('units', 'pixels', 'position', [LB+[0,         hgt*(h-1)+spc*(h-2)], 120 hgt],    'String', ' SAVE Current', 'callback', @saveCurCellPSTHPref, 'parent', pref_fig, 'style', 'pushbutton', 'tag', 'current');
       h = 3;                uicontrol('units', 'pixels', 'position', [LB+[0,         hgt*(h-1)+spc*(h-2)], 120, hgt], 'String', ' SAVE Auto', 'callback', @saveCurCellPSTHPref, 'parent', pref_fig, 'style', 'pushbutton', 'tag', 'auto');       
       h = 2;                uicontrol('units', 'pixels', 'position', [LB+[0,         hgt*(h-1)],            120, hgt], 'String', ' SAVE as NAN', 'callback', @saveCurCellPSTHPref, 'parent', pref_fig, 'style', 'pushbutton', 'tag', 'nan');       
       h = 1;                uicontrol('units', 'pixels', 'position', [LB+[0,         hgt*(h-1)],            120, hgt], 'String', ' SAVE as Zero', 'callback', @saveCurCellPSTHPref, 'parent', pref_fig, 'style', 'pushbutton', 'tag', 'empty');       
                     
    end
    tfToVis = @(tf) iff(tf, 'on', 'off');

%     new_nStim = 10;
%     new_L_bin = 5;
%     new_R_bin = 12;    
    

    function saveCurCellPSTHPref(hObject, eventdata, handles) %#ok<INUSD>
        switch get(gcbo,'Tag')
            case 'current', [lbin, rbin, intSize] = deal(L_bin, R_bin, intervalSize);
            case 'auto',    [lbin, rbin, intSize] = deal(L_bin_auto, R_bin_auto, intervalSize);
                    if (lbin==0), intSize = 0; end
            case 'nan',     [lbin, rbin, intSize] = deal(nan); % = p''''undecided
            case 'empty',   [lbin, rbin, intSize] = deal(0);   % = not reproducible
            
        end
%         cell_idx1 = find(allGidsX == Gid & allCellIdsX == cellId, 1);    
        cell_idx1 = find([allPrefs.Gid] == Gid & [allPrefs.cellId] == cellId, 1);            
        if isempty(cell_idx1)
            cell_idx1 = length(allPrefs)+1;
        end        
        allPrefStruct = struct('Gid', Gid, 'cellId', cellId, 'L_bin', lbin, 'R_bin', rbin, 'intSize', intSize);        
        allPrefs(cell_idx1) = allPrefStruct;
        
        updatePrefFields;
        save(psthWindowPrefsFile, 'allPrefs', '-v6');        
%         fprintf('Saved (now %d saved total)\n',length(allPrefs))
    end


    function [l, r, n] = loadCurCellPSTHPref(Gid_load, cellId_load)
        [l,r,n] = deal([]);
        cell_idx1 = find([allPrefs.Gid] == Gid_load & [allPrefs.cellId] == cellId_load, 1);    
        if ~isempty(cell_idx1)
            s = allPrefs(cell_idx1);
            [l,r,n] = deal(s.L_bin, s.R_bin, s.intSize);
        end
    end

    function updatePrefFields
        [l_pref, r_pref, intSize_pref] = loadCurCellPSTHPref(Gid, cellId);
        set(s_gc_labl, 'String', sprintf('Gid = %d. cellId = %d', Gid, cellId));
        
        set(s_L_val, 'String', sprintf('%d [[ %d ]]', cur_L_bin_wind, l_pref));
        set(s_R_val,  'String', sprintf('%d [[ %d ]]', cur_R_bin_wind, r_pref));
        set(s_nbin_val, 'String', sprintf('%d [[ %d ]]', intervalSize, intSize_pref));
        if (l_pref>0)
            set(h_l_margin_saved, 'xdata', psthWindow(1)+(binWidth * [l_pref-1])*[1;1], 'visible', 'on', 'linestyle', '-.');
            set(h_r_margin_saved, 'xdata', psthWindow(1)+(binWidth * [r_pref  ])*[1;1], 'visible', 'on', 'linestyle', '-.');                
        else
            set([h_l_margin_saved, h_r_margin_saved], 'visible', 'off');             
        end            
    end

    
    function updatePlot(new_cell_idx, newPsthMethod, newStatName, newIntervalSize, new_nStdTh, newObjTh, new_nStim, new_L, new_R, newTargetWindowMode, newAutoMethod, newPhaseCCType, showWindows, showLines)
   
        cellIdxChanged      = cell_idx ~= new_cell_idx;
        psthMethodChanged   = ~strcmp(psthMethod, newPsthMethod);
        statNameChanged     = ~strcmp(statName, newStatName);        
        windowChanged       = any([L_bin R_bin] ~= [new_L new_R]);
        intervalSizeChanged = intervalSize ~= newIntervalSize;
        nStdThChanged       = nStdTh ~= new_nStdTh;
        objThChanged        = objTh ~= newObjTh;
        nStimChanged        = cur_nStim ~= new_nStim ;        
%         targetWindowModeChanged = ~strcmp(targetWindowMode, newTargetWindowMode);
        autoMethodChanged   = ~strcmp(autoMethod, newAutoMethod);
        phaseCCTypeChanged  = ~strcmp(newPhaseCCType, phaseCCType);
                
        [    cell_idx,    psthMethod, statName,    L_bin, R_bin, intervalSize,    cur_nStim,         nStdTh,    targetWindowMode,    autoMethod,   objTh, phaseCCType] = deal(...
         new_cell_idx, newPsthMethod, newStatName, new_L, new_R, newIntervalSize, new_nStim,  new_nStdTh, newTargetWindowMode, newAutoMethod, newObjTh, newPhaseCCType);        

        cphSplitTime_ms = 70;    
     
        nrm = @(x) x/sum(x(:));
        cph_separateBeforeAfterTrials = true;
        tgtWind_osp_ph = 0;
%         ospFuncName = 'max';
        ospFuncName = psthStatsSettings.ospPhCompressFcn;
        
        ospFunc = switchh(ospFuncName, {'max', 'mean'}, {@(osp) max(osp, [], 3), @(osp) mean(osp, 3) } );
        
        [showCurWindow, showFullWindow, showAutoWindow, showSavedWindow] = dealV(showWindows);
        [showStatLines, showStatSmLines, showText] = dealV(showLines);
        
        if cellIdxChanged || psthMethodChanged  % load new cell data
            
            Gid = allGidsX(cell_idx);
            cellId = allCellIdsX(cell_idx);    
            grp = mcg(find(mcg_Gids == Gid, 1));
            stimType = strtok(grp.stimType, 'Movie:Flashed_Gratings:');
            frmLngth_ms = grp.frameLength_ms;
            nSpikes = grp.nSpikes( grp.cellIds == cellId);
            [uori, usp, uph] = dbGetUniqueOriSpPh('Gid', Gid);
            [nOri, nSp, nPh] = deal(length(uori), length(usp), length(uph));
            histOpts = struct('psthWindow_ms', psthWindow, 'psthMethod', psthMethod, 'trialGrouping', 'odd/even');
            [bins, stimPSTH_vals_oe, meanRate] = dbGetCellSpkStimHists(Gid, cellId, histOpts);    
            binWidth = diff(bins(1:2));
            nBins = length(bins);            
            if cph_separateBeforeAfterTrials && (size(stimPSTH_vals_oe, 3) == 4)
                bins_before =  ( bins <  cphSplitTime_ms );
                bins_after  =  ( bins >= cphSplitTime_ms );                
                stimPSTH_vals_oe = cat(1, stimPSTH_vals_oe(bins_before, :, [1 2], :), ...
                                          stimPSTH_vals_oe(bins_after,  :, [3 4], :) );
            end                  
            stimPSTH_vals = mean(stimPSTH_vals_oe, 3);            
            nBinsPerStim = round(frmLngth_ms/binWidth);

            idx0 = find(bins > 0, 1);
            stimWindow_osp = sum(stimPSTH_vals(idx0 + [0:nBinsPerStim-1],:),1);
            v = 1;%getValFor1Spk(stimWindow_osp);    
            nSpikesInPres = sum(stimWindow_osp(:)/v); 
            
            psthCycle = mean(stimPSTH_vals(1:nBinsPerStim,:), 2);
%             stimInducedF1oDC = getF1oDC( meanAllPsths );
            stimInducedF1oDC = (max(psthCycle)-min(psthCycle))/min(psthCycle);

%             psth_title_str =  sprintf('(( %s: %d )) Cell %d:%d : %s. t = %.0f ms. Nspk = %d [%d]. <r> = %.4f.', selectStimType, cell_idx, Gid, cellId, stimType, frmLngth_ms, nSpikesInPres, nSpikes, meanRate );
            
            psth_title_str =  sprintf('( %s# %d ) Cell %d:%d : %s. t = %.0f ms. <r> = %.4f.', selectStimType, cell_idx, Gid, cellId, stimType, frmLngth_ms, meanRate );
            
%             'F1/DC = %.2f', stimInducedF1oDC
            set(h_psth_tit, 'string', psth_title_str, 'fontsize', 11 );            
            
            fullPsthWindow = [20, 120];
            l_full = find( bins > fullPsthWindow(1), 1, 'first');
            r_full = find( bins < fullPsthWindow(2), 1, 'last');
                        
            fullPsth_bins =  l_full : r_full ;
            psth_means = mean(stimPSTH_vals(fullPsth_bins,:), 1);
%             psth_maxes = max(stimPSTH_vals(fullPsth_bins,:), [], 1);                
%             new_order = ord(psth_means .* psth_maxes, 'descend'); 
            new_order = ord(psth_means, 'descend'); 
            
            stimPSTH_vals_ordered = stimPSTH_vals(:,new_order);
            stimPSTH_vals_cum = cummean(stimPSTH_vals_ordered,2);             
%             nStimMax = 50; %find(stimPSTH_sums / stimPSTH_sums(1) <= nStimTh, 1);                    
            set(h_l_margin_full, 'xdata', psthWindow(1)+(binWidth * [l_full-1])*[1;1], 'visible', 'on', 'linestyle', '-.');
            set(h_r_margin_full, 'xdata', psthWindow(1)+(binWidth * [r_full  ])*[1;1], 'visible', 'on', 'linestyle', '-.');                
            
            
            doRvsNStim = false;
            if doRvsNStim
                nSums = 200;
                psth_sums = sum( stimPSTH_vals_ordered(:,1:nSums), 1);                
                set(h_rvsn, 'xdata', 1:nSums, 'ydata', psth_sums);                        
                Lx = cell(length(th)*2, 1);
                Lx(:) = {nan(1,2)}; Ly = Lx;            
                set(h_rvsn_lines, 'visible', 'off');
                ylims = get(h_rvsn_ax, 'ylim');
                for th_i = 1:length(th)
                    idx = find(psth_sums / psth_sums(1) <= th(th_i), 1);
                    if isempty(idx), break; end
                    Lx(2*th_i + [-1,0]) = {[0, idx], [idx, idx]};
                    Ly(2*th_i + [-1,0] ) = {psth_sums(idx) * [1,1], [psth_sums(idx), ylims(1)]};
                end                    
                set(h_rvsn_lines, {'xdata'}, Lx, {'ydata'}, Ly);
                set(h_rvsn_lines, 'visible', 'on');
            end            
                                                
%             osp_max = reshape(max(stimPSTH_vals,[], 1), [nOri, nSp, nPh]);
            
            allStats = getPSTHwindowData(Gid, cellId, allStatNames, psthWindow, psthMethod);                
            for j = 1:nstats
                if ~isempty(strfind(allStatNames{j}, '_p'))
%                     allStats.(allStatNames{j}) = pval2NegLogPval( allStats.(allStatNames{j}) );
                    allStats.(allStatNames{j}) = fixNegLogPval( allStats.(allStatNames{j}) );
                end                    
%                 if strcmp(allStatNames{j}, 'r_entropy')
%                     allStats.(allStatNames{j}) = -( allStats.(allStatNames{j}) );
%                 end                                    
            end
                        
%             [bins_allTrials, stimPSTH_vals_allTrials] = dbGetCellSpkStimHists(Gid, cellId, struct('trialGrouping', 'individual'));
%             binOffset = find( abs(bins - bins_allTrials(1)) < 1e-5, 1) - 1;
%             nTrials = size(stimPSTH_vals_allTrials,3);
%             stimPSTH_vals_allTrials = reshape(stimPSTH_vals_allTrials, [length(bins_allTrials), nOri, nSp, nPh, nTrials]);
            
            
            fullWindow_osp_ph = reshape(mean(stimPSTH_vals(fullPsth_bins,:),1), [nOri, nSp, nPh]);
            fullWindow_osp_ph = fullWindow_osp_ph * meanRate/mean(fullWindow_osp_ph(:));
            fullWindow_osp = ospFunc(fullWindow_osp_ph);    
            fullFracBins = nnz(fullPsth_bins)/nBinsPerStim; 
            updateOspFigure(h_fullWind_osp_ax, fullWindow_osp_ph, fullWindow_osp, 'Full PSTH', fullFracBins, nSpikesInPres )
            
            binE = binCent2edge(bins);            
            
            if doPsth1CurAll
                psth_1 = ext_rep( nrm (stimPSTH_vals_cum(:,1) ) );
                psth_end = ext_rep( nrm ( stimPSTH_vals_cum(:,end) ) );
                set(h_cum_psth_1,   'xdata', binE, 'ydata', psth_1);
                set(h_cum_psth_end, 'xdata', binE, 'ydata', psth_end );
                set(h_cum_psth_ax, 'xlim', binE([1, end]));                
            end
            
            set(h_cur_psth_ax, 'xlim', binE([1, end]));
            
%             if cur_nStim > 0
%                 cur_psth = stimPSTH_vals_cum(:,cur_nStim);
% %                 set(h_cum_psth_n,   'xdata', binE, 'ydata', ext_rep( cur_psth / sum(cur_psth) ) );            
%             end
            if showStatImages
                for stat_i = 1:nstats
                    data = allStats.(allStatNames{stat_i});
                    idx_gt0 = find(bins > 0);
                    binEdges = binCent2edge(bins);
                    xylims = [binEdges(idx_gt0(1)), binEdges(end)];
                    set(h_stat_im(stat_i), 'xdata', bins(idx_gt0), 'ydata', bins(idx_gt0), 'cdata', data(idx_gt0,idx_gt0));
                    set(h_stat_im_ax(stat_i), 'xlim', xylims, 'ylim', xylims);
%                     set(h_stat_im_tit(stat_i), 'string', sprintf('%s [Cur: %.2f]', allStatNames{stat_i}, data(R_bin, L_bin) ) );
                    set(h_stat_im_tit(stat_i), 'string', sprintf('%s', allStatNames{stat_i} ), 'interpreter', 'tex' );
                end    
            end
            
        end                
                
                
        if nStimChanged || cellIdxChanged || psthMethodChanged  
            if cur_nStim > 0
                nStim1 = min(cur_nStim, size(stimPSTH_vals_cum,2));
                cur_psth = stimPSTH_vals_cum(:,nStim1);
                cur_psth_ext = ext_rep( (cur_psth) );
                set(h_cur_psth, 'xdata', binE, 'ydata', cur_psth_ext );                            
            end
            
            h_lr = [h_l_margin, h_r_margin, h_l_margin_auto, h_r_margin_auto h_l_margin_saved, h_r_margin_saved, h_l_margin_full, h_r_margin_full];
            ylims1 = bestLimsFromData( cur_psth );
            ylims1 = [0 ylims1(2)];
            set(h_lr, 'ydata', ylims1);          
            set(h_cur_psth_ax(1), 'ylim', ylims1, 'ytickMode', 'auto')            
            if doWhichStim
                cnt = zeros(nOri, nSp, nPh);
                cnt(new_order(1:cur_nStim)) = 1;
                set(h_count_im, 'cdata', sum(cnt,3) );            
            end                
            
        end 
        
        L_bin_ms = psthWindow(1) + (binWidth * [L_bin-1]);
        R_bin_ms = psthWindow(1) + (binWidth * [R_bin  ]);
        if windowChanged                         
            set(h_l_margin, 'xdata', L_bin_ms*[1;1]);
            set(h_r_margin, 'xdata', R_bin_ms*[1;1]);            
            if showStatImages
                set(h_stat_im_cur, 'xdata', bins(L_bin), 'ydata', bins(R_bin), 'markersize', 3);            
            end
        end
        
        
        updatePsthStatsAndWindow = autoMethodChanged || cellIdxChanged || psthMethodChanged  || intervalSizeChanged || statNameChanged || nStdThChanged || objThChanged;
%         updatePsthStatStdLines = updatePsthStats || nStdThChanged;                
%         updateAutoWindow       = updatePsthStats;        
                

        if updatePsthStatsAndWindow
            [L_bin_auto, R_bin_auto, stat_bin_x, stat_val_y, statVals_m, statVals_s, tmp, comment, allStats] = ...
                getLRbin(Gid, cellId, bins, allStats, intervalSize, statName, stimType, nStdTh, objTh, meanRate);
            
            set(h_cur_stat, 'xdata', stat_bin_x, 'ydata', stat_val_y, 'marker', 'o', 'markersize', 3    );                        
            
            stat_sm = gaussSmooth(stat_val_y, 1);
            id_locmin_y = findLocalMinima(stat_val_y);
            id_locmin_sm = findLocalMinima(stat_sm, 1);
%             id_ba = union(id_locmin_sm, id_locmin_y);
%             id_b  = setdiff(id_locmin_sm, id_locmin_y);
            set(h_cur_stat_sm(1), 'xdata', stat_bin_x, 'ydata', stat_sm, 'marker', '.');
            set(h_cur_stat_sm(2), 'xdata', stat_bin_x(id_locmin_sm), 'ydata', stat_sm(id_locmin_sm), 'marker', 'd', 'markerfacecolor', 'r');
            set(h_cur_stat_sm(3), 'xdata', stat_bin_x(id_locmin_y), 'ydata', stat_val_y(id_locmin_y), 'marker', 'o', 'color', 'r', 'markerfacecolor', 'none');
            set(h_psth_ann, 'string', comment );                    

            if ~all(isnan(stat_val_y))
                ylims2 = bestLimsFromData(stat_val_y);
                set(h_cur_psth_ax(2), 'ylim', ylims2, 'ytickMode', 'auto')                        
                set(h_mean, 'xdata',   [bins(1) bins(end)], 'ydata', statVals_m * [1,1], 'color', 'g', 'linestyle', '-');                
        
                % update stdev lines of current statistic.    
                set(h_std_lo, 'xdata', psthWindow, 'ydata', (statVals_m - statVals_s*nStdTh) * [1,1], 'color', 'g', 'linestyle', ':');
                set(h_std_hi, 'xdata', psthWindow, 'ydata', (statVals_m + statVals_s*nStdTh) * [1,1], 'color', 'g', 'linestyle', ':');                                    
            end                
                
            if (L_bin_auto>0) 
                set(h_l_margin_auto, 'xdata', psthWindow(1) + (binWidth * [L_bin_auto-1])*[1;1], 'visible', 'on', 'color', 'm', 'linewidth', 3);
                set(h_r_margin_auto, 'xdata', psthWindow(1) + (binWidth * [R_bin_auto  ])*[1;1], 'visible', 'on', 'color', 'm', 'linewidth', 3);            
            else
                set([h_l_margin_auto, h_r_margin_auto], 'visible', 'off')
                [L_bin_auto, R_bin_auto] = deal(0);
            end                            
%             set(h_area_tit, 'string', sprintf('Area'));            
        end

        switch targetWindowMode
            case 'auto',  [L_bin_wind, R_bin_wind] = deal(L_bin_auto, R_bin_auto);
%                 if (L_bin_auto==0)
%                     [L_bin_wind, R_bin_wind] = deal(L_bin, R_bin); 
%                 end
            case 'manual',[L_bin_wind, R_bin_wind] = deal(L_bin, R_bin); 
            case 'saved', [L_bin_wind, R_bin_wind] = loadCurCellPSTHPref(Gid, cellId); 
        end
%         L_bin_wind_ms = psthWindow(1) + (binWidth * [L_bin_wind-1]);
%         R_bin_wind_ms = psthWindow(1) + (binWidth * [R_bin_wind  ]);
        
        
        targetWindowChanged = isempty(cur_L_bin_wind) || any([cur_L_bin_wind, cur_R_bin_wind] ~= [L_bin_wind, R_bin_wind]);
        [cur_L_bin_wind, cur_R_bin_wind] = deal(L_bin_wind, R_bin_wind);
        
        if (targetWindowChanged || cellIdxChanged || psthMethodChanged  || statNameChanged || nStdThChanged || phaseCCTypeChanged) 
            
            % update OSP (from targeted window)
            isRealWindow = (L_bin_wind <= R_bin_wind) && ~isnan(L_bin_wind) && (L_bin_wind>0);
            if isRealWindow
                psths_wind = stimPSTH_vals(L_bin_wind:R_bin_wind,:);    
                tgtWind_osp_ph = reshape( mean ( psths_wind , 1), [nOri, nSp, nPh]);
                tgtWind_osp_ph = tgtWind_osp_ph * meanRate/mean(tgtWind_osp_ph(:));
                tgtWind_osp = ospFunc(tgtWind_osp_ph);                                
                
            else
                tgtWind_osp_ph = [];
                tgtWind_osp    = nan(nOri, nSp);                
            end
            tgtFracBins = (R_bin_wind-L_bin_wind+1)/nBinsPerStim;
            updateOspFigure(h_os_wind_ax, tgtWind_osp_ph, tgtWind_osp, 'Tgt PSTH', tgtFracBins, nSpikesInPres, Gid, cellId, cell_idx )
            
            % Plot odd vs. even trials
            if isRealWindow
                [L_bin_wind_oe, R_bin_wind_oe] = deal(L_bin_wind, R_bin_wind); oe_label = 'Tgt Window';
            else
                [L_bin_wind_oe, R_bin_wind_oe] = deal(fullPsth_bins(1), fullPsth_bins(end)); oe_label = 'Full Window';
            end
                        
            switch phaseCCType
                case 'cc', ccStatName = 'phaseTC_CCs'; clims = [-1 1];
                case 'dot', ccStatName = 'phaseTC_Dots'; clims = [0 1];
            end
                        
            
            if true || isRealWindow
%                 cur_psth_rep = iff(cur_nStim > 0, cur_psth, []);
                shiftCent = false;                
                [curStat, r, r_odd, r_even] = getOspStatVsPsthBinning(Gid, bins, stimPSTH_vals_oe, L_bin_wind_oe, R_bin_wind_oe, shiftCent, [], {statName} );
%                 [tmp, rep_p1] = corr(r_odd(:), r_even(:), 'tail', 'gt');
%                 [tmp, rep_p2a] = regressionSlopeTtest(r_odd(:), r_even(:), .01, '+');
%                 [tmp, rep_p2b] = regressionSlopeTtest(r_even(:), r_odd(:), .01, '+');
%                 mxdiff = max(abs([ rep_p1-rep_p2a,  rep_p2a-rep_p2b, rep_p2a-rep_p1]));
%                 fprintf('%.4g,  %.4g, %.4g. [%.4g]\n', rep_p1, rep_p2a, rep_p2b, mxdiff);
                v = 1; %getValFor1Spk([r_odd(:); r_even(:)]);
%                 oe_stats = getOEstats(r_odd, r_even);
                [r_ph, r_odd_ph, r_even_ph] = calcOspForPsthWindow(Gid, {bins, stimPSTH_vals_oe}, L_bin_wind_oe, R_bin_wind_oe, shiftCent, [], 'phase');
                if showPhRepAllStim && ishandle(h_ph_rep);
                    
                    [ph_ccs, ph_cc_ps] = getAllPhaseTCccs(r_odd_ph, r_even_ph);
                    set(h_ph_rep, 'cdata', ph_ccs);
                    set(h_ph_rep_ax, 'clim', clims);
                end
                

                if showPhRepAllStim_sampled && ishandle(h_ph_rep_smp);
%                     'phaseTC_CCs', 'phaseTC_Dots'
                    ph_ccs_smp = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_wind_oe, R_bin_wind_oe, [], ccStatName);
%                     [ph_ccs_smp, ph_cc_ps_smp] = getAllPhaseTcCCs_sampled(stimPSTH_vals_allTrials, L_bin_wind_oe-binOffset, R_bin_wind_oe-binOffset);
                    set(h_ph_rep_smp, 'cdata', ph_ccs_smp);
                    
                    tgtWind_osp_tmp = tgtWind_osp;
                    if ~isRealWindow
                        tgtWind_osp_tmp = ospFunc(r_ph);
                    end
                    
                    tgtWind_osp_norm = tgtWind_osp_tmp / max(tgtWind_osp_tmp(:));
                    ph_ccs_wgt = ph_ccs_smp .* tgtWind_osp_norm;
                    set(h_ph_rep_wgt, 'cdata', ph_ccs_wgt);
                    cc_ptcc_osp = corr(ph_ccs_wgt(:), tgtWind_osp_norm(:)); 
                    set(h_ph_rep_wgt_x, 'string', sprintf('mean : %.2g.  cc: %.2g', mean(ph_ccs_wgt(:)), cc_ptcc_osp )  );                    
                    
                    set(h_ph_rep_smp_t, 'string', ['Sampled ' phaseCCType]);
                    set(h_ph_rep_wgt_t, 'string', {'Sampled, ', ['weighted ' phaseCCType ]});
                    set([h_ph_rep_smp_ax h_ph_rep_wgt_ax], 'clim', clims);
                end
                if showF1oDCallStim && ishandle(h_f1odc);
                    f1odcs = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_wind_oe, R_bin_wind_oe, [], 'stimF1oDCs');                    
                    f1odc_max = f1odcOfMax(f1odcs, tgtWind_osp_norm, 0.8);
                    
                    set(h_f1odc, 'cdata', f1odcs);
                    set(h_f1odc_x, 'string', sprintf('f1odc : %.2f', f1odc_max ) );                    
                    
%                     f1odcs_wgt = f1odcs .* tgtWind_osp_norm;
%                     set(h_f1odc_wgt, 'cdata', f1odcs_wgt);
                end
                
                
                
                curStat = iff(~isempty(strfind(statName, '_p')), 10.^(-curStat), curStat);
                rand('state', 0);
                        [r_odd, r_even] = nudgeOddEvenTrials(r_odd / v,r_even / v);
                set(h_rep, 'xdata', r_odd, 'ydata', r_even, 'marker', '.');                
                L = max([r_odd(:); r_even(:)])+1;
                set(h_oe_ax, 'xlim', [0 L], 'ylim', [0 L]);                
                statStr = sprintf('%s = %.3g', statName, curStat );
%                 oeNums = sprintf('z: %d, h: %d, f: %d', oe_stats);
                set(h_oe_title, 'string', statStr, 'interpreter', 'none');
%                 set(h_oe_xlab, 'string', oe_label);
                set(h_oe_xlab, 'string', 'odd # trials');
                set(h_oe_ylab, 'string', 'even # trials');
            end

                           
            
        end
            
        
        if (targetWindowChanged || cellIdxChanged || psthMethodChanged || statNameChanged || nStimChanged || nStdThChanged)                
        
            isRealWindow = (L_bin_wind <= R_bin_wind) && ~isnan(L_bin_wind) && (L_bin_wind>0);
            if isRealWindow                                
                targetedPSTH = cur_psth(L_bin_wind:R_bin_wind);
                targetedPSTH = targetedPSTH*length(targetedPSTH)/sum(targetedPSTH);
                isRealWindow = ~isnan(targetedPSTH(1));
            end            
            
            %update Targeted & Weighted OSP (& Phase Tuning Curve)            
            if isRealWindow                                                
                psth_wgt_wind_oe = bsxfun(@times, stimPSTH_vals_oe(L_bin_wind:R_bin_wind,:), targetedPSTH);
                wgtWind_osp_ph_oe = reshape( mean ( psth_wgt_wind_oe, 1), [nOri, nSp, nPh, 2]);    
                wgtWind_osp_ph_oe = wgtWind_osp_ph_oe * meanRate/mean(tgtWind_osp_ph(:));
                wgtWind_osp_ph = mean(wgtWind_osp_ph_oe, 4);
                wgtWind_osp    = ospFunc(wgtWind_osp_ph);
            else
                wgtWind_osp_ph = [];
                wgtWind_osp    = nan(nOri, nSp);                
            end                  
            if doTgtWgtWindowProfile
                updateOspFigure(h_os_wgt_wind_ax, wgtWind_osp_ph, wgtWind_osp, 'Wgt PSTH', tgtFracBins, nSpikesInPres )
            end
            

            if isRealWindow && doPhaseTC                    
                ph = linspace(0,360, nPh+1);
                [tmp, ori_sp_ind] = maxElement(max(wgtWind_osp,[], 3));
                phase_tc_odd = squeeze( wgtWind_osp_ph_oe(ori_sp_ind(1), ori_sp_ind(2), :, 1) );
                phase_tc_even = squeeze( wgtWind_osp_ph_oe(ori_sp_ind(1), ori_sp_ind(2), :, 2) );                    
%                     phase_tc = (phase_tc_odd + phase_tc_even) / 2;
                set(h_phase_tc(1), 'xdata', ph, 'ydata', wrp(phase_tc_odd), 'color', 'b');
                set(h_phase_tc(2), 'xdata', ph, 'ydata', wrp(phase_tc_even), 'color', 'm');
                [phase_tc_r, phase_tc_r_p] = corr(phase_tc_odd, phase_tc_even, 'type', 'pearson', 'tail', 'gt');
                [phase_tc_rho, phase_tc_rho_p] = corr(phase_tc_odd, phase_tc_even, 'type', 'spearman', 'tail', 'gt');                    
                set(h_phase_tc_tit, 'string', {sprintf('r = %.2f. \\rho = %.2f', phase_tc_r, phase_tc_rho), sprintf('p_r = %.3f. p_{\\rho} = %.3f', phase_tc_r_p, phase_tc_rho_p)} );
                set(h_phase_tc_xlab, 'string', sprintf('ori = %d, sp = %d', ori_sp_ind(1), ori_sp_ind(2)));
            end
            if isRealWindow && doOSP8 
                clims = [0, max(wgtWind_osp_ph(:))];
                for ph_i = 1:nPh
                    set(h_wgtWind_osp_im(ph_i), 'cdata', wgtWind_osp_ph(:,:,ph_i));
                    set(h_wgtWind_osp_ax(ph_i), 'clim', clims);
                end
                for ph_i = nPh+1:8;
                    set(h_wgtWind_osp_im(ph_i), 'cdata', zeros(nOri, nSp));
                end

            end
            
        end
      
        
        
        if recordPreferences
            updatePrefFields;
            
%             [l_saved, r_saved] = loadCurCellPSTHPref(Gid, cellId);
%             if ~isempty(l_saved)
%                 L_bin_wind_ms = psthWindow(1) + (binWidth * [L_bin_wind-1]);
%                 R_bin_wind_ms = psthWindow(1) + (binWidth * [R_bin_wind  ]);
%                 
%                 
%             else
%                 
%             end
            
        end

        
        
        if doCheckStimOrder && (L_bin <= R_bin)
            idx_relevant = (bins < 120);
            nStimMax = 50;
            psth_means = mean(stimPSTH_vals(idx_relevant,:), 1);
            psth_maxes = max(stimPSTH_vals(idx_relevant,:), [], 1);
            
            stim_wgt{1} = psth_means;
            stim_wgt{2} = psth_maxes;
            stim_wgt{3} = psth_means .* psth_maxes;
            
            stimPSTH_vals_ordered = stimPSTH_vals(:,new_order);
            stimPSTH_vals_cum = cummean(stimPSTH_vals_ordered,2);
            stim_th = .5;
            idx_top_stim = find(wgtWind_osp(:) > stim_th*max(wgtWind_osp(:)) );
            actual_resp = wgtWind_osp(idx_top_stim);
            
            stim_order = cell(1,nChks);
            for chk_i = 1:nChks
                stim_order{chk_i} = ord(stim_wgt{chk_i}, 'descend');
                stimPSTH_vals(:,stim_order{chk_i}(1:nStimMax));
                
                ydata = stim_wgt{chk_i}(idx_top_stim);
                r = pearsonR(actual_resp(:), ydata(:));
                p = spearmanRho(actual_resp(:), ydata(:));
                
                set(cso_L(chk_i), 'xdata', actual_resp, 'ydata', ydata);
                set(cso_tit(chk_i), 'string', sprintf('r = %.2f. \\rho = %.2f', r, p));
            end
        end

        set([h_l_margin, h_r_margin], 'visible', tfToVis(showCurWindow));
        set([h_l_margin_full, h_r_margin_full], 'visible', tfToVis(showFullWindow));        
        set([h_l_margin_auto, h_r_margin_auto], 'visible', tfToVis(showAutoWindow));        
        set([h_l_margin_saved, h_r_margin_saved], 'visible', tfToVis(showSavedWindow));
        set(h_cur_stat_sm, 'visible', tfToVis(showStatSmLines));
        set(h_cur_stat, 'visible', tfToVis(showStatLines));
        set(h_psth_ann, 'visible', tfToVis(showText));
        
        
        
%         if (statNameChanged || intervalSizeChanged || nStdThChanged || objThChanged) && ~isempty(statVarHandles)
%             manipulateSet(statVarHandles, {'statName', 'intSize', 'nStdTh', 'objTh'}, {statName, intSize, nStd, objTh});
%         end
        
    end
        
    tf = [true, false]; ft = [false, true];

    psthAutoSelectMethods = {'std'}; %{'std', 'fracArea'};
    intSize0 = 5;
    nStdTh0 = 3;
    objTh0 = 3;
    nStim0 = 20;
    
    toShowWind = {'Current', 'Full', 'Auto', 'Saved'};
    toShowWindVals = {tf, tf, tf, tf};
    toShowWindVals0 = [true, true, true, false];

    toShowLines = {'stat', 'stat_smoothed', 'text'}; 
    toShowLinesVals = {tf, tf, tf};
    toShowLinesVals0 = [true, false, false];
    
    updatePlot(cell_idx_start, allPsthMethods{1}, allStatNames{1}, intSize0, nStdTh0, objTh0, nStim0, 58, 65, 'auto', 'std', 'cc', toShowWindVals0, toShowLinesVals0);    
    
    args = { {'Cell Index', [1:nCellsTot], cell_idx}, ...             
             {'psthMethod', allPsthMethods, psthMethod}, ...
             {'stat', allStatNames}, ...
... %              {'statsShow', repmat(tf, nstats,1), [], [], allStatNames};
             {'intSize', [1:nBins-1], intervalSize}, ...
             {'nStdTh', [.1:.1:4], nStdTh}, ...
             {'objTh', [.1:.1:4], objTh}, ...
             {'nStim', [0:2880], cur_nStim}, ...
             {'L_margin', [1:nBins], L_bin}, ...
             {'R_margin', [1:nBins], R_bin}, ...             
             {'targetWindow', {'auto', 'manual', 'saved'}, targetWindowMode}, ...
             {'selectMethod', psthAutoSelectMethods}, ...
             {'phaseTC', {'cc', 'dot'}}, ...
             {'showWind', toShowWindVals, toShowWindVals0, [], toShowWind}, ...
             {'showLine', toShowLinesVals, toShowLinesVals0, [], toShowLines}, ...
             };

    psthVarHandles = manipulate(@updatePlot, args, 'FigId', 112);
    

    
end



function updateOspFigure(h_ax, osp_ph, osp, mainTitle, fracBins, nSpkTotal, Gid, cellId, cell_idx)
    h_title = get(h_ax, 'title');
    h_xlabel = get(h_ax, 'xlabel');
    h_im = get(h_ax, 'children');
    set(h_im, 'cdata', osp);
    if nargin >= 7
        gc_str = '';%sprintf('cell %d : %d (%d)', Gid, cellId, cell_idx);
    else
        gc_str = '';
    end
    if all(isnan(osp))
        [osp_noisiness, osp_nunique, osp_frac] = deal(nan);
    else
        osp_noisiness = ospNoisiness(osp);           
        osp_nunique   = length(unique(osp(:)));      
    %     osp_frac      = sum(osp(:))/sum(osp_ph(:));  
        if ~strcmp(mainTitle, 'Wgt PSTH')
            v = 1;%getValFor1Spk(osp_ph);    
%             osp_nspk      = sum(osp_ph(:)/v); 
%             osp_frac = (osp_nspk / nSpkTotal) / fracBins;
%             v2 = getValFor1Spk(osp);    
%             assert(abs(v-v2)<1e-5);
            osp_nspk2 = max(osp(:)/v); 
            osp_frac = (osp_nspk2);
        else
            osp_frac = 0;
        end
    end
    
    noise_str = sprintf('Ns: %.2f', osp_noisiness );
    nunique_str = sprintf('U: %.2f', osp_nunique);
    frac_str = sprintf('Frac: %.2f', osp_frac);
    
    title_str = iff(isempty(gc_str), mainTitle, {gc_str, mainTitle});
    set(h_title, 'string', title_str);
%     set(h_xlabel, 'string', {noise_str, nunique_str, frac_str}  );
    set(h_xlabel, 'string', '' );
end

% function y = ff(x)
%     
%     
%     
% 
% 
% end
% 


%                 window_auto_ms = getBestTimeWindowFromPSTH(bins, cur_psth);
%                 L_bin_auto = find(bins > window_auto_ms(1), 1, 'first');
%                 R_bin_auto = find(bins < window_auto_ms(2), 1, 'last');


% 24) 2346:5(10x1) 41ms
% 25) 2771:2(2x8)  16ms
% 26) 5226:4(16x1) 16ms


% figure(1); imagesc(s.rho); axis xy equal tight; colorbar; title('\rho')
% figure(2); imagesc(s.tau); axis xy equal tight; colorbar; title('\tau')
% figure(3); imagesc(s.rep_tstat); axis xy equal tight; colorbar; title('rep_t')
% figure(4); imagesc(s.rvar); axis xy equal tight; colorbar; title('var')
% figure(11); imagesc(pval2NegLogPval(s.rho_p)); axis xy equal tight; colorbar; title('\rho - pval')
% figure(12); imagesc(pval2NegLogPval(s.tau_p)); axis xy equal tight; colorbar; title('\tau - pval')
% figure(13); imagesc(pval2NegLogPval(s.rep_p)); axis xy equal tight; colorbar; title('rep - pval')
% figure(21); plot(s.rho(:), s.tau(:), '.'); axis equal; xlabel('rho'); ylabel('tau');
% figure(22); plot(pval2NegLogPval(s.rho_p(:)), pval2NegLogPval(s.tau_p(:)), '.'); axis equal; xlabel('rho_p'); ylabel('tau_p');
% figure(22); plot(pval2NegLogPval(s.tau(:)), pval2NegLogPval(s.tau_p(:)), '.'); axis equal; xlabel('rho_p'); ylabel('tau_p');


% rho_p vs H
% scatterplot: peak bin, or bin width, color coded by max rho_p


% Let r', h' be zero-mean versions of 
% r, h; 
% compute covariance matrix  C=( <r'^2>,  <r' h'>; <r' h'> <h'^2> ); 
% then in a Gaussian approximation, letting v be the vector (r;h) 
% and m the mean of this vector, we have p(r,h)=( 1/(2\pi sqrt(Det(C)) ) e^{ -.5 (v-m)' C^{-1} (v - m) )


% if doGetParsedSpikes
%     figure(8); clf;
%     [h, h_os_gps] = imageOSP(tmp, 'mean:ph');
%     h_os_tit = title('SpikesView: (#spk: )');
% end
% if doGetParsedSpikes
%     frameStimIds = getStimulusFrameSequence(Gid, 'OSP');
%     [uStimIds, stimIdsIdx] = uniqueList(frameStimIds);
% end
% 
% if doGetParsedSpikes
%     relContrOfFrameToSpike = getParsedSpikes('frame', Gid, cellId, [L_bin_wind_ms, R_bin_wind_ms]);
%     relContrOfFrameToSpike = [relContrOfFrameToSpike{:}];
%     nSpks = sum(relContrOfFrameToSpike);
%     R_gps = zeros([nOri, nSp, nPh]);
%     for stim_i = 1:length(uStimIds)
%         R_gps(stim_i) = sum(relContrOfFrameToSpike(stimIdsIdx{stim_i}));
%     end
%     set(h_os_gps, 'cdata', mean(R_gps,3));
%     set(h_os_tit, 'string', sprintf('SpikesView: (#spk: %.2f)', nSpks));
%     
% end

function oe_stats = getOEstats(r_odd, r_even)
    z = nnz(     r_odd == 0 & r_even == 0);
    h = nnz( xor(r_odd == 0, r_even == 0) );
    f = nnz(     r_odd > 0 & r_even > 0);
    oe_stats = [z h f];
end

function p = fixNegLogPval(p)    
    % we assigned p==0 a small amount of 1e-100. But this is a huge outlier
    % if the other p values are not so significant. so now we set those
    % pvalues to more intermediate values (close to the original)
    amt = 1e-100;
    defaultAmt_idx = (p == -log10(amt) );
    realMaxP = max(  p(~defaultAmt_idx(:))  );
    p(defaultAmt_idx) = ceil(realMaxP) + 1;
    
    % we also assigned 'logp'-values of -1 to p-values of 1. This can be too
    % negative- we want it to be only slightly below 0 compared with the
    % amounts above zero.
    r = .01;
    maxP = max(p(p ~= 0));
    if maxP < 0, maxP = -1; end
    p(p < 0) = -abs(maxP) * r;        
end


function [ccs, ph_cc_ps] = getAllPhaseTCccs(r1_ph, r2_ph)
    [ccs, ph_cc_ps] = deal( zeros(36,10) ); 
    for oi = 1:36
        for si = 1:10            
            [ccs(oi,si), ph_cc_ps(oi,si)] = doPearsonCorr(r1_ph(oi, si, :), r2_ph(oi, si, :));
        end
    end
end

function y = f1odcOfMax(f1odcs, osp, th)
    idx = (osp(:) > max(osp(:))*th);
    y = mean(f1odcs(idx));
end
