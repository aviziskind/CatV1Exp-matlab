function compareComparisonsAfterSmoothing
    
    gratingType = curGratingType('');    
    smFermi_b_hp = 100;
    
    global smoothPhases
                
    [locations, measures, OSPs, pairData] = deal([]);
    chkLocs = false;
    [nPairs, allNPhases, allPairTypes, pairTypeIdxs, nPhaseIdxs] = deal([]);
    
    
    function cmps = getCmpData(cmpAction, cmpN)
        if strcmp(cmpAction, 'noSmooth')
            smoothPhases = {};
            ext = '';
            dir_ext = '';
        else
            if strcmp(cmpAction, 'Fermi')
                cmpN = [cmpN, cmpN / smFermi_b_hp];
            end
            
            smoothPhases = {cmpAction, cmpN};
            ext = ['_' cmpAction '_' num2str(cmpN(1))];
            dir_ext = 'smooth_test\';            
        end
        
        data_dir = [CatV1Path dir_ext];        
        cmpDatafile = [data_dir gratingType 'GratingComparisonData_DB_tuning_Wcc' ext '.mat'];
                
        if ~exist(cmpDatafile, 'file') 
            gen3;
        end
        
        if isempty(OSPs)             
            ospDatafile_noSmooth = [data_dir gratingType 'GratingCells_DB.mat'];            
            S1 = load(ospDatafile_noSmooth);
            OSPs = S1.allCells;    
        end
        
        S2 = load(cmpDatafile);        
        if isempty(pairData)
            pairDataTmp = S2.pairData;
            oriSp_cmps = {pairDataTmp.oriSp_cmp}';
            pairData = struct('Gids', {pairDataTmp.Gids}', 'cellIds', {pairDataTmp.cellIds}', 'oriSp_cmp', oriSp_cmps);
        end
                        
        cmps = S2.allStatsC;     
        
        if isempty(locations)
            locations = S2.locations;
            measures = S2.measures;
            allNPhases = unique(  arrayfun(@(s) length(s.ph), OSPs) );
            allPairTypes = S2.pairTypes;
            
            Gids = cat(1, pairData.Gids);
            for i = 1:length(allPairTypes)                
                switch allPairTypes{i}
                    case 'Wcc', idxs = ( Gids(:,1) == Gids(:,2) );
                    case 'Bcc', idxs = ( Gids(:,1) ~= Gids(:,2) );
                end
                pairTypeIdxs.(allPairTypes{i}) = idxs;                
            end
            
            n_phases = [pairDataTmp.n_phases]';
            for i = 1:length(allNPhases)                                
                idxs = allNPhases(i) == n_phases;                
                nPhaseIdxs.(['ph' num2str(allNPhases(i))]) = idxs;                
            end            
            nPairs = length(pairDataTmp);
        end
        
    end


    fermiDOFtoK = @(k) k/2;

    if strcmp(gratingType, 'flashed')
        smGaussAmt_dof = [3:.1:7];
        smFermiAmt_dof = [3:.5:7];
        smAveragesAmt_dof = [3:.5:7];
        
        b = [12.9415   -1.3612    2.0189    0.7425    2.9712];
        gaussWtoDOF = @(b, w) max(min( [b(1)./((w-b(2)).^(w.*b(3) + b(4) ) ) + b(5)], 8*ones(size(w)) ), 3*ones(size(w))) ;                        
        gaussDOFtoW = @(dof) arrayfun( @(d) fminsearch(@(w) abs(gaussWtoDOF(b, w)-d), 1), dof);        
        
    elseif strcmp(gratingType, 'drifting')
        b = [31.1195   -0.3064    1.2191    1.5947];

        gaussDOFtoW = @(w)  (b(1)./(w-b(4))).^(1/b(3)) + b(2);
        
        smGaussAmt_dof = [3:16, 16.1:.1:18.5, 19:.5:19.5, 20:2:60, 65];
        smFermiAmt_dof = [2:20, 22:2:60, 65];
        smAveragesAmt_dof = [2:20,30,60];

%     smGaussAmt_dof = [6:10:50];
%     smFermiAmt_dof = [6:10:50];
%     smAveragesAmt_dof = [2:20,30,60];
    end
    
    smGaussAmt = gaussDOFtoW(smGaussAmt_dof);
    smFermiAmt = fermiDOFtoK(smFermiAmt_dof);
    smAveragesAmt = smAveragesAmt_dof;
        
    smoothTypes = {'Gauss', 'Fermi', 'Alias'};    
    smoothAmts = {smGaussAmt, smFermiAmt, smAveragesAmt};    
    smoothAmts_dof = {smGaussAmt_dof, smFermiAmt_dof, smAveragesAmt_dof};    
    smoothLabels = {'Gauss', 'Fermi', 'Alias'};

    smON = logical([1 0 0]);
    smoothTypes = smoothTypes(smON);
    smoothAmts  = smoothAmts(smON);
    smoothAmts_dof = smoothAmts_dof(smON);
    smoothLabels = smoothLabels(smON);    
    
    nSmTypes = nnz(smON);
    
    sm_N = cellfun(@length, smoothAmts);
    sm_cols = 'brg';

    noSm_Cmps = getCmpData('noSmooth', []);    
    
    sm_Cmps = cell(nSmTypes,1);
    for ti = 1:nSmTypes
        sm_Cmps{ti} = cell(1,sm_N(ti));
        for sm_i = 1:sm_N(ti)
            sm_Cmps{ti}{sm_i} = getCmpData(smoothTypes{ti}, smoothAmts{ti}(sm_i));
        end
    end

    function tf = allTheSame(x)
        tf = true;
        for i = 2:length(x)
            if x(i) ~= x(1)
                tf = false;
                return;
            end
        end        
   end
    smoothAmts_x = smoothAmts_dof;
    
    
    
    % check that all oriSp_cmps are the same
    if chkLocs 
        for ti = 1:nSmTypes        
            pd_sm_i = sm_pairData{ti};
            oriSp_cmps_cell = cellfun(@(pd) cat(1, pd.oriSp_cmp),pd_sm_i, 'un', 0);
            oriSp_cmps = cat(3, oriSp_cmps_cell{:}) ;
            [nPairs, nLocs, nSms] = size(oriSp_cmps);
            for p_i = 1:nPairs
                for loc_i = 1:nLocs
                    cur_oriSp_cmp = squeeze( oriSp_cmps(p_i, loc_i, :) )';
                    assert( allTheSame( cur_oriSp_cmp ) )
                end
            end
        end
        fprintf('all locations matched\n');
        return;
    end
    [curPairs,curCmps,curVals] = deal({});
    
      
    msId = 1;
    locId = 1;
    curMeasure = measures{msId};
    curPairIdx = [];
    curSmModeIdx = 1;
    curSmAmtIdx = 1;

    
	[vals_forSmType, vals_forSmType_sum, meanVals_forSmType, meanVals_forSmType_pos, meanVals_forSmType_neg, curPairVals_smType] = ...
        deal(cell(1,nSmTypes));
    [vals_forNoSmooth, meanVal_noSmooth, curPairVal_noSmooth] = deal([]);    
    amtIdx_smType = ones(1,nSmTypes);
        
    figure(1); clf; hold on;
    for ti = 1:nSmTypes        
        h1ValMeans_sum(ti) = plot(0, 0, [sm_cols(ti) '.-']);      %#ok<AGROW>
        h1ValMeans_pos(ti) = plot(0, 0, [sm_cols(ti) '.--']);     %#ok<AGROW>
        h1ValMeans_neg(ti) = plot(0, 0, [sm_cols(ti) '.--']);     %#ok<AGROW> 
        h1CurSmAmt(ti) = plot(0, 0, [sm_cols(ti) 's'], 'markerfacecolor', sm_cols(ti), 'visible', 'off'); %#ok<AGROW>
    end
    title('Means');
    xlab_mean = xlabel('');
    ylab_mean = ylabel('');    
    h_line_mean_noSmooth = line([1 1], [0 0], 'linestyle', ':', 'color', 'k');
    legendLabs = {'gauss', 'fermi', 'Alias'};
    legend(h1ValMeans_sum(1:nSmTypes), legendLabs(1:nSmTypes), 'location', 'NE');
    
        
    binEdges = [-1 1];
    binCent = binEdge2cent(binEdges);

    fmt = '%.2f';
    str_mn_std = @(x) sprintf( [fmt ' \\pm ' fmt ' (%d) '], mean(x), stderr(x), length(x));

    figure(2); clf;
    h2CC_ax = subplot(1,2,1);
    h2CC_pts = plot(0, 0, 'bo', 'markersize', 2); hold on;
    h2CurPoint = plot(0,0, 'ro', 'markerfacecolor', 'r', 'visible', 'off');
    h2XEqY = plot(binEdges, binEdges, 'k:');
        
    h_title = title('', 'interpreter', 'none');
    h_lab(1) = xlabel('', 'interpreter', 'none');
    h_lab(2) = ylabel('', 'interpreter', 'none');    
    h_hst_ax(1) = subplot(2,2,2); hCC_hst(1) = bar(0, 0, 1);  hCC_hst_tit(1) = title(''); hCC_hst_ylab(1) = ylabel('');
    h_hst_ax(2) = subplot(2,2,4); hCC_hst(2) = bar(0, 0, 1);  hCC_hst_tit(2) = title(''); hCC_hst_ylab(2) = ylabel(''); 
    set(h_hst_ax, 'activepositionproperty', 'outerposition')
    
    uicontrol('parent', 2, 'style', 'pushbutton', 'String', ' + ', 'units', 'pixels', 'position', [3 3, 20 20], 'callback', @updateTuningCurves, 'tag', 'Button')
    
    
    figure(3); clf;
    subplot(3,1,1);
    h3TC(1,1,1) = plot(0,0, 'bo-'); hold on;
    h3TC(1,2,1) = plot(0,0, 'go-');    
    h3TC(1,1,2) = plot(0,0, 'b:', 'visible', 'off');
    h3TC(1,2,2) = plot(0,0, 'g:', 'visible', 'off');
    h3TC_ylab(1) = ylabel('');
    h3TC_tit(1) = title('');    
    xlim([0 360]);
    
    subplot(3,1,2);
    h3TC(2,1,1) = plot(0,0, 'bo-'); hold on;
    h3TC(2,2,1) = plot(0,0, 'go-');
    h3TC(2,1,2) = plot(0,0, 'b:', 'visible', 'off');
    h3TC(2,2,2) = plot(0,0, 'g:', 'visible', 'off');
    h3TC_ylab(2) = ylabel('');
    h3TC_tit(2) = title('');
    xlim([0 360]);

    subplot(3,1,3); hold on;
    for ti = 1:nSmTypes
        h3CurPairVals(ti)  = plot(0, 0, [sm_cols(ti) 'o:']); hold on; %#ok<AGROW>
        h3CurPairCurVal(ti) = plot(0, 0, ['rs'], 'visible', 'off'); %#ok<AGROW>
    end    
    h_tit_curPair = title('Values for current pair');    
    h_xlab_curPair = xlabel('');
    h_line_curPair_noSmooth = line([1 1], [0 0], 'linestyle', ':', 'color', 'k');    
    

    allGids =    [OSPs.Gid];
    allCellIds = [OSPs.cellId];
    curIdxs = [];
    
    wrp = @(x) x([1:end, 1]);
    ext = @(x) [x(:); x(end)+diff(x(1:2))]';
    
    function c = computeCorr(tc1, tc2)
        switch lower(curMeasure)
            case 'dot',  c = normDotProd(tc1, tc2);
            case 'cc',  c = pearsonR(tc1, tc2);
            case 'rho', c = spearmanRho(tc1, tc2);
            case 'dphi', c = abs( deltaPhi(tc1, tc2, 'cross-correlation') );
            case 'dF1', c = deltaPhi(tc1, tc2, 'angle');
        end
    end
        
    
    function updateTC(h, ph, tc)
%         set(h(1), 'xdata', ext(ph), 'ydata', wrp(tc));        
        set(h(1), 'xdata', ph, 'ydata', tc);        
        if strcmp(curMeasure, 'dF1');
            [phi, f_cos1, t1] = getF1phase(deg2rad(ph), tc, 2*pi);            
            set(h(2), 'xdata', rad2deg(ext(t1)), 'ydata', wrp(f_cos1), 'visible', 'on');
        else
            set(h(2), 'visible', 'off');
        end
        
    end
    
    function updateTuningCurves(hObject, eventdata, handles) %#ok<INUSD>
        if strcmp(get(gcbo,'Tag'), 'Button')
            selectPairOfTuningCurves;        
        end        

        if isempty(curPairIdx)
            set(h2CurPoint, 'visible', 'off')        
            return;
        end
        
        xs = get(h2CC_pts, 'xdata'); ys = get(h2CC_pts, 'ydata');
        set(h2CurPoint, 'xdata', xs(curPairIdx), 'ydata', ys(curPairIdx), 'visible', 'on')        

        if ~isempty(curPairIdx)  % data for 3.3: cur pair vs smoothing amount.
            for smMode_i = 1:nSmTypes
                sm_curPairVals = cellfun(@(v) v(curPairIdx), vals_forSmType{smMode_i});
                set(h3CurPairVals(smMode_i), 'xdata', smoothAmts_x{smMode_i}, 'ydata', sm_curPairVals, 'visible', 'on');
                if smMode_i == curSmModeIdx
                    set(h3CurPairCurVal(smMode_i), 'xdata', smoothAmts_x{curSmModeIdx}(curSmAmtIdx), 'ydata', sm_curPairVals(curSmAmtIdx), 'visible', 'on');                    
                else
                    set(h3CurPairCurVal(smMode_i), 'visible', 'off');
                end
            end
            x_lims = [min(smoothAmts_x{smMode_i}), max(smoothAmts_x{smMode_i})];
            curPairVal_noSmooth = vals_forNoSmooth(curPairIdx);
            set(h_line_curPair_noSmooth, 'xdata', x_lims, 'ydata', [1;1]*curPairVal_noSmooth);
        end
        
        % update 2 tuning curves (before / after)
        for i = 1:2
        
            pd = pairData(curIdxs(curPairIdx));
            if i == 1
                oriSp_cmp = pd.oriSp_cmp(:,locId);                        
                [Gid_a, Gid_b] = dealV(pd.Gids);
                [cellId_a, cellId_b] = dealV(pd.cellIds);
                idx_a = find(Gid_a == allGids  &  cellId_a == allCellIds, 1);
                idx_b = find(Gid_b == allGids  &  cellId_b == allCellIds, 1);

                ph = OSPs(idx_a).ph;
                
                tc_a = squeeze( OSPs(idx_a).R(oriSp_cmp(1), oriSp_cmp(2), :) );
                tc_b = squeeze( OSPs(idx_b).R(oriSp_cmp(1), oriSp_cmp(2), :) ); 
                
                [tc_a1, tc_b1] = deal(tc_a, tc_b);
                
%                 c1 = computeCorr(tc_a, tc_b);
%                 assert(curVals{1}(curPairIdx) == c1);
                
            elseif i == 2
%                 assert(isequal([Gid_a, Gid_b], curPairs{i}(curPairIdx).Gids))
%                 assert(isequal([cellId_a, cellId_b], curPairs{i}(curPairIdx).cellIds))                                
%                 tc_a = squeeze( curOsps{i}(idx_a).R(oriSp_cmp(1), oriSp_cmp(2), :) );
%                 tc_b = squeeze( curOsps{i}(idx_b).R(oriSp_cmp(1), oriSp_cmp(2), :) ); 
                
                tc_a = compressOSP_Phs(tc_a, smoothTypes{curSmModeIdx}, smoothAmts{curSmModeIdx}(curSmAmtIdx));
                tc_b = compressOSP_Phs(tc_b, smoothTypes{curSmModeIdx}, smoothAmts{curSmModeIdx}(curSmAmtIdx));                
                
                [tc_a2, tc_b2] = deal(tc_a, tc_b);
%                 c2 = computeCorr(tc_a, tc_b);
%                 assert(curVals{2}(curPairIdx) == c2);
            end
            
            updateTC(h3TC(i,1,:), ph, tc_a);
            updateTC(h3TC(i,2,:), ph, tc_b);      
            cmp_str = [curMeasure ' = ' num2str(curVals{i}(curPairIdx), '%.3f')];            
            if i == 1, 
                title_str = {sprintf('cell1: [%d, %d]. cell2: [%d, %d]', Gid_a, cellId_a, Gid_b, cellId_b), cmp_str};
            else
                title_str = cmp_str;
            end
            set(h3TC_tit(i), 'string', title_str);
                        
        end        
                
    end
%     figure(3); clf;
%     hAx1 = subplot(2,1,1); 
%     hAx2 = subplot(2,1,2);
    
    function selectPairOfTuningCurves
        figure(2);
        axes(h2CC_ax);
%         
        xlims0 = xlim;
        ylims0 = ylim;    
        xlims = xlims0;
        ylims = ylims0;
        zoom_level = 1;  max_zoom_level = 4;
        zoom_factor = 1.5;
    
        xs = get(h2CC_pts, 'xdata');
        ys = get(h2CC_pts, 'ydata');
        
        [cur_x,cur_y, button] = ginput(1);
        curPairIdx = [];
        done = false;
        while ibetween(cur_x, xlims) && ibetween(cur_y, ylims) && (gca == h2CC_ax) && ~done

            if button == 1  % find closest data point, redo stats for that point
                da = get(gca, 'DataAspectRatio');
                da = da(2)/da(1);            
                curPairIdx = indmin( normV([xs(:)-cur_x, (ys(:)-cur_y)/da], 2) );
%                 set(h2CurPoint, 'xdata', xs(idx), 'ydata', ys(idx), 'markerEdgeColor', 'k' );            
                done = true;            

            elseif button > 1
                zoom_level = cycle(zoom_level, 1:max_zoom_level);
                if zoom_level > 1
                    x_window = (diff(xlims) / (zoom_factor ^ zoom_level))/2;
                    y_window = (diff(ylims) / (zoom_factor ^ zoom_level))/2;
                    xlim( cur_x + [-x_window, + x_window]);
                    ylim( cur_y + [-y_window, + y_window]);                
                else
                    xlim(xlims0);
                    ylim(ylims0);
                end
                xlims = xlim;
                ylims = ylim;
                [cur_x,cur_y, button] = ginput(1);
            end
            
        end

        
    end
    
    


    function updateAllFigures(varargin)
        
        [showPosNeg, showMethods, smType, pairType, nPhasesTF, x_axis,    locName, msName] = deal(varargin{[1:6, end-1:end]});
        curSmAmts = [varargin{7:7+nSmTypes-1}];   % curSmAmts = [sG_val, sF_val, av_val];
        
        [showSum, showPos, showNeg] = dealV(showPosNeg);  
        showMethods = showMethods(1:nSmTypes);
        activeSmModeIds = find(showMethods(:)');
        inActiveSmModeIds = find(~showMethods(:)');
        onOff = @(tf) iff(tf, 'on', 'off');
        if any(strcmpi(msName, {'dot', 'cc', 'rho'}))
            nullMean = 0;
        elseif any(strcmpi(msName, {'dphi', 'dF1'}))
            nullMean = 90;
        end
        
        nPhasesActive = allNPhases(nPhasesTF);
        ph_idx = false(nPairs,1);
        for nph = nPhasesActive(:)';            
            ph_idx( nPhaseIdxs.(['ph' num2str(nph) ]) ) = true;
        end        
        idxs = find(pairTypeIdxs.(pairType) & ph_idx);
        
        if isempty(curIdxs)
            curIdxs = idxs;
        end
        
        if ~isequal(curIdxs, idxs)
            curPairIdx = [];
        end
        
        curMeasure = msName;
        msId = find(strcmp(msName, measures));
        locId = find(strcmp(locName, locations));
        
        
        for t_i = activeSmModeIds
            vals_forSmType{t_i} = cellfun(@(cmp) cmp{locId,msId}.val(idxs), sm_Cmps{t_i}, 'un', 0);
            meanVals_forSmType{t_i} = cellfun(@(v) nanmean(v), vals_forSmType{t_i});
            meanVals_forSmType_pos{t_i} = cellfun(@(v) nanmean(v(v>nullMean)), vals_forSmType{t_i});
            meanVals_forSmType_neg{t_i} = cellfun(@(v) nanmean(v(v<nullMean)), vals_forSmType{t_i});
                
        end        
        vals_forNoSmooth = noSm_Cmps{locId,msId}.val(idxs);
        meanVal_noSmooth = nanmean(vals_forNoSmooth);  
        
        % 1) update figure 1: mean (& pos/neg) vs smooth amount
        
        
        switch x_axis
            case 'degrees of freedom', 
                smoothAmts_x = smoothAmts_dof;
            case 'frequency cutoff',
                smoothAmts_x{1} = 10*sqrt(2)./smoothAmts{1};
                smoothAmts_x{2} = smoothAmts{t_i};
                smoothAmts_x{3} = smoothAmts_dof{t_i};                
            case 'smooth parameter',
                smoothAmts_x = smoothAmts;
        end
        set([h_xlab_curPair, xlab_mean], 'string', x_axis);
                
        for t_i = activeSmModeIds
            set(h1ValMeans_sum(t_i), 'xdata', smoothAmts_x{t_i}, 'ydata', meanVals_forSmType{t_i}, 'visible', onOff(showSum));
            set(h1ValMeans_pos(t_i), 'xdata', smoothAmts_x{t_i}, 'ydata', meanVals_forSmType_pos{t_i}, 'visible', onOff(showPos));
            set(h1ValMeans_neg(t_i), 'xdata', smoothAmts_x{t_i}, 'ydata', meanVals_forSmType_neg{t_i}, 'visible', onOff(showNeg));
        end        
        for t_i = inActiveSmModeIds
            set([h1ValMeans_sum(t_i), h1ValMeans_pos(t_i), h1ValMeans_neg(t_i)], 'visible', 'off')
        end
        

        curSmModeIdx = find(strcmp(smType, smoothTypes));        
        otherSmModeIds = setdiff([1:nSmTypes], curSmModeIdx);        
        amtIdx_smType(curSmModeIdx) = indmin( abs( smoothAmts{curSmModeIdx} - curSmAmts(curSmModeIdx) ));
        curSmAmtIdx = amtIdx_smType(curSmModeIdx);
        
        set(h1CurSmAmt(curSmModeIdx), 'xdata', smoothAmts_x{curSmModeIdx}(curSmAmtIdx), 'ydata', meanVals_forSmType{curSmModeIdx}(curSmAmtIdx), 'visible', 'on');
        set(h1CurSmAmt(otherSmModeIds), 'visible', 'off');

        x_lims = [min(smoothAmts_x{curSmModeIdx}), max(smoothAmts_x{curSmModeIdx})];
        set(h_line_mean_noSmooth, 'xdata', x_lims, 'ydata', [1; 1]*meanVal_noSmooth);
                                
        if ~isempty(curPairIdx) % for individual pair curve on figure 3                         
            for smMode_i = 1:nSmTypes                
                curPairVals_smType{smMode_i} = cellfun(@(v) v(curPairIdx), vals_forSmType{smMode_i});                

                set(h3CurPairVals(smMode_i), 'xdata', smoothAmts_x{smMode_i}, 'ydata', curPairVals_smType{smMode_i})
                if smMode_i == curSmModeIdx
                    set(h3CurPairCurVal(smMode_i), 'xdata', smoothAmts_x{smMode_i}(curSmAmtIdx), 'ydata', curPairVals_smType{smMode_i}(curSmAmtIdx), 'visible', 'on');
                else
                    set(h3CurPairCurVal(smMode_i), 'visible', 'off');
                end
            end
            curPairVal_noSmooth = vals_forNoSmooth(curPairIdx);
            set(h_line_curPair_noSmooth, 'xdata', x_lims, 'ydata', [1;1]*curPairVal_noSmooth);
            
        end
        
        % 1. update plot #1 (with cc vs smoothing level)                                
                                
        val_strs = {'0'; num2str(curSmAmts(curSmModeIdx))};
                                
        
        % 2a. update data for figure 2  (with ccs 1 vs ccs 2)
        curCmps = {noSm_Cmps, sm_Cmps{curSmModeIdx}{curSmAmtIdx}};        
        curVals = {curCmps{1}{locId, msId}.val(idxs),  curCmps{2}{locId, msId}.val(idxs)};
        fmt = '%.3f';
        str_mn_std = @(x) sprintf( [fmt ' \\pm ' fmt ' (%d) '], mean(x), stderr(x), length(x));        
        set(h2CC_pts, 'xdata', curVals{1}, 'ydata', curVals{2});        
        set(h_title, 'string', [msName ' at ' locName]);
        binEdges = curCmps{1}{locId, msId}.binEdges;        
        binCent = binEdge2cent(binEdges);
        axis(h2CC_ax, [binEdges(1), binEdges(end), binEdges(1), binEdges(end)]);
        set(h_hst_ax, 'xlim', [binEdges(1), binEdges(end)])
        set(h2XEqY, 'xdata', binEdges, 'ydata', binEdges);
        set(ylab_mean, 'string', ['mean of ' msName ' at ' locName]);
        
        % 2b. update plot #2  (with ccs 1 vs ccs 2)
        for i = 1:2   % x and y axes of figure 2
            set(hCC_hst_tit(i), 'string', str_mn_std(curVals{i}(~isnan(curVals{i}))) );
            set(h_lab(i), 'string', [smType ' = ' val_strs{i}]);
            set([h3TC_ylab(i) hCC_hst_ylab(i)], 'string', [smType ' = ' val_strs{i}]);            
            ns = histcnt(curVals{i}, binEdges);
            set(hCC_hst(i), 'xdata', binCent, 'ydata', ns);
        end
        
        % 3. update plot #3  (with individual tuning curves)
        updateTuningCurves;
        
        curIdxs = idxs;
    end



    show1_0 = [true, false, false];
    showMeth_0 = [true, true, true];
    x_axisOptions = {'degrees of freedom', 'frequency cutoff', 'smooth parameter'};
    
    nPhaseArg = {'nPhases', repmat({[true, false]}, 1, length(allNPhases)), [], [], arrayfun(@(nph) ['ph' num2str(nph)], allNPhases, 'un', 0 ) };
         
    mth_vars = { {'Gauss'}, {'Fermi'}, {'Alias'} };
    mth_vars = mth_vars(smON);
    tf = [true, false];
    showLabels = {'sum', 'pos', 'neg'};
    methArgs = cellfun(@(lab, amt) {lab, amt}, smoothLabels, smoothAmts, 'un', 0);

    updateAllFigures(show1_0, showMeth_0, smoothTypes{1}, allPairTypes{1}, [true, false], x_axisOptions{1}, smGaussAmt_dof(1), smFermiAmt_dof(1), smAveragesAmt_dof(1), locations{1}, measures{1});
    
    
    args = { {'showPosNeg', {tf,tf,tf}, show1_0, [], showLabels}, ...        
             {'showSmMeth', repmat({tf}, 1, nSmTypes), showMeth_0(smON), [], mth_vars}, ...        
             {'method', smoothTypes, smoothTypes{1}, mth_vars}, ...     
             {'pairType', allPairTypes}, ...
             nPhaseArg, ...             
             {'x_axis', x_axisOptions}, ...        
             methArgs{:}, ...
             {'location', locations}, {'measure', measures}};
             
    manipulate(@updateAllFigures, args, 'FigId', 100);


end