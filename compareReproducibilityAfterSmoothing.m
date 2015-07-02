function compareReproducibilityAfterSmoothing

%     curGratingType('d');
    gratingType = curGratingType('');    
    smFermi_b_hp = 100;
    
    global smoothPhases showWorking
                
    showWorking = [];
    S1 = load([gratingType 'GratingCells_full.mat']);
    OSPs = S1.allOSPs;
    n_phases = arrayfun(@(s) length(s.ph), OSPs);
    allNPhases = unique( n_phases );
    nOSPs = length(OSPs);
    
    function stats = getStatsForSmoothing(cmpAction, cmpN)
        if strcmp(cmpAction, 'noSmooth')
            smoothPhases = {};
            ext = '';
        else
            if strcmp(cmpAction, 'Fermi')
                cmpN = [cmpN, cmpN / smFermi_b_hp];
            end
            smoothPhases = {cmpAction, cmpN};
            ext = ['_' cmpAction '_' num2str(cmpN(1))];
        end
        
        data_dir = [CatV1Path 'smooth_test\stats_test\'];        
        statsDatafile = [data_dir gratingType 'OSPstats' ext '.mat'];
        
        if ~exist(statsDatafile, 'file')
            recalculateStatsWithSmoothing(cmpAction, cmpN, statsDatafile)
        end
        
        S = load(statsDatafile);
        stats = S.allstats;        
    end

        
%     fermiKtoDOF = @(k) min(2*k, 30);
    fermiDOFtoK = @(k) k/2;
    
    if strcmp(gratingType, 'flashed')
        smGaussAmt_dof = [3:.1:7];
        smFermiAmt_dof = [3:.5:7];
        smAveragesAmt_dof = [3:.5:7];
        
        b = [12.9415   -1.3612    2.0189    0.7425    2.9712];
        gaussWtoDOF = @(b, w) max(min( [b(1)./((w-b(2)).^(w.*b(3) + b(4) ) ) + b(5)], 8*ones(size(w)) ), 3*ones(size(w))) ;                        
        gaussDOFtoW = @(dof) arrayfun( @(d) fminsearch(@(w) abs(gaussWtoDOF(b, w)-d), 1), dof);        
        
    elseif strcmp(gratingType, 'drifting')
        b = [87.2120   -0.6554   -1.1062    1.6746    2.0664];        
        gaussWtoDOF = @(b, w) max(min( [b(1)./((w-b(2)).^((w).^b(3) + b(4) ) ) + b(5)], 60*ones(size(w)) ), 3*ones(size(w))) ;
        gaussDOFtoW = @(dof) arrayfun( @(d) fminsearch(@(w) abs(gaussWtoDOF(b, w)-d), 1), dof);        
                
        smGaussAmt_dof = [3:2:60];
        smFermiAmt_dof = [2:20, 22:2:60];
        smAveragesAmt_dof = [2:20,30,60];

        smGaussAmt_dof = [5:5:55];
        
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
        
    sm_stats = cell(nSmTypes,1);
    for ti = 1:nSmTypes
        sm_stats{ti} = cell(1,sm_N(ti));
        for sm_i = 1:sm_N(ti)
            sm_stats{ti}{sm_i} = getStatsForSmoothing(smoothTypes{ti}, smoothAmts{ti}(sm_i));
        end
    end
    noSm_stats = getStatsForSmoothing('noSmooth', []);

    smoothAmts_x = smoothAmts_dof;
    
    
    [curStats,curVals] = deal({});
          
    curOspIdx = [];
    curSmModeIdx = 1;
    curSmAmtIdx = 1;

    
	[vals_forSmType, vals_forSmType_sum, meanVals_forSmType, curOspVals_smType] = ...
        deal(cell(1,nSmTypes));
    [vals_forNoSmooth, meanVal_noSmooth, curOspVal_noSmooth] = deal([]);    
    amtIdx_smType = ones(1,nSmTypes);
        
    figure(1); clf; hold on;
    for ti = 1:nSmTypes        
        h1ValMeans(ti) = plot(0, 0, [sm_cols(ti) '.-']);      %#ok<AGROW>
        h1CurSmAmt(ti) = plot(0, 0, [sm_cols(ti) 's'], 'markerfacecolor', sm_cols(ti), 'visible', 'off'); %#ok<AGROW>
    end
    title('Means');
    xlab_mean = xlabel('');
    ylab_mean = ylabel('');    
    h_line_mean_noSmooth = line([1 1], [0 0], 'linestyle', ':', 'color', 'k');
    legendLabs = {'gauss', 'fermi', 'Alias'};
    legend(h1ValMeans(1:nSmTypes), legendLabs(1:nSmTypes), 'location', 'SE');
        
    
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
%     h3LR(1,1,1) = plot(0,0, 'bo-'); hold on;
%     h3LR(1,2,1) = plot(0,0, 'go-');    
%     h3LR(1,1,2) = plot(0,0, 'b:', 'visible', 'off');
%     h3LR(1,2,2) = plot(0,0, 'g:', 'visible', 'off');
%     h3LR_ylab(1) = ylabel('');
%     h3LR_tit(1) = title('');    
%     xlim([0 360]);
    
    figure(4); clf;
%     subplot(3,1,2);
%     h3LR(2,1,1) = plot(0,0, 'bo-'); hold on;
%     h3LR(2,2,1) = plot(0,0, 'go-');
%     h3LR(2,1,2) = plot(0,0, 'b:', 'visible', 'off');
%     h3LR(2,2,2) = plot(0,0, 'g:', 'visible', 'off');
%     h3LR_ylab(2) = ylabel('');
%     h3LR_tit(2) = title('');
%     xlim([0 360]);

    figure(5); clf;
    
    h_rep_xy(1) = plot(0,0, 'bo'); hold on;
    h_rep_xy(2) = plot(0,0.5, 'r.');

    figure(6); clf;
%     subplot(3,1,3); hold on;
    for ti = 1:nSmTypes
        h3CurPairVals(ti)  = plot(0, 0, [sm_cols(ti) 'o:']); hold on; %#ok<AGROW>
        h3CurPairCurVal(ti) = plot(0, 0, [sm_cols(ti) 's'], 'markerfacecolor', sm_cols(ti), 'visible', 'off'); %#ok<AGROW>
    end    
    h_tit_curOsp = title('Values for current pair');    
    h_xlab_curOsp = xlabel('');
    h_line_curOsp_noSmooth = line([1 1], [0 0], 'linestyle', ':', 'color', 'k');    
    
%     smoothedPairData.Gids]
    allGids =    [OSPs.Gid];
    allCellIds = [OSPs.cellId];
    curIdxs = [];
    curField = '';
        
    
    function updateTuningCurves(hObject, eventdata, handles) %#ok<INUSD>
        if strcmp(get(gcbo,'Tag'), 'Button')
            selectOSP;        
        end        

        if isempty(curOspIdx)
            set(h2CurPoint, 'visible', 'off')        
            return;
        end
        
        xs = get(h2CC_pts, 'xdata'); ys = get(h2CC_pts, 'ydata');
        set(h2CurPoint, 'xdata', xs(curOspIdx), 'ydata', ys(curOspIdx), 'visible', 'on')        

        if ~isempty(curOspIdx)  % data for 3.3: cur pair vs smoothing amount.
            for smMode_i = 1:nSmTypes
                sm_curOspVals = cellfun(@(v) v(curOspIdx), vals_forSmType{smMode_i});
                set(h3CurPairVals(smMode_i), 'xdata', smoothAmts_x{smMode_i}, 'ydata', sm_curOspVals, 'visible', 'on');
                if smMode_i == curSmModeIdx
                    set(h3CurPairCurVal(smMode_i), 'xdata', smoothAmts_x{curSmModeIdx}(curSmAmtIdx), 'ydata', sm_curOspVals(curSmAmtIdx), 'visible', 'on');                    
                else
                    set(h3CurPairCurVal(smMode_i), 'visible', 'off');
                end
            end
            x_lims = [min(smoothAmts_x{smMode_i}), max(smoothAmts_x{smMode_i})];
            curOspVal_noSmooth = vals_forNoSmooth(curOspIdx);
            set(h_line_curOsp_noSmooth, 'xdata', x_lims, 'ydata', [1;1]*curOspVal_noSmooth);
        end
        
        % update 2 tuning curves (before / after)
        
        osp{1} = OSPs(curOspIdx);
        osp{2} = compressOSP_Phs(osp{1}, smoothTypes{curSmModeIdx}, smoothAmts{curSmModeIdx}(curSmAmtIdx));

        [x, y] = deal(cell(2,1));
        showField = strrep(curField, 'rep_', '');
        showField = strrep(showField, '_pval', '');
        showField = strrep(showField, '_str', '');
        
        for i = 1:2
                        
            [R, R_full, f] = deal(osp{i}.R, osp{i}.R_full, osp{i}.R_full_factor);
            figId = 2+i;
            showWorking = [showField num2str(figId)];
            stats = calcStatsFromOSPfull(R, R_full, f);
            h_ax = get(figId, 'children');
            h_line = get(h_ax, 'children');  
            x{i} = get(h_line(3), 'xdata');
            y{i} = get(h_line(3), 'ydata');
            set([h_ax; h_line], 'visible', 'off')
            
            if i == 2
                scl1 = mean([x{1}, y{1}]);
                scl2 = mean([x{2}, y{2}]);
                x{2} = x{2}*scl1/scl2;
                y{2} = y{2}*scl1/scl2;
            end
            
            title([ curField ' = ' num2str(stats.(curField), '%.3g') ]);
                        
            figure(5);
            set(h_rep_xy(i), 'xdata', x{i}, 'ydata', y{i});
            
            figure(10+i);
            imageOSP(osp{i}, 'pref:ori', 'SPO');
        end        
        
        3;
                
    end



    function selectOSP
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
        curOspIdx = [];
        done = false;
        while ibetween(cur_x, xlims) && ibetween(cur_y, ylims) && (gca == h2CC_ax) && ~done

            if button == 1  % find closest data point, redo stats for that point
                da = get(gca, 'DataAspectRatio');
                da = da(2)/da(1);            
                curOspIdx = indmin( normV([xs(:)-cur_x, (ys(:)-cur_y)/da], 2) );
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
        
        [showMethods, smType, nPhasesTF, repField, x_axis    ] = deal(varargin{[1:5]});
        curSmAmts = [varargin{6:6+nSmTypes-1}];   % curSmAmts = [sG_val, sF_val, av_val];
                
        showMethods = showMethods(1:nSmTypes);
        activeSmModeIds = find(showMethods(:)');
        inActiveSmModeIds = find(~showMethods(:)');
        onOff = @(tf) iff(tf, 'on', 'off');
        
        nPhasesActive = allNPhases(nPhasesTF);
        ph_idx = false(nOSPs,1);
        for nph = nPhasesActive(:)';            
            ph_idx( n_phases == nph ) = true;
        end        
        idxs = find(ph_idx);
        
        curField = repField;
        
        if ~isempty(strfind(curField, '_pval'));            
%             r = @(x)  x .* isfinite(x);
            f = @(x) pToLogP(x);
        else
            f = @(x) x;
        end
        
        if isempty(curIdxs)
            curIdxs = idxs;
        end
        
        if ~isequal(curIdxs, idxs)
            curOspIdx = [];
        end
                        
        for t_i = activeSmModeIds
            vals = cellfun(@(st) f([st(idxs).(repField)]), sm_stats{t_i}, 'un', 0);
            vals_forSmType{t_i} = vals;                
            meanVals_forSmType{t_i} = cellfun(@(v) nanmean(v), vals_forSmType{t_i});                
        end        
        vals_forNoSmooth = f([noSm_stats(idxs).(repField)]);
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
        set([h_xlab_curOsp, xlab_mean], 'string', x_axis);
                
        for t_i = activeSmModeIds
            set(h1ValMeans(t_i), 'xdata', smoothAmts_x{t_i}, 'ydata', meanVals_forSmType{t_i}, 'visible', 'on');
        end        
        for t_i = inActiveSmModeIds
            set([h1ValMeans(t_i)], 'visible', 'off')
        end
        

        curSmModeIdx = find(strcmp(smType, smoothTypes));        
        otherSmModeIds = setdiff([1:nSmTypes], curSmModeIdx);        
        amtIdx_smType(curSmModeIdx) = indmin( abs( smoothAmts{curSmModeIdx} - curSmAmts(curSmModeIdx) ));
        curSmAmtIdx = amtIdx_smType(curSmModeIdx);
        
        set(h1CurSmAmt(curSmModeIdx), 'xdata', smoothAmts_x{curSmModeIdx}(curSmAmtIdx), 'ydata', meanVals_forSmType{curSmModeIdx}(curSmAmtIdx), 'visible', 'on');
        set(h1CurSmAmt(otherSmModeIds), 'visible', 'off');

        x_lims = [min(smoothAmts_x{curSmModeIdx}), max(smoothAmts_x{curSmModeIdx})];
        set(h_line_mean_noSmooth, 'xdata', x_lims, 'ydata', [1; 1]*meanVal_noSmooth);
                                
        if ~isempty(curOspIdx) % for individual pair curve on figure 3                         
            for smMode_i = 1:nSmTypes                
                curOspVals_smType{smMode_i} = cellfun(@(v) v(curOspIdx), vals_forSmType{smMode_i});                

                set(h3CurPairVals(smMode_i), 'xdata', smoothAmts_x{smMode_i}, 'ydata', curOspVals_smType{smMode_i})
                if smMode_i == curSmModeIdx
                    set(h3CurPairCurVal(smMode_i), 'xdata', smoothAmts_x{smMode_i}(curSmAmtIdx), 'ydata', curOspVals_smType{smMode_i}(curSmAmtIdx), 'visible', 'on');
                else
                    set(h3CurPairCurVal(smMode_i), 'visible', 'off');
                end
            end
            curOspVal_noSmooth = vals_forNoSmooth(curOspIdx);
            set(h_line_curOsp_noSmooth, 'xdata', x_lims, 'ydata', [1;1]*curOspVal_noSmooth);
            
        end
        
        % 1. update plot #1 (with cc vs smoothing level)                                
                                
        val_strs = {'0'; num2str(curSmAmts(curSmModeIdx))};
                                
        
        % 2a. update data for figure 2  (with ccs 1 vs ccs 2)
        curStats = {noSm_stats, sm_stats{curSmModeIdx}{curSmAmtIdx}};        
        curVals = {f([curStats{1}(idxs).(repField)]),  f([curStats{2}(idxs).(repField)])};
%         curVals = cellfun(@(x) x(~isnan(x)), curVals, 'un', 0);
        fmt = '%.3f';
        str_mn_std = @(x) sprintf( [fmt ' \\pm ' fmt ' (%d) '], mean(x), stderr(x), length(x));        
        set(h2CC_pts, 'xdata', curVals{1}, 'ydata', curVals{2});        
%         set(h_title, 'string', '');
        binEdges = linspace(floor(min(curVals{1})), ceil(max(curVals{1})), 10);
        binCent = binEdge2cent(binEdges);
        axis(h2CC_ax, [binEdges(1), binEdges(end), binEdges(1), binEdges(end)]);
        set(h_hst_ax, 'xlim', [binEdges(1), binEdges(end)])
        set(h2XEqY, 'xdata', binEdges, 'ydata', binEdges);
%         set(ylab_mean, 'string', ['mean of ' msName ' at ' locName]);
        
        % 2b. update plot #2  (with ccs 1 vs ccs 2)        
        for i = 1:2   % x and y axes of figure 2            
            set(hCC_hst_tit(i), 'string', {str_mn_std(nonnans(curVals{i})), num2str(nnz(curVals{i}>2))  });
            set(h_lab(i), 'string', [smType ' = ' val_strs{i}]);
%             set([h3LR_ylab(i) hCC_hst_ylab(i)], 'string', [smType ' = ' val_strs{i}]);            
            ns = histcnt(curVals{i}, binEdges);
            set(hCC_hst(i), 'xdata', binCent, 'ydata', ns);
            
        end
        
        % 3. update plot #3  (with individual tuning curves)
        updateTuningCurves;
        
        curIdxs = idxs;
    end


    allRepFields = fieldnames( noSm_stats(1) );

    showMeth_0 = [true, true, true];
    x_axisOptions = {'degrees of freedom', 'frequency cutoff', 'smooth parameter'};
    
    nPhaseArg = {'nPhases', repmat({[true, false]}, 1, length(allNPhases)), [], [], arrayfun(@(nph) ['ph' num2str(nph)], allNPhases, 'un', 0 ) };
         
    mth_vars = { {'Gauss'}, {'Fermi'}, {'Alias'} };
    mth_vars = mth_vars(smON);
    tf = [true, false];
    methArgs = cellfun(@(lab, amt) {lab, amt}, smoothLabels, smoothAmts, 'un', 0);

    updateAllFigures(showMeth_0, smoothTypes{1}, [false, true], allRepFields{1}, x_axisOptions{1}, ...
        smGaussAmt_dof(1), smFermiAmt_dof(1), smAveragesAmt_dof(1) );
    
    
    args = { {'showSmMeth', repmat({tf}, 1, nSmTypes), showMeth_0(smON), [], mth_vars}, ...        
             {'method', smoothTypes, smoothTypes{1}, mth_vars}, ...                  
             nPhaseArg, ...  
             {'rep_field', allRepFields}, ...
             {'x_axis', x_axisOptions}, ...        
             methArgs{:}, ...
             };
             
    manipulate(@updateAllFigures, args, 'FigId', 100)



end


function recalculateStatsWithSmoothing(cmpAction, cmpN, filename)
    persistent allOSPs N;
    
    if isempty(allOSPs)        
        gratingType = curGratingType('');
        osp_file = [CatV1Path gratingType 'GratingCells_full.mat'];            
        if ~exist(osp_file, 'file')
            error('Need OSP_full file');
        end        
        S = load(osp_file);
        allOSPs = S.allOSPs;
        N = length(allOSPs);
%         allstats = [allOSPs.stats];
    end
        
    fprintf('computing stats : %s \n', filename);
    progressBar('init-', N);
    for oi = 1:N
        progressBar(oi);
        smOSP = allOSPs(oi);
        if ~strcmp(cmpAction, 'noSmooth')
            smOSP = compressOSP_Phs(allOSPs(oi), cmpAction, cmpN);
        end
        [R_full, R, f] = deal(smOSP.R_full, smOSP.R, smOSP.R_full_factor);   
        newstats = calcStatsFromOSPfull(R, R_full, f);
        if oi == 1
            allstats = repmat(newstats, N, 1);
        else        
            allstats(oi) = newstats;
        end
    end
    progressBar('done');
    
    save(filename, 'allstats');
end
    

function y = pToLogP(x)
    y = -log10(x);
    y(x == 0) = 10;
end