function viewStat1VsStat2_new(allOSPs)
%     if nargin < 1
%         S = load('flashedGratingOSPs.mat', 'allOSPs');
%         allOSPs = S.allOSPs;
%     end
    global showWorking;
    clear X Y nX nY dist ind_closestXY E_tot E_specific inds;

    gratingType = curGratingType('');  % FLASHED_GRATING = 1; DRIFTING_GRATING = 2;
    ospDatafile = [CatV1Path gratingType 'GratingOSPs.mat'];    
        
%     if ~exist('allOSPs', 'var')
        S = load(ospDatafile);
        allOSPs = S.allOSPs;
%     end
    ranks = [allOSPs.manualRank];
    allstats = [allOSPs.stats];
    showMarginalDistributions = true;

%     stat_name{1} = 'rep_ori_sp_avPhase_str';
%     stat_name{2} = 'rep_ori_sp_ph_str';
%     stat_name{2} = 'rep_ori_sp_maxPhase_str';

%     stat_name{1} = 'rep_ori_sp_avPhase_pval';
%     stat_name{2} = 'rep_ori_sp_ph_pval';
%     stat_name{2} = 'rep_ori_sp_maxPhase_pval';

%     stat_name{1} = 'rep_ori_sp_avPhase_pval';
%     stat_name{2} = 'rep_ori_sp_avPhase_str';

    stat_name{1} = 'rep_ori_sp_ph_pval';
    stat_name{2} = 'rep_ori_sp_ph_str';

% 'rep_ori_sp_avPhase_p', x = [allstats.rep_ori_sp_avPhase_p];
% 'rep_ori_sp_avPhase_sgn', x = [allstats.rep_ori_sp_avPhase_sgn];
% 'rep_ori_sp_avPhase_str',x = [allstats.rep_ori_sp_avPhase_str];
% 'rep_ori_sp_maxPhase_p',x = [allstats.rep_ori_sp_maxPhase_p];
% 'rep_ori_sp_maxPhase_sgn',x = [allstats.rep_ori_sp_maxPhase_sgn];
% 'rep_ori_sp_maxPhase_str',x = [allstats.rep_ori_sp_maxPhase_str];
% 'rep_ori_sp_ph_p',x = [allstats.rep_ori_sp_ph_p];
% 'rep_ori_sp_ph_sgn',x = [allstats.rep_ori_sp_ph_sgn];
% 'rep_ori_sp_ph_str',x = [allstats.rep_ori_sp_ph_str];
    

%     stat_name{1} = 'resp_size';
%     stat_name{2} = 'ori_sel_str';

%     ranksToPlot = [-1, 0 3];
%     rankColors = 'krb';
    ranksToPlot = [ 0 1 2 3];
    rankColors = 'rmgb';

    drawPerceptronBoundaries = false;
    
    relRisks = [1:7:50; 50:-7:1]';   
%     relRisks = [1 1];   
    
%     colorRef = {'r', [1 0 0]; 'g', [0 1 0]; 'b', [0 0 1]};
%     colorOrder = arrayfun(@(c) colorRef{ strfind([colorRef{:,1}], c)  ,2}', rankColors', 'un', 0);
%     colorOrder = [colorOrder{:}]';
%     statLabel = {};
    nR = length(ranksToPlot);    
    showP_0 = false;
    for i = 1:2
        switch stat_name{i}
            
            case 'rep_ori_sp_avPhase_pval', x = [allstats.rep_ori_sp_avPhase_pval];
            case 'rep_ori_sp_avPhase_sgn', x = [allstats.rep_ori_sp_avPhase_sgn];
            case 'rep_ori_sp_avPhase_str',x = [allstats.rep_ori_sp_avPhase_str];
            case 'rep_ori_sp_maxPhase_pval',x = [allstats.rep_ori_sp_maxPhase_pval];
            case 'rep_ori_sp_maxPhase_sgn',x = [allstats.rep_ori_sp_maxPhase_sgn];
            case 'rep_ori_sp_maxPhase_str',x = [allstats.rep_ori_sp_maxPhase_str];
            case 'rep_ori_sp_ph_pval',         x = [allstats.rep_ori_sp_ph_pval];
            case 'rep_ori_sp_ph_sgn',   x = [allstats.rep_ori_sp_ph_sgn];
            case 'rep_ori_sp_ph_str',x = [allstats.rep_ori_sp_ph_str];
                             
            otherwise, error('Unknown stat name');
        end
        
        if ~isempty(strfind(stat_name{i}, '_pval')) 
            idx0 = find(x == 0); idxInf = find(x == inf);
            mn = min(setdiff(x, 0));     x(idx0) = mn;
            mx = max(setdiff(x, inf));   x(idxInf) = mx;
            x = -log10(x); 
            x(idx0) = -log10(mn)+ 2;%randn(size(idx0))/3;
            if ~isempty(idx0)                
                showP_0 = true;
                p_0 = -log10(mn)+ 2;
            end
        end
        if ~isempty(strfind(stat_name{i}, 'resp_size'))
            ind_hi = (x > 20);   x(ind_hi) = 20;
        end        
        stat{i} = x;
        
    end
    showWorking = [stat_name{:}];
    
    
    mainFig = 1001;
    figure(mainFig); clf; 
    if showMarginalDistributions        
        [a,b] = deal(1,3);
        
        hM = subplot(b,b, subplotInds(b,b,[1 1], [b-a, b-a])); hold on; box on;
        hR = subplot(b,b, subplotInds(b,b,[1, b-a+1], [b-a, b])); 
        hD = subplot(b,b, subplotInds(b,b,[b-a+1, 1], [b, b-a])); 
%         set([hM, hR, hD], 'colorOrder', colorOrder);
    else
        hM = gca;
    end
    box on;
    
    % plot each set of points.
    for ri = 1:nR
        inds{ri} = find(ranks == ranksToPlot(ri)); %#ok<*AGROW>
        
        X{ri} = stat{1}(inds{ri});
        Y{ri} = stat{2}(inds{ri});
        
        plot(hM, X{ri}, Y{ri}, [rankColors(ri) 'o'], 'markersize', 2);
    end
    allX = [X{:}]; allY = [Y{:}];
    q = -log10(.05);
    fprintf('# of cells with cc > 0: %d / %d (%.3f)\n', nnz(allX>q), length(allX),nnz(allX>0)/length(allX) );

    nbins = 15;
    xedges = linspace(min(allX)-eps, max(allX)+eps, nbins+1);
    if any(~isfinite(xedges)),  error('bad values in stat1'); end
    yedges = linspace(min(allY)-eps, max(allY)+eps, nbins+1);
    if any(~isfinite(yedges)),  error('bad values in stat2'); end
    nX = zeros(nbins, nR);
    for ri = 1:nR
        nX(:, ri) = histcnt(X{ri}, xedges);
        nY(:, ri) = histcnt(Y{ri}, yedges);
    end
    hbR = barh(hR, binEdge2cent(yedges), nY, 'stacked');
    hbD = bar(hD, binEdge2cent(xedges), nX, 'stacked', 'barwidth', .85);
    for i = 1:nR
        set(hbR(i), 'facecolor', rankColors(i));
        set(hbD(i), 'facecolor', rankColors(i));
    end
    
    if showP_0
        axes(hM)
        text(double(p_0), .1, {'\uparrow', 'p = 0'}, 'horiz', 'center' )
%         text(double(p_0), -.3, {'$$p=0$$'}, 'horiz', 'center', 'interpreter', 'latex', 'fontsize', 14)
        3;
    end
    
    
    % draw Perceptron boundaries    
    if drawPerceptronBoundaries
%         relRisks = [1:7:50; 50:-7:1]';   
        relRisks = [1, 1];   
        
        axes(hM);
        x1s = linspace(min(allX)-1, max(allX)+1, 150);
        x2FromWt = @(wt, x1) (-wt(1)-wt(2)*x1)/wt(3);

        axis(axis);
        Tid = cellfun(@(X, id) ones(1, length(X))*id, X([1,end]), num2cell([1,2]), 'un', 0);
        X_m = [X{[1,end]}; Y{[1,end]}];  Tid = [Tid{:}];
        for si = 1:size(relRisks,1);
            [Wt, E_tot(si), E_specific(si,:)] = perceptron(X_m, Tid, relRisks(si,:) );

            x2s = x2FromWt(Wt, x1s);
            plot(x1s, x2s, 'k-');            
%             [w0, w1, w2] = dealV(Wt);
%             c = -w0/w2;
%             m = -w1/w2;
%             x_test = [1; 0; c+1];        
%             sgn = -sign( Wt' * x_test); % negative b/c treats 
    %         plot(x1s, m*x1s+c, 'g:');    
            fprintf('Relative risk : %d-%d: Err tot: %d (%.2f%%). E1 = %d, E2 = %d\n',...
                 relRisks(si,1), relRisks(si,2), E_tot(si), E_tot(si)/length(Tid), E_specific(si,1), E_specific(si,2))

            Wts(:,si) = perceptron(X_m, Tid, relRisks(si,:) );
        end
    end
    
    
    
    % setup main figure for interactive selection of points
    figure(mainFig);
    axes(hM);
%     fplot(@(x) x, xlim, 'k:');
    hCur = plot(0,0, 'wo');    
%     xlim([-20 30]);
%     title(statLabel{1}, 'interpreter', 'none');
%     ylabel(statLabel{2}, 'interpreter', 'none');
    title(stat_name{1}, 'interpreter', 'none');
    ylabel(stat_name{2}, 'interpreter', 'none');
%     axis([0 10 0 1]);
    xlims0 = xlim;
    ylims0 = ylim;    
    xlims = xlims0;
    ylims = ylims0;
    zoom_level = 1;  max_zoom_level = 4;
    zoom_factor = 1.5;
%     while true
%         fprintf('x = %.2f, y = %.2f\n', cur_x,cur_y)        
%     end

    set(dispFigId, 'windowButtonDownFcn', {@selectPointsInFigure, @updateStatsOfChosenCell_grp});
    function updateStatsOfChosenCell_grp(glob_id, grp_id, loc_id, x, y) %#ok<INUSL>

            r = indmin(dist);
            idx = inds{r}(ind_closestXY(r));

            set(hCur, 'xdata', X{r}(ind_closestXY(r)), 'ydata', Y{r}(ind_closestXY(r)), 'markerEdgeColor', 'k', 'markerFaceColor', rankColors(r) );

            [Gid, cellId] = deal(allOSPs(idx).GroupId, allOSPs(idx).cellId);
            figure(123); clf; 
            if strcmp(gratingType, 'flashed')
                imageOSP(allOSPs(idx), 'mean:ph');
            elseif if strcmp(gratingType, 'drifting')
                imageOSP(allOSPs(idx), 'mean:ph');
                figure(124); clf; 
                imageOSP(allOSPs(idx), 'subplots:horizontal', 'SPO');
            end
            
            title(sprintf('Gid: %d, cellId: %d', Gid, cellId));

            redoStats_id = [Gid, cellId]; redoStats;

            figure(mainFig);
    
    
        glob_id = glob_idx{grp_id}(loc_id);
        Gid = allGids(glob_id);
        cellId = allCellIds(glob_id);
        set(h_diff_vs_p_t, 'string', sprintf('(%d) Gid = %d, cellId = %d. [X = %.2f, Y = %.2f]', glob_id, Gid, cellId, x, y));
        cellIdx = find(allGidsX == Gid & allCellIdsX == cellId);
        manipulateSet(psthVarHandles, 'Cell Index', cellIdx );
    
    
    end

    
    
    
    
end


%     stat_name{1} = 'oriSpf_rep_p_sgn';
%     stat_name{2} = 'oriSpf_rep_p_sgn2';


%     stat_name{1} = 'oriSpf_rep_p_sgn';
%     stat_name{2} = 'oriSpf_rep_str';


%     stat_name{1} = 'ori_sel_p';
%     stat_name{2} = 'ori_sel_str';
%     stat_name{1} = 'cov_oo';
%     stat_name{2} = 'cov_ss';

%     stat_name{1} = 'ori_rep_p_sgn';
%     stat_name{2} = 'ori_rep_str';

%     stat_name{1} = 'ori_rep_p';
%     stat_name{2} = 'ori_rep_str';


%         switch stat_name{i}
%             case 'ori_sel_p',      x = [allstats.orientationSelectivePval];        statLabel{i} = 'Orientation Selectivity (-log(p-value))';
%              case 'ori_sel_str',   x = [allstats.orientationSelectivityStrength];  statLabel{i} = 'Orientation Selectivity (strength)';
%             case 'ori_rep_p',      x = [allstats.orientationReproduciblePval];     statLabel{i} = 'Orientation Reproducibility (-log(p-value))';
%              case 'ori_rep_p_sgn', x = [allstats.orientationReproduciblePval_sgn]; statLabel{i} = 'Orientation Reproducibility (-log(p-value))';
%              case 'ori_rep_str',   x = [allstats.orientationReproducibleStrength]; statLabel{i} = 'Orientation Reproducibility (regression slope)';
% 
%             case 'spf_sel_p',     x = [allstats.spatialfrequencySelectivePval];     statLabel{i} = 'Spatial Frequency Selectivity (-log(p-value))';
%              case 'spf_sel_str',  x = [allstats.spatialfrequencySelectiveStrength]; statLabel{i} = 'Spatial Frequency Selectivity (strength)';
%             case 'spf_rep_p',     x = [allstats.spatFreqReproduciblePval];          statLabel{i} = 'Spatial Frequency Reproducibility (-log(p-value))';
%              case 'spf_rep_p_sgn',x = [allstats.spatFreqReproduciblePval_sgn];      statLabel{i} = 'Spatial Frequency Reproducibility (-log(p-value))';
%              case 'spf_rep_str',  x = [allstats.spatFreqReproducibleStrength];      statLabel{i} = 'Spatial Frequency Reproducibility (regression slope)';
%                 
%             case 'oriSpf_rep_p',          x = [allstats.responseReproduciblePval];     statLabel{i} = 'Reproducibility (-log(p-value))';
%              case 'oriSpf_rep_p_sgn',     x = [allstats.responseReproduciblePval_sgn]; statLabel{i} = 'Reproducibility (-log(p-value))';
%              case 'oriSpf_rep_str',       x = [allstats.responseReproducibleStrength]; statLabel{i} = 'Reproducibility (regression slope)';
%             case 'oriSpf_rep_p_top',      x = [allstats.responseReproduciblePval_top]; statLabel{i} = 'Reproducibility (-log(p-value))';
%              case 'oriSpf_rep_p_sgn_top', x = [allstats.responseReproduciblePval_sgn_top]; statLabel{i} = 'Reproducibility (-log(p-value))';
%              case 'oriSpf_rep_str_top',   x = [allstats.responseReproducibleStrength_top]; statLabel{i} = 'Reproducibility (regression slope)';
%             case 'oriSpf_rep_p_smth',     x = [allstats.responseReproduciblePval_smoothed]; statLabel{i} = 'Reproducibility (-log(p-value))';
%              case 'oriSpf_rep_p_sgn_smth',x = [allstats.responseReproduciblePval_sgn_smoothed]; statLabel{i} = 'Reproducibility (-log(p-value))';
%              case 'oriSpf_rep_str_smth',  x = [allstats.responseReproducibleStrength_smoothed]; statLabel{i} = 'Reproducibility (regression slope)';
%                 
%             case 'resp_size', x = [allstats.responseSizeFrac]; statLabel{i} = 'Response size (#std dev)';
% 
%             case 'cov_oo',   x = [allstats.crossCov_ori_ori];  statLabel{i} = 'orientation variance';
%              case 'cov_ss',  x = [allstats.crossCov_spf_spf];  statLabel{i} = 'spatial frequency variance';
%              case 'cov_os',  x = [allstats.crossCov_ori_spf];  statLabel{i} = 'ori / spf covariance';
%                  
%             otherwise, error('Unknown stat name');
%         end

%{
    [cur_x,cur_y, button] = ginput(1);
    while ibetween(cur_x, xlims) && ibetween(cur_y, ylims) && (gca == hM)

        if button == 1  % find closest data point, redo stats for that point
            da = get(gca, 'DataAspectRatio');
            da = da(2)/da(1);
            for i = 1:nR;
                if ~isempty(X{i} > 0)
                    [dist(i), ind_closestXY(i)] = min( normV([(X{i}(:)-cur_x), (Y{i}(:)-cur_y)/da], 2) );
                else
                    dist(i) = inf;
                end
            end            
            r = indmin(dist);
            idx = inds{r}(ind_closestXY(r));

            set(hCur, 'xdata', X{r}(ind_closestXY(r)), 'ydata', Y{r}(ind_closestXY(r)), 'markerEdgeColor', 'k', 'markerFaceColor', rankColors(r) );

            [Gid, cellId] = deal(allOSPs(idx).GroupId, allOSPs(idx).cellId);
            figure(123); clf; 
            if strcmp(gratingType, 'flashed')
                imageOSP(allOSPs(idx), 'mean:ph');
            elseif if strcmp(gratingType, 'drifting')
                imageOSP(allOSPs(idx), 'mean:ph');
                figure(124); clf; 
                imageOSP(allOSPs(idx), 'subplots:horizontal', 'SPO');
            end
            
            title(sprintf('Gid: %d, cellId: %d', Gid, cellId));

            redoStats_id = [Gid, cellId]; redoStats;

            figure(mainFig);

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
        end
        
        [cur_x,cur_y, button] = ginput(1);
    end
%}