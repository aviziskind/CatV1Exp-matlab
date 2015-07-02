function plotCCvsElectrodeDistance

    %%%%%%%%%%%%%%%%%%%%%%%% SET PARAMETERS %%%%%%%%%%%%%%
    th_toRemove = []; %10^(-9.1)
    showRemovedPointsInRed = false;
    printHighestOverlaps = false;

%     minFracRs = [0, .5, .6, .7, .8, .9];
    minFracRs = [0.0];
    corrType = 'spearman';
%     corrType = 'pearson';

    % Parameters of sliding window:
        window_frac1 = 1/15; % fraction of data in each window for plotting mean +-std 
        window_frac2 = 0.5; % fraction of data in each sliding window.

    fontsize = 10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    [gratingType, gratingType_s] = curGratingType;  % FLASHED_GRATING = 1; DRIFTING_GRATING = 2;
%     [cmpType, cmpType_s] = curCmpType;    
    cmpType = curCmpType('');
    [pt_ids, pairTypes] = curPairTypes;
%     pairTypes_str = [pairTypes{:}];
    
    ospDatafile = getFileName('osps');    
    pairDatafile = getFileName('pairs');    
    cmpDatafile = getFileName('comparisons');
    statsDatafile = getFileName('controls');
%     ospDatafile = [CatV1Path gratingType_s 'GratingCells' curCellsType '_' cmpType_s '.mat'];    
%     pairDatafile= [CatV1Path gratingType_s 'GratingPairs' curCellsType '_' cmpType_s '.mat'];        
%     cmpDatafile = [CatV1Path gratingType_s 'GratingComparisonData' curCellsType '_' cmpType_s '_' pairTypes_str '.mat'];             
    
    
    S1 = load(ospDatafile);        
    allCells = S1.allCells;
    allGids = [allCells.Gid];
    allCellIds = [allCells.cellId];

%     allCellClosenessMeasures = {'negAmps_dist', 'negAmp_overlap', 'fullWaveform_cc', 'fullWaveform_ed', 'channelWvfm_meanCC'};    
    
    featureName = 'negAmps_overlap';
    
    switch featureName
        case 'negAmps_dist',       xscale = 'log';    cellDistanceLabel = 'Euclidean Distance between 4D Negative Amplitudes';              label_short = 'Amp ED';
        case 'negAmps_overlap',    xscale = 'linear'; cellDistanceLabel = '-log( Overlap between fitted 4D Negative Amplitude Gaussians )'; label_short = 'Amp Overlaps';
        case 'fullWaveform_cc',    xscale = 'linear'; cellDistanceLabel = '1 - Correlation coefficient between concatenated waveforms'; label_short = 'Waveform CC';
        case 'fullWaveform_ed',    xscale = 'linear'; cellDistanceLabel = 'Euclidean Distance between concatenated waveforms';      label_short = 'Waveform ED';
        case 'channelWvfm_meanCC', xscale = 'linear'; cellDistanceLabel = '1 - Mean correlation coefficient between individual channel waveforms'; label_short = 'Wvfm Channel CC';
        case 'diffFWHM',           xscale = 'linear'; cellDistanceLabel = 'Difference in FWHM (Full-width at Half Max)';  label_short = 'dFWHM';
        case 'diffPtPWidth',       xscale = 'linear'; cellDistanceLabel = 'Difference in Peak-toPeak Width'; label_short = 'Wvfm Channel CC';
        case 'minID',              xscale = 'linear'; cellDistanceLabel = 'Smaller IsolationDistance of the two cells';         label_short = 'MinID';
        case 'maxL_ratio',         xscale = 'log';    cellDistanceLabel = 'Larger L_{ratio} of the two cells';         label_short = 'MaxLr';
    
    end
                
    switch cmpType
        case 'phase',   
%             comparison_measure = 'cc'; null_value = 0;
            comparison_measure = 'dphi'; null_value = 90;
            location = 'maxMinFracR';
%             location = 'p50MinFracR';
                        
        case 'degree',  comparison_measure = 'Dw_spf'; %'D_spf_pref';

    end
    
    switch comparison_measure
        case 'cc', cmp_measure_short = 'cc'; cmp_measure_long = 'Correlation Coefficient';
        case 'dphi', cmp_measure_short = '\Delta\phi'; cmp_measure_long = 'Peak of cross correlation';
        
        case 'Dw_ori_glob_si', cmp_measure_short = 'Global Ori Width'; cmp_measure_long = 'Global Orientation Tuning Width';
        case 'Dw_ori_loc_si', cmp_measure_short = 'Local Ori Width'; cmp_measure_long = 'Local Orientation Tuning Width';
        case 'D_dsi_glob_si', cmp_measure_short = 'DSI'; cmp_measure_long = 'Direction Selectivity Index (DSI)';
        case 'Dw_spf', cmp_measure_short = 'Sp. Freq Width'; cmp_measure_long = 'Spatial Frequency Tuning Width';        
    end
    
    
%     'Dw_ori_glob_si', 'Dw_ori_loc_si', 'Dw_ori_glob_ss', 'Dw_ori_loc_ss', 'D_dsi_si', 'D_dsi_ss', 'D_ori_pref', 'D_dir_pref', 'Dw_spf', 'D_spf_pref'};
        

    S2 = load(pairDatafile);        
    [Wcc_pairIdxs, Wcm_pairIdxs, Bcc_pairIdxs, Bcm_pairIdxs, Wrcc_pairIdxs, Wrcm_pairIdxs] = ...
        deal(S2.Wcc_idxs, S2.Wcm_idxs, S2.Bcc_idxs, S2.Bcm_idxs, S2.Wrcc_idxs, S2.Wrcm_idxs);   
%     [pairIdxs, pairIdxList] = useSubsetOfPairIdxs(S2, pairTypes);    
            
%     flds_pd = fieldnames(pairData(1));    

    S3 = load(cmpDatafile); 
    [pairData, S, pairTypes, measures, locations] = deal(S3.pairData, S3.allStatsC, S3.pairTypes, S3.measures, S3.locations);
    [pairIdxs, pairIdxList, idxMtx] = useSubsetOfPairIdxs(S2, pairTypes, length(allCells));        
    measure_idx = find(strcmp(comparison_measure, measures), 1);
    
   
    switch cmpType
        case 'phase',  location_idx = find(strcmp(location, locations));
        case 'degree', location_idx = 1;
    end
    
    Wcc_idxs = idxMtx(Wcc_pairIdxs);
    M_all = S{location_idx, measure_idx}.val( Wcc_idxs );  
    D_all = [pairData.(featureName)(Wcc_idxs)];        
        
    idx_ok = ~isnan(M_all) & ~isnan(D_all);  
    
        
    M = M_all(idx_ok);
    D = D_all(idx_ok);    
    Wcc_idxs_ok = Wcc_idxs(idx_ok);
   
    if exist(statsDatafile, 'file')
        S_controls = load(statsDatafile);
        M_control = S_controls.(comparison_measure).vals_Bcc;
                
    end
    
%     nPairs = length(pairData);
    
    D = D(:);
                    

%     figure(500); clf;
%     h_ax(1) = subplot(1,2,1);
%     h_wvfm(1,1:4) = plot(ones(2,4),ones(2,4), 'o-'); hold on;
%     h_ax(2) = subplot(1,2,2);
%     h_wvfm(2,1:4) = plot(ones(2,4),ones(2,4), 'o-'); hold on;
%     curPairIdx = 100;
    
    
%     h_mn = plot(0,0, 'ko-');
        
%     minAbove0 = min(D(D>0));        
% %     D(idxTooSmall) = 10^(min_logx);
    
    M_orig = M;
    D_orig = D;
%     logD_orig = log10(D_orig);
    3;
    if strcmp(xscale, 'log');
        D_orig = log(D_orig);
    end
        
    minD = min(D_orig);
    maxD = max(D_orig);        

    plotPcc_vs_Dist = 1;
    plotPcc_MvsD    = 0;    
    
    for mf_i = 1:length(minFracRs)
        %%
        minFracR = minFracRs(mf_i);
        
        M = M_orig;
        D = D_orig;        
        

        if ~isempty(th_toRemove)        
            idx_rm = find( D > th_toRemove );    
        else
            idx_rm = [];    
        end

        if minFracR > 0
            idx_rm = [pairData.loc_minFracOfMaxes(Wcc_idxs_ok, 1)] < minFracR;
        else
            idx_rm = [];
        end
        
        if printHighestOverlaps
            idx_high = ord(D_orig, 'descend');
            for i = 1:30
                pd = pairData(idx_high(i));
                fprintf('%d) Gid: %d, cells %d, %d. (overlap = %.2g)\n', i, pd.Gids(1), pd.cellIds, D_orig(idx_high(i)) );
            end
            3;
        end


        % plot data
        figure(sum(pt_ids)+gratingType*10); clf; hold on; box on;
        plot(D, M, 'o'); 
        if showRemovedPointsInRed && ~isempty(idx_rm)
            M_rm = M_orig(idx_rm);
            D_rm = D_orig(idx_rm);
            plot(D_rm, M_rm, 'ro'); 
        end
        M(idx_rm) = [];
        D(idx_rm) = [];

        nWccPairs_i = length(M);
        
%         set(gca, 'xscale', xscale);
        drawHorizontalLine(0); % 

        % calc M.            
        X = D; Y = M;                 
        [rho,p_rho] = corr(X, Y, 'rows', 'complete', 'type', corrType);
        if minFracR == 0
            mfr_str = '';
        else            
            mfr_str = sprintf('minFracR > %.2f, ', minFracR);
        end
        grating_n_str = [titleCase(gratingType_s) ' Gratings (' mfr_str 'N = ' num2str(nWccPairs_i) ')'];
        
        
        % Add to plot #1: divide up into bins, and plot mean /stderr of each bin             
        
        nWindows = round(1/window_frac1);
        windowRanges = round( linspace(0, 1, nWindows+1)*nWccPairs_i );
        windowIdxs = arrayfun(@(i1, i2) i1:i2, windowRanges(1:end-1)+1, windowRanges(2:end), 'un', 0);

        D_binCenters = zeros(1, nWindows);
        M_mean = zeros(1, nWindows);
        M_std = zeros(1, nWindows);
                
        D_pct95 = prctile(D, 99);
        xlim([minD, D_pct95]);

        idx_ord = ord(D, 'ascend');
        for wi = 1:nWindows
            w_idx = idx_ord ( windowIdxs{wi} );
            D_binCenters(wi) = mean(D(w_idx));
            M_mean(wi) = nanmean(M(w_idx));
            M_std(wi)  = nanstderr(M(w_idx));
        end
        hold on; box on;
        errorbar(D_binCenters, M_mean, M_std, 'r.-');
%         plot(,, 'ro-', 'markerfacecolor', 'r' )
%         plot(M_inBins_m+, 10.^(ov_binsC), 'gs', 'markerfacecolor', 'g' )
        plot(D_binCenters, M_mean-M_std, 'gs', 'markerfacecolor', 'g', 'markersize', 3 );
        plot(D_binCenters, M_mean+M_std, 'gs', 'markerfacecolor', 'g', 'markersize', 3 )
        binSize = round(mean(cellfun(@length, windowIdxs)));
        binSize_str = sprintf('Bin size = %d points (%.0f%%)', binSize, binSize/nWccPairs_i*100);        
        h_cur = plot(0,0, 'rs', 'markerfacecolor', 'k');

%         
%         s_cc = sprintf('r = %.2f. p = %.2g', cc, p_cc);
        s_rho = sprintf('\\rho = %.2f. (p = %.2g)', rho, p_rho);        
%         if ~isempty(th_toRemove)
%     %         str_intro = iff(showRemovedPointsInRed, '(without red points):', '');
%     %         s3 = sprintf('%s r = %.2f. p = %.2g', str_intro, cc2, p_cc2);
%     %         s3 = sprintf('(without red points): r = %.2f. p = %.2g', cc2, p_cc2);
%             xyfig_tit_str = {cmp_measure_long, grating_n_str, s_rho};
%         else
%             tit_str = {s0_fig1, s_rho};
%         end
        xyfig_tit_str = {cmp_measure_long, grating_n_str, s_rho};
                
        xlabel(cellDistanceLabel, 'fontsize', fontsize);
        prefix_str = iff(strcmp(cmpType, 'degree'), 'Difference in ', '');
        ylabel( sprintf('%s%s', prefix_str, cmp_measure_short), 'interpreter', 'none', 'fontsize', fontsize ); 
        h = title(xyfig_tit_str, 'fontsize', fontsize);
%         set(h, 'fontsize', 12)


        %%% Plot #2 : show continuous relationship between M, and cc of
        %%% overlap vs phase-cc
        
        switch cmpType
            case 'phase',                
                switch comparison_measure 
                    case 'cc', sigTest = @(Wcc_vals)   pvalTtest(Wcc_vals, 0);
                    case 'dphi', 
                        sigTest = @(Wcc_vals) pvalTtest(Wcc_vals, 90);
%                         sigTest = @(Wcc_vals) pvalTtest(Wcc_vals, 90);
                        
                end
            case 'degree',                               
                sigTest = @(Wcc_vals)  pvalUtest(Wcc_vals, M_control);                                
                
        end
                
        nPerWindow = floor(nWccPairs_i*window_frac2);
        nWindows = nWccPairs_i-nPerWindow+1;
        windowIdxs = arrayfun(@(i1) [i1:i1+nPerWindow-1], 1:nWindows, 'un', 0);

        [cc_MvsD, pcc_MvsD, p_ccDist] = deal( zeros(1, nWindows) );        
        idx_ord = ord(D, 'ascend');
        for w_i = 1:nWindows
            w_idx = idx_ord ( windowIdxs{w_i} );
            
            [cc_MvsD(w_i), pcc_MvsD(w_i)] = corr(M_orig(w_idx), D_orig(w_idx),    'rows', 'complete', 'type', corrType);
%                     [h_tmp, p_ccDist(w_i)] = ttest(M_orig(w_idx));
            p_ccDist(w_i) = sigTest( M_orig(w_idx) );
        end
        figure(150+mf_i + gratingType*10); clf; hold on; box on;
        leg_strs = {};
        h_plot = [];
        if plotPcc_vs_Dist
            wind_fracs = [1:nWindows]/nWindows;
            h_plot(1) = plot(wind_fracs, -log10(p_ccDist), 'b-', 'linewidth', 2); 
            leg_strs = {  sprintf('-log(p[%s])', cmp_measure_short) };         
        end
            

        if plotPcc_MvsD
            idx_neg = cc_MvsD<0;                        
            idx_neg_C = continuousGroupings(find(idx_neg));

            h_MvsD(1) = plot(wind_fracs, -log10(pcc_MvsD), 'b-', 'linewidth', 2);
            h_plot(end+1) = h_MvsD(1);
            if nnz(idx_neg) > 1
                for i = 1:length(idx_neg_C)
                    h_MvsD(1+i) = plot(wind_fracs(idx_neg_C{i}), -log10(pcc_MvsD(idx_neg_C{i})), 'r', 'linewidth', 2); 
                end
            end
            corrLegends = {sprintf('-log(p[%s vs. %s])', cmp_measure_short, label_short)};
            leg_strs = [leg_strs, corrLegends];
        end
        nPairsLeg = {};
        h(mf_i) = gca;

        xlabel(sprintf('window range over values of %s', label_short))
%                 xlim([0 1]);              
        window_size_str = sprintf('windowSize = %d points (%.0f%%)', nPerWindow, window_frac2*100);
                        
        
%         legend(h_plot, leg_strs, 'location', 'NE', 'interpreter', 'none')
%         s2 = sprintf( 'minFracR > %.2f', minFracR);
        title({cmp_measure_long, window_size_str}, 'fontsize', fontsize);
        drawHorizontalLine(-log10(.01));
        drawHorizontalLine(-log10(.05), 'linestyle', ':', 'color', 'k');
        ylabel('-log_{10}( p_{cc} )', 'fontsize', fontsize);
        3;

    end
    matchAxes('Y', h);
   
%     set(1, 'WindowKeyPressFcn', @updateHelper);
    
    function updateHelper(src, evnt)
        if isfield(evnt, 'Modifier') && any(strcmp( evnt.Modifier , 'alt'))
            mult = 10;
        else
            mult = 1;
        end
        
        if isfield(evnt, 'Key') 
            switch evnt.Key
                case 'home',      curPairIdx = 1;
                case 'end',       curPairIdx = nWccPairs_i;
                case 'uparrow',   curPairIdx = curPairIdx-1*mult;
                case 'downarrow', curPairIdx = curPairIdx+1*mult;
                case 'pageup',    curPairIdx = curPairIdx-100*mult;
                case 'pagedown',  curPairIdx = curPairIdx+100*mult;
            end
        end
        
        curPairIdx = min(curPairIdx, nWccPairs_i);
        curPairIdx = max(curPairIdx, 1);
        updateWvfms;
    end
    
    
    function updateWvfms
        pairIdx = idx_ord(curPairIdx);
        
        set(h_cur, 'xdata', 10.^logD_orig(pairIdx), 'ydata', M(pairIdx));
        
        gids = pairData(pairIdx).Gids;
        cellids = pairData(pairIdx).cellIds;
                
        Gid = gids(1);
        
%         cell_idx2 = find(allGids == Gid & allCellIds == cellids(2), 2);
        
        
%         spk = getSpikes(Gid);
%         cellIds = spk(:,2);
%         wvfm_file = ['C:\ExperimentDB\Spikes\spiker\Group_' num2str(Gid) '_spkWaveforms_spkr.mat'];
%         S = load(wvfm_file);
%         allWvfms = S.spikeWaveforms;
%         
%         [M,C] = getChannelMeansAndCovariance(Gid, 1);
% %         [V,D] = eig(C);        
%         chnl_scl = [1./sqrt(diag(D))]';
        
%         wvfms = cell(1,2);
        for ci = 1:2
            cell_idx = find(allGids == Gid & allCellIds == cellids(ci), 1);
            wvfms_i = allCells(cell_idx).spkFeatures.meanWaveform;
%             wvfms_i = mean( allWvfms(:,:, idx_cell), 3);
%             wvfms_i = bsxfun(@minus, wvfms_i, M');
%             wvfms{ci} = bsxfun(@times, wvfms_i, chnl_scl);                        
            for chnl_i = 1:4
                set(h_wvfm(ci,chnl_i), 'xdata', [-.25:.05:1.3], 'ydata', wvfms_i(:,chnl_i));
            end            
            set(h_ax, 'xlim', [-.25, 1.3])
        end        
                
%         ylims = ylims + [-1, 1]*diff(ylims)/30;
%         set(h_wvfm_ax, 'ylim', ylims, 'xlim', spkWindow-spk_start);
%         set(h_horLine, 'xdata', spkWindow-spk_start);
%         set(h_spkLine(1), 'xdata', 0*[1; 1] , 'ydata', ylims);        
%         set(h_spkLine(2), 'xdata', double(spk_end-spk_start)*[1; 1], 'ydata', ylims, 'visible', 'on');        
%         set(h_tit, 'string', sprintf('th = %.2f. spk idx = %d', Xdata(spkIdx), idx));
        
    end
  
    3;

end

function p = pvalTtest(X, nullVal)
    [~, p] = ttest(X, nullVal);
end

function p = pvalUtest(X, Y)
    p = ranksum(X,Y);
    
end


function [pairIdxs, pairIdxList, idxMtx] = useSubsetOfPairIdxs(allPT, pairTypes, nUnits)
    function y = catIfCell(C)
        if iscell(C)
            y = [C{:}];
            y = y(:);
        else
            y = sort(C(:));
        end
    end
    if isstruct(allPT)
        [Wcc_pairIdxs, Wcm_pairIdxs, Bcc_pairIdxs, Bcm_pairIdxs, Wrcc_pairIdxs, Wrcm_pairIdxs] = ...
        deal(allPT.Wcc_idxs, allPT.Wcm_idxs, allPT.Bcc_idxs, allPT.Bcm_idxs, allPT.Wrcc_idxs, allPT.Wrcm_idxs);   
        allPairIdxs = {Wcc_pairIdxs, Wrcc_pairIdxs, Bcc_pairIdxs, Wcm_pairIdxs, Wrcm_pairIdxs, Bcm_pairIdxs};
    elseif iscell(allPT)
        assert(length(allPT) == 6);
    end
    allPairTypes = {'Wcc', 'Wrcc', 'Bcc',  'Wcm', 'Wrcm', 'Bcm'};
    
    % hack to save memory
%     if ~any(strcmp(pairTypes, 'Wrcc'))
%         allPairIdxs{ find( strcmp(allPairTypes, 'Wrcc'), 1) }  = [];        
%     end    
%     allPairIdxs_list = cellfun(@catIfCell, allPairIdxs, 'un', 0);
    
    % remove pairTypes that are not available
    pairTypesAvailable = cellfun(@(pr) any ( strcmp(allPairTypes, pr)), pairTypes);
    pairTypes = pairTypes(pairTypesAvailable);
    
    % reorder pairTypes to proper order
    pairTypes = pairTypes(  ord(cellfun(@(s) find(strcmp(s, allPairTypes)), pairTypes)) );
    
    % find which pairTypes are requested.
    pairTypesRequested = cellfun(@(pr) any ( strcmp(pairTypes, pr)), allPairTypes);
    
    pairIdxs = allPairIdxs(pairTypesRequested);
        
%     pairIdxs = cellfun( @(pr) allPairIdxs{ strcmp(allPairTypes, pr) },  pairTypes, 'un', false);

%     allPairIdxs_list = cellfun(@catIfCell, allPairIdxs, 'un', 0);
%     pairIdxList = unique(cat(1, allPairIdxs_list{pairTypesRequested}));

    pairIdxs_list = cellfun(@catIfCell, pairIdxs, 'un', 0);
    pairIdxList = unique(cat(1, pairIdxs_list{:}));
    
    
    idxMtx = zeros(nUnits, nUnits, 'uint32');
    idxMtx(pairIdxList) = 1:length(pairIdxList);

end






%{
switch slidingThreshType
            case 'threshold',                    
                ov_bins = linspace(minD, maxD, 30);
                ov_binsC = binEdge2cent(ov_bins);
                exp_ovBinsC = 10.^ov_binsC;
                M_inBins_m = zeros(1, length(ov_bins)-1);
                M_inBins_s = zeros(1, length(ov_bins)-1);
                for bi = 1:length(ov_bins)-1
                    idx = find(  ov_bins(bi) < logD & logD <= ov_bins(bi+1) );
                    M_inBins_m(bi) = mean(M(idx));
                    M_inBins_s(bi) = stderr(M(idx));
                end
                hold on;  box on;           
                errorbar(exp_ovBinsC, M_inBins_m, M_inBins_s, 'r.-');
        %         plot(,, 'ro-', 'markerfacecolor', 'r' )
        %         plot(M_inBins_m+, 10.^(ov_binsC), 'gs', 'markerfacecolor', 'g' )
                plot(exp_ovBinsC, M_inBins_m-M_inBins_s, 'gs', 'markerfacecolor', 'g', 'markersize', 3 )
                plot(exp_ovBinsC, M_inBins_m+M_inBins_s, 'gs', 'markerfacecolor', 'g', 'markersize', 3 )


        switch slidingThreshType
            case 'threshold', 
                ths = 0:-.1:-30;
                pcc_MvsD_log = zeros(1, length(ths));
                pcc_MvsD_lin = zeros(1, length(ths));
                p_ccDist = zeros(1, length(ths));
                nPairs_i   = zeros(1, length(ths));
                for th_i = 1:length(ths)
                    idx_ok = find(logD < ths(th_i));
                    nPairs_i(th_i) = length(idx_ok);
                    if nPairs_i(th_i) > 0
                        [~, pcc_MvsD_log(th_i)] = corr(M_orig(idx_ok), logD_orig(idx_ok), 'rows', 'complete', 'type', corrType);
                        [~, pcc_MvsD_lin(th_i)] = corr(M_orig(idx_ok), D_orig(idx_ok), 'rows', 'complete', 'type', corrType);
                                                
                        p_ccDist(th_i) = sigTest( M_orig(idx_ok) );
                    end

                end
                figure(150+mf_i); clf; hold on;
                if strcmp(corrType, 'pearson')
                    plot(ths, -log10(pcc_MvsD_log), 'b-');    
                    plot(ths, -log10(pcc_MvsD_lin), 'g-'); 
                elseif strcmp(corrType, 'spearman')
                    plot(ths, -log10(pcc_MvsD_lin), 'b-'); 
                end
                plot(ths, -log10(p_ccDist), 'r-'); 
                plot(ths,  log10(nWccPairs_i), 'k:'); 
                h(mf_i) = gca;

                xlabel('log_{10}threshold')
%                 xlim([-30 0]);
                nPairsLeg = {'log(Npairs)'};


%}