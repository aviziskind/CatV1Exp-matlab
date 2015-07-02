function fig_spatialFreqTuningCurves

    opt.spfSmoothW = 0.6;
    opt.maxPeakFitSpfTC = 1.5;
    opt.minSpfWparam = 0.09;
    opt.maxSpfSparam = 2;
    opt.useWeightsIfAllPositive = 0;
    
    xlabel_fsize = 12;
    ylabel_fsize = 12;
    title_fsize = 11;
    legend_fsize = 8;

%     label_fontSize = 13;
%         legend_fontsize = 10;
%         axis_fontSize = 11;
    
    %%
%     741, 4 **
%     486,1 **
%     Gids_plot = [438, 3405, 2877, 2905, 771];
%     cellIds_plot = [4, 1, 7, 2, 3];
    Gids_plot    = [5000, 703, 486, 741, 771];
    cellIds_plot = [1,    2,    1,   4,    3];
    
%     Gids_plot = [2034        2629        3045        3075        3089        4087        5166];
%     cellIds_plot = [2     1     2     8     6     2     2];
%      Gids_plot    = [701 714 714 716 830 1140 2034 2629 2739 2863 3045 3075 3089 4001 4001 4087 4504 4540 4624 5000 5050 5148 5148 5166 5172 5172 5272]; 
%      cellIds_plot = [ 2  3 6 1 2 1 2 1 3 2 2 8 6 1 3 2 3 4 1 1 1 1 6 2 1 3 2];

     % 5, 20
%      idx = [5, 20, 21, 22];
%      Gids_plot = Gids_plot(idx);
%      cellIds_plot = cellIds_plot(idx);
     
%%
    testMode = 0;
    S_d = load('driftingGratingCells_GLFcuw8_degree_SS_spf.mat');
    allGids_d = [S_d.allCells.Gid];
    allCellIds_d = [S_d.allCells.cellId];
    
    
    S_f = load('flashedGratingCells_GLFcuw8_degree_SS_spf.mat');
    allGids_f = [S_f.allCells.Gid];
    allCellIds_f = [S_f.allCells.cellId];
    
%%
    if 0
        %%
        allStats_d = nestedFields(S_d, 'allCells', 'tuningStats', 'spfStats_ss', 1);
        allStats_f = nestedFields(S_f, 'allCells', 'tuningStats', 'spfStats_ss', 1);
%         allStats = allStats_f;
        allStats = allStats_d;
        
        allSLNp = [allStats.SLNparams];
        maxStd = arrayfun(@(s) max(s.spf_tc ./ s.spfTC_stderr), allStats);
        rsqr = [allStats.rsqr];
        w_spf = [allStats.w_spf];
        s = [allSLNp.s];
        
        
    end
    
    subN = 2;
    subM = ceil(length(Gids_plot)/subN);
    
    if testMode
        idx_ok =  [    10    14    43    50    66    68    78   115   119   123   151   159   167   210   212   215   217   241   245   246   248   251];
        Gids_plot = allGids_f(idx_ok);
        cellIds_plot = allCellIds_f(idx_ok);
        
%         Gids_plot = allGids_d(idx_ok);
%         cellIds_plot = allCellIds_d(idx_ok);
        
%     allStats_f = nestedFields(S_d, 'allCells', 'tuningStats', 'spfStats_ss', -1);
        subN = 3;
        subM = ceil(length(Gids_plot)/subN);

    end
    
  %%      
%     subM = 4;
    fig_id = 2; 
    figure(fig_id); clf;
    
    m_spacing = [.03, .02, .03];
    n_spacing = [.03, .01, .05];
    %%
    showLegendOn = subN;
    doUnconstrainedFit = 1 && ~testMode;
    
    for i = 1:length(Gids_plot)
        doUnconstrainedFit_now = i == length(Gids_plot) && doUnconstrainedFit;
        showLegend = i == showLegendOn;
        %%
        Gid = Gids_plot(i);
        cellId = cellIds_plot(i);
        
        if any(allGids_d == Gid)
            gType = 'drifting';
            idx = find( allGids_d == Gid & allCellIds_d == cellId, 1);
            cellData = S_d.allCells(idx);
        else
            gType = 'flashed';
            idx = find( allGids_f == Gid & allCellIds_f == cellId, 1);
            cellData = S_f.allCells(idx);
        end        
        
        stats = cellData.tuningStats.spfStats_ss;
        sd = siteDataFor(Gid);

        spfs_cpd = 1./(sd.spPeriod_pix * sd.stimulusInfo.degreesPerBlock);                
        [spfs_cpd, idx_cpd] = sort(spfs_cpd(:));
        spfTC = stats.spf_tc;
        spfTC = spfTC(idx_cpd);
        spfTC = spfTC(:);        
        
        
        spfTC_std = zeros(size(spfTC));
        if isfield(stats, 'spfTC_stderr')
            spfTC_std = stats.spfTC_stderr;
            spfTC_std = spfTC_std(idx_cpd);
        end
            
%%
        logBase = 2;
        logBase_factor = log(logBase);        
        
        log_spfs_cpd = log(spfs_cpd);
        
        log_spfs_cpd_fine = linspace(log_spfs_cpd(1), log_spfs_cpd(end), 2000);
        
        logBase_spfs_cpd = log(spfs_cpd)/logBase_factor;
        logBase_spfs_cpd_fine = linspace(logBase_spfs_cpd(1), logBase_spfs_cpd(end), 2000);

%%        
        
        [result_txt, w_spf, f_peak, normBw, spfParams, spfParams_ci, fit_stats] = ...
            getSLNfit(spfs_cpd, spfTC, spfTC_std, opt, 1, 0);  
        cfun = fit_stats.cfun;
        spfTC_fine = feval(cfun, log_spfs_cpd_fine);
        
        if doUnconstrainedFit_now
            noConstraint = 1;    
            [result_txt, w_spf, f_peak, normBw, spfParams, spfParams_ci, fit_stats_uc] = ...
                getSLNfit(spfs_cpd, spfTC, spfTC_std, opt, 1, 0, noConstraint);   
            cfun_uc = fit_stats_uc.cfun_uc;
            spfTC_fine_uc = feval(cfun_uc, log_spfs_cpd_fine);
        end
        
          %%
%           figure(68); clf; hold on;
%         logSLN_func = fittype( 'skewLogNormal_exp(x, f_opt_log, r_max, r_bkg, w, s)' );
          
        ylims = lims([spfTC + spfTC_std; spfTC - spfTC_std; spfTC_fine], .05);

        xlims_log = lims(logBase_spfs_cpd, .02);
         
        lastPlotInCenter = i == length(Gids_plot) && mod(i,2) ~= 0;
        if lastPlotInCenter        
            n_spacing_last = [0.25 + n_spacing(1), 0, .25+n_spacing(end)];
            subplotGap(subM, 1, subM, [], m_spacing, n_spacing_last);
        else
            subplotGap(subM, subN, i, [], m_spacing, n_spacing);
            
        end
        
%         subplotGap(subM, subN, i, [], m_spacing, n_spacing);
%         plot(log_spfs_cpd, spfTC, 'bs'); hold on
        h_plot(1) = errorbar(logBase_spfs_cpd, spfTC, spfTC_std, 'ko');
        hold on;
%         plot(log_spfs_cpd, spfTC_sm, 'g.'); 
        xlim(xlims_log);
        ylim(ylims);
        %%
        xticks = ceil(xlims_log(1)) : 1 : floor(xlims_log(end));
        set(gca, 'xtick', xticks);  

        if ~testMode
            xlabel(sprintf('Log_{%d} sp. frequency (cyc/deg)', logBase), 'fontsize', xlabel_fsize);
            ylabel('Firing rate (Hz)', 'fontsize', ylabel_fsize);         
        else
            set(gca, 'xtick', [])
        end
        
%         set(gca, 'xtick', xticks, 'xticklabel', arrayfun(@(n) cellstr( num2str(2.^xticks(:)) ) );  
%         set(gca, 'xtick', ceil(xlims_log(1)) : 1 : floor(xlims_log(end)))
        
%         set(gca, 'XTickLabel',[])                      %# suppress current x-labels
%         yl = get(gca, 'YLim');
%         xt = get(gca, 'XTick');
%         str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
%         hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
%             'Interpreter','tex', ...                   %# specify tex interpreter
%             'VerticalAlignment','top', ...             %# v-align to be underneath
%             'HorizontalAlignment','center');           %# h-aligh to be centered
        %%
        
        spfTC_fine = feval(cfun, log_spfs_cpd_fine);
        
        if ~doUnconstrainedFit_now
            h_plot(2) = plot(logBase_spfs_cpd_fine, spfTC_fine, 'k');         
        
        else
            h_plot(2) = plot(logBase_spfs_cpd_fine, spfTC_fine, 'k');     
            h_plot(3) = plot(logBase_spfs_cpd_fine, spfTC_fine_uc, 'k:');     
            
        end
        
        log_f_lo = log(fit_stats.f_lo)/logBase_factor;
        log_f_hi = log(fit_stats.f_hi)/logBase_factor;
        
        r_peak = cfun.r_max;
        r_bckg = cfun.r_bkg;  
    
        r_halfMax = (r_peak-r_bckg)/2 + r_bckg;
        if ~doUnconstrainedFit_now
        
            h_plot(3) = plot([log_f_lo, log_f_hi], r_halfMax+[0, 0], 'k.:');
        end
        
        if ~testMode
    %         stat_str = sprintf('f_{opt} = %.2f cyc/deg.  w_{SF} = %.2f.   R^2 = %.2f',  f_peak, w_spf, fit_stats.rsquare);
            stat_str = sprintf('f_{pref} = %.2f cyc/deg.  w_{SF} = %.2f',  f_peak, w_spf);
            h_t = title({stat_str}, 'fontsize', title_fsize);
    %         if exist('log_f_lo', 'var')
    %             plot([log_f_lo, log_f_hi]/logBase_factor, r_halfMax+[0, 0], 'k.:');
    %         end

            %%
            if showLegend
                drawnow;
%                 legend(h_plot, {'a', 'b', 'c'});
                legend(h_plot, {'Tuning Curve', 'SLN Fit', 'Half-max'});
                drawnow;
                h_legend_main = legend(h_plot, {'Tuning Curve', 'SLN Fit', 'Half-max'}, 'location', 'best', 'fontsize', legend_fsize);
            end
        else
            title(sprintf('[%d,%d] s = %.2f. w = %.2f', Gid, cellId, stats.SLNparams.s, stats.w_spf), 'fontsize', 8);
        end
%         axis tight;
        
        if doUnconstrainedFit_now
            %%
%             coeffVals_uc = coeffvalues(cfun_uc);
%             coeff_uc_exp_C = num2cell(coeffVals_uc); coeff_uc_exp_C{1} = exp( coeff_uc_exp_C{1} );
%             cfun_uc_exp = cfit(SLN_func, coeff_uc_exp_C{:});    
            
%             h = [];
%             xlims = lims(spfs_cpd, .05, [], 1);
%             figure(16); clf; hold on; box on;
%             spfs_cpd_tmp = 10.^(log_spfs_cpd_tmp); % linspace(spfs_cpd(1), spfs_cpd(end), 100);        
%             h(1) = plot(spfs_cpd, spfTC, 'bs'); 
%             h(1) = errorbar(log10_spfs_cpd, spfTC, spfTC_std, 'bo');
%             h(2) = plot(log10_spfs_cpd_tmp, feval(cfun, log_spfs_cpd_tmp), 'r-');  
%             h(3) = plot(log10_spfs_cpd_tmp, feval(cfun_uc, log_spfs_cpd_tmp), ['-'], 'color', [0 .6 0] );              
            r_peak = cfun.r_max;
            set(gca, 'xlim', xlims_log);
            y1 = min(0, min(spfTC-spfTC_std*1));
            ylim([y1, r_peak*1.5]);
            
            drawnow;
%                 legend(h_plot, {'a', 'b', 'c'});
            legend(h_plot, {'Tuning Curve', 'SLN Fit', 'Half-max'});
            drawnow;
            h_legend_unc = legend(h_plot, {'Tuning Curve', 'Constrained Fit', 'Unconstrained Fit'}, 'location', 'best', 'fontsize', legend_fsize);                      
            title({stat_str}, 'fontsize', title_fsize);
            3;
            
        end
        
        3;
        if lastPlotInCenter
            h_let(i) = addSubplotLetter(subM, 1, subM, [], m_spacing, n_spacing_last, char('A'+i-1), [], 'fontname', 'Helvetica');
        else           
            h_let(i) = addSubplotLetter(subM, subN, i, [], m_spacing, n_spacing, char('A'+i-1), [], 'fontname', 'Helvetica');
        end

    end
    set(fig_id, 'windowStyle', 'normal')
    3;
    %%
    
    figureFolder = [CatV1Path 'Figures' filesep 'DegreePaper' filesep];
    fig_filename = sprintf('%sFigure2_spfExampleFits.pdf', figureFolder);
    set(fig_id, 'color', 'w', 'position', [1000, 250, 720, 700]);
    %%
    export_fig(fig_id, 'pdf', fig_filename);

    3;




end


%{
 good examples of standard (slightly narrow):
    448, 1
    703, 2 **

    2919, 22

 good examples of skew-left:
    1126, 5
    790,4
    741, 4 **


good examples of skew-right:
    486,1 **
    
    4870,1

%}

