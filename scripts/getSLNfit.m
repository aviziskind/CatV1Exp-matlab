function [result_txt, w_spf, f_peak, normBw, spfParams, spfParams_ci, fit_stats, stcProps] = ...
    getSLNfit(spfs_cpd, spfTC, spfTC_std, opt, cellOK, show_flag, noConstraint)    

    global GC

    show = 0;
%     show = 1;
    dbug = (nargin >= 6) && isequal(show_flag, 1);
    if nargin < 7
        noConstraint = 0;
    end
    
    err_str = '';
    
    constrainFpeak = 0;
    constrainW = 1;
    constrainS = 1;
    constrainRmax = 1;   
    
%     if noConstraint
%         [constrainFpeak, constrainW, constrainS, constrainRmax] = deal(0);        
%     end    
%    
        
    doUnconstrainedFit = noConstraint;
    
    if dbug
        show = 1;
        err_str = show_flag;
    end
%         Fitting stopped because the number of iterations or function evaluations exceeded the specified maximum.
    normBw = nan;
    w_spf = nan;        
    f_peak = nan;
    spfParams    = struct('f_opt', nan, 'r_max', nan, 'r_bckg', nan, 'w', nan, 's', nan);
    spfParams_ci = spfParams;
    fit_stats = struct('rsquare', nan);
    stcProps = struct('MaxAtLowest', nan, 'MaxAtHighest', nan, ...
        'PrefSpfLowerThanStim', nan, 'PrefSpfHigherThanStim', nan, ...
        'FracAtLowest_fit', nan, 'FracAtHighest_fit', nan, 'FracAtLowest_data', nan, 'FracAtHighest_data', nan, ...
        'FracAtZero', nan, 'FracAtInf', nan, 'PrefDistFromStim', nan, 'HalfMaxDistFromStim', nan, ...
        'nPeaks_itp', nan, 'nFits', nan);
    
    if ~cellOK
        result_txt = 'not reproducible';        
        return;
    end
    %%
    [spfs_cpd, idx_cpd] = sort(spfs_cpd(:));
    spfTC = spfTC(idx_cpd);         spfTC = spfTC(:);
    spfTC_std = spfTC_std(idx_cpd); spfTC_std = spfTC_std(:);
    
    [maxSpfTC, indMaxSpfTC] = max(spfTC);
    minSpfTC = min(spfTC);
    
    
    stcProps.MaxAtLowest  = indMaxSpfTC == 1;
    stcProps.MaxAtHighest = indMaxSpfTC == length(spfTC);
       
    
    % 1a. Estimate basic parameters (peak, background) of tuning curve
    log_spfs_cpd = log(spfs_cpd);
    dw = mean(diff(log_spfs_cpd));
    w_factor = opt.spfSmoothW;
    spfTC_sm = gaussSmooth_nu(log_spfs_cpd, spfTC, dw * w_factor);
    log_spfs_cpd_itp = linspace(log_spfs_cpd(1), log_spfs_cpd(end), 100);
    spfTC_sm_itp = interp1(log_spfs_cpd, spfTC_sm, log_spfs_cpd_itp, 'spline');
    
    if show
        figure(5); clf;
        plot(log_spfs_cpd, spfTC, 'bs'); hold on
        errorbar(log_spfs_cpd, spfTC, spfTC_std, 'bs');
%         plot(log_spfs_cpd, spfTC_sm, 'g.'); 
        plot(log_spfs_cpd_itp, spfTC_sm_itp, 'g');
        xlabel('Log Spatial Frequency (cyc/deg)');
        ylabel('Firing Rate (Hz)');
        
%         fplot(@(x) logSLNfunc(beta0, x), log_spfs_cpd([1, end]), 'm:');               
        
        figure(6); clf;
        plot(spfs_cpd, spfTC, 'bs'); hold on
        he = errorbar(spfs_cpd, spfTC, spfTC_std, 'bs-');
        xlabel('Spatial Frequency (cyc/deg)');
        ylabel('Firing Rate (Hz)');
%         plot(spfs_cpd, spfTC_sm, 'r.'); 
%         plot(exp(log_spfs_cpd_itp), spfTC_sm_itp, 'r');
%         fplot(@(x) SLNfunc(beta0_exp, x), spfs_cpd([1, end]), 'm:');               
        3;
    end
    
    % FIT (curve fitting)
    logSLN_func = fittype( 'skewLogNormal_exp(x, f_opt_log, r_max, r_bkg, w, s)' );
    SLN_func = fittype( 'skewLogNormal(x, f_opt, r_max, r_bkg, w, s)' );
    
    all_idxMax_sm = findLocalMaxima(spfTC_sm_itp, 1, [], [], 1);        

    [all_idxMax] = findLocalMaxima(spfTC, 1, [], [], 1);        
    idx_rm = [];
    if (all_idxMax_sm(1) == 1) && (all_idxMax(1) ~= 1)
        idx_rm = [1];
    end
    if (all_idxMax_sm(end) == length(log_spfs_cpd_itp)) && (all_idxMax(end) ~= length(log_spfs_cpd))
        idx_rm = [idx_rm, length(all_idxMax_sm)];
    end
    all_idxMax_sm(idx_rm) = [];

    nPeaks = length(all_idxMax_sm);
    all_cfun = cell(1, nPeaks);          
    
    for peak_i = 1:nPeaks
        idxMax_sm = all_idxMax_sm(peak_i);
%         maxVal_sm = all_maxVal_sm(peak_i);
        maxVal_sm = spfTC_sm_itp(idxMax_sm);

        log_f_peak_est = log_spfs_cpd_itp(idxMax_sm );    
        r_max_est = maxVal_sm;
        r_bkg_est = min(spfTC_sm_itp);

        % 1b. Estimate width / skewness of tuning curve
        halfMaxTh = (r_max_est-r_bkg_est)/2+r_bkg_est;    
        half_lo_est = log_spfs_cpd_itp(  find( spfTC_sm_itp(1:idxMax_sm) < halfMaxTh, 1, 'last' ) );    
        if isempty(half_lo_est)
            half_lo_est = log_spfs_cpd_itp(1);
        end
        half_hi_est = log_spfs_cpd_itp(  find( spfTC_sm_itp(idxMax_sm:end) < halfMaxTh, 1, 'first' )+idxMax_sm-1 );
        if isempty(half_hi_est)
            half_hi_est = log_spfs_cpd_itp(end);
        end

        w_est = (half_hi_est - half_lo_est)/(2);  % (empirical fit)
        s_est = (half_hi_est + half_lo_est - log_f_peak_est*2)/(1.25*w_est); % (empirical fit)
       

        max_Rmax = (maxSpfTC - r_bkg_est)*opt.maxPeakFitSpfTC + r_bkg_est;

%         max_Rmax_std = maxSpfTC + spfTC_std(indMaxSpfTC);
%         max_Rmax = min(max_Rmax, max_Rmax_std);

    %     logSLNfunc = @(b, f) skewLogNormal_exp(f, b(1), b(2), b(3), b(4), b(5));    
    %     SLNfunc = @(b, f) skewLogNormal(f, b(1), b(2), b(3), b(4), b(5));
%         logSLN_func_p = @(p, x) skewLogNormal_exp(x, p(1), p(2), p(3), p(4), p(5) );

        % 3. Curve-fitting (to get more accurate parameters)
        startPoint = [log_f_peak_est, r_bkg_est, r_max_est, s_est, w_est];    

        logf_peak_lims = iff(constrainFpeak, log_f_peak_est + [-dw, dw], [-50, 50]);
        r_max_lims = iff(constrainRmax,  [0, max_Rmax],   [0, maxSpfTC*8]);
        r_bkg_lims = [min(0, min(spfTC)), mean(spfTC)];
        w_lims = iff(constrainW,  [opt.minSpfWparam, inf],   [0 inf]);
        s_lims = iff(constrainS,  [opt.maxSpfSparam*[-1, 1]], [-inf, inf]);            

        lower_bounds = [logf_peak_lims(1), r_bkg_lims(1), r_max_lims(1), s_lims(1), w_lims(1)];
        upper_bounds = [logf_peak_lims(2), r_bkg_lims(2), r_max_lims(2), s_lims(2), w_lims(2)];    

        constraints_struct = struct('logf_peak', logf_peak_lims, 'r_max', r_max_lims, 'r_bkg', r_bkg_lims, 'w', w_lims, 's', s_lims);
        
        useWgts = opt.useWeightsIfAllPositive && all(spfTC_std > 0);
        if useWgts
            wgts = 1./spfTC_std;
        else
            wgts = [];
        end
        fit_opts = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', lower_bounds, 'Upper', upper_bounds, ...
            'StartPoint', startPoint, 'Weights', wgts);

        % I find it better to do the fitting in the logarithmic space, then we don't get any issues
        % with 0 or negative values of the f_peak.
        [cfun_i, all_fit_stats(peak_i), all_fit_results(peak_i)] = fit(log_spfs_cpd, spfTC, logSLN_func, fit_opts); %#ok<AGROW>
        
%         beta_exp = lsqcurvefit(logSLN_func_p, startPoint, spfs_cpd, spfTC, lower_bounds, upper_bounds, optimset('display', 'off'));
                
        all_cfun{peak_i} = cfun_i;  
        3;

    end
    all_sse = [all_fit_stats.sse];
    all_rsqr = [all_fit_stats.rsquare];
    bestFitIdx = indmin(all_sse);
%     assert( ~any(  all_sse(idx_rm) < all_sse(bestFitIdx)*.99 ));

    cfun = all_cfun{bestFitIdx};
    fit_stats = all_fit_stats(bestFitIdx);
    fit_results = all_fit_results(bestFitIdx);
    fit_stats.cfun = cfun;
        
    if doUnconstrainedFit
        fit_opts = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', coeffvalues(cfun), 'Weights', wgts);        
        [cfun_uc, fit_stats_uc, fit_results_uc] = fit(log_spfs_cpd, spfTC, logSLN_func, fit_opts); 
    end
        
    
    coeffNames = coeffnames(cfun);
    coeffVals = coeffvalues(cfun);
    coeff_ci = confint(cfun);

    log_f_peak = cfun.f_opt_log;
    f_peak = exp(log_f_peak);
    r_peak = cfun.r_max;
    r_bckg = cfun.r_bkg;  
    w_param = cfun.w;
    s_param = cfun.s;
    if r_peak < r_bckg
        r_peak = r_bckg+abs(r_peak-r_bckg); % this is necessary because we used abs(r_peak-r_bckg) in the SLN function
    end
    coeffVals(strcmp(coeffNames, 'r_max')) = r_peak;
            
    spfParams    = cell2struct( num2cell(coeffVals(:)), coeffNames(:), 1);
    spfParams_ci = cell2struct( num2cell(coeff_ci', 2), coeffNames(:), 1);
    
    if doUnconstrainedFit
        coeffVals_uc = coeffvalues(cfun_uc);
        i_rmax = find(strcmp(coeffNames, 'r_max'));   r_max_uc = coeffVals_uc(i_rmax);
        i_rbckg = find(strcmp(coeffNames, 'r_bkg')); r_bckg_uc = coeffVals_uc(i_rbckg); 
        if r_max_uc < r_bckg
            r_max_uc = r_bckg_uc+abs(r_max_uc-r_bckg_uc); % this is necessary because we used abs(r_peak-r_bckg) in the SLN function
            coeffVals_uc(i_rmax) = r_max_uc;
        end            
        
        for f_i = 1:length(coeffNames);
            spfParams.([coeffNames{f_i} '_uc']) = coeffVals_uc(f_i);
        end        
        fit_stats.cfun_uc = cfun_uc;
    end
    
    coeff_exp_C = num2cell(coeffVals); coeff_exp_C{1} = exp( coeff_exp_C{1} );
    cfun_exp = cfit(SLN_func, coeff_exp_C{:});    
        
    
    % plotting was here before
    3;
          
    
    r_halfMax = (r_peak-r_bckg)/2 + r_bckg;
    r_0   = feval(cfun, -100);    
    r_inf = feval(cfun,  100);    
    
    frac_of_peak_func = @(r) (r-r_bckg)/(r_peak-r_bckg);
    frac_of_peak_data = @(r) (r-minSpfTC)/(maxSpfTC-minSpfTC);
        
    stcProps.FracAtLowest_fit  = frac_of_peak_func(feval(cfun, log_spfs_cpd(1)  ) );
    stcProps.FracAtHighest_fit = frac_of_peak_func(feval(cfun, log_spfs_cpd(end)) );
    stcProps.FracAtLowest_data  = frac_of_peak_data(spfTC(1));
    stcProps.FracAtHighest_data = frac_of_peak_data(spfTC(end));
    
    stcProps.FracAtZero = frac_of_peak_func(r_0);
    stcProps.FracAtInf  = frac_of_peak_func(r_inf);
    stcProps.PrefSpfLowerThanStim = f_peak < spfs_cpd(1);
    stcProps.PrefSpfHigherThanStim = f_peak > spfs_cpd(end);
    stcProps.PrefDistFromStim = max( rectified( log_spfs_cpd(1) - log_f_peak ), ...
                                     rectified( log_f_peak - log_spfs_cpd(end)) );
    stcProps.nFits = nPeaks;

    [log_f_lo, log_f_hi] = getHalfHiLo(@(x) feval(cfun, x)-r_halfMax,  log_f_peak,  log_spfs_cpd([1, end]));    
    f_lo = exp(log_f_lo);
    f_hi = exp(log_f_hi);
    fit_stats.f_lo = f_lo;
    fit_stats.f_hi = f_hi;
    
    stcProps.HalfMaxDistFromStim = max( rectified( log_spfs_cpd(1) - log_f_lo ), ...
                                        rectified( log_f_hi - log_spfs_cpd(end)) );
    stcProps.constraints = constraints_struct;
    
    % calculate octave measure of spatial frequency tuning width 
    w_spf = log2(f_hi / f_lo);    
    
    % calculate "normalized bandwidth" (Tolhurst and Thompson 1981)
    normBw = (f_hi-f_lo)/f_peak;
    

    if show
        %%
        showAllFits = 0;
        log10_factor = log(10);        
        log_spfs_cpd_tmp = linspace(log_spfs_cpd(1), log_spfs_cpd(end), 2000);
        log10_spfs_cpd = log10(spfs_cpd);
        log10_spfs_cpd_tmp = linspace(log10_spfs_cpd(1), log10_spfs_cpd(end), 2000);

        xlims_log = lims(log10_spfs_cpd, .02);
        
        figure(5); clf; hold on; box on;        
%         plot(log_spfs_cpd, spfTC, 'bs'); hold on
        errorbar(log10_spfs_cpd, spfTC, spfTC_std, 'bo');
%         plot(log_spfs_cpd, spfTC_sm, 'g.'); 
        xlabel('Log_{10} Spatial Frequency (cyc/deg)');
        ylabel('Firing Rate (Hz)');        
        xlim(xlims_log);
        
        plot(log10_spfs_cpd_tmp, feval(cfun, log_spfs_cpd_tmp), 'r');     
        if ~isempty(GC)
            id_str = sprintf('Cell [%d, %d], %s', GC(1), GC(2));
        else
            id_str = '';
        end
        stat_str = sprintf('f_{opt} = %.2f cyc/deg.  B_l = %.2f.   R^2 = %.2f',  f_peak, w_spf, fit_stats.rsquare);
        title({stat_str});
        if exist('log_f_lo', 'var')
            plot([log_f_lo, log_f_hi]/log10_factor, r_halfMax+[0, 0], 'k.:');
        end
        if showAllFits
            for i = 1:nPeaks
               plot(log10_spfs_cpd_tmp, feval(all_cfun{i}, log_spfs_cpd_tmp), [color_s(i) ]);
                
            end
            
        end
        %%
        legend({'Tuning Curve', 'SLN Fit', 'Half-max'}, 'location', 'best')   
        
        if doUnconstrainedFit
            %%
            coeffVals_uc = coeffvalues(cfun_uc);
            coeff_uc_exp_C = num2cell(coeffVals_uc); coeff_uc_exp_C{1} = exp( coeff_uc_exp_C{1} );
            cfun_uc_exp = cfit(SLN_func, coeff_uc_exp_C{:});    
            
            h = [];
            xlims = lims(spfs_cpd, .05, [], 1);
            figure(16); clf; hold on; box on;
%             spfs_cpd_tmp = 10.^(log_spfs_cpd_tmp); % linspace(spfs_cpd(1), spfs_cpd(end), 100);        
%             h(1) = plot(spfs_cpd, spfTC, 'bs'); 
            h(1) = errorbar(log10_spfs_cpd, spfTC, spfTC_std, 'bo');
            h(2) = plot(log10_spfs_cpd_tmp, feval(cfun, log_spfs_cpd_tmp), 'r-');  
            h(3) = plot(log10_spfs_cpd_tmp, feval(cfun_uc, log_spfs_cpd_tmp), ['-'], 'color', [0 .6 0] );              
            
            xlabel('Spatial Frequency (cyc/deg)');
            ylabel('Firing Rate (Hz)');
            set(gca, 'xlim', xlims_log);
            y1 = min(0, min(spfTC-spfTC_std*1));
            ylim([y1, r_peak*2]);
            legend(h, {'Tuning Curve', 'Constrained Fit', 'Unconstrained Fit'}, 'location', 'best', 'fontsize', 9)                        
            title({stat_str});
            3;
            
        end
        
        
        
        
        %%
        figure(6);
        spfs_cpd_tmp = 10.^(log_spfs_cpd_tmp); % linspace(spfs_cpd(1), spfs_cpd(end), 100);        
        plot(spfs_cpd_tmp, feval(cfun_exp, spfs_cpd_tmp), 'r');  
        set(he, 'linestyle', 'none');
        title({id_str, stat_str});
        set(gca, 'xscale', 'log', 'xlim', lims(spfs_cpd, .05, [], 1) )
        3;       
        
    end         
    
    
    
    if show 
        figure(5);
        drawVerticalLine([log_f_lo, log_f_hi], 'linestyle', ':');
        drawHorizontalLine(r_halfMax, 'linestyle', ':');
        
        figure(6);
        drawVerticalLine([f_lo, f_hi], 'linestyle', ':');
        drawHorizontalLine(r_halfMax, 'linestyle', ':');        
        3;
    end
    3;
    result_txt = fit_results.message;        

end




function [lo, hi] = getHalfHiLo(f, x_peak, x_range)

    % verify that is peak
    %{
    [x_peak_search] = fminsearch(@(x) -f(x), x_peak);
    assert( abs(x_peak - x_peak_search) < 1e-3 );
    figure(55); fplot(f, x_range);
%}
    
    %extend x_range if necessary 
    x_range(1) = min(x_range(1), x_peak-1); % make sure range includes the peak
    x_range(2) = max(x_range(2), x_peak+1);
    
    while f(x_range(1)) > 0
        x_range(1) = x_range(1)-1;
    end
    while f(x_range(2)) > 0
        x_range(2) = x_range(2)+1;
    end    

%     [X, spfVals] = fplot(f, x_range);    
    X = linspace(x_range(1), x_range(2), 700);
    spfVals = f(X);
    
    crossings = find( diff(sign(spfVals)) ~= 0);
    
    dbug = 0;
    if dbug
       plot(X, spfVals, 'b.-')        
    end

    assert(length(crossings) == 2)
    
    lo_estimate = X( crossings(1)+[0,1] );
    hi_estimate = X( crossings(2)+[0,1] );    
    
    lo = fzero(@(x) f(x), lo_estimate);
    hi = fzero(@(x) f(x), hi_estimate);
    
    assert(lo < x_peak);
    assert(hi > x_peak);
    
end
