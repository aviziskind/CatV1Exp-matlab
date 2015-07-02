function allStats = calcDegreeOfTuningStats(R_full, bckgSamples, Gid, cellId)

    global doOddEvenTrials
    alpha_drifting = .01;
    alpha_flashed = .01;    
    opt.responseSizeTest = 'utest';
    opt.responseSizeNTop_flashed = 1;
    opt.responseSizeNTop_drifting = 1;

    opt.oriMaxSigma_deg = 60;
    opt.oriSmoothWidth_deg = 10;
    opt.oriWGlobal_includeOppDirection = 0;
    opt.oriFitNonPrefTuningCurve = 0;
    opt.oriFitUnconstrainedMeanTuningCurve = 1;

    opt.minRsqr = 0.4;
    opt.excludeOffPeakMaximumCases = 0;
    opt.localPeakMaxDiscrep_degrees = 15;
    opt.minPeakRatio = 0.75;
    
    opt.spfSmoothW = 0.6;
    opt.maxPeakFitSpfTC = 1.5;
    opt.minSpfWparam = 0.09;
    opt.maxSpfSparam = 2;
    opt.useWeightsIfAllPositive = 0;
    opt.includeTuningCurves = 1;
    
    opt.calculateAllTrialStats = 1;
    opt.calculateOddEvenTrialStats = 0;
        opt.oddEvenTrial_spontSubtracted = [true];  % = [false, true];
    opt.doJackKnifeStdErr = 1;
    opt.oriRep_smoothFilter = 'gaussian';
    opt.oriRep_gaussianSmoothWidth_degree = 5;
    opt.driftingOriSelectivityTest_combineOppositeDirections = 1;   
    
    opt.calcOriW_local = 1;  % temporary skip of ori local (used this once for fast calculation of global width for clusters)

    
    opt.applyStdErrorThresholds = true;
    opt.maxOriStdErr_deg = 5;
    opt.maxDSIStdErr     = 0.1;
    opt.maxSpfStdErr_oct = 0.5;
    opt.maxF1oDCStdErr   = 10;
    
    opt.doUnrectifiedVersionsForSpontSubtract = 0;
    doSpontIncluded = 1;
    doSpontSubtracted = 1;
        doSpontSubtracted_spf = 1;
    useResponseBckgSig = 1;    
    
    opt.doSLNcurveWithNoConstraints = 0;
    opt.useAlternateDefOfDSIforUnrectified = 1;
    
%     global show
    if ~isempty(doOddEvenTrials)
        opt.calculateOddEvenTrialStats = opt.calculateOddEvenTrialStats || doOddEvenTrials;
    end

    if isempty(bckgSamples)        

        
    end

    allStats = [];
    if iscell(R_full) % for multiple temporal frequency experiments
        allStats = cellfun(@(R1, R_full1, f1, uori1) ...
                calcDriftingGratingStats(R1, R_full1, f1, uori1, Gid, bckgSamples), R, R_full, f, 'un', 0);
        return;
    end
    
    stimType = getGratingStimType(Gid);
    sd = siteDataFor(Gid);
    
    if isstruct(R_full)
        R_full = decompress(R_full);
    end
    R_full = double(R_full);
        
    
    if isstruct(bckgSamples)
        bckgSamples = decompress(bckgSamples);      
    end
    bckgSamples = double(bckgSamples);

%     bckgRate = [mean(bckgSamples), std(bckgSamples)];
%     R_full_si = R_full;
%     R_full_ss = rectified(R_full - bckgRate(1));            
    
    switch stimType.gratingType
        case 'drifting',
            isOrientBatch = strcmp(stimType.driftingType, 'Orientation Batch');
            isSpatBatch   = strcmp(stimType.driftingType, 'Spatial Frequency Batch');  
            opt.alpha = alpha_drifting;            
        case 'flashed'
            isOrientBatch = 1;
            isSpatBatch = 1;
            opt.alpha = alpha_flashed;                    
    end

    
    
    
        
%     R_os_smoothed = mean(mean(R_full, 4),3);    
%     if isOrientBatch
%         dw_ori = 1;         
%         R_os_smoothed = gaussSmooth(R_os_smoothed, dw_ori, 1, 1);
%     end
%     if isSpatBatch        
%         log_spf = log(sd.spPeriod_pix);
%         dw_spf = mean(diff(log_spf));
%         R_os_smoothed = gaussSmooth_nu(log_spf, R_os_smoothed, dw_spf, 2);    
%     end        
% % %     R_os_smoothed = gaussSmooth_nu(log_spf, gaussSmooth(R_os, dw_ori, 1, circ_flag), dw_spf, 2);
% 
%     ori_pref_idx = indmax( sum(R_os_smoothed, 2) );
%     spf_pref_idx = indmax( R_os_smoothed(ori_pref_idx(1), :) );

    [idxs] = getIdxPrefStim(Gid, R_full);
    idx_sm = idxs.prefR_avP_sm;
    [ori_pref_idx, spf_pref_idx] = deal(idx_sm(1), idx_sm(2));

    % calculate F1/DC
    %%
    nPh = size(R_full, 3); phs = linspace(0, 360, nPh+1); phs = phs(1:nPh);
    nTrials = size(R_full, 4);
    phaseTC_atPref_trials = squeeze( R_full(ori_pref_idx, spf_pref_idx, :, :) );
    assert(isequal(size(phaseTC_atPref_trials), [nPh, nTrials]));
    phaseTC_atPref = mean( phaseTC_atPref_trials, 2);
        
    F1oDC = getF1oDC(phs, phaseTC_atPref(:), 360);
    opt.F1oDC_calculated = F1oDC;
    
%     phaseTC_atPref_jacks = arrayfun(@(i) mean(  phaseTC_atPref_trials(:, setdiff(1:nTrials, i) )  ,2), 1:nTrials, 'un', 0);
    phaseTC_atPref_jacks = jackknifeAverageTrials( phaseTC_atPref_trials, 2 );
    %%
    if opt.doJackKnifeStdErr
        F1oDC_jacks = deal(  zeros(1,nTrials) );
        for jack_i = 1:nTrials
            F1oDC_jacks(jack_i) = getF1oDC(phs, phaseTC_atPref_jacks{jack_i}(:), 360);
        end
        F1oDC_stderr_jack = jackknifeStdErr(F1oDC_jacks, F1oDC);
        opt.F1oDC_stderr_calculated = F1oDC_stderr_jack;
    end
    
    
    
    
    
    % Perform cell selection criteria tests that are done for both ori/spf batches:    
    if useResponseBckgSig
        responseSize_S = responseMagRelToBckg(R_full, bckgSamples, opt, stimType.gratingType, Gid, cellId);
    else
        responseSize_S = [];
    end

    
    % TESTS FOR ORIENTATION BATCHES;
    if opt.calculateAllTrialStats
        Ntrials = size(R_full, 4); 
        idxOdd  = 1:2:Ntrials; R_full_odd  = R_full(:,:,:,idxOdd);
        idxEven = 2:2:Ntrials; R_full_even = R_full(:,:,:,idxEven);                        
    end    
    
    if isOrientBatch        
        
        dirs_deg = sd.ori_deg(:);        
        
        doOrientationBatchTest_func = @(R_full_arg, subtractBkg_flag, oeStats) ...
            doOrientationBatchTests(R_full_arg, bckgSamples, dirs_deg, stimType.gratingType, opt, subtractBkg_flag, responseSize_S, ori_pref_idx, spf_pref_idx, oeStats);

        [oriStats_odd_even_si, oriStats_odd_even_ss] = deal([]);
        if opt.calculateOddEvenTrialStats
            if any(opt.oddEvenTrial_spontSubtracted == false)
                allStats.oriStats_si_even     = doOrientationBatchTest_func(R_full_even, 0, []);
                allStats.oriStats_si_odd      = doOrientationBatchTest_func(R_full_odd,  0, []);
                oriStats_odd_even_si = [allStats.oriStats_si_odd, allStats.oriStats_si_even];
            end
            if any(opt.oddEvenTrial_spontSubtracted == true)
                allStats.oriStats_ss_even     = doOrientationBatchTest_func(R_full_even, 1, []);
                allStats.oriStats_ss_odd      = doOrientationBatchTest_func(R_full_odd,  1, []);
                oriStats_odd_even_ss = [allStats.oriStats_ss_odd, allStats.oriStats_ss_even];
            end            
            
        end
        
        if opt.calculateAllTrialStats
            if doSpontIncluded
                allStats.oriStats_si = doOrientationBatchTest_func(R_full, 0, oriStats_odd_even_si);
            end
            if doSpontSubtracted
                allStats.oriStats_ss = doOrientationBatchTest_func(R_full, 1, oriStats_odd_even_ss);            
            end
        end


    end    
    
    % TESTS FOR SPATIAL FREQUENCY BATCHES;
    if isSpatBatch        
        spfs_cpd = 1./(sd.spPeriod_pix * sd.stimulusInfo.degreesPerBlock);

        if isOrientBatch
            ori_pref_idx_si = indmin(abs(allStats.oriStats_si.dir_pref_deg - dirs_deg));
            if doSpontSubtracted && doSpontSubtracted_spf
                ori_pref_idx_ss = indmin(abs(allStats.oriStats_ss.dir_pref_deg - dirs_deg));
            end
        else
            [ori_pref_idx_si, ori_pref_idx_ss] = deal([]);
        end
        
        
        [spfStats_odd_even_si, spfStats_odd_even_ss] = deal([]);
        if opt.calculateOddEvenTrialStats
            if any(opt.oddEvenTrial_spontSubtracted == false)
                allStats.spfStats_si_even = doSpatialFrequencyBatchTests(R_full_even, bckgSamples, spfs_cpd, stimType.gratingType, opt, 0, responseSize_S, ori_pref_idx_si, []);
                allStats.spfStats_si_odd  = doSpatialFrequencyBatchTests(R_full_odd,  bckgSamples, spfs_cpd, stimType.gratingType, opt, 0, responseSize_S, ori_pref_idx_si, []);
                spfStats_odd_even_si = [allStats.spfStats_si_odd, allStats.spfStats_si_even];
            end
            if any(opt.oddEvenTrial_spontSubtracted == true) && doSpontSubtracted_spf
                allStats.spfStats_ss_even = doSpatialFrequencyBatchTests(R_full_even, bckgSamples, spfs_cpd, stimType.gratingType, opt, 1, responseSize_S, ori_pref_idx_ss, []);
                allStats.spfStats_ss_odd  = doSpatialFrequencyBatchTests(R_full_odd,  bckgSamples, spfs_cpd, stimType.gratingType, opt, 1, responseSize_S, ori_pref_idx_ss, []);                                
                spfStats_odd_even_ss = [allStats.spfStats_ss_odd, allStats.spfStats_ss_even];
            end
            
        end        
        
        if opt.calculateAllTrialStats
            if doSpontIncluded
                allStats.spfStats_si = doSpatialFrequencyBatchTests(R_full, bckgSamples, spfs_cpd, stimType.gratingType, opt, 0, responseSize_S, ori_pref_idx_si, spfStats_odd_even_si);
            end
            if doSpontSubtracted && doSpontSubtracted_spf
                allStats.spfStats_ss = doSpatialFrequencyBatchTests(R_full, bckgSamples, spfs_cpd, stimType.gratingType, opt, 1, responseSize_S, ori_pref_idx_ss, spfStats_odd_even_ss);
            end            
        end
        
    end
%     profile viewer;
    3;
        
    
end
    
%%%%%% FUNCTIONS FOR ORIENTATION %%%%%%%-----------------------------------------------------------

function oriStats = doOrientationBatchTests(R_full, bckgSamples, dirs_deg, gratingType, opt, subtractBckgFlag, responseSize_S, ori_pref_idx, spf_pref_idx, oeStats)
    % note replace R_full with R_full_si - subtract & rectify after averaging

    [nDir, nSpf, nPh, nTrials] = size(R_full); 
    assert(nDir > 2);
    alpha = opt.alpha;
    minRsqr = opt.minRsqr;
    
    oriMax = switchh(gratingType, {'flashed', 'drifting'}, [180, 360]);            
    oriSampleBiasVector = sum([cos(deg2rad(2*dirs_deg(:))), sin(deg2rad(2*dirs_deg(:)))], 1);
    oriSampleBias = normV(oriSampleBiasVector,2);
    assert(oriSampleBias < .01);
    
    R_os = mean(mean(R_full,4),3);
    maxFiringRate = max(R_os(:));
    if maxFiringRate > 500;
        3;
    end
    F1oDC = opt.F1oDC_calculated;
    F1oDC_stderr_jack = opt.F1oDC_stderr_calculated;
    
    haveOEstats = exist('oeStats','var') && ~isempty(oeStats);
    
%     R_ori_sp = mean(mean(R_full, 3),4);
    bckgMeanRate = mean(bckgSamples);
    
 
    doUnrectifiedVersions = opt.doUnrectifiedVersionsForSpontSubtract && subtractBckgFlag;

%     R_full_avPh_avTr = mean(mean(R_full, 3),4);
    
    R_full_prefSpf = R_full(:,spf_pref_idx,:,:);    
%     R_full_meanSpf = mean(R_full,2);
                
    
    if subtractBckgFlag
        %1. r_k (Ori tuning curve (rk))  - for ori width, DSI,
        %2. R_full / (really: just r_k) for ori tuning test  -- should be just at pref spf for SS
        %3. R_full at pref_spf (with odd/even) for reproducibility

        % find the orientations that, when averaged over phases & trials, are below background (for pref spf) 
%         R_full_prefSpf_avPh_avTr = mean(mean(R_full_prefSpf,3),4);
        dims_average = [2, 3, 4]; % average over phases, and trials;
                
        R_full_prefSpf_unrectified = R_full_prefSpf;
        R_full_prefSpf = rectifyAverage(R_full_prefSpf - bckgMeanRate, dims_average, 1);
         
        R_full_allSpf = rectifyAverage(R_full - bckgMeanRate, dims_average, 1);
        
%         R_full_oriTuningTest = R_full_prefSpf_ss_rect;  % drifting: only 1 spf anyway. flashed: in case averaging takes below zero.
        R_full_oriTuningTest = R_full_allSpf;  % drifting: only 1 spf anyway. flashed: in case averaging takes below zero.
        
    else
        R_full_oriTuningTest = R_full;        

        R_full_prefSpf_unrectified = R_full_prefSpf;
    end
    
    
    R_ori_tr = reshape(mean( R_full_prefSpf, 3), [nDir, nTrials]);
    
    
    r_k_dir = mean(R_ori_tr, 2);    
    r_k_dir_std = std(R_ori_tr, [], 2);
    
    if opt.doJackKnifeStdErr
        r_k_dir_jacks = jackknifeAverageTrials( R_ori_tr, 2 );
%         r_k_dir_jacks2 = arrayfun(@(i) mean(  R_ori_tr(:, setdiff(1:nTrials, i) )  ,2), 1:nTrials, 'un', 0);
%         assert(isequal(r_k_dir_jacks, r_k_dir_jacks2));        
    end
    
    if doUnrectifiedVersions
        R_ori_tr_unrec = reshape(mean( R_full_prefSpf_unrectified, 3), [nDir, nTrials]);
        r_k_dir_unrec = mean(R_ori_tr_unrec, 2);
    end
        
    
    if subtractBckgFlag        
%         r_k_dir_ss = rectified(r_k_dir - bckgMeanRate);        
        
%         assert(~any(r_k_dir(idx_ori_zero_response_prefSpf)));
%         assert(~any(r_k_dir_std(idx_ori_zero_response_prefSpf)));
        3;
    end
    3;
%     wgts = 1
%     if nSpf == 1        
%         r_k_dir = R_ori_sp;  % r_k = ori/dir tuning curve (average response elicited by gratings of the kth orientation/direction.)
%     else        
%         r_k_dir = R_ori_sp(:, spf_pref_idx); % % use orientation tuning curve at preferred spatial frequency.                
%     end
    
%     r_k_dir_orig = r_k_dir;
    
    
    % 1. ORIENTATION SELECTIVITY test
    % test that, for each set of stimuli that includes one orientation each,
    % the components of the response along the preferred orientation form a distribution
    % with a mean significantly different from 0.
    
    combine_flag = opt.driftingOriSelectivityTest_combineOppositeDirections;
    [ori_pref_deg, ori_sel_pval_stats] = getPreferredOriFromOSP(R_full_oriTuningTest, dirs_deg, spf_pref_idx, gratingType, combine_flag);
    if doUnrectifiedVersions
        [ori_pref_deg_unrec, ori_sel_pval_stats_unrec] = getPreferredOriFromOSP(R_full_prefSpf_unrectified, dirs_deg, spf_pref_idx, gratingType, combine_flag);
    end
    [ori_pref_deg_spfPref] = getPreferredOriFromOSP(R_full_prefSpf, dirs_deg, spf_pref_idx, gratingType, combine_flag);
%     [ori_pref_deg_prefSpf, ori_sel_pval_stats_prefSpf] = getPreferredOriFromOSP(R_full_prefSpf, dirs_deg, spf_pref_idx, gratingType, combine_flag);
            
    ori_sel_pval = ori_sel_pval_stats.resultant_prob;
    if strcmp(gratingType, 'flashed') && ~isnan(ori_pref_deg)
        assert(ori_pref_deg <= 180);   % ori_pref_deg = mod(ori_pref_deg, 180);
    end
    
   

    

    % 2. ORIENTATION REPRODUCITIBILITY test
    % test that odd & even trials (at preferred spatial frequency) are correlated, (after smoothing).
    %%
    switch opt.oriRep_smoothFilter
        case 'gaussian', 
            gaussianSmoothWidth_oriSteps = opt.oriRep_gaussianSmoothWidth_degree / diff(dirs_deg(1:2));
            smoothingFilter = gaussian(-4:4, 0, gaussianSmoothWidth_oriSteps);
            smoothingFilter = smoothingFilter/sum(smoothingFilter);
        case 'triangular', 
            if nDir == 72
                smoothingFilter = [1, 2, 3, 2, 1]/9; %gaussian(-1:1, 0, 1);
            elseif any(nDir == [30, 36])
                smoothingFilter = [1, 2, 1]/4; %gaussian(-1:1, 0, 1);
            end
    end
    
    ori_rep_stats.pval_gs_spfPref = testOrientationRep(R_full_prefSpf, [], smoothingFilter);
%     ori_rep_stats.pval_gs_spfAv   = testOrientationRep(R_full_meanSpf, [],           gs_smoothingFilter);
%     ori_rep_stats.pval_tr_spfPref = testOrientationRep(R_full_prefSpf, [], tr_smoothingFilter);
%     ori_rep_stats.pval_tr_spfAv   = testOrientationRep(R_full_meanSpf, [],           tr_smoothingFilter);
        
    ori_rep_pval = ori_rep_stats.pval_gs_spfPref;    
    %%
    
    % 3. RESPONSE SIZE SIGNIFICANCE (compared to background); this is done above        
    if subtractBckgFlag 
        response_size_pval = responseSize_S.pval;
    else
        response_size_pval = 0;  % ie. automatically pass the response size test 
    end

    
    
    % DIRECTION SELECTIVITY test (for drifting gratings)
    alternateDSI = opt.useAlternateDefOfDSIforUnrectified;
    switch gratingType
        case 'drifting',
            [dir_pref_deg, DSI_global] = getPreferredDirection(r_k_dir, ori_pref_deg, dirs_deg);            
            if doUnrectifiedVersions
                [dir_pref_deg_unrec, DSI_global_unrec] = getPreferredDirection(r_k_dir_unrec, ori_pref_deg_unrec, dirs_deg, alternateDSI);
            end
                        
        case 'flashed',  
            [dir_pref_deg, DSI_global] = deal(ori_pref_deg, nan);
            if doUnrectifiedVersions
                 [dir_pref_deg_unrec, DSI_global_unrec] = deal(ori_pref_deg_unrec, nan);
            end
    end
    
    cellOK_soFar = (ori_sel_pval(1) < alpha) && (ori_rep_pval < alpha) && (response_size_pval < alpha);  % instead of (nStdAboveBckg > 3)
    
    % Orientation tuning widths (global & local):   
    if opt.oriWGlobal_includeOppDirection
        ori_diff_max = 180;
    else 
        ori_diff_max = 90;
    end
    
    idx_ori_use = getOriIndexInRange(dirs_deg, dir_pref_deg, oriMax, ori_diff_max);
    
    if isnan(dir_pref_deg)
        [w_ori_global, ori_OSI, ori_circVar] = deal(nan);
    else        
        w_ori_global = calcOriGlobalWidth(r_k_dir, dirs_deg, dir_pref_deg, gratingType, ori_diff_max, idx_ori_use); % no fitting required for this - so can always get answer
        ori_OSI = getOriSelectivityIndex(r_k_dir, dirs_deg, idx_ori_use); % no fitting required for this - so can always get answer
        ori_circVar = getCircularVariance(r_k_dir, dirs_deg, idx_ori_use);    
        
     
    end
    
    if opt.calcOriW_local
        cellOK_soFar_arg = cellOK_soFar;    
    else  % skip ori tuning curve fitting
        cellOK_soFar_arg = 0;
    end
    
    [w_ori_local, R_null, R_pref, DSI_local, result_txt, oriParams, oriParams_ci, otc_stats, r_k_dir_cent, ori_peak_unconstrained_mean] = ...
        getOriLocalWidth(r_k_dir, r_k_dir_std, dir_pref_deg, dirs_deg, oriMax, cellOK_soFar_arg, opt, 0); 


     

    if doUnrectifiedVersions
%     [w_ori_local_unrec, R_null_unrec, R_pref_unrec, DSI_local_unrec, result_txt_unrec, oriParams_unrec, oriParams_ci_unrec, otc_stats_unrec, r_k_dir_cent_unrec] = ...
%         getOriLocalWidth(r_k_dir_unrec, r_k_dir_std, dir_pref_deg_unrec, dirs_deg, oriMax, cellOK_soFar_arg, opt, 0); 
    end
%     if DSI_global_unrec < 0 || DSI_local_unrec < 0; % ~isnan(DSI_local) && isnan(DSI_local_unrec)
%         3;
%     end
        
    
     stderr_jack = struct; 
     if opt.doJackKnifeStdErr
        %%
        [ori_pref_deg_jacks, dir_pref_deg_jacks, DSI_global_jacks, w_ori_global_jacks, w_ori_local_jacks] = deal(  zeros(1,nTrials) );
         
        for jack_i = 1:nTrials
            R_full_ori_use = R_full_oriTuningTest(:,:,:, setdiff(1:nTrials, jack_i));
            [ori_pref_deg_jacks(jack_i)] = getPreferredOriFromOSP(R_full_ori_use, dirs_deg, spf_pref_idx, gratingType, combine_flag);            
            if strcmp(gratingType, 'drifting')
                [dir_pref_deg_jacks(jack_i), DSI_global_jacks(jack_i)] = getPreferredDirection(r_k_dir_jacks{jack_i}, ori_pref_deg_jacks(jack_i), dirs_deg);
            else
                dir_pref_deg_jacks(jack_i) = ori_pref_deg_jacks(jack_i);
            end
            if ~isnan(dir_pref_deg)
                w_ori_global_jacks(jack_i) = calcOriGlobalWidth(r_k_dir_jacks{jack_i}, dirs_deg, dir_pref_deg_jacks(jack_i), gratingType, ori_diff_max, idx_ori_use);
                w_ori_local_jacks(jack_i) = getOriLocalWidth(r_k_dir_jacks{jack_i}, r_k_dir_std, dir_pref_deg_jacks(jack_i), dirs_deg, oriMax, cellOK_soFar_arg, opt, 0);
            end
        end
        
        ori_pref_stderr_jack = jackknifeStdErr(ori_pref_deg_jacks, ori_pref_deg, 'circ', 180);
        DSI_global_stderr_jack = jackknifeStdErr(DSI_global_jacks, DSI_global);
        w_ori_global_stderr_jack = jackknifeStdErr(w_ori_global_jacks, w_ori_global);
        w_ori_local_stderr_jack = jackknifeStdErr(w_ori_local_jacks, w_ori_local);
        stderr_jack = struct('ori_pref', ori_pref_stderr_jack, 'DSI_global', DSI_global_stderr_jack, ...
                             'w_ori_global', w_ori_global_stderr_jack, 'w_ori_local', w_ori_local_stderr_jack, 'F1oDC', F1oDC_stderr_jack);
     end
    
     
     

     
     
    otc_fit_rsqr = otc_stats.gof.rsquare;    

    if ~isnan(R_null)
        assert(R_null < 200);
    end
    
    if opt.calcOriW_local
        cellOK = cellOK_soFar && strncmp(result_txt, 'Success', 7) && (otc_fit_rsqr > minRsqr);               
    else
        cellOK = cellOK_soFar;
    end        
    
    
    [w_ori_global_err, w_ori_global_err_rel, w_ori_local_err, w_ori_local_err_rel, DSI_global_err, DSI_global_err_rel, ori_pref_deg_err] = deal(nan);
    if haveOEstats
        %%
        [w_ori_global_err, w_ori_global_err_rel] = diffOverAv(oeStats(1).w_ori_global, oeStats(2).w_ori_global);
        [w_ori_local_err, w_ori_local_err_rel]   = diffOverAv(oeStats(1).w_ori_local, oeStats(2).w_ori_local);
        [DSI_global_err, DSI_global_err_rel]     = diffOverAv(oeStats(1).DSI_global, oeStats(2).DSI_global);
        ori_pref_deg_err = circDist(oeStats(1).ori_pref_deg, oeStats(2).ori_pref_deg, 180);

    end
    error_est_OE = struct('w_ori_global_err', w_ori_global_err, 'w_ori_global_err_rel', w_ori_global_err_rel, ...
        'w_ori_local_err', w_ori_local_err, 'w_ori_local_err_rel', w_ori_local_err_rel, ... 
        'DSI_global_err', DSI_global_err, 'DSI_global_err_rel', DSI_global_err_rel, ...
        'ori_pref_deg_err', ori_pref_deg_err);
        
    
    
%     cellOK = strcmp(result_txt, 'ok');
%     if cellOK_soFar && w_ori_local(1) > 50 
%     if DSI_local > .99 && DSI_global < .2
%         [w_ori_local, R_null, R_pref, DSI_local, result_txt, oriParams, oriParams_ci, otc_stats] = ...
%             getOriLocalWidth(r_k_dir, r_k_dir_std, dir_pref_deg, dirs_deg, oriMax, cellOK_soFar, opt, 1); 
%     end
    if doUnrectifiedVersions
        unrec_fields = {'DSI_global_unrec', DSI_global_unrec, 'DSI_local_unrec', DSI_local_unrec, 'dir_pref_deg_unrec', dir_pref_deg_unrec};         
    else
        unrec_fields = {};
    end

    if w_ori_global > 100
        3;
    end
    
    origValues = struct('ori_pref_deg', ori_pref_deg, 'dir_pref_deg', dir_pref_deg, 'DSI_global', DSI_global, 'w_ori_global', w_ori_global, 'w_ori_local', w_ori_local, 'F1oDC', F1oDC);
    
    if opt.applyStdErrorThresholds
        if ori_pref_stderr_jack > opt.maxOriStdErr_deg
            [ori_pref_deg, dir_pref_deg] = deal(nan);
        end
        if DSI_global_stderr_jack > opt.maxDSIStdErr
            DSI_global = nan;
        end
        if w_ori_global_stderr_jack > opt.maxOriStdErr_deg
            w_ori_global = nan;
        end
        if w_ori_local_stderr_jack > opt.maxOriStdErr_deg
            w_ori_local = nan;
        end        
        if F1oDC_stderr_jack > opt.maxF1oDCStdErr 
            F1oDC = nan;
        end        
    end
    
    oriStats = struct('cellOK', cellOK, 'cellOK_asMU', cellOK_soFar, ...
        'ori_sel_pval',ori_sel_pval, 'ori_rep_pval', ori_rep_pval, 'ori_rep_stats', ori_rep_stats, ...
        'ori_sel_pval_stats', ori_sel_pval_stats, 'rsqr', otc_fit_rsqr, 'response_size_pval', response_size_pval, ...
        'w_ori_global', w_ori_global, 'w_ori_local', w_ori_local, 'w_local_fit_txt', result_txt, ...
        'DSI_global', DSI_global, 'DSI_local', DSI_local, 'OSI', ori_OSI, 'circVar', ori_circVar, ...        
        'ori_pref_deg', ori_pref_deg, 'dir_pref_deg', dir_pref_deg, ... ...
        'ori_pref_deg_spfPref', ori_pref_deg_spfPref, ...
        'ori_pref_deg_unc', ori_peak_unconstrained_mean, ...
        unrec_fields{:}, ...
        'error_OE', error_est_OE, ...
        'error_jack', stderr_jack, ...
        'ori_pref_idx', ori_pref_idx, 'spf_pref_idx', spf_pref_idx, ...
        'oriParams', oriParams, 'oriParams_ci', oriParams_ci, ...
            'response_sig', responseSize_S, 'otc_stats', otc_stats, ...
            'R_spont_abs', bckgMeanRate,        'R90_total_abs', R_null,        'R90_stim_abs', (R_null-bckgMeanRate), ...
            'R_spont_rel', bckgMeanRate/R_pref, 'R90_total_rel', R_null/R_pref, 'R90_stim_rel', (R_null-bckgMeanRate)/R_pref, ...
            'F1oDC', F1oDC, ...
            'r_k_dir', single(r_k_dir), ...
            'r_k_dir_cent', single(r_k_dir_cent), ...
            'orig', origValues ...
        );        
    3;
end


function [d_abs, d_rel] = diffOverAv(a,b)
    d_abs = abs(a-b);
    d_rel = d_abs / abs( (a + b) / 2);
end


function R = rectifyAverage(R, dims_av, dims_idx)
    R_av = R;
    for i = 1:length(dims_av)
        R_av = mean(R_av, dims_av(i));
    end
    idx_neg = R_av < 0;
    
    assert(dims_idx == 1);
    R(idx_neg,:,:,:) = 0; % if dims_idx == 2, R(:,idx_neg,:,:) = 0;
    
end

function ori_rep_pval = testOrientationRep(R_full, spf_pref_idx, smoothingFilterCoefs)            	
    [nOri, nSpf, nPh, nReps] = size(R_full);
    
    % (1) select ori's at preferred spf, averaging over phases, (2) smooth, and (3) compare odd & even trials    
    if ~isempty(spf_pref_idx)
        R_av_ph = reshape(mean(R_full(:, spf_pref_idx, :, :), 3), [nOri, nReps]);
    else
        R_av_ph = reshape(mean(mean(R_full(:, :, :, :), 3),2), [nOri, nReps]);
    end

    odd_idx  = 1:2:nReps;
    even_idx = 2:2:nReps;    
    
    R_odd  = mean(R_av_ph(:,odd_idx ), 2);
    R_even = mean(R_av_ph(:,even_idx), 2);
    
    R_odd_sm  = smoothedR(R_odd,  smoothingFilterCoefs);
    R_even_sm = smoothedR(R_even, smoothingFilterCoefs);

    [~, ori_rep_pval] = corr(R_odd_sm, R_even_sm, 'tail', 'right'); % null: cc<=0. alternative hypothesis: cc>0.
    
%     [~, ori_rep_pval] = regressionSlopeTtest(x, y, alpha, '+', showWorkingFlag);    

end

function idx_ori_use = getOriIndexInRange(oris_deg, ori_pref_deg, oriMax, ori_diff_max)
    if isnan(ori_pref_deg)
        idx_ori_use = [];
        return;
    end

    ori_diff_max = min(ori_diff_max, oriMax/2);

    oris_deg_wrapped = [oris_deg - oriMax, oris_deg, oris_deg + oriMax];            
    diff_fromPref = oris_deg_wrapped - ori_pref_deg;
    all_inRange = -ori_diff_max < diff_fromPref & diff_fromPref <= ori_diff_max;        
    idx_ori_use = [find(all_inRange(:,1)); find(all_inRange(:,2));  find(all_inRange(:,3)) ];

    
    if (oriMax == 180)
        assert(length(idx_ori_use) == length(oris_deg) );
    elseif (oriMax == 360) && (ori_diff_max < 180)  % drifting gratings
        assert(length(idx_ori_use)/length(oris_deg) == 1/2 );        
    end
    3;

    
%     % make sure differences are all 
%     idx_tooHi = (ori_diffs > oriMax/2);
%     idx_tooLo = (ori_diffs < -oriMax/2);
%             
%     ori_diffs(idx_tooHi) = ori_diffs(idx_tooHi)-oriMax;
%     ori_diffs(idx_tooLo) = ori_diffs(idx_tooLo)+oriMax;
%     
%     assert( isequal(abs(ori_diffs), circDist(ori_diffs_orig, oriMax)) );    
end


function osi = getOriSelectivityIndex(r_k, oris_deg, idx_ori_use)        
    % for orientations: F1 / DC
%     if oriMax == 180
%         f_harmonic = 1;        
%     elseif oriMax == 360
%         f_harmonic = 2;        
%     end
%     [fn, f0] = getF1oDC(oris_deg(idx_ori_use), r_k(idx_ori_use), oriMax, f_harmonic);
%     osi = fn/f0;
    if isempty(idx_ori_use)
        osi = nan;
        return;
    end
    
    osi = getF1oDC( r_k(idx_ori_use) );
    
    % for directions: F2 / DC
        
end


function CV = getCircularVariance(r_k_full, dirs_deg, idx_ori_use)
    %%    
    if isempty(idx_ori_use)
        CV = nan;
        return;
    end
    
    oris_deg = mod(dirs_deg, 180);    
    theta = deg2rad( oris_deg(idx_ori_use) );
    r_k   = r_k_full(idx_ori_use);
    
    cos_term = sum(r_k .* cos(2.*theta));
    sin_term = sum(r_k .* sin(2.*theta));
    
    CV = (1- sqrt(cos_term^2 + sin_term^2)/sum(r_k));
    %%
    if CV < 0
        assert(abs(CV) < 1e-10);
        CV = 0;
    end
    assert(ibetween(CV, 0, 1));
end


function [w_ori_local, R_null, R_pref, DSI_local, result_txt, oriParams, oriParams_ci, otc_stats, r_k_plot, ori_peak_unconstrained] = ...
    getOriLocalWidth(r_k, r_k_std, dir_pref_deg, dirs_deg, oriMax, cellOK, opt, show_flag)
    % only include responses in half of the directions (within 90 of preferred ori)
    global GC
    w_ori_local = nan;
    R_null = nan;
    R_pref = nan;
    r_k_plot = nan;
    DSI_local = nan;
    ori_peak_unconstrained = nan;
    oriParams = struct;
    oriParams_ci = struct;
    otc_stats = struct('peakRatio', nan, 'peakDist', nan, 'gof', struct('rsquare', nan) );
    if ~cellOK
        result_txt = 'cell not ok'; % skip fitting if cell is not ori-tuned or reproducible
        return;
    end            
            
    show = exist('show_flag', 'var') && ~isempty(show_flag) && (show_flag == true);
    
%     show = 1;
    % find responses within [-90, 90] window
    dirs_deg_wrapped = [dirs_deg - oriMax, dirs_deg, dirs_deg + oriMax];        
    dir_range_max = 90;
    all_inRange = abs([dirs_deg_wrapped - dir_pref_deg]) < dir_range_max;        
    idx_inRange = [find(all_inRange(:,1)); find(all_inRange(:,2));  find(all_inRange(:,3)) ];                    
    
    dir_inRange = dirs_deg_wrapped(all_inRange);
    r_k_inRange = r_k(idx_inRange);
    r_k_std_inRange = r_k_std(idx_inRange);

    fitNonPrefTuningCurve = opt.oriFitNonPrefTuningCurve && oriMax == 360; % ie. drifting gratings.
    fitUnconstrainedMeanTuningCurve = opt.oriFitUnconstrainedMeanTuningCurve;
    
    if fitNonPrefTuningCurve % ie. drifting gratings
        dir_opp_deg = mod(dir_pref_deg+180, 360);
        all_oppRange = abs([dirs_deg_wrapped - dir_opp_deg]) < dir_range_max;        
        idx_oppRange = [find(all_oppRange(:,1)); find(all_oppRange(:,2)); find(all_oppRange(:,3))];
    
        dir_oppRange = dirs_deg_wrapped(all_oppRange);
        r_k_oppRange = r_k(idx_oppRange);
        r_k_std_oppRange = r_k_std(idx_oppRange);        
    end        
    
        
    % Estimate B
    r_k_smoothed_w1 = gaussSmooth(r_k, 1, [], 1);    
    B_min = min(r_k_smoothed_w1);  
    B_est = B_min;     

    % Find smoothed, interpolated ori-tuning curve, estimate A and sigma:
    oriSmoothW = opt.oriSmoothWidth_deg/diff(dirs_deg(1:2));

    r_k_smoothed = gaussSmooth(r_k, oriSmoothW, [], 1);
    r_k_inRange_smoothed = r_k_smoothed(idx_inRange);    
    
    itpFactor = 10;
    dir_inRange_itp = linspace(dir_inRange(1), dir_inRange(end), (length(dir_inRange)-1)*itpFactor+1);
    r_k_inRange_smoothed_itp = interp1(dir_inRange, r_k_inRange_smoothed, dir_inRange_itp, 'spline');    
    [r_k_smooth_itp_max, idx_max_sm_itp] = max(r_k_inRange_smoothed_itp); 
    
    ind_pref_inRange_itp = indmin(abs(dir_pref_deg - dir_inRange_itp));
        
    peak_est = r_k_inRange_smoothed_itp(ind_pref_inRange_itp);        
    half_height = (r_k_smooth_itp_max - B_min)/2 + B_min;
    
    idxAboveHalf = continuousGroupings(find(r_k_inRange_smoothed_itp > half_height));
    if length(idxAboveHalf) > 1
        idx_withMax = cellfun(@(i) any(i == idx_max_sm_itp), idxAboveHalf);
        idxAboveHalf = idxAboveHalf(idx_withMax);
    end
    peak_width = diff(dir_inRange_itp(1:2))*length(idxAboveHalf{1});
    
    sig_est = peak_width/2;
    sig_est = min(sig_est, 45);        
    A_est = peak_est-B_est;
    
    % calculate peak ratio, peak distance
    if ~any(idx_max_sm_itp == [1, length(r_k_inRange_smoothed_itp)])
        idx_aroundMax_itp = idx_max_sm_itp + [-1, 0, 1];    
        [dir_max_qitp, r_k_smooth_qitp_max] = quadInterp(dir_inRange_itp(idx_aroundMax_itp), r_k_inRange_smoothed_itp(idx_aroundMax_itp));
    else
        dir_max_qitp = dir_inRange_itp(idx_max_sm_itp);
        r_k_smooth_qitp_max = r_k_inRange_smoothed_itp(idx_max_sm_itp);
    end
    dir_pref_local = dir_max_qitp;
    peakRatio = (peak_est-B_min)/(r_k_smooth_qitp_max-B_min);        
    localPeakDiscrep_degrees = abs(dir_pref_local - dir_pref_deg);    
    otc_stats.peakRatio = peakRatio;
    otc_stats.peakDist = localPeakDiscrep_degrees;    
    
    if opt.excludeOffPeakMaximumCases  && ...
        ((localPeakDiscrep_degrees > opt.localPeakMaxDiscrep_degrees) || (peakRatio < opt.minPeakRatio))
        result_txt = 'off-peak maximum';
        return;
    end
            
    
%     show = 1;

        showShiftedTo0 = 1;
        showErrorBars = 0;
        showSmoothed = 0;
        
        if showShiftedTo0
%             dir_range_max_ext = 180;
            dir_range_max_ext = oriMax/2;
            all_inRange_ext = abs([dirs_deg_wrapped - dir_pref_deg]) < dir_range_max_ext;    
            idx_inRange_ext = [find(all_inRange_ext(:,1)); find(all_inRange_ext(:,2));  find(all_inRange_ext(:,3)) ];            
            dir_inRange_ext = dirs_deg_wrapped(all_inRange_ext);
            r_k_inRange_ext = r_k(idx_inRange_ext);
            r_k_std_inRange_ext = r_k_std(idx_inRange_ext);                        
            
            dirs_plot = dir_inRange_ext-dir_pref_deg;
            r_k_plot = r_k_inRange_ext;
            r_k_std_plot = r_k_std_inRange_ext;
            d_cent = 0;
            xlims = dir_inRange_ext([1, end])-dir_pref_deg;
            xticks = [-180:45:180];
        else
            dirs_plot = dirs_deg;
            r_k_plot = r_k;
            r_k_std_plot = r_k_std;                        
            d_cent = dir_pref_deg;
            xlims = lims(dirs_plot);
        end
        
    if show
        %%
                
        figure(1); clf; 
        plot(dirs_plot, r_k_plot, 'b.-'); hold on;        
        if showErrorBars
            errorbar(dir_inRange_ext-dir_pref_deg, r_k_inRange_ext, r_k_std_inRange_ext);
        end
%         plot(dir_inRange-dir_pref_deg, r_k_inRange, 'go-'); 
        if showSmoothed
            plot(dir_inRange, r_k_inRange_smoothed, 'm:')
            plot(dir_inRange-dir_pref_deg, r_k_inRange_smoothed, 'ms')
            plot(dir_inRange_itp-dir_pref_deg, r_k_inRange_smoothed_itp, 'm-');
        end
        xlim(xlims);
        set(gca, 'xtick', xticks);
        drawVerticalLine(d_cent);
        title(sprintf('Gid = %d. cellId = %d', GC(1), GC(2) ));
%         plot(dir_max_qitp-dir_pref_deg, r_k_smooth_qitp_max, 'kp');

        figure(2); clf;
%         plot(dirs_deg, r_k, 'b.-'); hold on;        
%         errorbar(dirs_deg, r_k, r_k_std);
%         plot(dir_inRange-dir_pref_deg, r_k_inRange, 'go-'); 
        polar(deg2rad(dirs_deg), r_k + 5)
        
        
%         fplot(@(x) gaussFunc(beta0, x+dir_pref_deg)+B_min, xlims, 'k:')        
       
        if fitNonPrefTuningCurve && 0;
            all_oppRange_ext = abs([dirs_deg_wrapped - dir_opp_deg]) < dir_range_max_ext;    
            idx_oppRange_ext = [find(all_oppRange_ext(:,1)); find(all_oppRange_ext(:,2));  find(all_oppRange_ext(:,3)) ];            
            dir_oppRange_ext = dirs_deg_wrapped(all_oppRange_ext);
            r_k_oppRange_ext = r_k(idx_oppRange_ext);
            r_k_std_oppRange_ext = r_k_std(idx_oppRange_ext);                        
            
            figure(11);
            plot(dir_oppRange_ext-dir_opp_deg, r_k_oppRange_ext, 'b.-'); hold on;        
            errorbar(dir_oppRange_ext-dir_opp_deg, r_k_oppRange_ext, r_k_std_oppRange_ext);
            plot(dir_oppRange-dir_opp_deg, r_k_oppRange, 'go-'); 
            3;
            figure(12);
            polar(deg2rad([dir_oppRange_ext([1:end, 1])-dir_opp_deg]), r_k_oppRange_ext([1:end, 1]))
            
        end


        figure(2); clf;
        polar(deg2rad([dir_inRange_ext([1:end, 1])-dir_pref_deg]), r_k_inRange_ext([1:end, 1]))
        
        3;
    end


    % FIT to Gaussian with center at ori_pref
    gaussFun       = fittype( 'gaussOri(A, sigma, B, x)' );
    gaussFun_unconstrainedMean = fittype( 'gaussOri_withMean(A, mu, sigma, B, x)' );
%     doubleGaussFun = fittype( 'gaussOri(A1, A2, sigma, B, x)' );

    startPoint   = [A_est B_est sig_est];

    A_lims   = [0, inf];
    B_lims   = [B_min, inf];
    mu_lims = [-20, 20];
    sig_lims = [0, opt.oriMaxSigma_deg];        
    lower_bounds = [A_lims(1), B_lims(1), sig_lims(1)];
    upper_bounds = [A_lims(2), B_lims(2), sig_lims(2)];
           

    
%     useWgts = all(r_k_std_inRange > 0);
    useWgts = opt.useWeightsIfAllPositive && all(r_k_std > 0);
    wgts = [];
    wgts_opp = [];
    if useWgts
        wgts = 1./r_k_std_inRange;
        if fitNonPrefTuningCurve
            wgts_opp = 1./r_k_std_oppRange;
        end
    end
    
    fit_opts = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', lower_bounds, 'Upper', upper_bounds, 'StartPoint', startPoint, 'Weights', wgts);
    [cfun, gof, fit_results] = fit(dir_inRange-dir_pref_deg, r_k_inRange, gaussFun, fit_opts);
    
    if fitNonPrefTuningCurve
        fit_opts_opp = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', lower_bounds, 'Upper', upper_bounds, 'StartPoint', startPoint, 'Weights', wgts_opp);
        [cfun_opp, gof_opp, fit_results_opp] = fit(dir_oppRange-dir_opp_deg, r_k_oppRange, gaussFun, fit_opts_opp);
    end
    if fitUnconstrainedMeanTuningCurve
        startPoint_uncm   = [A_est 0, B_est sig_est];
        lower_bounds_uncm = [A_lims(1), mu_lims(1), B_lims(1), sig_lims(1)];
        upper_bounds_uncm = [A_lims(2), mu_lims(2), B_lims(2), sig_lims(2)];

        fit_opts_uncm = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', lower_bounds_uncm, 'Upper', upper_bounds_uncm, 'StartPoint', startPoint_uncm, 'Weights', wgts);
        [cfun_uncm] = fit(dir_inRange-dir_pref_deg, r_k_inRange, gaussFun_unconstrainedMean, fit_opts_uncm);
        ori_peak_unconstrained = mod(dir_pref_deg + cfun_uncm.mu, 360);        
        3;
        
    else
        ori_peak_unconstrained = nan;
    end
%     fit_opts_nowgt = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', lower_bounds, 'Upper', upper_bounds, 'StartPoint', startPoint, 'Weights', []);
    %     [cfun_nowgt, gof] = fit(dir_inRange-dir_pref_deg, r_k_inRange, gaussFun, fit_opts_nowgt);

%     prob_fit = gammainc(resnorm/2, df/2, 'upper');
                
    coefNames = coeffnames(cfun);
    coefVals  = coeffvalues(cfun);
    coef_ci   = confint(cfun);
                        
        
    if show
        %%
        figure(1); hold on;
        dir_plot_tmp = linspace(dir_inRange(1), dir_inRange(end), 360)-dir_pref_deg;
        h = plot(dir_plot_tmp, feval(cfun, dir_plot_tmp), 'r');            
        hx = xlabel('\phi - \phi_{pref}'); hy = ylabel('firing rate (Hz)');
        
        tit_str = sprintf('Gid = %d. cellId = %d. \\sigma = %.1f\\circ', GC(1), GC(2), cfun.sigma);
        ht = title(tit_str);
        set([hy, ht], 'fontsize', 11)
        set([hx], 'fontsize', 13)
%         title(sprintf('Gid = %d. cellId = %d. Ratio = %.2f. Dist = %.1f\\circ. R2 = %.2f', GC(1), GC(2), peakRatio, localPeakDiscrep_degrees, gof.rsquare));
%         title(sprintf('Gid = %d. cellId = %d. Ratio = %.2f. Dist = %.1f\\circ. R2 = %.2f', GC(1), GC(2), peakRatio, localPeakDiscrep_degrees, gof.rsquare));
        
        figure(2); hold on;
        thetas = dir_inRange_ext-dir_pref_deg; 
        rhos = feval(cfun, thetas);
        polar(deg2rad(thetas), rhos, 'r');
        title(tit_str);
    
        3;

    end
    
    w_ori_local = cfun.sigma; 
            
    oriParams    = cell2struct( num2cell(coefVals(:)),  coefNames(:), 1);
    oriParams_ci = cell2struct( num2cell(coef_ci', 2),  coefNames(:), 1);
                
    R_null = feval(cfun, 90);
    R_pref = feval(cfun, 0);
    
    if fitNonPrefTuningCurve
        A_pref = feval(cfun, 0);
        A_opp = feval(cfun_opp, 0);
        if A_opp > A_pref
            [A_pref, A_opp] = deal(A_opp, A_pref);
        end
        if ~opt.useAlternateDefOfDSIforUnrectified
            DSI_local = (A_pref-A_opp)/(A_pref+A_opp);
        else
            DSI_local = (A_pref-A_opp)/(A_pref);
        end
        if DSI_local > .999
            3;
        end
    else
        DSI_local = nan;
    end    
    

    
    result_txt = fit_results.message; % 'ok';
    otc_stats.gof = gof;
    if fitNonPrefTuningCurve
        otc_stats.sigma_opp = cfun_opp.sigma;
    end
    
end

function [dir_pref_deg, DSI] = getPreferredDirection(r_k, ori_pref_deg, dirs_deg, useAlternateDef_flag)

    dir_pref_deg = ori_pref_deg;
    dir_opp_deg = mod(ori_pref_deg + 180, 360);

    dirs_deg_wrapped = [dirs_deg - 360, dirs_deg, dirs_deg + 360];        
    dir_range_max = 90;
    
    idx_inRange_pref = any( abs([dirs_deg_wrapped - dir_pref_deg]) < dir_range_max, 2) ;  % use < (not <=) to exclude 90 deg directions.
    idx_inRange_opp  = any( abs([dirs_deg_wrapped - dir_opp_deg ]) < dir_range_max, 2) ;    
    
    N_pref = sum( r_k(idx_inRange_pref) );
    N_opp  = sum( r_k(idx_inRange_opp) );
    
    if N_opp > N_pref
        dir_pref_deg = dir_opp_deg;
        [N_opp, N_pref] = deal(N_pref, N_opp);
    end
    
    useAlternateDef = exist('useAlternateDef_flag','var') && isequal(useAlternateDef_flag, 1);
    if ~useAlternateDef
        DSI = (N_pref - N_opp)/(N_opp + N_pref);
    else
        DSI = abs(N_pref - N_opp)/abs(N_pref);
    end
    
    
%     show = any( abs(DSI - [.17, .40, .71]) < .01);
%     show = DSI > .85;
%         show = 1;
    show = 0;

    if show
        figure(93); clf;
        dirs_deg_plot = [dirs_deg(:); dirs_deg(1)];
        r_k_plot = [r_k(:); r_k(1)];
        polar(deg2rad(dirs_deg_plot-dir_pref_deg), r_k_plot, 'b.-');        
        dir_perps = 0+[-90, 90];
%         polar(deg2rad(dirs_deg), r_k, 'b.-');        
%         dir_perps = ori_pref_deg+[-90, 90];
        
        r_max = max(r_k);
%         r_max = max(axis)*.9;
        hold on;
        h(1) = polar(deg2rad(dir_perps(1))*[1,1], [0, r_max], 'r-');        
        h(2) = polar(deg2rad(dir_perps(2))*[1,1], [0, r_max], 'r-');        
        set(h, 'linewidth', 2)
        title(sprintf('DSI = %.2f', DSI), 'fontsize', 13);
        3;
        
    end

end


function S = responseMagRelToBckg(R_full, bckgSamples, opt, gratingType, Gid, cellId) 
    % Test RESPONSE SIZE SIGNIFICANCE
    % test that the cell's response to the best stimulus is significantly higher than would be
    % expected by chance from the spontaneous activity.
%     global GC
    [nOri, nSp, nPh, nTrials] = size(R_full);
    bckgNSamples = 1000;

    % for drifting gratings, bckg bins are already the size of an entire cycle, so no need to average. 
    % for flashed gratings, bckg bin size is the same size as stimulus bin size.
    nBckgAv  = switchh(gratingType, {'flashed', 'drifting'}, [nPh, 1]);
    nTopStim = switchh(gratingType, {'flashed', 'drifting'}, [opt.responseSizeNTop_flashed, opt.responseSizeNTop_drifting]);
    
    
%     [~, idx_bestOriSpf] = maxElement(mean(mean(R_full, 4),3));
    
    show = 0;    
    % 1. (a) average over phases (b) used smoothed R (av phase, av trials, smoothed across ori&spf)
    %       to find best N ori&sp. (c) use all trials at best ori/sp
    doAllNTop = 1;
    
    
    all_nTop_stim = 1:5;
    
    R_avP = mean(R_full,3);
    R_avP_os = reshape(R_avP, [nOri*nSp, nTrials]);
    R_avP_avT_os = mean(R_avP_os, 2);
%     R_avP_avT_os = mean(R_avP_os,2);
    
%     Gid = GC(1); cellId = GC(2);
    
    idxMaxOS = indmax_n(R_avP_avT_os, max(all_nTop_stim));

    if doAllNTop
        B_avP_trialsAtBestOS_all = dbGetSavedSampledBckgDistrib(Gid, cellId, bckgSamples, nOri*nSp, nBckgAv, nTrials, max(all_nTop_stim), bckgNSamples);
        % dimensions of B_avP_trialsAtBestOS_all are : nTrials, nTopStim, nSamples]
        all_trials_idx = 1:nTrials;
        for i = 1:length(all_nTop_stim)   
            %%
            nTop = all_nTop_stim(i);
            R_avP_trialsAtBestOS_i = R_avP_os(idxMaxOS(1:nTop), all_trials_idx);            
            B_avP_trialsAtBestOS_i = B_avP_trialsAtBestOS_all(:, 1:nTop, :);
            if show
                %%
                figure(10+i); clf; hnd = hist2({R_avP_trialsAtBestOS_i(:), B_avP_trialsAtBestOS_i(:)}, 40, 'stairs', 'norm');
                set(hnd, 'linewidth', 2); set(hnd(2), 'color', 'r'); title(num2str(i));
            end
            S.pval_avP_trialsAtBestOS_U(i) = myRanksum(R_avP_trialsAtBestOS_i(:), B_avP_trialsAtBestOS_i(:), 'right');
            [h, S.pval_avP_trialsAtBestOS_t(i)] = ttest2(R_avP_trialsAtBestOS_i(:), B_avP_trialsAtBestOS_i(:), .01, 'right');
            3;
        end
    end
            
    if nTrials > 4        
        idx_use = find(nTopStim == all_nTop_stim, 1);
    else
        idx_use = find(nTopStim*2 == all_nTop_stim, 1);
    end
        
%     R_avP_trialsAtBestOS_main = R_avP_os(idxMaxOS(1:nTopStim), :);    
%     B_avP_trialsAtBestOS_main = dbGetSavedSampledBckgDistrib(Gid, cellId, bckgSamples, nBckgAv, nOri*nSp, nTopStim, bckgNSamples);
        
    switch lower(opt.responseSizeTest)
        case 'utest', S.pval = S.pval_avP_trialsAtBestOS_U(idx_use) ; %myRanksum(R_avP_trialsAtBestOS_main(:), B_avP_trialsAtBestOS_main(:), 'right');
        case 'ttest', S.pval = S.pval_avP_trialsAtBestOS_t(idx_use); % ttest2(R_avP_trialsAtBestOS_main(:), B_avP_trialsAtBestOS_main(:), .01, 'right');
    end
    
%     title(sprintf('Gid = %d. cellId = %d. Means (stim:%.2f, spont:%.2f). Medians(stim%.2f, spont:%.2f)', 2911, 3, ...
%         mean(R_avPhaseTrials(:)), mean(bckgNullDistrib(:)), median(R_avPhaseTrials(:)), median(bckgNullDistrib(:)) ));     
        
%     R = mean(R_full, 4);
%     S.pvalR_full          = myRanksum(R_full(:), bckgSamples(:), 'right');
    
    
    % old test: test that the response at the preferred orientation/sp is >3 * std.dev + spontaneous firing rate.    

    bckgMean = mean(bckgSamples);
    bckgStd  = std(bckgSamples);
    respToBestOriSp = mean(R_avP_os(idxMaxOS(1), :)); % response to best ori/sp (averaged over trials)
    S.nStdAboveBckg = (respToBestOriSp-bckgMean) / bckgStd;    
            
end




%%%%%% FUNCTIONS FOR SPATIAL FREQUENCY  %%%%%%%-----------------------------------------------------

        
function spfStats = doSpatialFrequencyBatchTests(R_full, bckgSamples, spfs_cpd, gratingType, opt, subtractBckgFlag, responseSize_S, ori_pref_idx, oeStats)

    [nOri, nSpf, nPh, nTrials] = size(R_full); 
    assert(nSpf > 1);
    
    alpha = opt.alpha;
    minRsqr = opt.minRsqr;
    
    bckgMeanRate = mean(bckgSamples);
    R_ori_sp = mean(mean(R_full, 3),4); 
    
    if isempty(ori_pref_idx)  % for drifting gratings spf batches, don't have all orientations,
%     
        ori_pref_idx = indmax(sum(R_ori_sp,2));
%         idx_bestSpf = indmax(R_ori_sp(idx_bestOri,:));                
    end
%     idxBestOriSpf = [idx_bestOri, idx_bestSpf];
    
    % get spatial frequency tuning curves for each trial. (average over phases) 
%     spfTCs = mean(mean(  R_full, 3),4); % average over phases, over trials
        
    % select orientation with strongest response
%     [tmp, ind_max] = maxElement(spfTCs);
%     ori_best = ind_max(1);        
        
%     if strcmp(gratingType, 'flashed')
% %         ori_pref_idx_ext = mod(ori_pref_idx + [-1, 0, 1], nOri);
%         ori_pref_idx_ext = ori_pref_idx;
%     else
%         ori_pref_idx_ext = ori_pref_idx;
%     end
%     end

    F1oDC = opt.F1oDC_calculated;
    F1oDC_stderr_jack = opt.F1oDC_stderr_calculated;

    % Spatial Frequency Reproducibilty    
    doUnrectifiedVersions = opt.doUnrectifiedVersionsForSpontSubtract && subtractBckgFlag;
    
    haveOEstats = exist('oeStats','var') && ~isempty(oeStats);

%     R_full_avPh_avTr = mean(mean(R_full, 3),4);
    %%
    R_full_prefOri = R_full(ori_pref_idx,:,:,:);
    R_full_prefOri_unrec = R_full_prefOri;
%     R_full_meanSpf = mean(R_full,2);
                    
    if subtractBckgFlag
        %1. r_k (Ori tuning curve (rk))  - for ori width, DSI,
        %2. R_full / (really: just r_k) for ori tuning test  -- should be just at pref spf for SS
        %3. R_full at pref_spf (with odd/even) for reproducibility

        % find the spatial frequencies, that, when averaged over phases & trials, are below background (for pref ori)
        R_full_prefOri_avPh_avTr = mean(mean(R_full_prefOri,3),4);
        idx_spf_zero_response_prefOri = R_full_prefOri_avPh_avTr < bckgMeanRate;                        
        
        R_full_prefOri = R_full_prefOri - bckgMeanRate;
        R_full_prefOri_unrec = R_full_prefOri;
        R_full_prefOri(:,idx_spf_zero_response_prefOri,:,:) = 0; % rectify based on rectified ori tuning curve;                
    end
    
    
    spfTC_1 = reshape(mean(mean(R_full_prefOri,3),4), [nSpf, 1]);
    spfTC_unrec = reshape(mean(mean(R_full_prefOri_unrec,3),4), [nSpf, 1]);
    
    %%
    if subtractBckgFlag        
%         r_k_dir_ss = rectified(r_k_dir - bckgMeanRate);        
        
        assert(~any(spfTC_1(idx_spf_zero_response_prefOri)));        
        3;
    end    
    
    
%     [spfTC, spfTC_std, spf_rep_pval] = getSpfFreqReproducibility(R_full, ori_pref_idx, gratingType);
    [spfTC, spfTC_std, spf_rep_pval] = getSpfFreqReproducibility(R_full_prefOri, 1, gratingType);
    assert(isequalToPrecision(spfTC_1, spfTC, 1e-8))
    
    
    if opt.doJackKnifeStdErr
        %%
        spfTC_trials = reshape( mean(R_full_prefOri,3), [nSpf, nTrials]);
        assert(isequalToPrecision(mean(spfTC_trials,2), spfTC, 1e-10));
        
        spfTC_jacks = jackknifeAverageTrials( spfTC_trials, 2 );
%         spfTC_jacks2 = arrayfun(@(i) mean(  spfTC_trials(:, setdiff(1:nTrials, i) )  ,2), 1:nTrials, 'un', 0);
%         assert(isequal(spfTC_jacks, spfTC_jacks2))
    end

    
    

    spfTC_stderr = spfTC_std / sqrt(nTrials);
    
    assert(max(abs(spfTC_1-spfTC)) < 1e-10);

    if subtractBckgFlag                
        assert(~any(spfTC(idx_spf_zero_response_prefOri)));
        assert(~any(spfTC_std(idx_spf_zero_response_prefOri)));
        3;
    end    
    
    
    % Magnitude of best response compared to background (spontaneous) activity.
%     responseSize_S = responseMagRelToBckg(R_full, bckgSamples, stimType.gratingType);    
    if ~isempty(responseSize_S)
        response_size_pval = responseSize_S.pval;
    else
        response_size_pval = 0;
    end
    
    cellOK_soFar = (spf_rep_pval(1) < alpha) && (response_size_pval < alpha);
        
    
    [result_txt, w_spf, f_opt, normBw, spfParams, spfParams_ci, fit_stats, stcProps] = ...
        getSLNfit(spfs_cpd, spfTC, spfTC_stderr, opt, cellOK_soFar, 0, opt.doSLNcurveWithNoConstraints);  % Measure spatial frequency tuning width (if reproducible)      
    spf_fit_rsqr = fit_stats.rsquare;

%     if opt.doSLNcurveWithNoConstraints
%         
%         [result_txt2, w_spf2, f_opt2, normBw2, spfParams2, spfParams_ci2, fit_stats2, stcProps2] = ...
%             getSLNfit(spfs_cpd, spfTC, spfTC_stderr, opt, cellOK_soFar, 1, noConstraint);  % Measure spatial frequency tuning width (if reproducible)      
%         
%         spfParams.paramsNoConstraint = spfParams2;
%     end
    
    f_opt_distFromStimRange = stcProps.PrefDistFromStim;
    GoodFitToSLNcurve = (spf_fit_rsqr > minRsqr);
    
    cellOK = cellOK_soFar && GoodFitToSLNcurve && (f_opt_distFromStimRange == 0);

%     if ~strncmp(result_txt, 'Success', 7)  % for debugging
%         [result_txt, w_spf, f_opt, normBw, spfParams, spfParams_ci, fit_stats] = ...
%             getSLNfit(spfs_cpd, spfTC, spfTC_stderr, opt, cellOK_soFar, 'err');  % Measure spatial frequency tuning width (if reproducible)                      
%     end
            
    % calculate F1/DC
%     spf_pref_idx = indmax(spfTC);
%     phs = linspace(0, 360, nPh+1); phs = phs(1:nPh);
%     phaseTC = mean( R_full(ori_pref_idx, spf_pref_idx, :, :), 4);
%     F1oDC = getF1oDC(phs, phaseTC(:), 360);
    if isnan(f_opt)
        spf_peak_idx = nan;
    else
        spf_peak_idx = indmin( abs(log2(spfs_cpd) - log2(f_opt) ) );
    end


    [w_spf_err, w_spf_err_rel, f_opt_err, f_opt_err_rel] = deal(nan);
    if haveOEstats
        %%
        [w_spf_err, w_spf_err_rel] = diffOverAv(oeStats(1).w_spf, oeStats(2).w_spf);
        [f_opt_err, f_opt_err_rel] = diffOverAv(log2(oeStats(1).f_opt), log2(oeStats(2).f_opt));

    end
    error_est_OE = struct('w_spf_err', w_spf_err, 'w_spf_err_rel', w_spf_err_rel, ...
                          'f_opt_err', f_opt_err, 'f_opt_err_rel', f_opt_err_rel);

                    
    stderr_jack = struct;
    if opt.doJackKnifeStdErr
        %%
        [w_spf_jacks, f_opt_jacks] = deal(  zeros(1,nTrials) );
        
        for jack_i = 1:nTrials
            [~, w_spf_jacks(jack_i), f_opt_jacks(jack_i)] = ...
                getSLNfit(spfs_cpd, spfTC_jacks{jack_i}, spfTC_stderr, opt, cellOK_soFar, 0, opt.doSLNcurveWithNoConstraints);  % Measure spatial frequency tuning width (if reproducible)      
        end
        
        w_spf_stderr_jack = jackknifeStdErr(w_spf_jacks, w_spf);
        f_opt_stderr_jack = jackknifeStdErr(f_opt_jacks, f_opt, 'log');
        stderr_jack = struct('w_spf', w_spf_stderr_jack, 'f_opt', f_opt_stderr_jack, 'F1oDC', F1oDC_stderr_jack);
    end

    origValues = struct('w_spf', w_spf, 'f_opt', f_opt, 'F1oDC', F1oDC);

    if opt.applyStdErrorThresholds
        if w_spf_stderr_jack > opt.maxSpfStdErr_oct
            w_spf = nan;
        end
        if f_opt_stderr_jack > opt.maxSpfStdErr_oct
            f_opt = nan;
        end        
        if F1oDC_stderr_jack > opt.maxF1oDCStdErr 
            F1oDC = nan;
        end        
    end

    spfStats = struct('cellOK', cellOK, 'cellOK_asMU', cellOK_soFar,...
                      'spf_rep_pval', spf_rep_pval, 'response_size_pval', response_size_pval, 'rsqr', spf_fit_rsqr, ...
                      'GoodFitToSLNcurve', GoodFitToSLNcurve, 'f_opt_distFromStimRange', f_opt_distFromStimRange, ...                      
                      'response_sig', responseSize_S, 'SLNfitResult', result_txt, 'stcProps', stcProps, ...    
                      'w_spf', w_spf, 'f_opt', f_opt, 'Bn', normBw, 'SLNparams', spfParams, 'spfParams_ci', spfParams_ci, ...
                      'F1oDC', F1oDC, ...  % use instead F1/DC calculated at smoothed OSP
                      'ori_pref_idx', ori_pref_idx, 'spf_peak_idx', spf_peak_idx, ...
                      'error_OE', error_est_OE, ...
                      'error_jack', stderr_jack, ...
                      'spf_tc', spfTC, ...
                      'spfTC_stderr', spfTC_stderr, ...
                      'orig', origValues ...
                      );
    
end

function [spfTC, spfTC_std, spf_rep_pval] = getSpfFreqReproducibility(R_full, ori_pref_idx, gratingType)
    [nOri, nSpf, nPh, nReps] = size(R_full);
%     nOri_test = length(ori_pref_idx)
%     spfTC_trials = reshape( mean(R_full(ori_pref_idx,:,:,:),3), [nSpf, nReps]); % average over phases, but not trials.
    
    phaseAction = switchh(gratingType, {'flashed', 'drifting'}, {'counterphase', 'average'});
    
    spfTC_ph_trials = reshape( R_full(ori_pref_idx,:,:,:), [nSpf, nPh, nReps]);     
    
    % average all phases
    spfTC_avPh = reshape( mean(spfTC_ph_trials, 2), [nSpf, nReps]);
    spf_rep_pval(1) = pval_odd_even(spfTC_avPh);
    
    % odd/even phases
    idx_oe = {1:2:nPh, 2:2:nPh};
    spfTC_oe = cellfun(@(i) mean(spfTC_ph_trials(:,i,:), 2), idx_oe, 'un', 0);
    spfTC_oe = reshape( cat(1, spfTC_oe{:}), [nSpf*2, nReps ]);
    spf_rep_pval(2) = pval_odd_even(spfTC_oe);
    
    % all phases individually
    spfTC_indivPh = reshape( spfTC_ph_trials, [nSpf*nPh, nReps]);
    spf_rep_pval(3) = pval_odd_even(spfTC_indivPh);
    
    
    % 'counterphase'  % average counter-phase phases
    idx_cph = arrayfun(@(i) [i, i+nPh/2], 1:nPh/2, 'un', 0);                
    spfTC_cph = cellfun(@(i) mean(spfTC_ph_trials(:,i,:), 2), idx_cph, 'un', 0);
    spfTC_cph = reshape( cat(1, spfTC_cph{:}), [nSpf*nPh/2, nReps ]);    
    spf_rep_pval(4) = pval_odd_even(spfTC_cph);

    % consecutive phases
    idx_consec = arrayfun(@(i) [i, i+1], 1:2:nPh, 'un', 0);
    spfTC_consec = cellfun(@(i) mean(spfTC_ph_trials(:,i,:), 2), idx_consec, 'un', 0);
    spfTC_consec = reshape( cat(1, spfTC_consec{:}), [nSpf*nPh/2, nReps ]);
    spf_rep_pval(5) = pval_odd_even(spfTC_consec);
    
        
%     n_k_ij_cPh = cellfun(@(i) mean(n_k_ij_nPhXnRep(:, i, :),2), idx_cph, 'un', 0);
%     n_k_ij_cPh = reshape( [n_k_ij_cPh{:}], [nOri, nReps*nPh/2]);
%     n_k_ij = n_k_ij_cPh;
    
%     idx_odd  = 1:2:nReps;
%     idx_even = 2:2:nReps;
%     spfTC_odd  = double( mean(spfTC_trials(:,idx_odd),2) );
%     spfTC_even = double( mean(spfTC_trials(:,idx_even),2) );
%     
%     [tmp, spf_rep_pval] = pearsonR(spfTC_odd, spfTC_even);        
    
    spfTC_trials = mean(spfTC_ph_trials, 3);
    spfTC = mean(spfTC_trials,2);
    spfTC_std = std(spfTC_trials, [], 2);
    
    
end

function pval = pval_odd_even(spfTC_cph_trials)
    [n1, nReps] = size(spfTC_cph_trials);
    idx_odd  = 1:2:nReps;
    idx_even = 2:2:nReps;
    spfTC_odd  = double( mean(spfTC_cph_trials(:,idx_odd), 2) );
    spfTC_even = double( mean(spfTC_cph_trials(:,idx_even),2) );
    
    [~, pval] = corr(spfTC_odd, spfTC_even, 'tail', 'right');    
end






function ori_tc_sm = smoothedR(ori_tc, coefs)

    n = length(ori_tc);

    coefs = coefs(:)'/sum(coefs);
    nCoefs = length(coefs);
    if ~odd(nCoefs), 
        error('Must have odd number'); 
    end;
    idx = [1:nCoefs]-ceil(nCoefs/2);
        
    all_idx = mod(bsxfun(@plus, idx(:), 0:n-1), n)+1;
        
    ori_tc_sm = zeros(size(ori_tc));
    for i = 1:n
        ori_tc_sm(i) = coefs * ori_tc(all_idx(:,i));
    end    

end
    
    
function [rep, rep_pval, rep_pval_sgn, rep_str] = signRegressionTest(x, y, alpha, showWorkingFlag)


    [rep, rep_pval, rep_slope] = regressionSlopeTtest(x, y, alpha, '+', showWorkingFlag);
    if rep_slope < 0
        [tmp1, rep_pval2] = regressionSlopeTtest(x, y, alpha, '-');
        rep_pval_sgn = 1/rep_pval2;
    elseif rep_slope >= 0
        rep_pval_sgn = rep_pval;
    end

    rep_str = rep_slope;
    if rep_slope < 0
        rep_str = max(rep_slope, 1/rep_slope);
    elseif rep_slope > 0
        rep_str = min(rep_slope, 1/rep_slope);
    end                            

%     if ~isempty(showWorkingFlag) && showWorkingFlag
%         xlims = xlim; ylims = ylim; L = max(xlims(2), ylims(2));
% %         L = roundToNearest(max([x(:); y(:)]), 5, 'up');        
%         axis equal;
%         axis([0 L 0 L]); 
%     end    
end

function [x_min, y_min] = quadInterp(x, y)
    % find the x coordinate of the min (or max) of a parabola given 3 points on the parabola.
    
    % for stability/accuracy, shift the x-axis to the center of the 3 points.
    x_shft = -x(2);
    x = x+x_shft;
    
    Dx = [x(2)^2 - x(1)^2,  x(2) - x(1); 
          x(3)^2 - x(1)^2,  x(3) - x(1)];
    dy = [y(2) - y(1); 
          y(3) - y(1)];
    coef = Dx\dy;
    a = coef(1); b = coef(2); 
    c = y(1) - a*x(1)^2 - b*x(1); 
    
    x_min = -b/(2*a);    
    
    if x_min < x(1) || x_min > x(3)
        keyboard;
        plot(x, 10*y, 'ko')
    end
    y_min = a*x_min.^2 + b*x_min + c;
    
    x_min = x_min-x_shft;
end



%{
    minNforCmp = (nOri*nSpf)* (.3);
    if length(indOriSpfSig) > minNforCmp % just take top 20%
        inds = ord(R_ori_spf(:), 'descend');
        indOriSpfSig = inds(1:minNforCmp);        
    end
        
    if length(indOriSpfSig) > minNforCmp
        responses  = mean( R_full(indOriSig, indSpfSig, :,:) ,3); % at peak ori, sum across phases
        oriSpfInds = oriSpfTuningCurve > max(oriSpfTuningCurve)*fracOfMax;
        if (nnz(oriSpfInds) < minNforCmp) % maybe was a little too strict.
            if (nSpf > minNforCmp)  
                oriSpfInds = oriSpfTuningCurve > max(oriSpfTuningCurve)*fracOfMax/2;  % relax requirements a little
                if (nnz(oriSpfInds) < minNforCmp) && (length(oriSpfInds) >= minNforCmp)  % if still too few data points
                    [tmp, ind] = sort(oriSpfTuningCurve, 'descend');
                    oriSpfInds = ind(1:minNforCmp);  % just take most top N most responsive.
                end
            else
                oriSpfInds = 1:nSpf;
            end
        end
        
    end
  %}  

%{
    f_lo = fzero(@(f) cell_SLNfunc(f)-r_halfMax, half_lo_est);
    if isnan(f_lo)
        result_txt = 'nan f_lo';
        return;
    end

    f_hi = fzero(@(f) cell_SLNfunc(f)-r_halfMax, half_hi_est);    
    ii = 0;    
    while f_hi == f_lo
        f_hi = fzero(@(f) cell_SLNfunc(f)-r_halfMax, f_peak + ii*(f_peak-f_lo));    
        ii = ii+1;
        if ii > 20
            error('something is wrong')
        end
    end


%     r_end = logSLNfunc(beta, spfs_cpd(end));
%     if r_end > r_halfMax
%         result_txt = 'Above r_peak/2 at 4';
%         return;
%     end
        
%     idx_f_afterPeak = find(spfs_cpd > f_peak, 1);
%     idx_f_belowHalf = find(spfTC(idx_f_afterPeak:end) < r_peak/2, 1) + idx_f_afterPeak-1;

%}

%{ 
difference between Preferred direction and local peak in orientation
    Gid, CellId, Difference
    1170, 6       40
    1528, 2       33
    2318, 3       36
    2330, 2       24
    3909, 5       28 (weird cell) 

    2474, 1       20.96, and cell looks ok.
    So a threshold of 23 degrees seems reasonable
%}

%{
    R_avP_maxOS = max(R_avP_os, [], 1);      
    B_avP_maxOS = sampleBckgDistrib_av_max(bckgSamples, nBckgAv, nOri*nSp, 1000);

    figure(11); clf; hist2({R_avP_maxOS, B_avP_maxOS}, 40, 'stairs', 'norm');
    S.pval_avP_maxOS_U = myRanksum(R_avP_maxOS(:), B_avP_maxOS(:), 'right');
    [h, S.pval_avP_maxOS_t] = ttest2(R_avP_maxOS(:), B_avP_maxOS(:), .01, 'right');
    


    function bckgSmp = sampleBckgDistrib_av_max(bckgSamples, nAv, nMax, nSamples)
    %     for i = 1:100
    %         [uVal, valCount] = uniqueCount(bckgSamples);
    %         valProb = valCount/sum(valCount);
    %         cumProb = [0, cumsum(valProb)];        
    %         allVals = uVal(binarySearch( cumProb, rand(nTogether, nSamples), [], -1));
    %     end
        allVals = reshape(randsample(bckgSamples, nAv*nMax*nSamples, true), [nAv, nMax, nSamples]);

        bckgSmp = max(mean(allVals, 1), [], 2);    
        bckgSmp = bckgSmp(:);
    end
%}


%{
function i = indmax_n(x, n)
    if n == 1
        i = indmax(x);
    else
        [~, idx] = sort(x, 'descend');
        i = idx(1:n);
    end
end


function bckgSmp = sampleBckgDistrib_av_allAtMax(bckgSamples, nAv, nMax, all_nTop_stim, nSamples)
    rand('state', 0);
    allVals = reshape(randsample(bckgSamples, nAv*nMax*nSamples, true), [nAv, nMax, nSamples]);
    
    bckgSmp_avPh = mean(allVals, 1);    
    bckgSmp_avPh_avSamp = mean(bckgSmp_avPh, 3);
    
    bckgSmp = cell(1,length(all_nTop_stim));
    for i = 1:length(all_nTop_stim)
        bestStim_idxs = indmax_n(bckgSmp_avPh_avSamp, all_nTop_stim(i));        
        bckgSmps = bckgSmp_avPh(1, bestStim_idxs, :);        
        bckgSmp{i} = bckgSmps(:);
    end
    if length(all_nTop_stim) == 1
        bckgSmp = bckgSmp{1};
    end
        
end
%}

%{
%     ind_max = indmin(abs(dir_pref_deg - dirs_deg));
%     ind_near_max = ind_max + [-1, 0, 1];
%     ind_near_max = mod(ind_near_max-1, length(dirs_deg))+1;
    ind_90degLeft_from_max = indmin(abs( mod(dir_pref_deg - 90, oriMax) - dirs_deg)) + [-5:5];        
    ind_90degLeft_from_max = mod(ind_90degLeft_from_max-1, length(dirs_deg))+1;

    ind_90degRight_from_max = indmin(abs( mod(dir_pref_deg + 90, oriMax) - dirs_deg)) + [-5:5];        
    ind_90degRight_from_max = mod(ind_90degRight_from_max-1, length(dirs_deg))+1;
%}



%{

%     B_min = min(r_k);  
%     B_est = mean(r_k([ind_90degLeft_from_max, ind_90degRight_from_max]))-B_min;
    B_est = B_min; %mean(r_k([ind_90degLeft_from_max, ind_90degRight_from_max]))-B_min;
    




% f(theta) = A( (theta - dir_pref)^2/ 2sigma^2) + B
%     gauss = @(A, sigma, B, theta) A .* exp( -(theta - dir_pref_deg).^2 ./(2*(sigma).^2)) + abs(B)+B_min;
    gauss = @(A, sigma, B, theta) A .* exp( -(theta - dir_pref_deg).^2 ./(2*(sigma).^2)) + B;    
    gaussFunc = @(beta, th) gauss(beta(1), beta(2), beta(3), th);        

    [beta,resnorm,residual,~,~,~,jacobian] = ...
        lsqcurvefit(gaussFunc, beta0, dir_inRange, r_k_inRange, lower_bnds, upper_bnds, optimset('display', 'off'));    
    beta_ci = nlparci(beta, residual, 'jacobian', jacobian);



%     beta = nlinfit(dir_inRange, r_k_inRange-B_min, gaussFunc, beta0);        


    msgstr = lastwarn;
    if ~isempty(msgstr)
        % try one more time, starting off with estimated coefficients, and adding points on end
%         lastwarn('');
        
%         dir_inRange_plus = [dir_inRange(1)-90 + [0; 5; 10]; dir_inRange; dir_inRange(end)+90- [0; 5; 10]];
%         r_k_inRange_plus = [B_min*[1;1;1]; r_k_inRange; B_min*[1;1;1]];
%         dir_inRange_plus = [dir_inRange(1)-90; dir_inRange; dir_inRange(end)+90];
%         r_k_inRange_plus = [B_min; r_k_inRange; B_min];
        
%         r_k_inRange_smoothed_plus = [B_min*[1;1;1]; r_k_inRange_smoothed; B_min*[1;1;1]];   
        
        
%         beta = nlinfit(dir_inRange, r_k_inRange-B_min, gaussFunc, beta);        
%         beta = nlinfit(dir_inRange_plus, r_k_inRange_smoothed_plus, gaussFunc, beta);        
    end                

%     for i = 1:length(errorIfGetWarnings)
%         warning('on', errorIfGetWarnings{i});
%     end
    
%     [msgstr, msgid] = lastwarn;
% %     if any(strcmp(msgid, errorIfGetWarnings))       
%         result_txt = msgid;
%         return;
%     end


        if show
%             beta_exp = nlinfit(spfs_cpd,     spfTC, SLNfunc,    beta0_exp);
%             lower_bnds(1) = 0;
%             beta_exp = lsqcurvefit(SLNfunc,  beta0_exp, spfs_cpd, spfTC, lower_bnds, upper_bnds, optimset('display', 'off'));
        end

        for i = 1:length(suppressWarnings)
            warning('on', suppressWarnings{i});
        end

        [msgstr, msgid] = lastwarn;
        if any(strcmp(msgid, errorIfGetWarnings))
            throw(MException(msgid, msgstr));
        end



%         figure(11); clf;
%         plot(dir_inRange_ext, r_k_inRange_ext, 'b.-');
%         drawVerticalLine(dir_pref_deg); xlim(dir_inRange_ext([1, end]))
%         set(gca, 'xtick', [-180:45:180]);
%         
%         figure(12); clf;        
%         polar(deg2rad([dirs_deg; 0]), [r_k; r_k(1)]); hold on;
%         polar(deg2rad([dir_pref_deg, dir_pref_deg]), [0 max(r_k)], 'r-')

%}


%{
%     errorIfGetWarnings = {'stats:nlinfit:IterationLimitExceeded', 'stats:nlinfit:IllConditionedJacobian', 'MATLAB:rankDeficientMatrix', 'MATLAB:nearlySingularMatrix', 'stats:nlinfit:ModelConstantWRTParam'};
%     acceptIfGetWarnings = {'MATLAB:singularMatrix'};
%     suppressWarnings = [errorIfGetWarnings, acceptIfGetWarnings];
        
%     try
%         for i = 1:length(suppressWarnings)
%             warning('off', suppressWarnings{i});
%         end
%         lastwarn('');

        
%         [beta_nlin, r,J,COVB,mse] = nlinfit(spfs_cpd,  spfTC, logSLNfunc, beta0);
%         ci1 = nlparci(beta_nlin, r, 'jacobian', J);
%         w_lims = [0.2, inf];
%         lower_bnds = [logf_peak_lims(1), r_max_lims(1), r_bkg_est(1), w_lims(1), s_lims(1)];


%         beta = nlinfit(log_spfs_cpd, spfTC, logSLNfunc, beta0 );
%         [beta,resnorm,residual,exitflag,output,lambda,jacobian] = ...
%             lsqcurvefit(logSLNfunc, beta0, log_spfs_cpd, spfTC, lower_bnds, upper_bnds, optimset('display', 'off'));
        
%         beta_ci = nlparci(beta, residual, 'jacobian', jacobian);
%}

%     all_idxMax = findLocalMaxima(spfTC_sm, 1, [], [], 1);
%     % remove first & last maxima in smoothed/interpolated curve if not a max in original curve:
%     if (length(all_idxMax_sm) > 1) && (all_idxMax_sm(end) == length(spfTC_sm_itp)) && ~any(all_idxMax == length(spfTC) )       
%         all_idxMax_sm = all_idxMax_sm(1:end-1);
%     end
%     if (length(all_idxMax_sm) > 1) && (all_idxMax_sm(1) == 1) && ~any(all_idxMax == 1)
%         all_idxMax_sm = all_idxMax_sm(2:end);
%     end            
%     all_maxVal_sm = spfTC_sm_itp(all_idxMax_sm);


%         if stcProps.MaxAtLowest
%             result_txt = 'Peak spf at 0 cyc/deg';
%             return;
%         end    
%         if stcProps.MaxAtHighest
%             result_txt = 'Peak at highest spf';
%             return;
%         end    
%         if spfTC(end) >= maxSpfTC * 0.9
%             result_txt = 'Close to max at highest spf';
%             return;
%         end    


%     if (f_peak < spfs_cpd(1))
%         result_txt = 'Preferred spf is lower than lowest stimulus spf';
%         return;
%     end
%     if (f_peak > spfs_cpd(end))
%         result_txt = 'Preferred spf is higher than highest stimulus spf';
%         return;
%     end
%     
%     if r_0 > r_halfMax
%         result_txt = 'Above r_peak/2 at 0';
%         return;
%     end
%     if r_inf > r_halfMax
%         result_txt = 'Above r_peak/2 at inf';
%         return;
%     end    
%     if f_lo > spfs_cpd(1)
%         result_txt = 'f_lo_half is below measured spf (~Above r_peak/2 at 0)';  % these tend to overestimate the spf width.
%         return;   % about 90% of tuning curves have this!!!
%     end    


%{
%     hi = x_peak;
%     count1 = 1;
%     count2 = 1;
%     exitflag = 0;
%     nTries = 20;

    
        
    % try to find the half-max (lo) point at the best     
    
    % estimate f_lo
    if nargin >= 3
        offsets = lo_est - x_peak;
    else
        offsets = [];
    end    
    offsets = [offsets,  -.1*(-1.2).^[0:nTries-1]];
    
    lo = x_peak;    
    while (x_peak <= lo) && (count1 < nTries)
        [lo, ~, exitflag] = fzero(@(x) f(x), lo_estimate);
        if exitflag ~= 1
            lo = x_peak;
        end
        count1 = count1+1;
    end
    
    % estimate f_hi
    if nargin >= 4
        offset = hi_est - x_peak;
    else
        offset = .1;
    end
    
    while (x_peak >= hi) && (count1 < 20)
        [hi, ~, exitflag] = fzero(@(x) f(x), x_peak + offset);
        if exitflag ~= 1
            hi = x_peak;
        end        
        offset = offset * -1.5;
        count1 = count1+1;
    end
    
    if (count1 >= 15) || (count2 >= 15)
        3;
    end
%     if lo > hi
%         [lo, hi] = deal(hi, lo);
%     end
%}


%{
  changes made 11/2
    orientation selectivity tests --  now use resultant vectors (ignore trial variations)
    flashed gratings:
        estimating preferred orientation - use all spatial frequencies instead of preferred 
        use calculated ori_pref to select preferred ori for spf tuning curve.
    

%}

%{
    if strcmp(gratingType, 'drifting');
        [ori_pref_deg, ori_sel_pval_stats] = getPreferredOriFromOSP(R_full, dirs_deg, spf_pref_idx, gratingType, [], combine_flag);
%         [ori_pref_deg2, ori_sel_pval_stats2] = getPreferredOriFromOSP(R_full, dirs_deg, spf_pref_idx, gratingType, [], 0);
%         ori_sel_pval_stats.stats_no_combine = ori_sel_pval_stats2;
%         assert(abs(ori_pref_deg-ori_pref_deg2)<1e-7 );        
    else
        [ori_pref_deg, ori_sel_pval_stats] = getPreferredOriFromOSP(R_full, dirs_deg, spf_pref_idx, gratingType, [], combine_flag);
    end
%}

%{
%         R_full_meanSpf_avPh_avTr = mean(mean(R_full_meanSpf,3),4);
%         idx_ori_zero_response_meanSpf = R_full_meanSpf_avPh_avTr < bckgMeanRate;                        
%         
%         R_full_meanSpf = R_full_meanSpf - bckgMeanRate;
%         R_full_meanSpf(idx_ori_zero_response_meanSpf,:,:,:) = 0; % rectify based on rectified ori tuning curve;



    [ori_pref_deg_meanSpf, ori_sel_pval_stats_meanSpf] = getPreferredOriFromOSP(R_full_oriTuningTest, dirs_deg, spf_pref_idx, gratingType, combine_flag);
%     [ori_pref_deg_prefSpf, ori_sel_pval_stats_prefSpf] = getPreferredOriFromOSP(R_full_prefSpf, dirs_deg, spf_pref_idx, gratingType, combine_flag);
    
    ori_pref_deg = ori_pref_deg_meanSpf;
    ori_sel_pval_stats = ori_sel_pval_stats_meanSpf;
    ori_sel_pval_stats.ori_sel_pval_stats_prefSpf = ori_sel_pval_stats_prefSpf;
    ori_sel_pval_stats.ori_pref_deg_prefSpf = ori_pref_deg_prefSpf;


%}

%{
    Notes:
    Spf tuning curve examples:

    2905,2
    2877,9 - example of unconstrained fit.



%}