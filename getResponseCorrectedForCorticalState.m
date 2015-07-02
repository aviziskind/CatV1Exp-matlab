function [r_corrected, r_corrected_full, stats, corticalGain] = getResponseCorrectedForCorticalState(Gid, cellId, R_full_arg, trialsDivision, inputGainFunction)

    doSmoothingForDriftingGratings = 0;
%     error('tmp', 'tmp!')
%     Gid = 2288;
%     cellId = 3;
%     stateFitType = 'sin';
    
    if nargin < 3 || isempty(R_full_arg)
        if flashedOrDrifting(Gid) == 1
            prevRespType = curResponseType;
            curResponseType('raw');
            PSTHdata = getPSTHforCell(Gid, cellId);
            LR_bins = PSTHdata.timeWindow_bins;
            [L_bin, R_bin, windowProfile] = deal(LR_bins(1), LR_bins(2), PSTHdata.windowProfile);
            
            curResponseType(prevRespType);
        else
            [L_bin, R_bin, windowProfile] = deal(1, 1, []);
        end
        [allBins, psthVals, meanRate] = dbGetCellSpkStimHists(Gid, cellId);
        R_full = calcOspForPsthWindow(Gid, {allBins, psthVals}, L_bin, R_bin, false, windowProfile, 'osp_full', meanRate);
    else
        R_full = R_full_arg;
    end
    R_full = double(R_full);
    
    [nOri, nSpf, nPh, nTrials] = size(R_full);
    nStim = nOri*nSpf*nPh;
    if doSmoothingForDriftingGratings && flashedOrDrifting(Gid) == 2
        [smoothPhase_Width, smoothPhase_Method] = curPhaseSmoothing;
        R_full = smoothOSP_Phs(R_full, smoothPhase_Method, smoothPhase_Width);
    end 
    R_full_orig = R_full;
    
    show = 1;

    
    
    if nargin < 4 || isempty(trialsDivision)
        trialsDivision = 'all';
    end
    trialsDivision = switchh(trialsDivision, {'osp_ph', 'osp_ph_oe', 'osp_ph_hoe'}, {'all', 'oe', 'hoe', trialsDivision});
    
    stateFitType = 'sin';
    nStateTerms = 8;
    
    haveInputGain = nargin == 5 && ~isempty(inputGainFunction);

    recur_t = getStimulusRecurrenceTime(Gid);
    decor_t = getStimulusDecorrelationTime(Gid);
    if isequalToPrecision(recur_t, decor_t, 1e-3) % stimuli are ordered - hard to estimate cortical state properly --> skip this one.
        fprintf('Gid = %d : stimulus is ordered\n', Gid);
        r_corrected_full = [];
        r_corrected = [];
        stats = [];
%         return
    end
    if flashedOrDrifting(Gid) == 1
        [allGids, allCellIds] = getAllGids('f');
    else
        [allGids, allCellIds] = getAllGids('d');
    end
    cell_glob_idx = find(allGids(:) == Gid & allCellIds(:) == cellId);
    fprintf('\n\nCell %d/%d: Gid = %d, cellId = %d\n', cell_glob_idx, length(allCellIds), Gid, cellId);    
    
    adjustMeanState = 1;  % constantly adjust between iterations so that mean gain of firing rate = 1
%         minCorticalState = 0.05; % rectify so that doesn't go below this value.
        opt.minCorticalState = 0.10; % rectify so that doesn't go below this value.
        
        opt.rectifyType= 'softplus';
%         opt.rectifyType = 'hard';
        opt.softPlusScale = 50;
        rectifyFunc = @(x) rectifyToValue(x, opt.minCorticalState, opt.rectifyType, opt.softPlusScale);
    show_opt_details = 'off'; % default: 'final';

    % Get times of each stimulus.
    R_full_t = double( dbGetStimulusTimes(Gid) );
    R_full_t_orig = R_full_t;
    assert(isequal(size(R_full_t), [nOri, nSpf, nPh, nTrials])); 
    
%     if ~haveInputGain
        % Get cortical gain/state function        
        [state_func, state_S] = getCorticalState(Gid, cellId, stateFitType, nStateTerms);
%     end    
    
    %% 
    
    nTrials_orig = nTrials;
    nStim_orig = nStim;    
    
    if strcmp(trialsDivision, 'all')
        reshapeToVecs = @(R) reshape(R, [nStim, nTrials]);
        reshapeBack   = @(R) reshape(R, [nOri, nSpf, nPh, nTrials]);

        oe_str = '';
%         odd_even_trials_factor = 1;
    else
%         assert(~odd(nTrials));
        %%

        if strcmp(trialsDivision, 'oe')
            idx_1 = 1:2:nTrials; % odd
            idx_2 = 2:2:nTrials; % even
            oe_str = 'odd/even separate';

        elseif strcmp(trialsDivision, 'hoe')
            idx_halfOdd = odd (  floor([1:nTrials]/2)  );
            idx_1 = find( idx_halfOdd );
            idx_2 = find( ~idx_halfOdd );
            oe_str = 'half odd/even separate';
        end
        
        %%
%         if ~odd(nTrials)
        reshapeToVecs = @(R) reshapeToStimVectors(R, [nOri, nSpf, nPh, nTrials], idx_1, idx_2);
        reshapeBack   = @(R) reshapeBackToOrigSize(R, [nOri, nSpf, nPh, nTrials], idx_1, idx_2);
            
%         odd_even_trials_factor = 0.5;
        nTrials = ceil(nTrials / 2); 
        nStim = nStim * 2;
        
    end
    assert(isequaln(R_full_t, reshapeBack(reshapeToVecs(R_full_t))));
    assert(isequaln(R_full,   reshapeBack(reshapeToVecs(R_full))));
    
    R_full_t = reshapeToVecs(R_full_t);
    R_full = reshapeToVecs(R_full);

    c_state_all = state_func(R_full_t);    
    meanRate = nanmean(R_full(:));

    
    [cc_R_gain_full_init, p_R_gain_full_init] = corr(c_state_all(:), R_full(:), 'rows', 'complete');
    fprintf('Corr between fitted state and R_mn : cc = %.2f, p = %.3g\n', cc_R_gain_full_init, p_R_gain_full_init);

    stats.cc_R_gain_full_init = cc_R_gain_full_init;
    stats.p_R_gain_full_init = p_R_gain_full_init;
        
%     if (cc_p > 0.05)
%         warning('Cortical state is not a good fit to the data for this cell')
%         r_corrected = [];
%         corticalGain = [];
%         return
%     end
%     
    R = nanmean(R_full, 2);
    [~, idx_bestStim] = sort(R(:), 'descend');
    [~, idx_worstStim] = sort(R(:), 'ascend');

    %%
    showInitial = 0;
    if showInitial
        mean_state = nanmean(c_state_all, 1);    
        mean_r = nanmean(R_full, 1);

        n_top_stim = min(500, nStim);
        stim_start = 1;

        top_r = nanmean( R_full(idx_bestStim(stim_start+[1:n_top_stim]-1), :), 1);
        bot_r = nanmean( R_full(idx_worstStim(stim_start+[1:n_top_stim]-1), :), 1);

        figure(1); clf; hold on;
        [poly_coef, S2] = polyfit(mean_state, mean_r, 1);
        line1 = @(beta, x) beta(1)*x;        
        [m, r1] = nlinfit(mean_state, mean_r, line1, poly_coef(1));    

        S1.normr = sqrt( sum(r1.^2) );
        S1.df = length(mean_r)-1;
        pval = nestedFtest(S1, S2);

    %         if pval < .05 % likely that model 2 is correct
            fit_fun = @(x) polyval(poly_coef, x);
            fun_s = sprintf('y = %.1f x + %.1f', poly_coef(1), poly_coef(2));
    %         else
    %             fit_fun = @(x) line1(m, x);
    %             fun_s = sprintf('y = %.1f x', m);
    %         end

        idx = ord(mean_state, 'ascend');
        plot(mean_state(idx), top_r(idx), ['ro-']); 
        box on;
        hold on;
        plot(mean_state(idx), mean_r(idx), ['bo-']); 
        plot(mean_state(idx), bot_r(idx), ['ko-']); 
        xlabel('g(t)'); ylabel('Mean firing rate');

        xlims = [0, max(mean_state)];
        fplot(fit_fun, xlims, ['b:']);
    %         title( sprintf('Group %d, cell %d. (F-pval : %.2g) %s\n', Gid, cellId, pval, fun_s) );
        title( sprintf('Group %d, cell %d.', Gid, cellId) );
        hold on;
        plot(0,0, 'w')
        drawHorizontalLine(0);
        ylims = ylim; ylim([-.5, ylims(2)]);
        if 0
           %%
           legend({sprintf('Top %d stimuli', n_top_stim), 'All Stimuli', sprintf('Bottom %d stimuli', n_top_stim)}, 'location', 'best' );

        end
    end

    corticalStateFactor = @(c_states_in) getCorticalScaleFactor(c_states_in, rectifyFunc);

    if haveInputGain
        As = arrayfun(@(i) inputGainFunction.(sprintf('a%d', i)), 1:nStateTerms)';
        bs = arrayfun(@(i) inputGainFunction.(sprintf('b%d', i)), 1:nStateTerms)';
        cs = arrayfun(@(i) inputGainFunction.(sprintf('c%d', i)), 1:nStateTerms)';
      
    else
        S_gainParams = cell2struct(num2cell(coeffvalues(state_func))', coeffnames(state_func), 1);
        As = arrayfun(@(i) S_gainParams.(sprintf('a%d', i)), 1:nStateTerms)';
        bs = arrayfun(@(i) S_gainParams.(sprintf('b%d', i)), 1:nStateTerms)';
        cs = arrayfun(@(i) S_gainParams.(sprintf('c%d', i)), 1:nStateTerms)';
    end
    
    useAdditiveConstant = 1 && 0;
    optimizeB = false;
    optimizeC = false;
    
        
    if adjustMeanState && ~haveInputGain
        init_scale_factor = corticalStateFactor(c_state_all);

        c_state_all = rectifyFunc( c_state_all * init_scale_factor );
        assert( isequalToPrecision (nanmean(c_state_all(:)), 1, 1e-5) );

%         c_state_all = c_state_all * init_scale_factor;
%         c_state_all(c_state_all < minCorticalState) = minCorticalState;
%         scale_f2 = corticalStateFactor(c_state_all);
%         assert( abs(scale_f2 - 1) < 1e-3 )

        As = As * init_scale_factor;    
    end

    A0 = [As];

    if optimizeB
        A0 = [A0; bs];
    end
    if optimizeC
        A0 = [A0; cs];
    end
    
    const0 = 0;
    if useAdditiveConstant
        A0 = [const0; A0]; %#ok<*AGROW>
    end
    R0 = R;
    R_A0 = [R0(:); A0];

    c_state_v = c_state_all(:);
    R_full_t_rep = repmat(R_full_t(:), [1, nStateTerms]);


%%      
    switch stateFitType
        case 'sin'
            bxplusC = bsxfun(@plus, bsxfun(@times, bs(:)', R_full_t_rep ), cs(:)') ;
            corticalStateSines = sin(bxplusC);
            corticalStateTerms = corticalStateSines;
        case 'gauss'
            corticalStateGausses = exp( -(bsxfun(@rdivide, bsxfun(@minus, R_full_t_rep, bs(:)'), cs(:)') .^2)); % a1*exp(-((x-b1)/c1)^2)
            corticalStateTerms = corticalStateGausses;                
    end

    opt.useAdditiveConstant = useAdditiveConstant;
    opt.bs = bs;
    opt.cs = cs;
    opt.optimizeB = optimizeB;
    opt.optimizeC = optimizeC;
    opt.nStateTerms = nStateTerms;
    opt.cost_loglikelihood = true;
    opt.cost_sumStimFirst = false;
    opt.cost_sumTrialsFirst = false;
    opt.softPlus = strcmp(opt.rectifyType, 'softplus');
    opt.fixGradients = 0;
    opt_params = opt; % save these options that just contain parameters

    
    opt.R_full = R_full;
    opt.R_full_t = R_full_t;
    opt.corticalStateTerms = corticalStateTerms;
    opt.corticalStates = c_state_all; 
    opt.r_true = R0;
    opt.const = const0;
    opt_forFirstIteration = opt;

    %%
    doGradientChecks = 0 && 1;
    %%
    
    % test cost function - test Ai gradient    
    opt_A = opt; opt_A.optimizeB = false; opt_A.optimizeC = false; params_A = [A0];
    opt_A.r_true = R;
    C2 = costFunction_response_CorticalState(params_A, 'state', opt_A);
    A_cost_func = @(A_in) costFunction_response_CorticalState(A_in, 'state', opt_A);
    if doGradientChecks
        testCostFunctionGradient(A_cost_func, params_A);
    end
    
    % test cost function - test response gradient
    if 0
        opt.corticalStates = c_state_all; opt.const = const0;
        opt_tmp = opt; opt_tmp.fixGradients = true;
        R_cost_func = @(R_in) costFunction_response_CorticalState(R_in, 'response', opt_tmp);
        C1 = R_cost_func(R0);
        if doGradientChecks
            testCostFunctionGradient(R_cost_func, R);
        end
    end
    %%
    
%     assert(isequalToPrecision( C1, C2, 1e-7));

    %%
    if flashedOrDrifting(Gid) == 2
        %%
        [~, ~, ~, nCycles, nRep] = dbGetUniqueOriSpPh('Gid', Gid, 'length');
        nCycles_use = nCycles-1;
        assert(nTrials_orig == nCycles_use * nRep);
        nRep_use = min(nRep, nTrials); % normally use nRep. but in case that only have 1 cycle, just use nTrials
%         if nTrials < nTrials_orig
%             nRep_use = nRep / 2;
%         end
            seg_idxs = segIdxs(nTrials, nRep_use );
    end
    mean_state_eachTrial_cycAv = [];
    mean_r_corrected_eachTrial_cycAv = [];
    coeff_of_var_cycAv = [];
    
    
    delta_Ai = 1;  
    delta_R = 1;   
    th = 1e-4;
    th2 = .01;
    th3 = .05;
    
    allCosts = [];
    const_i = const0;
    c_state_all_i = c_state_all;
    R_i = R0;
    A_i = A0;
    R_i_prev = R_i;
    A_i_prev = A_i;
%         Ac_i = Ac0;
    count = 1;
%         nStateTerms = 4;
    stim_idx = ord(R_full_t(:), 'ascend');

    Ac_a_idx = 1:nStateTerms;
    Ac_c_idx = 1+nStateTerms:nStateTerms*2;
    A_a_idx = 1:nStateTerms;
    if useAdditiveConstant
        Ac_a_idx = Ac_a_idx + 1;
        Ac_c_idx = Ac_c_idx + 1;
        A_a_idx = A_a_idx + 1;
    end

    %%
    maxIter = iff(haveInputGain, 1000, 400);
    fmin_options_A = optimset('GradObj', 'on', 'Display', show_opt_details);
    fmin_options_R = optimset('GradObj', 'on', 'HessPattern', speye( numel(R_i) ), 'Display', show_opt_details, 'MaxIter', maxIter);
%         figure(10); clf; hold on; % ongoing estimate of g(t)
%         plot( R_full_t(stim_idx), c_state_all(stim_idx)/mean(c_state_all(:)), ['k-' ]);
    finished = false;
    %%
    while ~finished
        % 1. optimize Ai  (given a fixed response)
%         optimizeAi_now = ~haveInputGain;
        if ~haveInputGain
            fprintf('Iteration #%2d :  Optimizing Ai ... ', count); tic;

            % 1a. first update cost function for Ai
            opt.r_true = R_i;
            A_cost_f = @(A_in) costFunction_response_CorticalState(A_in, 'state', opt);

            % 1b. optimize
            [A_i, cost1_i, ~, optA_output] = fminunc(A_cost_f, A_i, fmin_options_A); 
            Ai = A_i(A_a_idx);
            c_state_all_i = corticalStateTerms * Ai; 
            

            % Adjust ai so that mean of cortical states == 1;

            if adjustMeanState
                scale_factor = corticalStateFactor(c_state_all_i);

                c_state_all_i = rectifyFunc( c_state_all_i * scale_factor );
                assert( isequalToPrecision (nanmean(c_state_all_i(:)), 1, 1e-5) );
    %             c_state_all_i(c_state_all_i < opt.minCorticalState) = opt.minCorticalState;

                Ai = Ai * scale_factor;
                A_i(A_a_idx) = Ai;
            end

            fprintf(' (it=%d) [%s]', optA_output.iterations, sec2hms(toc));
        else
            cost1_i = nan;
        end

%%
        % 2. optimize R_i  (given a fixed cortical state)
        fprintf(' Optimizing Ri ...'); tic;

        % 2a. first update cost function for R.
        opt.const = const_i;
        opt.corticalStates = c_state_all_i;
        R_cost_func = @(R_in) costFunction_response_CorticalState(R_in, 'response', opt);

        % 2b. optimize
        [R_i, cost2_i, ~, optR_output] = fminunc(R_cost_func, R_i, fmin_options_R); 
        neg_str = iff(any(R_i(:)<0), '!', '');
        R_i(R_i<0) = 0;
        fprintf(' (it=%d) %s [%s].  LL = %.10g;  ', optR_output.iterations, neg_str, sec2hms(toc), cost2_i);

        %%
        allCosts(:, count) = [cost1_i; cost2_i];   
        count = count + 1;


        delta_R = max( abs((R_i(:) - R_i_prev(:))./R_i_prev(:)) );
        delta_Ai = max( abs((A_i - A_i_prev)./A_i_prev) );

        R_i_prev = R_i;
        A_i_prev = A_i;

        converged  = (delta_Ai < th) && (delta_R < th);
        goodEnough = (count > 10 && (delta_Ai < th2) && (delta_R < th2));
        well_ok    = (count > 30 && (delta_Ai < th3) && (delta_R < th3));
        err        = count > 100;

        %%
        corticalStates_stim = reshape(c_state_all_i, [nStim, nTrials]);
        R_full_corrected = R_full ./ corticalStates_stim;
        mean_r_corrected_eachTrial = nanmean(R_full_corrected, 1);
        mean_state_eachTrial = nanmean(corticalStates_stim, 1);
        
        mean_r_eachTrial = nanmean(R_full, 1);
        
        [cc_R_av_gain, p_R_av_gain]         = corr(mean_r_eachTrial(:), mean_state_eachTrial(:));  

        [cc_Rcorr_av_gain, p_Rcorr_av_gain] = corr(mean_r_corrected_eachTrial(:), mean_state_eachTrial(:), 'rows', 'complete');  
        [cc_Rcorr_all_gain, p_Rcorr_all_gain] = corr(R_full_corrected(:), corticalStates_stim(:), 'rows', 'complete');  
        
        idx_trials_use = mean_r_corrected_eachTrial > 0; %< mean_r_corrected_eachTrial (can ignore trials with no response at all)
        coeff_of_var = std(mean_r_corrected_eachTrial(idx_trials_use))/mean(mean_r_corrected_eachTrial(idx_trials_use));
        3;
%         coeff_of_var_all = nanstd(R_full_corrected(:))/nanmean(R_full_corrected(:));
        if flashedOrDrifting(Gid) == 2
            %%
            mean_state_eachTrial_cycAv = cellfun(@(j) mean(mean_state_eachTrial(j)), seg_idxs);
            mean_r_corrected_eachTrial_cycAv = cellfun(@(j) mean(mean_r_corrected_eachTrial(j)), seg_idxs);
            coeff_of_var_cycAv = std(mean_r_corrected_eachTrial_cycAv)/mean(mean_r_corrected_eachTrial_cycAv);
            coeff_of_var_cycAv_str = sprintf(' [ca:%.2f]', coeff_of_var_cycAv);
        else
            coeff_of_var_cycAv_str = '';
        end
        
        
        if isnan(p_Rcorr_av_gain) && isequalToPrecision(abs(cc_Rcorr_av_gain), 1, 1e-4)
            p_Rcorr_av_gain = 0;
        end

%         [cc_Rcorr_nz_gain, p_Rcorr_nz_gain] = corr(R_full_corrected(idx_responded), corticalStates_stim(idx_responded));  
     
%         [hh, slope_pval] = regressionSlopeTtest(mean_r_corrected_eachTrial(:), mean_state_eachTrial(:), .05);

        %%
        %%%%%
        
        finished = converged || goodEnough || well_ok || err || haveInputGain;
        
        done_str = iff(finished, ' - done!', '');
        fprintf('cc = %.1f. p_cc = %.3g. c.o.v. = %.2f%s; dA = %.4f, dR = %.4f %s\n', cc_Rcorr_av_gain, p_Rcorr_av_gain, coeff_of_var, coeff_of_var_cycAv_str, delta_Ai,delta_R, done_str);
        if err
%             error('Too long!')
        end
        3;
        %%%%%%%%%%%%%%%%
        %%
        showProgress = show;
        showAtEnd = show;
       
        if finished;
            3;
        end
%         if finished && ~haveInputGain
%             R_cost_func0 = @(R_in) costFunction_response_CorticalState(R_in, 'response', opt_forFirstIteration);            
%             [R_i_0, cost2_i0, ~, optR_output0] = fminunc(R_cost_func0, R_i_0, fmin_options_R); 
%             cc_R_0 = corr(R_i(:), R_i_0(:));
%             fprintf('[From scratch: cost = %.10f; it=%d. cc_R = %.4f]\n', cost2_i0, optR_output0.iterations, cc_R_0);
%             3;
%         end
        

        if showProgress || (showAtEnd && finished)


            %%
            R_i_rep = R_i(:,ones(1,nTrials));
            corticalStates_stim = reshape(c_state_all_i, [nStim, nTrials]);
            
            rm_gt = R_i_rep .* corticalStates_stim;
            allErrors = R_full - (const_i + rm_gt);
            allErrors_scl = allErrors ./ (sqrt(corticalStates_stim) + .001);
%             allErrors_scl = allErrors ./ (sqrt(rm_gt) + .01);
%             allErrors_scl = allErrors ./ sqrt(corticalStates_stim);
%             allErrors_scl = allErrors; % ./ (sqrt(rm_gt) + .01);
            
            if ~useAdditiveConstant
                assert(const_i == 0);
            end
            R_full_corr_add = const_i + R_i_rep + allErrors_scl;
            R_full_corr = R_full ./ corticalStates_stim; % R_i_rep;
%             R_full_opt = res

            R_full_use = R_full;

            useRawFiringRates = 0;
            if useRawFiringRates
                if ~exist('allHistVals', 'var')
                    [bins, allHistVals] = dbGetCellSpkStimHists(Gid, cellId);
                end
                wind = [0, 100];
                l_bin = indmin(abs(bins-wind(1)));
                r_bin = indmin(abs(bins-wind(2)));
                %
                R_full_raw = mean(allHistVals(l_bin:r_bin,:,:),1);
                R_full_use = R_full_raw;
            end
                
%             R_full_use = R_full_2;
%             R_full_use = R_full;
            figure(11); clf; hold on; box on; % current snapshot of R_full_corr (vs orig)
            t_lims = lims(R_full_t(:), .01);
            binE = linspace(t_lims(1), t_lims(2), 100); dt = diff(binE(1:2));
            binC = binEdge2cent(binE);
            spkRates_orig = histcnt(R_full_t(:), binE, R_full_use(:)) / (dt * nanmean(R_full(:)));
            spkRates_orig_sm = gaussSmooth(spkRates_orig, 10);
            spkRates_corr = histcnt(R_full_t(:), binE, R_full_corr(:)) / (dt * nanmean(R_full(:)));
            spkRates_corr_add = histcnt(R_full_t(:), binE, R_full_corr_add(:)) / (dt * nanmean(R_full(:)));
            h_corr = plot(binC, spkRates_corr, 'rs-');
            h_corr_add = plot(binC, spkRates_corr_add, '^-', 'color', [0 .6 0]);
            h_orig = plot(binC, spkRates_orig, 'bo-', 'linewidth', 1);

%             [relContrOfFrameToSpike, meanStimFiringRate] = getParsedSpikes('frame',  Gid, cellId, [0 100])

            mean_spkRates_orig = mean(spkRates_orig(:));
            mean_spkRates_corr = mean(spkRates_corr_add(:));
            mean_spkRates_corr2 = mean(spkRates_corr(:));
            h_state_orig = plot( R_full_t(stim_idx), c_state_all(  stim_idx)/nanmean(c_state_all(  :))*mean_spkRates_orig, ['b-'], 'linewidth', 2);
            h_state_fit = plot( R_full_t(stim_idx), c_state_all_i(stim_idx)/nanmean(c_state_all_i(:))*mean_spkRates_orig, 'r-', 'linewidth', 2);
%             plot( R_full_t(stim_idx), c_state_all_i(stim_idx)/mean(c_state_all_i(:))*mean_spkRates_corr2, '-', 'linewidth', 2, 'color', [0 .7 0]);

            xlabel('time (s)');
            ylabel('Firing rate (Hz)');
            ylims = ylim;
            xlim(t_lims);
            ylim([ylims(1), ylims(2)*1]);
            legend([h_orig, h_corr, h_corr_add, h_state_orig, h_state_fit], ...
                {'Raw firing rates', 'R_{mn} / g(t)', 'r_m + (R_{mn} - r_m*g(t))',  'Initial cortical state (scl)', 'Final cortical state (scl)'}, 'location', 'best', 'fontsize', 8);
            title( sprintf('Group %d, cell %d. %s', Gid, cellId, oe_str) );


%           
            %%
%                 plot( R_full_t(stim_idx), c_state_all_i(stim_idx), ['-' color_s(count)]);

            %%
            figure(12); clf; hold on;

%                 R_full_use = R_full ;   r_label = 'Observed mean firing rate (Hz)'; showFitEqn = false;  doRMS = false;
%                 R_full_use = R_full_corr_add;  r_label = 'Corrected firing rate (Hz)'; showFitEqn = true; doRMS = false;
                R_full_use = R_full_corr;  r_label = 'Corrected firing rate (Hz)'; showFitEqn = true; doRMS = false;
%                 R_full_use = abs(allErrors).^2;  r_label = 'Errors : abs(observed - predicted) (Hz)';  showFitEqn = true ; doRMS = false;


            c_state_use = c_state_all_i;

            c_state_all_i_stim = reshape(c_state_all_i, [nStim, nTrials]);
            mean_state = nanmean(c_state_all_i_stim, 1);
            xlims = [0, max(mean_state)];

            n_top_stim = min(500, nStim);
            stim_start = 1;
            if doRMS
                mean_r = nanrms(R_full_use, 1);
                top_r = nanrms( R_full_use(idx_bestStim(stim_start+[1:n_top_stim]-1), :), 1);
                bot_r = nanrms( R_full_use(idx_worstStim(stim_start+[1:n_top_stim]-1), :), 1);

            else
                mean_r = nanmean(R_full_use, 1);
                top_r = nanmean( R_full_use(idx_bestStim(stim_start+[1:n_top_stim]-1), :), 1);
                bot_r = nanmean( R_full_use(idx_worstStim(stim_start+[1:n_top_stim]-1), :), 1);                            
            end

            [poly_coef, polyfit_S] = polyfit(mean_state, mean_r, 1);
            [cc, cc_p] = corr(mean_state(:), mean_r(:));
%             [Y,DELTA] = polyconf(p,X,S) takes outputs p and S from polyfit and generates 95% prediction intervals Y ± DELTA for new observations at the values in X.
        
            
            fit_fun = @(x) polyval(poly_coef, x);
            fun_s = sprintf('y = %.1f x + %.1f', poly_coef(1), poly_coef(2));


            idx_state = ord(mean_state, 'ascend');
            h_top = plot(mean_state(idx_state), top_r(idx_state), ['ro-'], 'linewidth', 2);
            box on;
            h_bot = plot(mean_state(idx_state), bot_r(idx_state), ['ko-'], 'linewidth', 2);
            h_av = plot(mean_state(idx_state), mean_r(idx_state), ['bo-'], 'linewidth', 2);
            hold on;
            xlabel('g(t)'); ylabel(r_label);
            
            if flashedOrDrifting(Gid) == 2
                [~, idx_state] = sort(mean_state_eachTrial_cycAv);
                plot(mean_state_eachTrial_cycAv(idx_state), mean_r_corrected_eachTrial_cycAv(idx_state), 'ks-', 'linewidth', 2);
                
%                 mean_state_eachTrial_cycles_C = cellfun(@(j) (mean_state_eachTrial(j)), seg_idxs, 'un', 0);
%                 mean_r_corrected_eachTrial_cycles_C = cellfun(@(j) (mean_r_corrected_eachTrial(j)), seg_idxs, 'un', 0);
%                 for seg_i = 1:length(seg_idxs)
%                     for cyc_j = 1: min( length(seg_idxs{seg_i}), 6);
%                         plot(mean_state_eachTrial_cycles_C{seg_i}(cyc_j), mean_r_corrected_eachTrial_cycles_C{seg_i}(cyc_j), [color_s(cyc_j+1) 's'], 'markersize', 12, 'linewidth', 3)
%                     end
%                 end
            end
               %% 
                
            

            logFactorial = @(N) gammaln(N+1);   % John D'Errico one line solution


            [x_fit, y_fit] = fplot(fit_fun, xlims);
            h_fit = plot(x_fit, y_fit, ['b:'], 'linewidth', 3);
            xlim(xlims);
            title( sprintf('Group %d, cell %d.', Gid, cellId) );
            hold on;
            plot(0,0, 'w')
            drawHorizontalLine(0)
            leg_labels = {sprintf('Top %d stimuli', n_top_stim), 'All Stimuli', sprintf('Bottom %d stimuli', n_top_stim)};
            plot_h = [h_top, h_av, h_bot];
            if showFitEqn
                leg_labels = [leg_labels, fun_s];
                plot_h = [plot_h, h_fit];

            end
%             legend(plot_h, leg_labels, 'location', 'best');
    
            %% figure 
            figure(13); clf; hold on; box on;

            plot(c_state_all_i(:), R_full(:), 'o')
            c_state_all_stim = reshape(c_state_all_i, [nStim, nTrials]);

%             R = mean(R_full, 2);
%             [~, idx_bestStim] = sort(R(:), 'descend');
%             [~, idx_worstStim] = sort(R(:), 'ascend');

            idx_show = idx_bestStim(1);
            c_states_best = c_state_all_stim(idx_show, :);
            r_obs = R_full(idx_show, :);

            r_true = R(idx_show);
            r_pred = const_i + r_true * c_states_best;

%                         diff_sqr = (r_obs  - r_pred).^2;
%                         cost_str = sprintf('Cost = %.1f', round(sum(diff_sqr)));

            ll = sum( -r_true .* c_states_best + r_obs .* log (r_true .* c_states_best) - logFactorial(r_obs) );
            ll_str = sprintf('Log-likelihood = %.1f', ll);

            h_obs = plot(c_states_best, r_obs, 'ro', 'linewidth', 2, 'markersize', 10);
            h_pred = plot(c_states_best, r_pred, 's', 'color', [0 .8 0], 'linewidth', 2, 'markersize', 10);
            legend([h_obs, h_pred], {'Observed', 'Predicted'});
            title(ll_str);

            %%

            figure(14); clf; hold on; box on;
            plot(corticalStates_stim(:), R_full_corrected(:), '.'); xlabel('gain'); ylabel('R-corr');
            plot(corticalStates_stim(:), R_full_corr_add(:), 'r.'); xlabel('gain'); ylabel('R-corr');
            
            title(sprintf('ALL: cc = %.1f., p = %.3g\n average: cc = %.1f., p = %.3g \n. ', cc_Rcorr_all_gain, p_Rcorr_all_gain, cc_Rcorr_av_gain, p_Rcorr_av_gain))
                        
            3;
            %%
            
%                                                 plot(meanR, meanErr, '.');  

            if 0
                %%
                absErrors = abs(allErrors(:));
                meanErr = nanmean( abs(allErrors), 2);
                meanR   = nanmean( R_full_corr, 2);
                
%                 plotVs = 'state';   nStateBins = 10;   nResponseBins = 7;
                plotVs = 'response';  nStateBins = 6;   nResponseBins = 10;
                    nLines = 6;

                dt_sec = median( diff(sort(R_full_t(:))));
                    
%                 response_use = R_i_rep(:); response_name = 'mean response (r_m)';
                response_use = R_i_rep(:) .* corticalStates_stim(:) * dt; response_name = 'expected response (r_m * g(t) )';
%                 r_use = R_i_rep(:);
                
                
                state_lims = lims(corticalStates_stim(:), .01); state_binE = linspace(state_lims(1), state_lims(2), nStateBins+1);
                r_lims = lims(response_use(:), .0); 
                r_lims(2) = 80;                
                r_binE = linspace(r_lims(1), r_lims(2), nResponseBins+1);
                [~, state_binIds] = histcnt(corticalStates_stim(:), state_binE);
                [~, r_binIds] = histcnt(response_use(:), r_binE);
                r_binC = binEdge2cent(r_binE); state_binC = binEdge2cent(state_binE);
                Z_err = zeros(length(r_binC), length(state_binC));
                for i = 1:length(r_binC)
                    r_binM(i) = mean( response_use (r_binIds == i) );
                    for j = 1:length(state_binC)
                        s_binM(j) = nanmean( corticalStates_stim (state_binIds == j) );
                        
                        idx_ij = r_binIds == i & state_binIds == j;
%                         Z_err(i,j) = nanmean( absErrors (idx_ij) ./ R_i_rep(idx_ij) ); ylab = 'mean relative error (err / r_m)'
%                         Z_err(i,j) = nanmean( absErrors (idx_ij) );  ylab = 'mean absolute error'
                        Z_err(i,j) = nanrms( absErrors (idx_ij) ) ;  ylab = 'RMS error';
%                         Z_err(i,j) = nanmean( absErrors(idx_ij).^2 );  ylab = 'variance (error^2)';
                        
%                         Z_err(i,j) = nanmean( absErrors(idx_ij).^2 ./ (R_i_rep(idx_ij).*corticalStates_stim(idx_ij)) ); ylab = 'variance / expected mean [err^2 / r_m * gain]';
%                         Z_err(i,j) = nanmean( R_i_rep(idx_ij) / rms(absErrors (idx_ij))   ); ylab = 'mean / RMS error'
                    end
                end


%                 figure(14); clf; hold on; box on;
%                 hist2dPlot(Z_err', state_binE, r_binE);
%                  ylabel('mean response');
%                  xlabel('cortical state');
%                 zlabel('mean Error');
%                 %%
%                 [s_g, r_g] = meshgrid(state_binC, r_binC);
% %                 plot3(s_g(:), r_g(:), Z_err(:), 'ro')
%                 
%                 idx_use = ~isnan(Z_err(:));
%                 p = fit2dPolySVD(s_g(idx_use), r_g(idx_use), Z_err(idx_use), 2)
%                 
% %                 Z_fit = reshape( eval2dPoly( s_g(:), r_g(:), p), size(Z_err));
%                 x = s_g(:); y = r_g(:);
%                 Z_fit = eval2dPoly( s_g(:), r_g(:), p);
% %                 Z_fit2 = p(1) + p(2)*y + p(3)*x;
%                 plot3(s_g(:), r_g(:), Z_fit(:), 'rs');
  %}              
                
%                 Z_fit = 
                

                figure(15); clf; hold on; box on;
                cols = get(gca, 'colorOrder');
                if strcmp(plotVs, 'state')
                    usePolyFit = 1;
                    Z_errs_plot = Z_err(1:nLines,:)';
                    x = state_binC;
                    xE = state_binE;
                    x_name = 'cortical state';
                    z_short = 'r_m';
                    z = r_binC; 
                    zE = r_binE;
                    zM = r_binM;
                    z_suff = ' Hz';
                    wid = 0;
                elseif strcmp(plotVs, 'response')
                    usePolyFit = 0;
                    Z_errs_plot = Z_err(:, 1:nLines);
                    x = r_binC;
%                     x = sqrt(r_binC);
                    
                    xE = r_binE;
                    z = state_binC;
                    zE = state_binE;
                    zM = s_binM;
                    x_name = response_name;
                    z_short = 's';
                    z_suff = '';
                    wid = 1;
                end
                plot(x, Z_errs_plot, '-o', 'linewidth', 2);
                xlabel(x_name);
                
                ylabel(ylab, 'interpreter', 'none');
                 
                
                 
                 fit_strs = cell(1,nLines);
                 for j = 1:nLines
                     idx_use = ~isnan(Z_errs_plot(:,j));
                     if usePolyFit
                         p = polyfit (x(idx_use), Z_errs_plot(idx_use,j), 1);
                         [x_fit,y_fit] = fplot(@(x) polyval(p, x), xlim);
                         fit_strs{j} = sprintf('y = %.1f x + %.1f', p(1), p(2));
                     else
                         func = @(b, x) b(1) * x.^b(2);
                         b0 = [1 1];
                         b = nlinfit(x(idx_use), Z_errs_plot(idx_use,j), func, [1 1]);
                         [x_fit,y_fit] = fplot(@(x) func(b, x), xlim);
                         fit_strs{j} = sprintf('y = %.1f x ^ %.1f', b(1), b(2));
                     end
                     
                     plot(x_fit,y_fit, ':', 'color', cols(j,:));
                     
                 end
                 
                 if usePolyFit
                      ylims = ylim;
                     ylim([0, ylims(2)*1.3]);
                 else
                     ylims = ylim;
                     ylim([0, ylims(2)*1.5]);
                 end
                 r1 = num2cell(zE(1:nLines));
                 r2 = num2cell(zE(2:nLines+1));
                 rM = num2cell(zM(1:nLines));
                 leg_strs = cellfun(@(rM, r1, r2,s) sprintf('<%s> = %.1f %s, [%.*f, %.*f], %s', z_short, rM, z_suff, wid, r1, wid, r2, s), rM, r1, r2, fit_strs, 'un', 0);
                 legend(leg_strs, 'location', 'NW', 'fontsize', 8, 'interpreter', 'none')

                 
            end
            
        end                        
                    %%


    end            


    if strcmp(trialsDivision, 'all')
        r_corrected = reshape(R_i, [nOri, nSpf, nPh]);
    elseif any(strcmp(trialsDivision, {'oe', 'hoe'}))
        r_corrected = reshape(R_i, [nOri, nSpf, nPh, 2]);        
    end
    
    r_corrected_full = reshapeBack( R_full_corr_add );
    
    %%
    stats.log_likelihood = -allCosts(end);
    
    stats.cc_Rcorr_av_gain = cc_Rcorr_av_gain;
    stats.p_Rcorr_av_gain = p_Rcorr_av_gain;

    if isnan(p_Rcorr_av_gain)
        beep;
        keyboard;
    end
    
    stats.cc_R_av_gain = cc_R_av_gain;
    stats.p_R_av_gain  = p_R_av_gain;
    
    OSP_prev = nanmean(R_full_orig,4);
    OSP_now  = nanmean(r_corrected,4);
    
    OS_prev = nanmean(OSP_prev,3);
    OS_now = nanmean(OSP_now,3);
    
    R_full_corr_reshaped = reshapeBack(R_full_corr);
    [cc_OSPt_before_after, p_OSPt_before_after] = corr(R_full_orig(:), R_full_corr_reshaped(:), 'rows', 'complete');
    [cc_OSP_before_after, p_OSP_before_after] = corr(OSP_prev(:), OSP_now(:), 'rows', 'complete');
    [cc_OS_before_after, p_OS_before_after] = corr(OS_prev(:), OS_now(:), 'rows', 'complete');
    
    stats.cc_OSPt_before_after = cc_OSPt_before_after;
    stats.p_OSPt_before_after = p_OSPt_before_after;
    stats.cc_OSP_before_after = cc_OSP_before_after;
    stats.p_OSP_before_after = p_OSP_before_after;
    stats.cc_OS_before_after = cc_OS_before_after;
    stats.p_OS_before_after = p_OS_before_after;
    
    stats.cc_Rcorr_gain = cc_Rcorr_all_gain; 
    stats.p_Rcorr_gain = p_Rcorr_all_gain;
    stats.opt = opt_params;
    
    stats.coeff_of_var = coeff_of_var;
    stats.r_corr_trial = mean_r_corrected_eachTrial(:);
    stats.state_trial = mean_state_eachTrial(:);
    
    stats.coeff_of_var_cycAv = coeff_of_var_cycAv;
    stats.r_corr_trial_cycAv = mean_r_corrected_eachTrial_cycAv(:);
    stats.state_trial_cycAv = mean_state_eachTrial_cycAv(:);

    
    for j = 1:nStateTerms
        S_gainParams.(sprintf('a%d', j)) = A_i(j);
        S_gainParams.(sprintf('b%d', j)) = bs(j);
        S_gainParams.(sprintf('c%d', j)) = cs(j);
    end
    
    stats.corticalGainParams = S_gainParams;
    %%
    stateFit_S = state_S.opt;
    stateFit_S.min_unrectfitValue = state_S.min_unrectfitValue;
    stateFit_S.fracRectified = state_S.fracRectified;
    
    stats.stateFit = stateFit_S;
%%    
    if nargout >= 3 % get cortical gain
        
        new_state_func = getCorticalStateFunction(state_func, A_i);
        corticalGain = new_state_func(R_full_t_orig);
       
        [R_full_t_chron, idx_origstim2chron] = sort(R_full_t_orig(:));
        assert(isequal(R_full_t_chron,  R_full_t_orig(idx_origstim2chron)))

        R_full_t_back = reshapeBack( R_full_t );
        idx_chron2stim = binarySearch(R_full_t_chron, R_full_t_back(:));
        assert(isequal(R_full_t_chron(idx_chron2stim),  R_full_t_back(:)))
        
        corticalGain_chron = corticalGain(idx_origstim2chron);

        corticalGain = corticalGain_chron;
        %%
%         binarySearch(R_full_t
%         if ~oddEven_differentStimuli
%             assert(isequal(idx_origStimOrder(:)', 1:nStim*nTrials));
%         end
        
%         corticalGain = reshape(c_state_all_i(idx_origStimOrder), size(R_full_orig) );
    
    end
    
    
    showDifferentOSPestimates = show;
    if showDifferentOSPestimates
        %%
        figure(67); clf;
        OSP_prev = mean(R_full_orig,4);
        OSP_now  = mean(r_corrected,4);
        
        OS_prev = mean(OSP_prev,3);
        OS_now = mean(OSP_now,3);
        if isvector(OS_prev)
            if size(OS_prev,1) > 1 % orienations:
                oris = linspace(0, 360, size(OS_prev,1)+1);
                oris = oris(1:end-1);                            
                plot(oris, OS_prev, 'bo-', oris, OS_now, 'rs-');
                xlim([0, 360]);
                set(gca, 'xtick', [0:90:360]);
            else
                spfs = 1:size(OS_prev,2);
                plot(spfs, OS_prev(:), 'bo-', spfs, OS_now(:), 'rs-');
            end
        else
            subplot(1,2,1);
            oris = 0:5:180;
            spfs = 1:10;
            imagesc(spfs, oris, OS_prev); 
            title('Original');
            subplot(1,2,2);
            imagesc(spfs, oris, OS_now);
            title('Corrected');
            
        end
        [~,idx1] = max( OS_prev(:) .* OS_now(:));
        [ori_i, sp_i] = ind2sub(size(OS_prev), idx1);
        ptc_prev = squeeze(OSP_prev(ori_i, sp_i, :));
        ptc_now  = squeeze(OSP_now(ori_i, sp_i, :));
        
        figure(68); subplot(2,1,1);
        plot(OSP_prev(:), OSP_now(:), '.');
        xlabel('Original'); ylabel('Corrected');
        L = max([OSP_prev(:); OSP_now(:)]);
        axis tight; ax = axis; axis([ax(1) ax(2)*1.02, ax(3) ax(4)*1.02])
        title( sprintf('Group %d, cell %d.', Gid, cellId) );
        
        subplot(2,1,2);
        ph = linspace(0, 360, nPh+1); ph = ph(1:end-1);
        if nPh > 10;
            ptc_prev = gaussSmooth(ptc_prev, 2);
            ptc_now = gaussSmooth(ptc_now, 2);
        end
        
        plot(ph, ptc_prev, 'o-', ph, ptc_now, 'rs-');
        xlim([0, 360]); set(gca, 'xtick', [0:90:360]);
        xlabel('Phase');
        M = max([ptc_now; ptc_prev]);
        ylim([0, M+1]);
        
        
        3;
        
    end
  
    3;
end

function R_resized = reshapeToStimVectors(R, size_orig, idx_1, idx_2)
     [nOri, nSpf, nPh, nTrials] = dealV(size_orig);
     nStim = nOri * nSpf *nPh ;
     if ~odd(nTrials)
            R_resized = [reshape(R(:,:,:,idx_1), [nStim, nTrials/2]);
                         reshape(R(:,:,:,idx_2), [nStim, nTrials/2]) ];
                              
     else
         assert( abs(length(idx_1) - length(idx_2)) <= 1);
         npad1 = iff(length(idx_1) < length(idx_2), 1, 0);
         npad2 = iff(length(idx_2) < length(idx_1), 1, 0);
         R1 = [reshape(R(:,:,:,idx_1), [nStim, length(idx_1)]),  nan(nStim,npad1)];
         R2 = [reshape(R(:,:,:,idx_2), [nStim, length(idx_2)]),  nan(nStim,npad2)];
         R_resized = [R1; R2];
                 
     end
     
end



function R_resized = reshapeBackToOrigSize(R, size_orig, idx_1, idx_2)
%%    
    [nOri, nSpf, nPh, nTrials] = dealV(size_orig);
    nStim = nOri * nSpf *nPh ;
    
    size_vec = [nStim*2, length(idx_1)];
    if isvector(R)
        R = reshape(R, size_vec);
    end 
%     if ~odd(nTrials)
%         R_resized = reshape(R, [nOri, nSpf, nPh, nTrials]);
%             
%     else
        R_resized = zeros(size_orig);
        R_resized(:,:,:,idx_1) = reshape( R(1:nStim,          1:length(idx_1)), [nOri, nSpf, nPh, length(idx_1)]);
        R_resized(:,:,:,idx_2) = reshape( R([1:nStim]+nStim,  1:length(idx_2)), [nOri, nSpf, nPh, length(idx_2)]);
%     end
    
end

function scale_factor = getCorticalScaleFactor(c_states, rectFunc)

%     scale_factor = 1;
%     getFactor = @(c_states_in) meanRate / mean(c_states_in(:) .* R_full(:));
%     getFactor = @(c_states_in) 1 / nanmean(c_states_in(:)); % .* R_full(:));
%     getFactor = @(c_states_in) 1 / mean(R_full(:) ./ c_states_in(:)); % .* R_full(:));
%     scale_factor_i = 1000;
    if 0
        %%
        allScl = [.001:.001:.1];
        m = zeros(size(allScl));
        for i = 1:length(allScl)
            m(i) = nanmean(rectFunc(c_states(:) .* allScl(i)));
        end
        
    end

    3;
    scale_factor = fzero(@(scl) nanmean(rectFunc(c_states(:) .* scl)) - 1, [1e-10, 10]);
    
%     while abs(scale_factor_i - 1) > 1e-5
% 
%         scale_factor_i = getFactor(c_states);
%         scale_factor = scale_factor * scale_factor_i;
% 
%         c_states = c_states * scale_factor_i;
%         c_states = rectFunc(c_states);
% %         c_states(c_states < minCorticalState) = minCorticalState;
%     end
    
end

function f = getCorticalStateFunction(state_func, newA)

    [a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4,a5,b5,c5,a6,b6,c6,a7,b7,c7,a8,b8,c8] = dealV(coeffvalues(state_func)); %#ok<ASGLU>
    [A1,A2,A3,A4,A5,A6,A7,A8] = dealV(newA);
    
    f_sin_func = @(A1,b1,c1,A2,b2,c2,A3,b3,c3,A4,b4,c4,A5,b5,c5,A6,b6,c6,A7,b7,c7,A8,b8,c8,x) ...
        (A1.*sin(b1.*x+c1) + A2.*sin(b2.*x+c2) + A3.*sin(b3.*x+c3) + A4.*sin(b4.*x+c4) + A5.*sin(b5.*x+c5) + A6.*sin(b6.*x+c6) + A7.*sin(b7.*x+c7) + A8.*sin(b8.*x+c8));
        
    f = @(x) f_sin_func(A1,b1,c1, A2,b2,c2, A3,b3,c3, A4,b4,c4, A5,b5,c5, A6,b6,c6, A7,b7,c7, A8,b8,c8, x);
        

end


function y = rectifyToValue(x, rectValue, rectType, softPlusScale)
    if strcmp(rectType, 'hard')
        y = x;
        y(x < rectValue) = rectValue;
    elseif strcmp(rectType, 'softplus')
        y = softplus(x, rectValue, softPlusScale);
%         y = (log(1+exp(  ((x - rectValue)*softPlusScale)    ))) /softPlusScale + rectValue;
    end


end


%                 

%{

        figure(cell_i); clf; hold on;
        [poly_coef, S2] = polyfit(mean_state, mean_r, 1);
        line1 = @(beta, x) beta(1)*x;        
        [m, r1] = nlinfit(mean_state, mean_r, line1, poly_coef(1));    
        
        S1.normr = sqrt( sum(r1.^2) );
        S1.df = length(mean_r)-1;
        pval = nestedFtest(S1, S2);
        
%         if pval < .05 % likely that model 2 is correct
            fit_fun = @(x) polyval(poly_coef, x);
            fun_s = sprintf('y = %.1f x + %.1f', poly_coef(1), poly_coef(2));
%         else
%             fit_fun = @(x) line1(m, x);
%             fun_s = sprintf('y = %.1f x', m);
%         end

%}


%{

    %%      
    do1StepOptimize = false;
    if do1StepOptimize
        %%
        nParamsTot = length(R_A0);
        HessPattern_combined = eye(nParamsTot); idx_stateParams = nParamsTot - (nStateTerms + useAdditiveConstant)+1 : nParamsTot;
        HessPattern_combined(idx_stateParams, idx_stateParams) = 1;
        options_combined = optimset('GradObj', 'on', 'HessPattern', HessPattern_combined, 'Display', show_opt_details);
        tic;

        opt_RA = opt_A; params_RA = [R_A0];
        RA_cost_func = @(A_in) costFunction_response_CorticalState(A_in, 'response_state', opt_A);
        C3 = costFunction_response_CorticalState(params_RA, 'response_state', opt_RA);
        if doGradientChecks
            testCostFunctionGradient(RA_cost_func, params_RA);
        end
        [RA_best, cost_1step] = fminunc(RA_cost_func, R_A0, options_combined); 
        toc;
    end        
    %%



                if do1StepOptimize 
                    %%

                    RA_cost_func = @(RA_in) costFunction_response_CorticalState(RA_in, 'response_state', opt);
%                     C3 = costFunction_response_CorticalState([R_i(:); Ai(:)], 'response_state', opt_RA);
                    [RA_best_i, cost_1step_i] = fminunc(RA_cost_func, [R_i(:); Ai(:)], options_combined);


                    nStim = nOri*nSpf*nPh;
                    R_best = reshape(RA_best_i(1:nStim), [nOri, nSpf, nPh]);
                    Ai_best = RA_best_i(nStim+1+useAdditiveConstant:end);
                    const_i_best = 0;

                    c_state_all_best = corticalStateTerms * Ai_best(:);
                    scale_factor_best = corticalStateFactor(c_state_all_best);

%                     Ai_best = Ai_best * scale_factor_best;
                    c_state_all_best = c_state_all_best * scale_factor_best;

                    R_i_rep_best = R_best(:,:,:,ones(1,nTrials));

                    corticalStates_stim_best = reshape(c_state_all_best, size(R_full));
                    allErrors_best = R_full - (const_i_best + R_i_rep_best .* corticalStates_stim_best);

                    if ~useAdditiveConstant
                        assert(const_i_best == 0);
                    end
                    R_full_corr_best = R_i_rep_best + allErrors_best;

                end

                if do1StepOptimize
                    spkRates_corr_best = histcnt(R_full_t(:), binE, R_full_corr_best(:)) / (dt * mean(R_full(:)));
                    plot(binC, spkRates_corr_best, 'g^:');
                    mean_spkRates_corr_best = mean(spkRates_corr_best(:));
                    plot( R_full_t(stim_idx), c_state_all_best(stim_idx)/mean(c_state_all_best(:))*mean_spkRates_corr_best, 'g:', 'linewidth', 2);
                end                

%}


%{

    % test "response" cost function 
    opt_tmp = opt; opt_tmp.fixGradients = true;
    R_cost_func = @(R_in) costFunction_response_given_corticalState(R_in, const0, c_state_all, opt_tmp);
    C1 = R_cost_func( R0 );

    if doGradientChecks
        testCostFunctionGradient(R_cost_func, R);
    end


    %% test "cortical state parameters" cost function 
    
    tmp_const = iff(useAdditiveConstant, const0, []);

    % test "cortical state parameters" cost function - just A
    opt_A = opt; opt_A.optimizeB = false; opt_A.optimizeC = false; params_A = [tmp_const; As];
    [C2, c2g] = costFunction_corticalState_given_response( params_A, R, opt_A);
    A_cost_func = @(A_in) costFunction_corticalState_given_response(A_in, R, opt_A);
    if doGradientChecks
        testCostFunctionGradient(A_cost_func, params_A);
    end

    
    if false && 0 %~opt.cost_sumStimFirst
        % test "cortical state parameters" cost function - A and c
        opt_Ac = opt; opt_Ac.optimizeB = false; opt_Ac.optimizeC = true; params_Ac = [tmp_const; As; cs];
        C3 = costFunction_corticalState_given_response( params_Ac, R, opt_Ac);
        Ac_cost_func = @(A_in) costFunction_corticalState_given_response(A_in, R, opt_Ac);
        testCostFunctionGradient(Ac_cost_func, params_Ac);

        % test "cortical state parameters" cost function - A and b and c
        opt_Abc = opt; opt_Abc.optimizeB = true; opt_Abc.optimizeC = true; params_Abc = [tmp_const; As; bs; cs];
        C4 = costFunction_corticalState_given_response( params_Abc, R, opt_Abc);
        Ac_cost_func = @(A_in) costFunction_corticalState_given_response(A_in, R, opt_Abc);
        testCostFunctionGradient(Ac_cost_func, params_Abc);

        assert(isequal(C1, C2, C3, C4))
    else
        assert( all(abs(diff([C1, C1b, C2, C2b]))<1e-7) )
    end
    3;

    % test "complete" cost function
%         opt_RA = opt; opt_A.optimizeB = false; opt_A.optimizeC = false; params_RA = [tmp_const; As; R0];
%         C2 = costFunction_corticalState_given_response( params_A, R, opt_A);
%         RA_cost_func = @(A_in) costFunction_response_CorticalStateA corticalState_given_response(A_in, R, opt_A);
%         testCostFunctionGradient(A_cost_func, params_A);


%}


        %{
        if finished
%             opt.corticalStates = c_state_all_i;
            R_full_v = R_full(:);
            opt_full = opt;
            opt_full.R_full = R_full_v;
            opt_full.fixGradients = 0;
            
            Rf_cost_func = @(Rf_in) costFunction_response_CorticalState(Rf_in, 'response', opt_full);
            C1 = Rf_cost_func(R_full_v);
            if doGradientChecks
%                 testCostFunctionGradient(Rf_cost_func, R_full_v);
            end
    
            % 2b. optimize
            fmin_options_Rfull = optimset('GradObj', 'on', 'HessPattern', speye( numel(R_full_v) ) , 'Display', show_opt_details  );
            [R_full_v_opt, cost2_i] = fminunc(Rf_cost_func, R_full_v, fmin_options_Rfull); 
            3;
        end
        %}