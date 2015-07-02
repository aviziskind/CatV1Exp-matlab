function plotFiringRateVsCorticalState(Gid, cellIds)

%     Gid = 2288;
%     cellId = 3;
%     stateFitType = 'gauss';
    stateFitType = 'sin';

    adjustMeanState = 1;
        minCorticalState = 0.05;
    show_opt_details = 'off'; % default: 'final';
    
    oddEven_differentStimuli = false;
    
    R_full_t = double( dbGetStimulusTimes(Gid) );
    [nOri, nSpf, nPh, nTrials] = size(R_full_t);
    state_func = getCorticalState(Gid, cellIds);
    c_state_all = reshape(feval(state_func, R_full_t), size(R_full_t));
    nT = size(R_full_t, 4);
    nStim = nOri*nSpf*nPh;
    
    
    nStateTerms = length(coeffvalues(state_func))/3;
    
    if nargin < 2 || isempty(cellIds)
        sd = siteDataFor('Gid', Gid, 1);
        cellIds = sd.cellIds;
        cellIds = cellIds(cellIds> 0);
    end
    
    mean_state = zeros(1,nT);
    for i = 1:nT
        ti = c_state_all(:,:,:,i);
        mean_state(i) = mean(ti(:));                
    end
    xlims = [0, max(mean_state)];
%%
   
    for cell_i = 1:length(cellIds)
        
        cellId = cellIds(cell_i);
    
        s = calculatePSTH_STAs_OSP_ForOneCell(Gid, cellId);
        R_full = double( decompress( s.OSP.R_full ) );
        R_full_ref = double( decompress( s.OSP.R_full ) );
        meanRate = mean(R_full(:));
        

        mean_r = zeros(1,nT);
        for i = 1:nT
            resp = R_full(:,:,:,i);
            mean_r(i) = mean(resp(:));
        end

        R = mean(R_full, 4);
        [~, idx_bestStim] = sort(R(:), 'descend');
        [~, idx_worstStim] = sort(R(:), 'ascend');
        
        %%
        
        n_top_stim = 500;
        stim_start = 1;
        [ori_i_top, sp_i_top, ph_i_top] = deal(zeros(1, n_top_stim));
        [ori_i_bot, sp_i_bot, ph_i_bot] = deal(zeros(1, n_top_stim));
        for j = 1:n_top_stim
            [ori_i_top(j), sp_i_top(j), ph_i_top(j)] = ind2sub(size(R), idx_bestStim(j+stim_start-1));
            [ori_i_bot(j), sp_i_bot(j), ph_i_bot(j)] = ind2sub(size(R), idx_worstStim(j+stim_start-1));
        end
        
        [top_r, bot_r] = deal(zeros(1,nT));
        for i = 1:nT
            resp_top = zeros(1, n_top_stim);
            resp_bot = zeros(1, n_top_stim);
            for j = 1:n_top_stim
                resp_top(j) = R_full(ori_i_top(j), sp_i_top(j), ph_i_top(j), i);
                resp_bot(j) = R_full(ori_i_bot(j), sp_i_bot(j), ph_i_bot(j), i);
            end
            top_r(i) = mean(resp_top(:));
            bot_r(i) = mean(resp_bot(:));
        end
        
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
   
        idx = ord(mean_state, 'ascend');
        plot(mean_state(idx), top_r(idx), ['ro-']); 
        box on;
        hold on;
        plot(mean_state(idx), mean_r(idx), [color_s(cell_i) 'o-']); 
        plot(mean_state(idx), bot_r(idx), ['ko-']); 
        xlabel('g(t)'); ylabel('Mean firing rate');
    

        fplot(fit_fun, xlims, [color_s(cell_i) ':']);
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
        
        
        %%
        
       
        
%         linFunc = @(beta, x) beta(1) + beta(2:end) .* x(:);
        % fix g(t) . maximize a_i;
%         beta_fit = nlinfit(c_state_all(:), R_full(:), linFunc, beta0);
        
                                     
                             
        %%
        corticalStateFactor = @(c_states_in) getCorticalScaleFactor(c_states_in, R_full, meanRate, minCorticalState);
        
%         b1 = state_func.b1;
%         state_func_copy = state_func
        
%         R_full_t
        S = cell2struct(num2cell(coeffvalues(state_func))', coeffnames(state_func), 1);
        
        As = arrayfun(@(i) S.(sprintf('a%d', i)), 1:nStateTerms)';
        bs = arrayfun(@(i) S.(sprintf('b%d', i)), 1:nStateTerms)';
        cs = arrayfun(@(i) S.(sprintf('c%d', i)), 1:nStateTerms)';
        useAdditiveConstant = 1 && 0;
        optimizePhase = false;
        optimizeFreqs = false;
        
        if adjustMeanState
            init_scale_factor = corticalStateFactor(c_state_all);
            
            c_state_all = c_state_all * init_scale_factor;
            c_state_all(c_state_all < minCorticalState) = minCorticalState;
            scale_f2 = corticalStateFactor(c_state_all);
            assert( abs(scale_f2 - 1) < 1e-5 )
            
            As = As * init_scale_factor;    
        end
        
        A0 = [As];
        
        if optimizeFreqs
            A0 = [A0; bs];
        end
        if optimizePhase
            A0 = [A0; cs];
        end
        if useAdditiveConstant
            A0 = [const0; A0]; %#ok<*AGROW>
        else
            const0 = 0;
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
        opt.R_full = R_full;
        opt.corticalStateTerms = corticalStateTerms;
        opt.R_full_t = R_full_t;
        opt.bs = bs;
        opt.cs = cs;
        opt.optimizeB = optimizeFreqs;
        opt.optimizeC = optimizePhase;
        opt.nStateTerms = nStateTerms;
        opt.cost_loglikelihood = true;
        opt.cost_sumStimFirst = false;
        opt.cost_sumTrialsFirst = false;
        opt.minCorticalState = minCorticalState;
        opt.fixGradients = 0;
        
        %%
        doGradientChecks = 1 && 0;
        % test "response" cost function 
        
        R_cost_func = @(R_in) costFunction_response_given_corticalState(R_in, const0, c_state_all, opt);
        C1 = R_cost_func( R0 );
        
        if doGradientChecks
            testCostFunctionGradient(R_cost_func, R);
        end
        %%
        % test 'general' cost function
        opt.corticalStates = c_state_all; opt.const = const0;
        R_cost_func2 = @(R_in) costFunction_response_CorticalState(R_in, 'response', opt);
        C1b = R_cost_func2(R0);
        if doGradientChecks
            testCostFunctionGradient(R_cost_func2, R);
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
%%
        opt_A = opt; opt_A.optimizeB = false; opt_A.optimizeC = false; params_A = [tmp_const; As];
        opt_A.r_true = R;
        C2b = costFunction_response_CorticalState(params_A, 'state', opt_A);
        A_cost_func2 = @(A_in) costFunction_response_CorticalState(A_in, 'state', opt_A);
        if doGradientChecks
            testCostFunctionGradient(A_cost_func2, params_A);
        end
3;
        %%
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
        
        %%
        
        
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
        
        options = optimset('GradObj', 'on', 'Display', show_opt_details);
        delta_Ai = 1;  
        delta_R = 1;   
        th = 1e-4;
        th2 = .01;
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
%         figure(10); clf; hold on; % ongoing estimate of g(t)
%         plot( R_full_t(stim_idx), c_state_all(stim_idx)/mean(c_state_all(:)), ['k-' ]);
        finished = false;
        while ~finished
            % 2. optimize Ai  (given a fixed response)
                % 2a. first update cost function for Ai
                fprintf('Iteration #%d :  Optimizing Ai ... ', count);
                tic;
                opt.r_true = R_i;
                
                if optimizePhase
                    
                    Ac_cost_f = @(A_in) costFunction_response_CorticalState(Ac_in, 'state', opt);

                    % 2b. optimize
                    [Ac_i, cost1_i] = fminunc(Ac_cost_f, Ac_i, options); 
                    if useAdditiveConstant
                        const_i = Ac_i(1);
                    end
                    Ai = Ac_i(Ac_a_idx);
                    ci = Ac_i(Ac_c_idx);
                    
                    bxplusC = bsxfun(@plus, bsxfun(@times, bs(:)', repmat(R_full_t(:), [1, nStateTerms])), ci(:)') ;
                    corticalStateTerms = sin(bxplusC);
                    c_state_all_i = corticalStateTerms * Ai(:);
                    
                else
                    
%                     A_cost_f = @(A_in) costFunction_corticalState_given_response(A_in, R_i, opt_A);
                    
                    A_cost_f = @(A_in) costFunction_response_CorticalState(A_in, 'state', opt);

                    % 2b. optimize
                   
                    [A_i, cost1_i] = fminunc(A_cost_f, A_i, options); 
                    
                    if useAdditiveConstant
                        const_i = A_i(1);
                    end
                    Ai = A_i(A_a_idx);
                    
                    c_state_all_i = corticalStateTerms * Ai; 
                end
                
                % Adjust ai so that mean of cortical states == 1;
                
                if adjustMeanState
                  
                    scale_factor = corticalStateFactor(c_state_all_i);
                    c_state_all_i = c_state_all_i * scale_factor;
                    c_state_all_i(c_state_all_i < opt.minCorticalState) = opt.minCorticalState;
                    
                    Ai = Ai * scale_factor;
                    A_i(A_a_idx) = Ai;
                    Ac_i(Ac_a_idx) = Ai;
                end
                
                fprintf(' [%s]', sec2hms(toc));
                
                
%                 c_state_all_i = corticalStateSines * Ai; 
%                 tic;
%                 [Aic_i3, cost2_i3] = fmincg(Ai_cost_f, Aic_i, options); 
%                 toc;
        
%%
                fprintf(' Optimizing Ri ...'); tic;

            % 1. optimize R_i  (given a fixed cortical state)
                % 1a. first update cost function for R.
%                 Ai = Aic_i(2:end);
%                 const_i = Aic_i(1);
                
                opt.const = const_i;
                opt.corticalStates = c_state_all_i;
%                 R_cost_func = @(R_in) costFunction_response_given_corticalState(R_in, const_i, c_state_all_i, opt);
                R_cost_func = @(R_in) costFunction_response_CorticalState(R_in, 'response', opt);
                
                % 1b. optimize
%                 options2 = op

                nParams = numel(R_i);
                options_R = optimset('GradObj', 'on', 'HessPattern', eye( nParams), 'Display', show_opt_details  );
                tic;
                [R_i, cost2_i] = fminunc(R_cost_func, R_i, options_R); 
                
                fprintf(' [%s].  LL = %.5g \n', sec2hms(toc), cost2_i);
                
%                 tic;
%                 [R_i2, cost1_i] = fmincg(R_cost_f, R_i, options2); 
%                 toc;
                allCosts(:, count) = [cost1_i; cost2_i];   
                count = count + 1;

                
                delta_R = max( abs((R_i(:) - R_i_prev(:))./R_i_prev(:)) );
                delta_Ai = max( abs((A_i - A_i_prev)./A_i_prev) );
                              
                R_i_prev = R_i;
                A_i_prev = A_i;
                
                converged = (delta_Ai < th) && (delta_R < th);
                goodEnough = (count > 4 && (delta_Ai < th2) && (delta_R < th2));

                finished = converged || goodEnough;
                3;
                %%%%%%%%%%%%%%%%
                %%
                showProgress = 1;
                
                if showProgress

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

                    %%
                    R_i_rep = R_i(:,:,:,ones(1,nTrials));

                    corticalStates_stim = reshape(c_state_all_i, size(R_full));
                    allErrors = R_full - (const_i + R_i_rep .* corticalStates_stim);

                    if ~useAdditiveConstant
                        assert(const_i == 0);
                    end
                    R_full_corr_add = const_i + R_i_rep + allErrors;
                    R_full_corr = R_full ./ corticalStates_stim; % R_i_rep;


    %                 figure(10);
    %                 clf; hold on; % ongoing estimate of g(t)
    %                 plot( R_full_t(stim_idx), c_state_all(stim_idx)/mean(c_state_all(:)), ['k-' ]);
    %                 plot( R_full_t(stim_idx), c_state_all_i(stim_idx), ['-' color_s(count)]);
    %                 if do1StepOptimize
    %                     plot( R_full_t(stim_idx), c_state_all_best(stim_idx), ['g:']);
    %                 end

%                     figure(10); clf; hold on;
%                     plot(c_state_all(:), R_full(:), 'o')
%                     
%                     [ori_i_top1, sp_i_top1, ph_i_top1] = ind2sub(size(R), idx_bestStim(3));
%                     c_states_best = squeeze(c_state_all(ori_i_top1, sp_i_top1, ph_i_top1, :));
%                     r_best = squeeze(R_full(ori_i_top1, sp_i_top1, ph_i_top1, :));
%                     idx = ord(c_states_best);
%                     plot(c_states_best(idx), r_best(idx), 'ro', 'linewidth', 2, 'markersize', 10);


                    figure(11); clf; hold on; box on; % current snapshot of R_full_corr (vs orig)

                    t_lims = lims(R_full_t(:), .01);
                    binE = linspace(t_lims(1), t_lims(2), 100); dt = diff(binE(1:2));
                    binC = binEdge2cent(binE);
                    spkRates_orig = histcnt(R_full_t(:), binE, R_full(:)) / (dt * mean(R_full(:)));
                    spkRates_corr = histcnt(R_full_t(:), binE, R_full_corr(:)) / (dt * mean(R_full(:)));
                    spkRates_corr_add = histcnt(R_full_t(:), binE, R_full_corr_add(:)) / (dt * mean(R_full(:)));
                    h_orig = plot(binC, spkRates_orig, 'bo-');
                    h_corr = plot(binC, spkRates_corr, 'rs-');
                    h_corr_add = plot(binC, spkRates_corr_add, '^-', 'color', [0 .6 0]);
                    
                    

                    mean_spkRates_orig = mean(spkRates_orig(:));
                    mean_spkRates_corr = mean(spkRates_corr_add(:));
                    mean_spkRates_corr2 = mean(spkRates_corr(:));
                    h_state_orig = plot( R_full_t(stim_idx), c_state_all_i(stim_idx)/mean(c_state_all_i(:))*mean_spkRates_orig, 'r-', 'linewidth', 2);
                    h_state_fit = plot( R_full_t(stim_idx), c_state_all(  stim_idx)/mean(c_state_all(  :))*mean_spkRates_orig, ['b-'], 'linewidth', 2);
%                     plot( R_full_t(stim_idx), c_state_all_i(stim_idx)/mean(c_state_all_i(:))*mean_spkRates_corr2, '-', 'linewidth', 2, 'color', [0 .7 0]);

                    xlabel('time (s)');
                    ylabel('Firing rate (Hz)');
                    ylims = ylim;
                    xlim(t_lims);
                    ylim([ylims(1), ylims(2)*1]);
                    legend([h_orig, h_corr, h_corr_add, h_state_orig, h_state_fit], ...
                        {'Raw firing rates', 'R_{mn} / g(t)', 'r_m + (R_{mn} - r_m*g(t))',  'Initial cortical state (scl)', 'Final cortical state (scl)'}, 'location', 'best', 'fontsize', 8);
                    title( sprintf('Group %d, cell %d.', Gid, cellId) );

                    if do1StepOptimize
                        spkRates_corr_best = histcnt(R_full_t(:), binE, R_full_corr_best(:)) / (dt * mean(R_full(:)));
                        plot(binC, spkRates_corr_best, 'g^:');
                        mean_spkRates_corr_best = mean(spkRates_corr_best(:));
                        plot( R_full_t(stim_idx), c_state_all_best(stim_idx)/mean(c_state_all_best(:))*mean_spkRates_corr_best, 'g:', 'linewidth', 2);
                    end                

                    %%
    %                 plot( R_full_t(stim_idx), c_state_all_i(stim_idx), ['-' color_s(count)]);

                    %%
                    figure(12); clf; hold on;


%                             R_full_use = R_full ;   r_label = 'Observed mean firing rate (Hz)'; showFitEqn = false;
%                             R_full_use = R_full_corr_add;  r_label = 'Corrected firing rate (Hz)'; showFitEqn = true; 
                            R_full_use = R_full_corr;  r_label = 'Corrected firing rate (Hz)'; showFitEqn = true; 
%                             R_full_use = abs(allErrors);  r_label = 'Errors : abs(observed - predicted) (Hz)';  showFitEqn = true ;
% 
                            
                            c_state_use = c_state_all_i;

                            c_state_all_i_stim = reshape(c_state_all_i, size(R_full));
                            mean_state = zeros(1,nT);
                            for i = 1:nT
                                ti = c_state_all_i_stim(:,:,:,i);
                                mean_state(i) = mean(ti(:));
                            end
                            xlims = [0, max(mean_state)];


                            mean_r = zeros(1,nT);
                            for i = 1:nT
                                resp = R_full_use(:,:,:,i);
                                mean_r(i) = mean(resp(:));
                            end

                            n_top_stim = 500;
                            stim_start = 1;
                            [ori_i_top, sp_i_top, ph_i_top] = deal(zeros(1, n_top_stim));
                            [ori_i_bot, sp_i_bot, ph_i_bot] = deal(zeros(1, n_top_stim));
                            for j = 1:n_top_stim
                                [ori_i_top(j), sp_i_top(j), ph_i_top(j)] = ind2sub(size(R), idx_bestStim(j+stim_start-1));
                                [ori_i_bot(j), sp_i_bot(j), ph_i_bot(j)] = ind2sub(size(R), idx_worstStim(j+stim_start-1));
                            end


                            [top_r, bot_r] = deal(zeros(1,nT));
                            for i = 1:nT
                                resp_top = zeros(1, n_top_stim);
                                resp_bot = zeros(1, n_top_stim);
                                for j = 1:n_top_stim
                                    resp_top(j) = R_full_use(ori_i_top(j), sp_i_top(j), ph_i_top(j), i);
                                    resp_bot(j) = R_full_use(ori_i_bot(j), sp_i_bot(j), ph_i_bot(j), i);
                                end
                                top_r(i) = mean(resp_top(:));
                                bot_r(i) = mean(resp_bot(:));
                            end

                            [poly_coef, S2] = polyfit(mean_state, mean_r, 1);
                            line1 = @(beta, x) beta(1)*x;
                            [m, r1] = nlinfit(mean_state, mean_r, line1, poly_coef(1));

    %                         S1.normr = sqrt( sum(r1.^2) );
    %                         S1.df = length(mean_r)-1;
    %                         pval = nestedFtest(S1, S2);
    % 
    %                         if pval < .05 % likely that model 2 is correct
                                fit_fun = @(x) polyval(poly_coef, x);
                                fun_s = sprintf('y = %.1f x + %.1f', poly_coef(1), poly_coef(2));
%                                 fun_s = '';
                                
    %                         else
    %                             fit_fun = @(x) line1(m, x);
    %                             fun_s = sprintf('y = %.1f x', m);
    %                         end


                            idx_state = ord(mean_state, 'ascend');
                            h_av = plot(mean_state(idx_state), mean_r(idx_state), [color_s(cell_i) 'o-'], 'linewidth', 2);
                            box on;
                            hold on;
                            h_top = plot(mean_state(idx_state), top_r(idx_state), ['ro-'], 'linewidth', 2);
                            h_bot = plot(mean_state(idx_state), bot_r(idx_state), ['ko-'], 'linewidth', 2);
                            xlabel('g(t)'); ylabel(r_label);

                            logFactorial = @(N) gammaln(N+1);   % John D'Errico one line solution
                            
                            [x_fit, y_fit] = fplot(fit_fun, xlims);
                            h_fit = plot(x_fit, y_fit, [color_s(cell_i) ':'], 'linewidth', 3);
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
                            legend(plot_h, leg_labels, 'location', 'best');

                            %% figure 
                            figure(13); clf; hold on; box on;

                            plot(c_state_all_i(:), R_full(:), 'o')
                            c_state_all_stim = reshape(c_state_all_i, size(R_full));

                            [ori_i_top1, sp_i_top1, ph_i_top1] = ind2sub(size(R), idx_bestStim(3));
                            c_states_best = squeeze(c_state_all_stim(ori_i_top1, sp_i_top1, ph_i_top1, :));
                            r_obs = squeeze(R_full(ori_i_top1, sp_i_top1, ph_i_top1, :));

                             idx = ord(c_states_best);

                            
                            r_true = R(ori_i_top1, sp_i_top1, ph_i_top1);
                            r_pred = const_i + r_true * c_states_best;

                            diff_sqr = (r_obs  - r_pred).^2;
                            cost_str = sprintf('Cost = %.1f', round(sum(diff_sqr)));
                            
                            ll = sum( -r_true .* c_states_best + r_obs .* log (r_true .* c_states_best) - logFactorial(r_obs) );
                            ll_str = sprintf('Log-likelihood = %.1f', ll);
                            
                            h_obs = plot(c_states_best(idx), r_obs(idx), 'ro', 'linewidth', 2, 'markersize', 10);
                            h_pred = plot(c_states_best(idx), r_pred(idx), 's', 'color', [0 .8 0], 'linewidth', 2, 'markersize', 10);
                            legend([h_obs, h_pred], {'Observed', 'Predicted'});
                            title(ll_str);

                end                        
                        %%
                        
        
                
                3;
                
        end            
        
        
        
    end
    
    
    keyboard;
    3;


    %  Set options for fminunc

        
        

    
    
    3;
end


function scale_factor = getCorticalScaleFactor(c_states, R_full, meanRate, minCorticalState)

    scale_factor = 1;
    getFactor = @(c_states_in) meanRate / mean(c_states_in(:) .* R_full(:));
    scale_factor_i = 1000;
    
    while abs(scale_factor_i - 1) > 1e-3

        scale_factor_i = getFactor(c_states);
        scale_factor = scale_factor * scale_factor_i;

        c_states = c_states * scale_factor_i;
        c_states(c_states < minCorticalState) = minCorticalState;
    end
    
end
%                 