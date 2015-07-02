function [J, grad] = costFunction_response_CorticalState(params, params_name, opt)
    % calculate the cost function & gradient for a particular set of cortical state
    % parameters, given a fixed response profile
    useAdditiveConstant = opt.useAdditiveConstant;
    R_full = opt.R_full;
    R_full_t = opt.R_full_t;
    
%     assert(ismatrix(R_full));
%     assert(ismatrix(R_full_t));

    [nStim, nTrials] = size(R_full);
%     [nOri, nSpf, nPh, nTrials] = size(R_full);      
            
    params = params(:);
    nParams = numel(params);
    nSines = opt.nStateTerms;
    softPlusRect = opt.softPlus;
    
    
    switch params_name
        case 'response',       calcResponseGrad = true;  calcStateGrad = false; 
        case 'state',          calcResponseGrad = false; calcStateGrad = true;  
        case 'response_state', calcResponseGrad = true;  calcStateGrad = true;
        otherwise error('Specify "response", "state", or "response_state"')
    end
    
    if calcResponseGrad
        % cortical state parameters are part of 'params' 
        idx_stim = 1:nStim;
        r_true = params(idx_stim);
        
    elseif ~calcResponseGrad
        % cortical state parameters are given (in opt)
        idx_stim = 0;
        r_true = opt.r_true;
        
    end
    
    r_true = r_true(:);
    r_true_rep = r_true(:,ones(1,nTrials));

    
    if calcStateGrad
        % cortical state parameters are part of 'params'
        
        idx_a = idx_stim(end) + [1 : nSines];
        idx_b = idx_a + nSines;
        idx_c = idx_a + nSines * (opt.optimizeB + opt.optimizeC);
        
        %%
        if useAdditiveConstant
            assert(~opt.cost_loglikelihood)
            idx_const = idx_a(1);
            
            c = params(idx_const);
            idx_a = idx_a + 1;
            idx_b = idx_b + 1;
            idx_c = idx_c + 1;
        else
            c = 0;
        end    

        Ai = params(idx_a); % 2-5
        assert(length(Ai) == nSines);
        
        if opt.optimizeB || opt.optimizeC
            if opt.optimizeB
                bi = params(idx_b); % 6-9
            else
                bi = opt.bs;
            end

            if opt.optimizeC
                ci = params(idx_c); % 10-13
            else
                ci = opt.cs;
            end
            bxplusC = bsxfun(@plus, bsxfun(@times, bi(:)', R_full_t(:)), ci(:)') ;

            corticalStateSines = sin(bxplusC);
            corticalStateCosines = cos(bxplusC);

        else
            corticalStateSines = opt.corticalStateTerms;
        end
        corticalStatesUnrect = corticalStateSines * Ai(:);
        
    elseif ~calcStateGrad
         % cortical state parameters are part of givens (in opt)
         corticalStatesUnrect = opt.corticalStates;
         c = opt.const; 
    end
    
    
    nParams_expected = calcResponseGrad*nStim +  calcStateGrad * ( (1+opt.optimizeB + opt.optimizeC)*nSines + useAdditiveConstant );
    assert(nParams_expected == nParams);
    
    %%
%     corticalStates_stim = reshape(corticalStates, [nStim, nTrials]);
    corticalStatesUnrect_stim = reshape(corticalStatesUnrect, [nStim, nTrials]);
    if softPlusRect
        corticalStatesUnrect_stim_offsetScl = (corticalStatesUnrect_stim - opt.minCorticalState) * opt.softPlusScale;
        
        corticalStates_stim = log(1 + exp( corticalStatesUnrect_stim_offsetScl  ) ) / opt.softPlusScale + opt.minCorticalState;
        if any(isinf(corticalStates_stim(:)))
            idx_inf = isinf(corticalStates_stim);
            corticalStates_stim(idx_inf) = corticalStatesUnrect_stim(idx_inf);
%             keyboard;
        end
        corticalStates_stim2 = softplus(corticalStatesUnrect_stim, opt.minCorticalState, opt.softPlusScale);
        nonnans_unrect = ~isnan(corticalStatesUnrect_stim);
        assert(isequal(corticalStates_stim2(nonnans_unrect), corticalStates_stim(nonnans_unrect)));
        
    else
        corticalStates_stim = corticalStatesUnrect_stim;
        corticalStates_stim(corticalStates_stim < opt.minCorticalState) = opt.minCorticalState;
        
    end
    

    allDiffs = (c + r_true_rep .* corticalStates_stim) - R_full;
    
    logFactorial = @(N) gammaln(N+1);   % John D'Errico one line solution
    
    if opt.cost_loglikelihood
        R_mn = R_full;
        r_m = r_true_rep;
        %%
        idx_use = r_m(:) > 0 & corticalStates_stim(:) > 0;
        LL = -r_m .* corticalStates_stim + R_mn .* log ( r_m .* corticalStates_stim) - logFactorial(R_mn);
        J = sum(-LL(idx_use)); % take negative because are minimizing this function (but want to maximize LogLikelihood)
        
    elseif opt.cost_sumStimFirst
        % first sum over all stimuli, before squaring, then summing over trials
        allDiffs_sumOverStim = sum(allDiffs,1);
        
        J = sum( allDiffs_sumOverStim .^2 );  
    elseif opt.cost_sumTrialsFirst
        % first sum over all stimuli, before squaring, then summing over trials
        allDiffs_sumOverTrials = sum(allDiffs,2);
        
        J = sum( allDiffs_sumOverTrials(:) .^2 );
    else
        % square first, then sum over all stimuli & trials
        J = sum( allDiffs(:) .^2); 
    end
    
    
    if nargout > 1
        % calculate GRADIENT
        grad = zeros(nParams, 1);
        %%
        if calcResponseGrad
            
            if opt.cost_loglikelihood
                R_full_over_r_true = R_full ./ r_true_rep;
                if opt.fixGradients
                    R_full_over_r_true(isnan(R_full_over_r_true)) = 0; %% makes gradient match finite-differences, but causes fminunc to fail  
                end
                r_grad_eachStim = - sum(-corticalStates_stim + R_full_over_r_true, 2);
%                 r_grad_eachStim = - sum(-corticalStates_stim + R_full ./ r_true_rep, 4);
                r_grad = r_grad_eachStim(:);
                r_grad(isnan(r_grad)) = 0;

            elseif opt.cost_sumStimFirst
                allDiffs_sumOverStim_rep = allDiffs_sumOverStim(ones(1, nStim), :);    

                r_grad_eachTrial = 2 * allDiffs_sumOverStim_rep .* corticalStates_stim;
                r_grad = sum(r_grad_eachTrial,2);

            elseif opt.cost_sumTrialsFirst
    %             allDiffs_sumOverTrial_rep = allDiffs_sumOverTrials(:,:,:, ones(1, nTrials));    

                corticalStates_sumOverTrials = sum(corticalStates_stim, 2);

                r_grad_eachStim = 2 * allDiffs_sumOverTrials .* corticalStates_sumOverTrials;
                r_grad = r_grad_eachStim(:);

            else
                r_grad_eachTrial = 2 * allDiffs .* corticalStates_stim;
                r_grad = sum(r_grad_eachTrial,2);
            end

            assert(length(r_grad) == nStim)
            grad(idx_stim) = r_grad(:);
            
        end
        
        if calcStateGrad
        
            if useAdditiveConstant

                if opt.cost_sumStimFirst
                    grad(idx_const) = sum( 2 * allDiffs_sumOverStim  );
                elseif opt.cost_sumTrialsFirst
                    grad(idx_const) = sum( 2 * allDiffs_sumOverTrials  );
                else
                    grad(idx_const) = sum( 2 * allDiffs(:)  );
                end
            end

            corticalStateSines_stim = reshape(corticalStateSines, [nStim, nTrials, nSines]);
            
            for j = 1:nSines
                % gradient for Ai
                if opt.cost_loglikelihood

                    if softPlusRect
                        dCorticalState_dAi = 1 ./ (1 + exp(-corticalStatesUnrect_stim_offsetScl)) .* corticalStateSines_stim(:,:,j);
                    else
                        dCorticalState_dAi = corticalStateSines; %1 ./ (1 + exp(-corticalStatesUnrect));
                    end
                        
                    if softPlusRect
  
%                         corticalStateSines_stim = reshape(corticalStateSines, [nStim, nTrials, nSines]);
%                         r_m_times_sines = r_m .* corticalStateSines_stim(:,:,j);
                        grad_allTrials = (R_full  ./ (corticalStates_stim) - r_m) .* dCorticalState_dAi;

                        grad(idx_a(j)) = nansum( -grad_allTrials(:) ); % take negative because are minimizing this function (but want to maximize LogLikelihood)
                        
                        
                    else

                        r_m_times_sines = r_m .* corticalStateSines_stim(:,:,j);
                        grad_allTrials = -r_m_times_sines + R_full .* r_m_times_sines ./ (r_m .* corticalStates_stim);

                        grad(idx_a(j)) = nansum( -grad_allTrials(:) ); % take negative because are minimizing this function (but want to maximize LogLikelihood)
                    end
                elseif opt.cost_sumStimFirst
                    

                    corticalState_times_r_eachTrial = sum( corticalStateSines_stim(:,:,j) .* r_true_rep, 1);

                    grad(idx_a(j)) = sum( 2 * allDiffs_sumOverStim .* corticalState_times_r_eachTrial );

                    % gradient for bi
                    if opt.optimizeB 
    %                     grad(idx_b(j)) = sum( 2 * allDiffs .* corticalStateCosines(:,j) .* Ai(j) .* r_true_rep(:) .* R_full_t(:) ) ;
                        error('not implemented yet');
                    end
                    % gradient for ci;
                    if opt.optimizeC
    %                     grad(idx_c(j)) = sum( 2 * allDiffs .* corticalStateCosines(:,j) .* Ai(j) .* r_true_rep(:) ) ;
                        error('not implemented yet');
                    end

                elseif opt.cost_sumTrialsFirst

                    corticalStateSines_stim = reshape(corticalStateSines, [nStim, nTrials, nSines]);                
                    corticalState_times_r_sumOverTrials = sum( corticalStateSines_stim(:,:,j) .* r_true_rep, 2);
                    grad(idx_a(j)) = sum( 2 * allDiffs_sumOverTrials(:) .* corticalState_times_r_sumOverTrials(:) );

                    % gradient for bi
                    if opt.optimizeB
    %                     grad(idx_b(j)) = sum( 2 * allDiffs .* corticalStateCosines(:,j) .* Ai(j) .* r_true_rep(:) .* R_full_t(:) ) ;
                        error('not implemented yet');
                    end

                    % gradient for ci;
                    if opt.optimizeC
    %                     grad(idx_c(j)) = sum( 2 * allDiffs .* corticalStateCosines(:,j) .* Ai(j) .* r_true_rep(:) ) ;
                        error('not implemented yet');
                    end

                else
                    %%
                    grad(idx_a(j)) = sum( 2 * allDiffs(:) .* corticalStateSines(:,j) .* r_true_rep(:) );

                    % gradient for bi
                    if opt.optimizeB
                        grad(idx_b(j)) = sum( 2 * allDiffs(:) .* corticalStateCosines(:,j) .* Ai(j) .* r_true_rep(:) .* R_full_t(:) ) ;
                    end

                    % gradient for ci;
                    if opt.optimizeC
                        grad(idx_c(j)) = sum( 2 * allDiffs(:) .* corticalStateCosines(:,j) .* Ai(j) .* r_true_rep(:) ) ;
                    end                

                end
            end % for j = 1:nSines


        end % if calcStateGrad
        
        
        
    end    

end






%{
function [J, grad] = costFunction_response_CorticalState(Ac, R_full_t, r_true, R_full, bi, useAdditiveConstant)
    % calculate the cost function & gradient for a particular set of cortical state
    % parameters, given a fixed response profile
        
    nTrials = size(R_full, 4);
    r_true_rep = r_true(:,:,:,ones(1,nTrials));
%     nSines = size(corticalStateSines, 2);
    nParams = length(Ac);
    nSines = (nParams-useAdditiveConstant)/2;

    idx_a = 1          : nSines;
    idx_c = 1 + nSines : 2*nSines;

    if useAdditiveConstant
        c = Ac(1);
        idx_a = idx_a + 1;
        idx_c = idx_c + 1;
    else
        c = 0;
    end    
    
    Ai = Ac(idx_a); % 2-5
    ci = Ac(idx_c); % 10-13
    
    bxplusC = bsxfun(@plus, bsxfun(@times, bi(:)', repmat(R_full_t(:), [1, nSines])), ci(:)') ;
    
    corticalStateSines = sin(bxplusC);
    corticalStateCosines = cos(bxplusC);
    
    assert(length(Ai) == nSines);
    
    corticalStates = corticalStateSines * Ai(:);

    allDiffs = (c + r_true_rep(:) .* corticalStates(:))  -  R_full(:);

    J = sum( allDiffs .^2);
    if nargout > 1
        grad = zeros(nParams, 1);
        
        if useAdditiveConstant
            grad(1) = sum( 2 * allDiffs  );
        end
        for j = 1:nSines
            % gradient for Ai
            grad(idx_a(j)) = sum( 2 * allDiffs .* corticalStateSines(:,j) .* r_true_rep(:) );
            
%             % gradient for bi
%             grad(idx_b(j)) = sum( 2 * allDiffs .* corticalStateCosines(:,j) .* Ai(j) .* r_true_rep(:) .* R_full_t(:) ) ;
            
            % gradient for ci;
            grad(idx_c(j)) = sum( 2 * allDiffs .* corticalStateCosines(:,j) .* Ai(j) .* r_true_rep(:) ) ;
            
        end
    end    

end



%}