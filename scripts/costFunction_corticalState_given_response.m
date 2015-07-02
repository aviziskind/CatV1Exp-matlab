function [J, grad] = costFunction_corticalState_given_response(params, r_true, opt)
    % calculate the cost function & gradient for a particular set of cortical state
    % parameters, given a fixed response profile
    %%
    useAdditiveConstant = opt.useAdditiveConstant;
    R_full = opt.R_full;
    R_full_t = opt.R_full_t;
    
    [nStim, nTrials] = size(R_full);
    r_true_rep = r_true(:,ones(1,nTrials));
%     nSines = size(corticalStateSines, 2);
    nParams = length(params);
    nSines = opt.nStateTerms;
   
    nParams_expected = (1+opt.optimizeB + opt.optimizeC)*nSines + useAdditiveConstant;
    assert(nParams_expected == nParams);
    
    idx_a = 1 : nSines;
    idx_b = idx_a + nSines;
    idx_c = idx_a + nSines * (opt.optimizeB + opt.optimizeC);
    
    if useAdditiveConstant
        assert(~opt.cost_loglikelihood)
        c = params(1);
        idx_a = idx_a + 1;
        idx_b = idx_b + 1;
        idx_c = idx_c + 1;
    else
        c = 0;
    end    
    
    Ai = params(idx_a); % 2-5
    
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
    
    assert(length(Ai) == nSines);
    
    corticalStates = corticalStateSines * Ai(:);
    corticalStates(corticalStates < opt.minCorticalState) = opt.minCorticalState;
    
    corticalStates_stim = reshape(corticalStates, size(R_full));

    allDiffs = (c + r_true_rep .* corticalStates_stim) - R_full;
    
    logFactorial = @(N) gammaln(N+1);   % John D'Errico one line solution

    if opt.cost_loglikelihood
        R_mn = R_full;
        r_m = r_true_rep;
        
        idx_use = r_m(:) > 0 & corticalStates_stim(:) > 0;
        LL = -r_m .* corticalStates_stim + R_mn .* log ( r_m .* corticalStates_stim) - logFactorial(R_mn);
        J = sum(-LL(idx_use)); % take negative because are minimizing this function (but want to maximize LogLikelihood)
        
    elseif opt.cost_sumStimFirst
        % first sum over all stimuli, before squaring, then summing over trials
        allDiffs_sumOverStim = sum(sum(sum(allDiffs,1),2),3);
        
        J = sum( allDiffs_sumOverStim .^2 );  
    elseif opt.cost_sumTrialsFirst
        % first sum over all stimuli, before squaring, then summing over trials
        allDiffs_sumOverTrials = sum(allDiffs,4);
        
        J = sum( allDiffs_sumOverTrials(:) .^2 );
    else
        % square first, then sum over all stimuli & trials
        J = sum( allDiffs(:) .^2); 
    end
    
    
    if nargout > 1
        grad = zeros(nParams, 1);
        %%
        if useAdditiveConstant
            
            if opt.cost_sumStimFirst
                grad(1) = sum( 2 * allDiffs_sumOverStim  );
            elseif opt.cost_sumTrialsFirst
                grad(1) = sum( 2 * allDiffs_sumOverTrials  );
            else
                grad(1) = sum( 2 * allDiffs(:)  );
            end
        end
        
        for j = 1:nSines
            % gradient for Ai
            if opt.cost_loglikelihood
                                
                corticalStateSines_stim = reshape(corticalStateSines, [size(R_full), nSines]);
                
%                 grad_allTrials = -r_m_times_sines + R_full .* r_m_times_sines ./ (r_m .* corticalStates_stim);
%%
                r_m_times_sines = r_m .* corticalStateSines_stim(:,:,j);
                r_m_times_sines_over_rmCort = r_m_times_sines ./ (r_m .* corticalStates_stim);
                r_m_times_sines_over_rmCort(isnan(r_m_times_sines_over_rmCort)) = 0;
                grad_allTrials = -r_m_times_sines + R_full .* r_m_times_sines_over_rmCort;
                
                grad(idx_a(j)) = sum( -grad_allTrials(:) ); % take negative because are minimizing this function (but want to maximize LogLikelihood)
 3;
%                 grad(idx_
                
            elseif opt.cost_sumStimFirst
                corticalStateSines_stim = reshape(corticalStateSines, [size(R_full), nSines]);
                
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
                
                corticalStateSines_stim = reshape(corticalStateSines, [size(R_full), nSines]);                
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
        end
        
    end    

end


