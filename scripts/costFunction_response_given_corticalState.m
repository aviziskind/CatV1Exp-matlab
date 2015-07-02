function [J, grad, hess] = costFunction_response_given_corticalState(r_true, c, corticalStates, opt)
    % calculate the cost function & gradient for a particular response profile
    % given a fixed cortical state function
    useAdditiveConstant = opt.useAdditiveConstant;
    R_full = opt.R_full;
    
    corticalStates(corticalStates < opt.minCorticalState) = opt.minCorticalState;
    corticalStates_stim = reshape(corticalStates, size(R_full));
    
    if ~useAdditiveConstant
        c = 0;
    end
    
    [nStim, nTrials] = size(R_full);
       
    r_true = r_true(:);
    r_true_rep = r_true(:,ones(1,nTrials));

    allDiffs = (c + r_true_rep .* corticalStates_stim)  -  R_full;
    

    logFactorial = @(N) gammaln(N+1);   % John D'Errico one line solution
    
    if opt.cost_loglikelihood
        R_mn = R_full;
        r_m = r_true_rep;
        
        idx_use = r_m(:) > 0 & corticalStates_stim(:) > 0;
        LL = -r_m .* corticalStates_stim + R_mn .* log ( r_m .* corticalStates_stim ) - logFactorial(R_mn);
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
    
%     fprintf('r');
    if nargout > 1
%        fprintf('R');
        %%
        if opt.cost_loglikelihood
            
            R_full_over_r_true = R_full ./ r_true_rep;
            if opt.fixGradients
                R_full_over_r_true(isnan(R_full_over_r_true)) = 0; %% makes gradient match finite-differences, but causes fminunc to fail            
            end
            r_grad_eachStim = - sum(-corticalStates_stim + R_full_over_r_true, 2);  %r_grad_eachStim = - sum(-corticalStates_stim + R_full ./ r_true_rep, 4);

            grad = r_grad_eachStim(:);
            grad(isnan(grad)) = 0;
           
        elseif opt.cost_sumStimFirst
            allDiffs_sumOverStim_rep = allDiffs_sumOverStim(ones(1, nStim), :);    
            
            grad_eachTrial = 2 * allDiffs_sumOverStim_rep .* corticalStates_stim;
            grad = sum(grad_eachTrial,2);
            
        elseif opt.cost_sumTrialsFirst
%             allDiffs_sumOverTrial_rep = allDiffs_sumOverTrials(:,:,:, ones(1, nTrials));    
            
            corticalStates_sumOverTrials = sum(corticalStates_stim, 2);
            
            grad_eachStim = 2 * allDiffs_sumOverTrials .* corticalStates_sumOverTrials;
            grad = grad_eachStim(:);
            
        else
            grad_eachTrial = 2 * allDiffs .* corticalStates_stim;
            grad = sum(grad_eachTrial,2);
        end
        
        grad = grad(:);

        assert(length(grad) == length(r_true) )
    end
    
    % calculate hessian:
    if nargout > 2
        
%         corticalStates_trials = reshape(corticalStates, [nOri*nSpf*nPh, nTrials]);
        
        corticalStates_sumSqrs = sum( corticalStates_stim.^2, 2);        
        hess = diag(  2 .* corticalStates_sumSqrs(:)  );        
        
    end

end
