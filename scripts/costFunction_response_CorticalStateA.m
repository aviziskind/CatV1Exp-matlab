function [J, grad] = costFunction_response_CorticalStateA(params, opt)
%                                                     (Ac, R_full_t, r_true, R_full, bi, useAdditiveConstant)
%%
    % calculate the cost function & gradient for a particular set of cortical state
    % parameters, given a fixed response profile
    useAdditiveConstant = opt.useAdditiveConstant;
    corticalStateSines = opt.corticalStateSines;
    R_full = opt.R_full;    
    
    [nOri, nSpf, nPh, nTrials] = size(R_full);
    nStim = nOri*nSpf*nPh;
    
    r_true = reshape(params(1:nStim), [nOri, nSpf, nPh]);    
    r_true_rep = r_true(:,:,:,ones(1,nTrials));
%     nSines = size(corticalStateSines, 2);

%     nStateParams = length(Ac);
    nSines = size(corticalStateSines, 2);

    idx_a = (nStim + 1) : nStim + nSines;
%     idx_c = nStim + 1 + nSines : 2*nSines;

    if useAdditiveConstant
        c = params(nStim+1);
        idx_a = idx_a + 1;
%         idx_c = idx_c + 1;
    else
        c = 0;
    end    
    
    Ai = params(idx_a); % 2-5
%     ci = Ac(idx_c); % 10-13
    
%     bxplusC = bsxfun(@plus, bsxfun(@times, bi(:)', repmat(R_full_t(:), [1, nSines])), ci(:)') ;
    
%     corticalStateSines = sin(bxplusC);
%     corticalStateCosines = cos(bxplusC);
    
%     assert(length(Ai) == nSines);
    
    corticalStates = corticalStateSines * Ai(:);

    allDiffs = (c + r_true_rep(:) .* corticalStates(:))  -  R_full(:);

    J = sum( allDiffs .^2);
    
    if nargout > 1
        grad = zeros(length(params), 1);
  
        % gradient for response 
        stim_grad_eachTrial = 2 * allDiffs .* corticalStates(:);
        stim_grad = sum(reshape(stim_grad_eachTrial, size(R_full)),4);
        grad(1:nStim) = stim_grad(:);
        
        % gradient for cortical state parameters
        if useAdditiveConstant
            grad(nStim + 1) = sum( 2 * allDiffs  );
        end
        
        for j = 1:nSines
            % gradient for Ai
            grad(idx_a(j)) = sum( 2 * allDiffs .* corticalStateSines(:,j) .* r_true_rep(:) );
            
%             % gradient for bi
%             grad(idx_b(j)) = sum( 2 * allDiffs .* corticalStateCosines(:,j) .* Ai(j) .* r_true_rep(:) .* R_full_t(:) ) ;
            
            % gradient for ci;
%             grad(idx_c(j)) = sum( 2 * allDiffs .* corticalStateCosines(:,j) .* Ai(j) .* r_true_rep(:) ) ;
            
        end
    end    

end


