function testCostFunctionGradient(func, p0, alsoHessianFlag)
%%
    nParams = numel(p0);
    eps = 1e-8;
    grad_finite = zeros(nParams, 1);
    
    doHessian = exist('alsoHessianFlag', 'var') && isequal(alsoHessianFlag, 1);
    doHessianShort = 1;
    
    if doHessian
        hessian_finite = zeros(nParams, nParams);
        [C0, grad, hess] = func(p0);
    else
        [C0, grad] = func(p0);
    end
    
    nJ_max = nParams;
    if doHessianShort && doHessian 
        nJ_max = min(40, nParams);
    end
    
    progressBar('init-', nParams, 40)
    for i = 1:nJ_max
        
        p_i_plus = p0;  p_i_plus(i)  = p0(i)+eps;
%         p_i_minus = p0; p_i_minus(i) = p0(i)-eps;
        
        C_i_plus  = func(p_i_plus);
%         C_minus = func(p_i_minus);
%         grad_finite(i) = (C_plus - C_minus)/(2*eps);
        grad_i = (C_i_plus - C0)/(eps);
        grad_finite(i) = grad_i;
        
        if doHessian

            
            for j = i:nJ_max
                p_i_plus_j_plus = p_i_plus;  p_i_plus_j_plus(j)  = p_i_plus(j)+eps;
                p_j_plus = p0;               p_j_plus(j)         = p_j_plus(j)+eps;
%                 p_j_minus = p0; p_i_minus(i) = p0(i)-eps;
                
                C_j_plus  = func(p_j_plus);

                C_i_plus_j_plus  = func(p_i_plus_j_plus);
                
                hessian_finite(i,j) = (C_i_plus_j_plus - C_i_plus - C_j_plus + C0)/(eps.^2);
                hessian_finite(j,i) = hessian_finite(i,j);
            end            
            
        end
        progressBar(i);
    end
    progressBar('done');
    
    %%
    d_grad = grad - grad_finite;
    d_grad_rel = (d_grad(:))./(abs(grad)+eps);
    idx_use = abs(grad) > max( abs(grad)) *.01;
    max_d_grad_rel = max(abs(d_grad_rel( idx_use )));
    
    lgrad = log(abs(grad(:))); lgrad_fd = log(abs(grad_finite(:)));
    dlgrad = abs(lgrad - lgrad_fd);
    d_log_grad_rel = dlgrad ./ (abs(lgrad) + eps);
    max_d_log_grad_rel = max(abs(d_log_grad_rel));
    %%
    cc = corr(grad, grad_finite);
    if (max_d_grad_rel > sqrt(eps)*10) && max_d_log_grad_rel > .05   && abs(1-cc) > 1e-3
        warning('Gradient Not good');
        keyboard;
        if 0
            %%
            figure(59); clf;
            plot(grad, grad_finite, 'o'); hold on
            fplot(@(x) x, xlim, 'r-')
            
            
            
        end
        
        
    end

    %%
    if doHessian
        if doHessianShort
            %%
            n = nJ_max;

            diag_hess = diag(hess(1:n, 1:n));
            diag_hess_fd = diag(hessian_finite(1:n, 1:n));
            
            hess_short = hess(1:n, 1:n);
            hess_fd_short = hessian_finite(1:n, 1:n);
            d_hess_short = hess_short - hess_fd_short;
            d_hess_rel = (d_hess_short(:))./(hess_fd_short(:)+1e-6);
        else
            d_hess = hess - hessian_finite;
            d_hess_rel = (d_hess(:))./(hess(:)+1e-6);
        end
        
        
        if max(abs(d_hess_rel)) > sqrt(eps)
            error('Hessian Not good');
        end

    end
end