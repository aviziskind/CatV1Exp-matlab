function [bestTh, bestErr, fp_rate, fn_rate] = optimalLinearThreshold(x1, x2, eps)
    % Determine the optimal threshold boundary between some data points (x1) and some noise points
    % (x2). eps determines the weight to put on type I (false positive) errors (including points
    % from x2) vs type I (false negative) errors (missing points from x1).

    if nargin < 3
        eps = 0.5;
    end

    Nx1 = length(x1);
    Nx2 = length(x2);
    Ntot = Nx1+Nx2;
    
    allVals = [x1(:); x2(:)];
    id = [ones(Nx1, 1); 2*ones(Nx2, 1)];
    [sortedVals, idx] = sort(allVals);
    sorted_id = id(idx);

    cumFrac_x1 = cumsum(sorted_id == 1) / Nx1;
    cumFrac_x2 = cumsum(sorted_id == 2) / Nx2;
%     edges = binCent2edge(allVals);
    3;
    %%
    if median(x1) < median(x2)  % x1 | x2
        falsePos_rate = cumFrac_x2;
        falseNeg_rate = 1-cumFrac_x1;
        x1_side = -1;
    else                        % x2 | x1
        falsePos_rate = 1-cumFrac_x2;
        falseNeg_rate = cumFrac_x1;
        x1_side = +1;
    end
    
    wgtErr = falsePos_rate * eps  + falseNeg_rate * (1-eps);
    [bestErr, bestTh_idx] = min( wgtErr(1:end-1) );
    
    if bestErr > 0.5
        error('should be less than 0.5'); % then could switch direction.
    end
        

    bestTh = mean(sortedVals([bestTh_idx, bestTh_idx+1]));
    
    show = 0;
    if show
        %%
        figure(644); clf;
        hist2({x1, x2}, 50, 'stairs', 'norm');
        drawVerticalLine(bestTh)
        figure(645); clf;        
        plot(sortedVals, cumFrac_x1, 'b-', sortedVals, cumFrac_x2, 'g-')
        drawVerticalLine(bestTh)
        figure(646); clf;        
        plot(sortedVals, falsePos_rate, 'g:', sortedVals, falseNeg_rate, 'r:', sortedVals, wgtErr, 'b')
        drawVerticalLine(bestTh);
        legend('false pos', 'false neg', 'wgt err')
            3;
    end
    
    doubleCheck = 0;
    if doubleCheck
        %%
        if x1_side == -1
            fp_rate = nnz(x1 > bestTh) / Nx1;
            fn_rate = nnz(x2 < bestTh) / Nx2;
        
            errRate = fp_rate * eps + fn_rate * (1-eps);
        elseif x1_side == 1;
            fp_rate = nnz(x1 < bestTh) / Nx1;
            fn_rate = nnz(x2 > bestTh) / Nx2;        
            
            errRate = fp_rate * eps + fn_rate * (1-eps);
        end
   
        assert(abs(errRate-bestErr) < 1e-15);
        
    end

    3;
    % double
    
    
%     for i = 1:Ntot+1
%         
%         
%         
%     end
%     
    





end