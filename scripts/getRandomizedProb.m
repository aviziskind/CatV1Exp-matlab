function [pval, two_sided_sgn, L, nPermutes] = getRandomizedProb(val_wcc, val_permute, tail)
                    
    nPermutes = max(length(val_permute), length(val_wcc));
    two_sided_sgn = 0;
    useCorrectProbs = 1;
    if useCorrectProbs
        probFunc = @(L,N) iff(L>0, (L+1)/(N+2),  0);
    else
        probFunc = @(L,N) L/N;
    end

    if nargin < 3
        tail = 'both';
    end
    
    if ischar(tail)
        tail = switchh(tail, {'left', 'right', 'both'}, [-1, 1, 0]);
    end
    
    if length(val_wcc) == 1
        switch tail
            case -1,  %left-tailed test: count # rand perms less than val_wcc. 
                % (null hypothesis: could get rand value as low as val_wcc by chance.)
                % e.g.  mean/medians/ks-stats of differences (smaller if spatial clustering)                
                L = nnz( val_permute <= val_wcc); % left-tailed test: detect if wcc is on the 
            case 1, % right-tailed test: count # rand perms greater than val_wcc. 
                % (null hypothesis: could get rand value as high as val_wcc by chance.)
                % e.g.  test of clustering indices: (larger if spatial clustering)
                L = nnz( val_permute >= val_wcc); 
            case 0,   
                % two-sided test: count # rand perms greater dist from M as val_wcc.
                % (null hypothesis: could get rand value deviating from M as much as val_wcc by chance.)
                % e.g.  test of correlation coefficient distribution                                
                M = mean(val_permute);
                L = nnz( abs(val_permute-M) >= abs(val_wcc-M) );
                if nnz(val_wcc>M)
                    two_sided_sgn = 1;
                else
                    two_sided_sgn = -1;
                end                    
        end
    else
        val_permute2 = val_wcc;
        if length(val_permute2) ~= length(val_permute)
            error('Must be same length');
        end

        switch tail
            case -1,  %e.g.  mean/medians/ks-stats
                L = nnz( val_permute2 <= val_permute); % left-tailed test
            case 1,
                L = nnz( val_permute2 >= val_permute2); % right-tailed test
            case 0,   %e.g.  cc
                L1 = nnz( val_permute2 <= val_permute); % left-tailed test
                L2 = nnz( val_permute2 >= val_permute); % right-tailed test
                L = min(L1, L2);
        end
        
        
        
    end
    
    pval = probFunc(L, nPermutes);
                    
end
