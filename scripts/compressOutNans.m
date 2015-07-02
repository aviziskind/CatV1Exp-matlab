function Y = compressOutNans(X)

    persistent ok_idx_forDiags L R;
    if isempty(ok_idx_forDiags)
        ok_idx_forDiags = {};
    end
    
    nBins = size(X,1);    
    assert(size(X,3) == 1);    
    
    assert(nBins == 120);
    
    % find first diagonal which is not all nans
    diag_n = nBins;
    diag_vals = X(nBins,1);
    while all(isnan(diag_vals)) && (diag_n >= 0)
        diag_n = diag_n-1;
        diag_vals = diag(X, -diag_n+1);    
    end
    if diag_n < 0
        error('No non-nans in array');
    end
        
    if (length(ok_idx_forDiags) < diag_n) || isempty(ok_idx_forDiags{diag_n})                    
        [L,R] = meshgrid(1:nBins, 1:nBins);
        ok_idx = arrayfun(@(l,r) r-l >= 0 && r-l < diag_n, L, R);            

        ok_idx_forDiags{diag_n} = ok_idx;
    else
        ok_idx = ok_idx_forDiags{diag_n};
    end
        
%     Y = X(ok_idx);
%     Y([end+1, end+2]) = [diag_n, nBins];

    Y = sparse(R(ok_idx), L(ok_idx), double(X(ok_idx)) );        
    
end

% function n = calcNumElements(L, nDiag)
%     n = L*(L+1)/2 - nDiag*(nDiag-1)/2;
% end


%     psthWindowDataFile = [CatV1Path 'MatLabDB_avi' filesep 'psthWindowData_tmp.mat'];        
%     allPsthWindowData_S = load(psthWindowDataFile);
%     fn = fieldnames(allPsthWindowData_S);
%     progressBar('init-', length(fn))
%     for i = 1:length(fn)
%         progressBar(i);
%         v = allPsthWindowData_S.(fn{i});
%         v_decomp = structfun(@decompressOutNans, v, 'un', 0);
%         v_recomp = structfun(@compressOutNans, v_decomp, 'un', 0);        
%         fn2 = fieldnames(v);
%         for j = 1:length(fn2)
%             d = v.(fn2{j})-v_recomp.(fn2{j});
%             assert( max(d) == 0);
%             assert( min(d) == 0);
%         end
%     end
%     progressBar('done')
%     save(psthWindowDataFile, '-struct', 'allPsthWindowData_S', '-v6');    


%{
function Y = compressOutNans(X)

    persistent nBins_cur diag_n_cur ok_idx_cur L R;

    nBins = size(X,1);
    assert(size(X,3) == 1);    
    
    % find first diagonal which is not all nans
    diag_n = nBins;
    diag_vals = X(nBins,1);
    while all(isnan(diag_vals)) && (diag_n >= 0)
        diag_n = diag_n-1;
        diag_vals = diag(X, -diag_n+1);    
    end
    if diag_n < 0
        error('No non-nans in array');
    end
        
    if isempty(ok_idx_cur) || (nBins_cur ~= nBins) || (diag_n_cur ~= diag_n)        
        [L,R] = meshgrid(1:nBins, 1:nBins);
        ok_idx = arrayfun(@(l,r) r-l >= 0 && r-l < diag_n, L, R);            
        [ok_idx_cur, nBins_cur, diag_n_cur] = deal(ok_idx, nBins, diag_n);
    else
        ok_idx = ok_idx_cur;
    end
        
%     Y = X(ok_idx);
%     Y([end+1, end+2]) = [diag_n, nBins];

    Y = sparse(R(ok_idx), L(ok_idx), double(X(ok_idx)) );    
    
    
end
%}