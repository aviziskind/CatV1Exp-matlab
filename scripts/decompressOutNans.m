function Y = decompressOutNans(X)

    persistent ok_idx_forDiags;
    
    assert( size(X,2) == 1);
    [diag_n, nBins] = dealV(X([end-1, end]));
    X([end-1, end]) = [];
    
    if isempty(ok_idx_forDiags)
        ok_idx_forDiags = {};
    end
    if (length(ok_idx_forDiags) < diag_n) || isempty(ok_idx_forDiags{diag_n})    
        [L,R] = meshgrid(1:nBins, 1:nBins);
        ok_idx = arrayfun(@(l,r) r-l >= 0 && r-l < diag_n, L, R);            
        
        ok_idx_forDiags{diag_n} = ok_idx;
    else
        ok_idx = ok_idx_forDiags{diag_n};
    end
    
    
%     [L,R] = meshgrid(1:nBins, 1:nBins);
% %     ok_idx = arrayfun(@(l,r) r-l >= 0 , L, R);                
%     ok_idx = arrayfun(@(l,r) r-l >= 0 && r-l < diag_n, L, R);            
    
        
%     if nStims == 1
        Y = nan(nBins, nBins, class(X));
        Y(ok_idx) = X;
        
%     else       
%         [ok_r,ok_c] = find(ok_idx);
% 
%         Y = nan(nBins, nBins, nStims, class(X) );
%         for n = 1:nStims
%             idx = sub2ind([nBins nBins nStims], ok_r, ok_c, ones(N,1)*n);
%             Y(idx) = X(:,n);
%         end        
%     end
end







% function Y = decompressOutNans(X)
% 
%     persistent nBins_cur diag_n_cur ok_idx_cur;
% %     persistent nBins_cur ok_idx_cur;
% 
%     N = size(X,1);
%     nBins = (-1 + sqrt(8*N+1))/2;
% 
%     nStims = size(X,2);
% %     
%     if isempty(ok_idx_cur) || (nBins_cur ~= nBins) || (diag_n_cur ~= diag_n)
%         [L,R] = meshgrid(1:nBins, 1:nBins);
% %         ok_idx = arrayfun(@(l,r) r-l >= 0 , L, R);                
%         ok_idx = arrayfun(@(l,r) r-l >= 0 && r-l < diag_n, L, R);            
%         [ok_idx_cur, nBins_cur, diag_n_cur] = deal(ok_idx, nBins, diag_n);
%     else
%         ok_idx = ok_idx_cur;
%     end
% %     [L,R] = meshgrid(1:nBins, 1:nBins);
% %     ok_idx = arrayfun(@(l,r) r >= l, L, R);            
%     
%     
%     if nStims == 1
%         Y = nan(nBins, nBins, class(X));
%         Y(ok_idx) = X;
%         
%     else       
%         [ok_r,ok_c] = find(ok_idx);
% 
%         Y = nan(nBins, nBins, nStims, class(X) );
%         for n = 1:nStims
%             idx = sub2ind([nBins nBins nStims], ok_r, ok_c, ones(N,1)*n);
%             Y(idx) = X(:,n);
%         end        
%     end
% end
