function idx_matches = matchUpCellsFromDifferentGroups(grpOsps)

    nGrps = length(grpOsps);

    nCells = cellfun(@length, grpOsps);
    
    idx_primary = indmax(nCells);
    idx_others  = setdiff(1:nGrps, idx_primary);
    
    idx_matches = arrayfun(@(n) zeros(1,n), nCells, 'un', 0);
    

    idx_matches = arrayfun(@(n) 1:n, nCells, 'un', 0);
    return;
    
    
%     figure(4);
    
    idx_matches{1} = [1:nCells(idx_primary)];
    I = idx_primary;
    for j = 1:length(idx_others)
        J = idx_others(j);
        
        r = zeros(nCells(I), nCells(J));
    
        for idx_i = 1:nCells(I)
            for idx_j = 1:nCells(J)            
                osp1 = mean(grpOsps{I}{idx_i},3); osp2 = mean(grpOsps{J}{idx_j},3);
                r(idx_i, idx_j) = spearmanRho(osp1(:), osp2(:));
            end
        end
        idx = indmax(r, [], 2);
        
        idx_matches{J} = idx;        
    end
%     figure(234); clf;
%     for i = 1:nCells(1)
%         subplot(1, nCells(1), i); imageOSP(grpOsps{1}{i}, 'mean:ph', 'OSP', 'nolabels', 'noTicks');
%     end    
%     for i = 1:nCells(2)
%         subplot(2, nCells(1), i+nCells(1)); imageOSP(grpOsps{2}{i}, 'mean:ph', 'OSP', 'nolabels', 'noTicks');
%     end
    
    
3;
   
end

