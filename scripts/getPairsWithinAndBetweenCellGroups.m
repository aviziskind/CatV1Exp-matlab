function [pairsWithinCellGroups, pairsBetweenCellGroups] = getPairsWithinAndBetweenCellGroups(cellGroupIds, grpF1oDCs, nOris, nSps, nPhs )
    
    % cellGroupIds is a vector [1..nCells], containing the group ids of each cell.
	% grpF1oDCs and nPhs are [1..nGroups] vectors, with the F1/DC and nPh of each group


    % check which groups have the same nori/nsp/nph
    uniqueGroupIds = unique(cellGroupIds);
    okPairs = false( length(uniqueGroupIds) );
    
    uniqueCombinations = unique([nOris(:), nSps(:), nPhs(:)], 'rows');
    disp(uniqueCombinations);
    
    for i = 1:length(uniqueGroupIds)
       for j = i+1:length(uniqueGroupIds)
            okPairs(i,j) = all([nOris(i) nSps(i) nPhs(i)] == [nOris(j) nSps(j) nPhs(j)]);
            okPairs(j,i) = okPairs(i,j);
       end
    end    
    
    [pairsWithinCellGroups, pairsBetweenCellGroups] = getPairsWithinAndBetweenGroups(groupIds, okPairs, maxPairs);
    nW = size(pairsWithinCellGroups,1);
    nB = size(pairsBetweenCellGroups,1);

    % create a 3rd column which contains the smaller F1/DC of the pair.
    for i = 1:nW
        inds = pairsWithinCellGroups(i,[1:2]);
        pairsWithinCellGroups(i,3) = min(grpF1oDCs(inds));        
    end
    for i = 1:nB
        inds = pairsBetweenCellGroups(i,[1:2]);
        pairsBetweenCellGroups(i,3) = min(grpF1oDCs(inds));        
    end

end