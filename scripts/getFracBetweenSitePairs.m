function Wcc_Bcc_ratio = getFracBetweenSitePairs(nCellsEachSite)

    nPairsForCells = @(n) n*(n-1)/2;
    nCellsTot = sum(nCellsEachSite);
    
    
    nPairsTot = nPairsForCells(nCellsTot);
    nWccPairsEachSite = arrayfun(nPairsForCells, nCellsEachSite);
    
    nWccTot = sum(nWccPairsEachSite);
    nBccTot = nPairsTot - nWccTot;
    Wcc_Bcc_ratio = nBccTot / nPairsTot;
end