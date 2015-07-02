function [Wrcc_pairIdxs, nSimp_Wrcc,   Brcc_pairIdxs, nSimp_Brcc] = getRandomizedPairs(randInfo, reset_rand_seed, nResamples)
    nUnits = randInfo.nUnits;
    idxMtx = randInfo.idxMtx;
    separateSimpleComplexCells = isfield(randInfo, 'unit_simple_tf')  && ~isempty(randInfo.unit_simple_tf);
    getBrccPairs = isfield(randInfo, 'getBrccPairs') && isequal(randInfo.getBrccPairs, 1);

    if ~reset_rand_seed
        randInfo.randseed = [];
    end

    if getBrccPairs
        [Wrcc_pairIdxs_M, ~, Brcc_pairIdxs_M] = generateListOfCellPairs(randInfo,  nResamples,  randInfo.subRanges);
    else
        [Wrcc_pairIdxs_M]                     = generateListOfCellPairs(randInfo,  nResamples,  randInfo.subRanges);
    end
    Wrcc_pairIdxs =   cellfun(@(idxs_M) idxMtx( idxs_M ), Wrcc_pairIdxs_M, 'un', 0);
    Wrcc_pairs = cellfun(@(idxs) ind2subV([nUnits, nUnits], idxs), Wrcc_pairIdxs_M, 'un', 0);     

    nSimp_Wrcc = [];
    if separateSimpleComplexCells
        unit_simple_tf = randInfo.unit_simple_tf;
        nSimp_Wrcc = cellfun( @(idxs) sum( unit_simple_tf(idxs), 2), Wrcc_pairs, 'un', 0);
    end

    Brcc_pairIdxs = [];
    nSimp_Brcc = [];
    if getBrccPairs
        Brcc_pairIdxs =   cellfun(@(idxs_M) idxMtx( idxs_M ), Brcc_pairIdxs_M, 'un', 0);
        Brcc_pairs = cellfun(@(idxs) ind2subV([nUnits, nUnits], idxs), Brcc_pairIdxs_M, 'un', 0);     
        if separateSimpleComplexCells
            nSimp_Brcc = cellfun( @(idxs) sum( unit_simple_tf(idxs), 2), Brcc_pairs, 'un', 0);
        end
    end

end
