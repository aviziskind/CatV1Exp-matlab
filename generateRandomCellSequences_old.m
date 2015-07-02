function [cellCellPairs, cellMUPairs] = generateRandomCellSequences_old(cellinfo, nSamples, subRanges, origFlag)
    
    %{ 
    Inputs:
    * cellinfo can either be :
        - a cell array containing vectors of cellIds {[0 1 2 3], [0 1 2], ...}
        - or a struct, containing two fields that are each [1..nGroups] vectors:
           nUnitsEachSite and  sitesWithMU
    
    * subRanges: if you want only certain cell types mixed with certain others, specify
    the subRanges variable as a cell array containing the lists of cells of each types
     eg: {[1 2 3 4], [5 6 7 8], ...}
    
    Outputs:
    the outputs are 1xnSamples cell arrays (if nSamples is defined), or
    just simple matrices (if nSamples is not defined).
    
    %}
    
    dbug = true;
    
    if iscell(cellinfo)
        cellinfo = cellinfo(:);
        nSites = length(cellinfo);                
        nUnitsEachSite = cellfun(@length, cellinfo);
        sitesWithMU = cellfun(@(cellIds) any(cellIds == 0), cellinfo);
        allCellIds = [cellinfo{:}];
        
    elseif isstruct(cellinfo)
        nUnitsEachSite = cellinfo.nUnitsEachSite;
        sitesWithMU    = cellinfo.sitesWithMU;
        nSites = length(nUnitsEachSite);
        allCellIds = arrayfun(@(tf, n) [iff(tf, 0, []) , 1:n], sitesWithMU, nUnitsEachSite, 'un', false);
    end

    if nargin < 2
        nSamples = 1;
    end

    randomizeCells = 0;% ~exist('origFlag', 'var') || isempty(origFlag);
    
    idxOfCells = find([allCellIds{:}] ~= 0);
    idxOfMUs = find([allCellIds{:}] == 0);

    cumsumUnits = cumsum(nUnitsEachSite);    
    nCellsEachSite = nUnitsEachSite - sitesWithMU;
    nMaxCellsPerSite = max(nCellsEachSite);    
    
    cellPairIds = cell(1,nMaxCellsPerSite);
    for n = 2:nMaxCellsPerSite
        cellPairIds{n} = nchoosek(1:n, 2);
    end
    nPairsForNCells = cellfun(@(x) size(x,1), cellPairIds);
    
    cumsumMUs = cumsum(sitesWithMU);
    cumsumMUs0 = [0, cumsumMUs];
    cumsumCells = cumsum(nCellsEachSite);
    cumsumCells0 = [0, cumsumCells];
    
    nUnits = cumsumUnits(end);
    nMUs = nnz(sitesWithMU);
    nCells = nUnits - nMUs;
    
    nCellCellPairs = sum(nPairsForNCells(nCellsEachSite(nCellsEachSite>1)));
    nCellMUPairs = sum(nCellsEachSite(sitesWithMU));
    
    cellCellPairs = cell(1,nSamples);
    cellMUPairs   = cell(1,nSamples);
    generateOneSample;
    
    function [cc_pairs, cm_pairs] = generateOneSample
        
        cellCellPairs_i = zeros(nCellCellPairs, 2);
        cellMUPairs_i = zeros(nCellMUPairs, 2);
        
        if randomizeCells
            newCell_order = randperm(nCells);
            newMU_order   = randperm(nMUs);
        else
            newCell_order = 1:nCells;
            newMU_order   = 1:nMUs;
        end
        cell_order = idxOfCells(newCell_order);
        mu_order   = idxOfMUs(newMU_order);
        
        
        if dbug
            disp('*** Separate (cell/mu) Indexing : ***')
            for si = 1:nSites
                cellIdxs = [cumsumCells0(si)+1 : cumsumCells0(si+1) ];
                muIdx = cumsumMUs0(si)+1;
                
                cellIdxs = newCell_order([cellIdxs]);
                muIdx = newMU_order(muIdx);
                
                disp(['Site ' num2str(si) ' : [ ' iff(sitesWithMU(si), num2str(muIdx), ' ') ' ] ' num2str(cellIdxs) ]);
            end
            
            disp('*** Grouped Indexing : *** ')
            for si = 1:nSites
                cellIdxs = [cumsumCells0(si)+1 : cumsumCells0(si+1) ];
                muIdx = cumsumMUs0(si)+1;
                
                cellIdxs = cell_order([cellIdxs]);
                muIdx = mu_order(muIdx);
                disp(['Site ' num2str(si) ' : [ ' iff(sitesWithMU(si), num2str(muIdx), ' ') ' ] ' num2str(cellIdxs) ]);
            end
        end
        
        cm_idx = 1;
        cc_idx = 1;
        
        for si = 1:nSites
            nU = nUnitsEachSite(si);
            isMU = sitesWithMU(si); 
            nC = nU - isMU;
            cellIdxs = [cumsumCells0(si)+1 : cumsumCells0(si+1) ];
            muIdx = cumsumMUs0(si)+1; 
            cellIdxs = cell_order([cellIdxs]);
            muIdx = mu_order(muIdx);
            
            if nC >= 2
                cc_pairs = cellIdxs(cellPairIds{nC});
                cellCellPairs_i(cc_idx : cc_idx+nPairsForNCells(nC)-1, :) = cc_pairs;
                cc_idx = cc_idx + nPairsForNCells(nC);
            end
            if isMU && (nC > 0)
                cm_pairs = [ones(nC,1) * muIdx,  cellIdxs(:)];
                cellMUPairs_i(cm_idx :cm_idx+nC-1,:) = cm_pairs;
                cm_idx = cm_idx + nC;
            end
            
        end
        %     assert(cc_idx-1 == nCellCellPairs);
        %     assert(cm_idx-1 == nCellMUPairs);
        
        % sort pairs in ascending order
        cellCellPairs_i = sort(cellCellPairs_i, 2);
        cellMUPairs_i   = sort(cellMUPairs_i, 2);
        
%         cellCellPairs{sample_i} = cellCellPairs_i;
%         cellMUPairs{sample_i} = cellMUPairs_i;


    end

    if nargin < 2
        cellCellPairs = cellCellPairs{1};
        cellMUPairs = cellMUPairs{1};
    end
    
end