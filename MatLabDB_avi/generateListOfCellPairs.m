function [WccPairs, WcmPairs, BccPairs, BcmPairs, BmmPairs] = generateListOfCellPairs(cellinfo, nSamples, subRangesIncl, origFlag)
    
    %{ 
    Inputs:
    * cellinfo can either be :
        - a vector - nUnitsEachSite - a list of the number of cells at each site. (assumes no MU)
        - a cell array containing vectors of cellIds {[0 1 2 3], [0 1 2], ...}
            what's important is that the multiunit is 0 and the other cells are non-zero.
        - or a struct, containing two fields that are each [1..nGroups] vectors:
           nUnitsEachSite and  sitesWithMU (optional. if left out, assumes there are no multi-units)
        
    
    * subRangesIncl: for between site comparison: if you want only certain groups
     compared with specific other groups, specify the subRangesIncl variable as a 
     cell array containing the lists of sites of each types
     eg: {[1 2 3 4], [5 6 7 8], ...}
        (then, cells in group 1 will only be compared with groups 2, 3, 4, etc)
    Typically, these groups correspond to # of oris / # spfs / # phases - you don't
    want to compare groups with (for example) other groups with a different # of phases.

    * subpopulations : for random permuting, this allows a cell array of cell indices which
    indicates list of subpopulations which should be permuted separately. (preserving the # of each
    subpopulation at each site).        
        
    
    Outputs:
    the outputs are 1xnSamples cell arrays (if nSamples is defined), or
    just simple matrices (if nSamples is not defined).
    
    %}
    outputClass = 'uint32';
    showWorking = false;
    convertToIndex = true;
    handleSubpopulations = 0;    
    sortOutputIdxs = 1;
    
    % input 1: cellinfo
    if iscell(cellinfo)
        groupedCellIds = cellinfo(:);        
        nUnitsEachSite = cellfun(@length, groupedCellIds);
        sitesWithMU = cellfun(@(cellIds) any(cellIds == 0), groupedCellIds);        
        nCellsEachSite = nUnitsEachSite-sitesWithMU;
    elseif isstruct(cellinfo)
        nUnitsEachSite = cellinfo.nUnitsEachSite(:);
        if isfield(cellinfo, 'sitesWithMU')
            sitesWithMU    = cellinfo.sitesWithMU(:);
        else
            sitesWithMU    = false(size(nUnitsEachSite));
        end
        nCellsEachSite = nUnitsEachSite-sitesWithMU;
        if isfield(cellinfo, 'subpopulations')
            subpopulations = cellinfo.subpopulations;
            handleSubpopulations = 1;
        end
        if isfield(cellinfo, 'randseed') && ~isempty(cellinfo.randseed);
            seed = cellinfo.randseed;
            rand('state', seed); %#ok<RAND>
        end            
        
        groupedCellIds = arrayfun(@(tf, n) [iff(tf, 0, []) , 1:n], sitesWithMU, nCellsEachSite, 'un', false);
    elseif isnumeric(cellinfo)
        nUnitsEachSite = cellinfo(:);
        sitesWithMU    = false(size(nUnitsEachSite));
        nCellsEachSite = nUnitsEachSite-sitesWithMU;
        
        groupedCellIds = arrayfun(@(n) [1:n], nCellsEachSite, 'un', false);        
    end
    
    nSites = length(nUnitsEachSite);
    
    idx = 1;
    groupedCellIndices = cell(1,nSites);
    for i = 1:nSites
        groupedCellIndices{i} = idx : idx+nUnitsEachSite(i)-1;
        idx = idx+nUnitsEachSite(i);
    end
%     assert(isequal([groupedCellIndices{:}], 1:sum(nUnitsEachSite)));
    doProgressBar = 1;
    if isfield(cellinfo, 'doProgressBar') 
        doProgressBar = cellinfo.doProgressBar;
    end
    
    % input 2: nSamples
    if nargin < 2
        nSamples = 1;
    end
    if doProgressBar
        progressBar('init', nSamples, 20);
    end
    
    % input 3: subRangesIncl
    if exist('subRangesIncl', 'var') && ~isempty(subRangesIncl)
        nRanges = length(subRangesIncl);        
        subRangesIncl = cellfun(@(x) x(:)', subRangesIncl, 'un', false); % make sure are horizontal vectors
        v = [subRangesIncl{:}];
        if length(v) ~= length(unique(v))
            error('Cannot have duplicate indices in the subRangesIncl cell array');
        end        
    else
        nRanges = 1;
        subRangesIncl = {[1:nSites]};
    end    
    
%     nCellsInEachRange = cellfun(@(x) length(x), subRangesIncl);
%     cumsumRanges0 = [0 cumsum(nCellsInEachRange)];
    
    % input 4: 'original order' flag
    randomizeCells = ~exist('origFlag', 'var') || isempty(origFlag);    

    % whether output between-cell pairs as well:
    doBetweenSitePairs = (nargout >= 3);
    
 
            
    nMaxCellsPerSite = max(nCellsEachSite);        
    cellPairIds = cell(1,nMaxCellsPerSite);
    for n = 2:nMaxCellsPerSite
        cellPairIds{n} = nchoosek(1:n, 2);
    end    
    nPairsForNCells = cellfun(@(x) size(x,1), cellPairIds);
    
    [rng_nCellsEachSite, rng_sitesWithMU,  rng_idxOfCells, rng_idxOfMUs, rng_unitsIndices] = deal(cell(1, nRanges)); 
    [rng_nWccPairs, rng_nWcmPairs, rng_nBccPairs, rng_nBcmPairs, rng_nBmmPairs] = deal(zeros(1, nRanges)); 
    
    % calculate # of pairs for each type of pair
    for ri = 1:nRanges
        
        rng_nCellsEachSite{ri} = nCellsEachSite(subRangesIncl{ri}); 
        rng_sitesWithMU{ri}    = sitesWithMU(subRangesIncl{ri}); 

        groupCellIdsForThisRange = groupedCellIds(subRangesIncl{ri});
        
        cellIdList = [groupCellIdsForThisRange{:}];                
        rng_unitsIndices{ri} = [groupedCellIndices{subRangesIncl{ri}}];
        
        rng_idxOfCells{ri} = find(cellIdList ~= 0);
        rng_idxOfMUs{ri}   = find(cellIdList == 0);
                
        N = rng_nCellsEachSite{ri}(:);     
        mu_sites = rng_sitesWithMU{ri};
        rng_nWccPairs(ri) = sum(nPairsForNCells(N(N>1)));
        rng_nWcmPairs(ri) = sum(N(mu_sites));        
        if doBetweenSitePairs
            rng_nBccPairs(ri) = ( sum(sum( N * N')) - sum(N.^2) )/2;
            sub = sum( arrayfun(@(nm) sum( N(setdiff(1:length(N), nm))), find(~mu_sites)))  ;
            
            rng_nBcmPairs(ri) = sum(N)*(length(N)-1) - sub;
            rng_nBmmPairs(ri) = nnz(mu_sites)*(nnz(mu_sites)-1)/2;
        end        
    end        
    cumsumWccRanges0 = [0, cumsum(rng_nWccPairs)];
    cumsumWcmRanges0 = [0, cumsum(rng_nWcmPairs)];
    cumsumBccRanges0 = [0, cumsum(rng_nBccPairs)];
    cumsumBcmRanges0 = [0, cumsum(rng_nBcmPairs)];
    cumsumBmmRanges0 = [0, cumsum(rng_nBmmPairs)];

%     cumsumMUs0   = [0; cumsum(sitesWithMU)];
%     cumsumCells0 = [0; cumsum(nCellsEachSite)];
    
    nUnits = sum(nUnitsEachSite);
    nMUsTot = nnz(sitesWithMU);
    nCellsTot = nUnits - nMUsTot;
       
    % preallocate memory:
    WccPairs = cell(1,nSamples);
    WcmPairs = cell(1,nSamples);
    if doBetweenSitePairs
        BccPairs = cell(1,nSamples);
        BcmPairs = cell(1,nSamples);
        BmmPairs = cell(1,nSamples);
    end
    
    
    function [Wcc_Pairs_i, Wcm_Pairs_i, Bcc_Pairs_i, Bcm_Pairs_i, Bmm_Pairs_i] = generateOneSample(...
        nCellsEachSite, sitesWithMU,  idxOfCells, idxOfMUs,   nPairs,  cellPairIds, subpop_cells, subpop_MUs)
        
        if ~doBetweenSitePairs            
            nWccPairs = nPairs(1); nWcmPairs = nPairs(2);
        else
            nWccPairs = nPairs(1); nWcmPairs = nPairs(2); 
            nBccPairs = nPairs(3); nBcmPairs = nPairs(4); nBmmPairs = nPairs(5);
        end
        
        nUnitsEachSite = nCellsEachSite+sitesWithMU;
        nSites = length(nCellsEachSite);
        nCells = length(idxOfCells);
        nMUs = length(idxOfMUs);
        cumsumMUs0   = [0; cumsum(sitesWithMU)];
        cumsumCells0 = [0; cumsum(nCellsEachSite)];        
%         cumsumUnits0 = [0; cumsum(nUnitsEachSite)];
    
        Wcc_Pairs_i = zeros(nWccPairs, 2);
        Wcm_Pairs_i = zeros(nWcmPairs, 2);
        if doBetweenSitePairs
            Bcc_Pairs_i = zeros(nBccPairs, 2);
            Bcm_Pairs_i = zeros(nBcmPairs, 2);
            Bmm_Pairs_i = zeros(nBmmPairs, 2);
        end
        
        if randomizeCells
            if ~handleSubpopulations            
                newCell_order = randperm(nCells);
                newMU_order   = randperm(nMUs);
            else                
                newCell_order = randperm_subpop(subpop_cells);
                newMU_order   = randperm_subpop(subpop_MUs);                                                
            end
            
            
        else
            newCell_order = 1:nCells;
            newMU_order   = 1:nMUs;
        end
        cell_positions = idxOfCells(newCell_order);
        mu_positions   = idxOfMUs(newMU_order);
        
        
        if showWorking
            %% 
            disp('*** Separate (cell/mu) Indexing : ***')
            for si = 1:nSites
                cellIdxs = [cumsumCells0(si)+1 : cumsumCells0(si+1) ];
                muIdx = cumsumMUs0(si)+1;
                
                cellIdxs = newCell_order([cellIdxs]);                
                muIdx = iff(sitesWithMU(si), @() newMU_order(muIdx), @() 0);
                
                disp(['Site ' num2str(si) ' : [ ' iff(sitesWithMU(si), num2str(muIdx()), ' ') ' ] ' num2str(cellIdxs) ]);
            end
            
            disp('*** Grouped Indexing : *** ')
            for si = 1:nSites
                cellIdxs = [cumsumCells0(si)+1 : cumsumCells0(si+1) ];
                muIdx = cumsumMUs0(si)+1;
                
                cellIdxs2 = cell_positions([cellIdxs]);
                muIdx2 = iff(sitesWithMU(si), @() mu_positions(muIdx), @() 0);
                disp(['Site ' num2str(si) ' : [ ' iff(sitesWithMU(si), num2str(muIdx2()), ' ') ' ] ' num2str(cellIdxs2) ]);
            end
        end
        
        Wcm_idx = 1;
        Wcc_idx = 1;
        if doBetweenSitePairs
            Bcc_idx = 1;
            Bcm_idx = 1;
            Bmm_idx = 1;
        end
        
        for si = 1:nSites
            nU = nUnitsEachSite(si);
            isMU = sitesWithMU(si); 
            nC = nU - isMU;
            cellIdxs = [cumsumCells0(si)+1 : cumsumCells0(si+1) ];            
            cellIdxs = cell_positions([cellIdxs]);
            
            % Within-site cell-cell pairs
            if nC >= 2
                Wcc_pairs = cellIdxs(cellPairIds{nC});
                Wcc_Pairs_i(Wcc_idx : Wcc_idx+nPairsForNCells(nC)-1, :) = Wcc_pairs;
                Wcc_idx = Wcc_idx + nPairsForNCells(nC);
            end
            
            % Within-site cell-multi-unit pairs
            if isMU 
                muIdx = cumsumMUs0(si)+1;
                muIdx = mu_positions(muIdx);
                if (nC > 0)
%                     Wcm_pairs = [ones(nC,1) * muIdx,  cellIdxs(:)];
                    Wcm_pairs = [muIdx( ones(nC,1) ),  cellIdxs(:)];
                    Wcm_Pairs_i(Wcm_idx :Wcm_idx+nC-1,:) = Wcm_pairs;
                    Wcm_idx = Wcm_idx + nC;
                end
            end
            
            if doBetweenSitePairs                                
                if si > 1
                    cellIdxsOfEarlierSites_raw = cumsumCells0(1)+1:cumsumCells0(si);
                    cellIdxsOfEarlierSites = cell_positions(cellIdxsOfEarlierSites_raw)';
                    
                    muIdxsOfEarlierSites_raw = cumsumMUs0(1)+1:cumsumMUs0(si);
                    muIdxsOfEarlierSites = mu_positions( muIdxsOfEarlierSites_raw );
                else
                    cellIdxsOfEarlierSites = [];
                    muIdxsOfEarlierSites = [];
                end
                
                if si < nSites
                    cellIdxsOfLaterSites_raw = cumsumCells0(si+1)+1:cumsumCells0(end);
                    cellIdxsOfLaterSites = cell_positions(cellIdxsOfLaterSites_raw)';
                    
                    muIdxsOfLaterSites_raw = cumsumMUs0(si+1)+1:cumsumMUs0(end);
                    muIdxsOfLaterSites = mu_positions(muIdxsOfLaterSites_raw)';

                else
                    cellIdxsOfLaterSites = [];
                    muIdxsOfLaterSites = [];
                end                                

%                 cellIdxsOfOtherSites = [cellIdxsOfEarlierSites; cellIdxsOfLaterSites];                                
                nCellsEarlier = length(cellIdxsOfEarlierSites);
                nCellsLater   = length(cellIdxsOfLaterSites);  
                
                nMUsEarlier = length(muIdxsOfEarlierSites);
                nMUsLater   = length(muIdxsOfLaterSites);  
                
%                 nMUsEarlier = length
                
                % Between-site cell-cell pairs
                if (nC >= 1) && (nCellsLater > 0)                   
                   nPairs = nC*nCellsLater;
                   
                   Bcc_pairs = zeros(nPairs,2);
                   for j = 1:nC
                       Bcc_pairs( (j-1)*nCellsLater+1:j*nCellsLater , :) = [ones(nCellsLater,1)*cellIdxs(j), cellIdxsOfLaterSites];
                   end      
                   Bcc_Pairs_i(Bcc_idx : Bcc_idx+nPairs-1, :) = Bcc_pairs;
                   Bcc_idx = Bcc_idx + nPairs;                   
                end
                
                % Between-site cell-multiunit pairs
                if isMU && (nCellsEarlier+nCellsLater > 0)
                    nPairs = nCellsEarlier+nCellsLater;
                    Bcm_pairs_earlier = [cellIdxsOfEarlierSites, ones(nCellsEarlier,1)* muIdx];
                    Bcm_pairs_later   = [ones(nCellsLater,1)* muIdx,  cellIdxsOfLaterSites   ];
                    Bcm_Pairs_i(Bcm_idx : Bcm_idx+nPairs-1, :) = [Bcm_pairs_earlier; Bcm_pairs_later];
                    Bcm_idx = Bcm_idx + nPairs;
                end
                
                % Between-site multiunit-multiunit pairs
                if isMU && (nMUsEarlier+nMUsLater > 0)
                    
                    nPairs = nMUsLater;
                    Bmm_pairs_later   = [ones(nMUsLater,1)* muIdx,  muIdxsOfLaterSites  ];
                    Bmm_Pairs_i(Bmm_idx : Bmm_idx+nPairs-1, :) = [Bmm_pairs_later];
                    Bmm_idx = Bmm_idx + nPairs;
                end
                
                
                
            end
        end
        
        assert(Wcc_idx-1 == nWccPairs);
        assert(Wcm_idx-1 == nWcmPairs);
        if doBetweenSitePairs
            assert(Bcc_idx-1 == nBccPairs);
            assert(Bcm_idx-1 == nBcmPairs);
            assert(Bmm_idx-1 == nBmmPairs);
        end
%         Bcm_Pairs_i(Bcm_idx:end,:) = [];
        
        % sort pairs in ascending order
        if ~randomizeCells
            assert(all(all( Wcc_Pairs_i == sort(Wcc_Pairs_i, 2))));
            assert(all(all( Wcm_Pairs_i == sort(Wcm_Pairs_i, 2))));
            if doBetweenSitePairs
                assert(all(all( Bcc_Pairs_i == sort(Bcc_Pairs_i, 2))));
                assert(all(all( Bcm_Pairs_i == sort(Bcm_Pairs_i, 2) )));
                assert(all(all( Bmm_Pairs_i == sort(Bmm_Pairs_i, 2) )));
            end
        
%             if doBetweenSitePairs
%                 assert(all(all( Bcc_Pairs_i == sort(Bcc_Pairs_i, 2))));
% %                 Bcc_Pairs_i = sort(Bcc_Pairs_i, 2);
%                 Bcm_Pairs_i = sort(Bcm_Pairs_i, 2);
%             end
            
        else
            Wcc_Pairs_i = sort(Wcc_Pairs_i, 2);
            Wcm_Pairs_i = sort(Wcm_Pairs_i, 2);
            if doBetweenSitePairs
                Bcc_Pairs_i = sort(Bcc_Pairs_i, 2);
                Bcm_Pairs_i = sort(Bcm_Pairs_i, 2);
                Bmm_Pairs_i = sort(Bmm_Pairs_i, 2);
            end
        end            
    end

    if convertToIndex
        f = @(v) mtxSub2ind(v, nUnits);
    else
        f = @(v) v;
    end

    %
    
    [subpop_cells, subpop_MUs] = deal(cell(1, nRanges));
    if handleSubpopulations
        assert(isequal(sort([subpopulations{:}]), 1:nUnits))
        for ri = 1:nRanges
            [~, subpop_cells{ri}] = cellfun(@(i) intersect(rng_idxOfCells{ri}, i), subpopulations, 'un', 0);
            [~, subpop_MUs{ri}]   = cellfun(@(i) intersect(rng_idxOfMUs{ri},   i), subpopulations, 'un', 0);
        end
    end
        
    doProgressBar = doProgressBar && nSamples > 10;
    for sample_i = 1:nSamples
        
        nCols = iff(convertToIndex, 1, 2);
        WccPairs_sample = zeros( sum(rng_nWccPairs), nCols );
        WcmPairs_sample = zeros( sum(rng_nWcmPairs), nCols );
        if doBetweenSitePairs
            BccPairs_sample = zeros( sum(rng_nBccPairs), nCols );
            BcmPairs_sample = zeros( sum(rng_nBcmPairs), nCols );
            BmmPairs_sample = zeros( sum(rng_nBmmPairs), nCols );
        end
                
        for ri = 1:nRanges
            Wcc_rng_i_idx = [cumsumWccRanges0(ri)+1:cumsumWccRanges0(ri+1)];
            Wcm_rng_i_idx = [cumsumWcmRanges0(ri)+1:cumsumWcmRanges0(ri+1)];
            if doBetweenSitePairs
                Bcc_rng_i_idx = [cumsumBccRanges0(ri)+1:cumsumBccRanges0(ri+1)];
                Bcm_rng_i_idx = [cumsumBcmRanges0(ri)+1:cumsumBcmRanges0(ri+1)];
                Bmm_rng_i_idx = [cumsumBmmRanges0(ri)+1:cumsumBmmRanges0(ri+1)];
            end            

            if doBetweenSitePairs
                nPairs = [rng_nWccPairs(ri), rng_nWcmPairs(ri), rng_nBccPairs(ri), rng_nBcmPairs(ri), rng_nBmmPairs(ri) ];
                [WccPairs_si_ri, WcmPairs_si_ri, BccPairs_si_ri, BcmPairs_si_ri, BmmPairs_si_ri] = generateOneSample(...
                rng_nCellsEachSite{ri}, rng_sitesWithMU{ri}, rng_idxOfCells{ri}, rng_idxOfMUs{ri},  nPairs,  cellPairIds, subpop_cells{ri}, subpop_MUs{ri});                
            else
                nPairs = [rng_nWccPairs(ri), rng_nWcmPairs(ri)];
                [WccPairs_si_ri, WcmPairs_si_ri] = generateOneSample(...
                rng_nCellsEachSite{ri}, rng_sitesWithMU{ri}, rng_idxOfCells{ri}, rng_idxOfMUs{ri},  nPairs,  cellPairIds, subpop_cells{ri}, subpop_MUs{ri});
            end
                            
            % convert from relative subscripts to subscripts of full array.
            WccPairs_si_ri_sub = rng_unitsIndices{ri}(WccPairs_si_ri);
            WcmPairs_si_ri_sub = rng_unitsIndices{ri}(WcmPairs_si_ri);
            if doBetweenSitePairs
                BccPairs_si_ri_sub = rng_unitsIndices{ri}(BccPairs_si_ri);
                BcmPairs_si_ri_sub = rng_unitsIndices{ri}(BcmPairs_si_ri);
                BmmPairs_si_ri_sub = rng_unitsIndices{ri}(BmmPairs_si_ri);
            end
                        
            % convert from subscripts to indices, and put into full sample vector 
            WccPairs_sample(Wcc_rng_i_idx, :) = f(WccPairs_si_ri_sub);
            WcmPairs_sample(Wcm_rng_i_idx, :) = f(WcmPairs_si_ri_sub);            
            if doBetweenSitePairs
                BccPairs_sample(Bcc_rng_i_idx, :) = f(BccPairs_si_ri_sub);
                BcmPairs_sample(Bcm_rng_i_idx, :) = f(BcmPairs_si_ri_sub);
                BmmPairs_sample(Bmm_rng_i_idx, :) = f(BmmPairs_si_ri_sub);
            end
        end
                
        WccPairs{sample_i} = WccPairs_sample;
        WcmPairs{sample_i} = WcmPairs_sample;
        if doBetweenSitePairs
            BccPairs{sample_i} = BccPairs_sample;            
            BcmPairs{sample_i} = BcmPairs_sample;
            BmmPairs{sample_i} = BmmPairs_sample;
        end
        
        if doProgressBar
            progressBar(sample_i);
        end
    end
    if doProgressBar
%         progressBar('done');
    end
    
    if (nargin < 2) || (nSamples == 1)
        WccPairs = WccPairs{1};
        WcmPairs = WcmPairs{1};
        if doBetweenSitePairs
            BccPairs = BccPairs{1};
            BcmPairs = BcmPairs{1};
            BmmPairs = BmmPairs{1};
        end
    end
    
    if ~strcmp(outputClass, 'double');
        if ~iscell(WccPairs)
            WccPairs = cast(WccPairs, outputClass);
            WcmPairs = cast(WcmPairs, outputClass);
        else
            WccPairs = cellfun(@(x) cast(x, outputClass), WccPairs, 'un', 0);
            WcmPairs = cellfun(@(x) cast(x, outputClass), WcmPairs, 'un', 0);            
        end            
        if doBetweenSitePairs
            %%
            if ~iscell(BccPairs)
                BccPairs = cast(BccPairs, outputClass);
                BcmPairs = cast(BcmPairs, outputClass);
                BmmPairs = cast(BmmPairs, outputClass);
            else
                BccPairs = cellfun(@(x) cast(x, outputClass), BccPairs, 'un', 0);
                BcmPairs = cellfun(@(x) cast(x, outputClass), BcmPairs, 'un', 0);            
                BmmPairs = cellfun(@(x) cast(x, outputClass), BmmPairs, 'un', 0);            
            end
        end
    end
    
    if sortOutputIdxs
                    
        if ~iscell(WccPairs)
            WccPairs = sort(WccPairs);
            WcmPairs = sort(WcmPairs);
        else
            WccPairs = cellfun(@sort, WccPairs, 'un', 0);
            WcmPairs = cellfun(@sort, WcmPairs, 'un', 0);            
        end            
        if doBetweenSitePairs
            if ~iscell(BccPairs)
                BccPairs = sort(BccPairs);
                BcmPairs = sort(BcmPairs);
                BmmPairs = sort(BmmPairs);
            else                
                BccPairs = cellfun(@sort, BccPairs, 'un', 0);
                BcmPairs = cellfun(@sort, BcmPairs, 'un', 0);            
                BmmPairs = cellfun(@sort, BmmPairs, 'un', 0);            
            end
        end
    end
        
        
end


% sub = sum( arrayfun(@(nm) sum( N(setdiff(1:length(N), nm))), find(~mu_sites)))  ;
%   is equivalent to : 
% no_mus = find(~mu_sites);
% sub = 0;
% for i = 1:length(no_mus)
%     N_tmp = N;
%     N_tmp(no_mus(i)) = 0;
%     sub = sub + sum(N_tmp);
% end



    %{
    % input 4: subRangesExcl
    if exist('subRangesExcl', 'var') && ~isempty(subRangesExcl)
        
        subRangesExcl = cellfun(@(x) x(:)', subRangesExcl, 'un', false); % make sure are horizontal vectors
        v = [subRangesExcl{:}];
        if length(v) ~= length(unique(v))
            error('Cannot have duplicate indices in the subRangesIncl cell array');
        end        
        
        relevantExclForEachSubRange = cell(1,nRanges);        
        for ri = 1:nRanges
            inRng = cellfun(@(idx) ~isempty(intersect(subRangesIncl{ri}, idx)), subRangesExcl);   
            relevantExclForEachSubRange{ri} = subRangesExcl(inRng);
        end
    else
        relevantExclForEachSubRange = cell(1,nRanges);
    end    
    %}

% %     * subRangesExcl: if, within cells in particular subGroup, you want to *exclude*
% %      certain group parings, specify this variable as a cell array containing the lists 
% %     of sites of each types
% %      eg: {[1 2 3], [4], [5 6], [7 8], ...}
% %      then cells in group 1 will not be compared with cells in group 2 or 3, but only with
% %      cells in group 4.
% %     I use this for making sure between-site comparisons are always done between sites
% %     that are in different animals, (ie. never within the same animal).

