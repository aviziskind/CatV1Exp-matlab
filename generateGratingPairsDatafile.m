function generateGratingPairsDatafile(cmpOriSpfType)

    tic;
    if nargin < 1
        cmpOriSpfType = '';
    end

    randSeed = 1234;
%     nResamples = 10000;
%     nResamples = 20;
    visualizeCellPairs = false;
    randomizeCellPairOrder = 0;

%     if curMatchDB
%         nResamples = 10;
%     end
    
    gratingType = curGratingType('');
    cmpType = curCmpType('');        
    
    if strcmp(cmpType, 'degree')
        nResamples = 10000;
%         nResamples = 100;
    elseif strcmp(cmpType, 'phase')        
%         nResamples = 10; % we use phase-shuffling, not cell-shuffling for control
        nResamples = 10; % tmp
    end
        
    
    
    F1oDC_field = 'F1oDC_maxR_avP_sm';
    showOriSpPhCombinations = false;
    preserveNSimpleComplexAtASite = strcmp(cmpType, 'degree') && curPreserveSimpleComplex;
    preserveNAlignedAtASite       = strcmp(cmpType, 'degree') && curPreserveAligned && strcmp(gratingType, 'drifting') && strcmp(cmpOriSpfType, 'ori');
    
    permuteCellsToSitesWithSameNCells = 0;
%     assert( ~ (curSamePen && curPreserveSimpleComplex) );
    
    maxNPairs = 1e10;
%     maxNPairs = 10000;
    maxBcmPairs = 100;

    bccType = curBccType;
    bccWP = strcmp(bccType, 'samePen'); % to only compare spatial frequencies at similar eccentricities.
    bccWA = strcmp(bccType, 'sameAnimal'); 
    bccBP_WA = strcmp(bccType, 'diffPen-sameAnimal'); 
    
    mixDiffFlashedGratingSpfSequences = true;
    mixDiffDriftingGratingSpfSequences = true;
    
    addCellMUPairs = strcmp(cmpType, 'phase') && 1;
    subtractSpont = curSubtractSpont;
        
    if strcmp(cmpType, 'degree') 
        if nargin == 0
            generateGratingPairsDatafile('ori');    
            if ~curPreserveAligned
                generateGratingPairsDatafile('spf'); % aligned/anti-aligned - only relevant for ori (direction)
            end
            return;
        end
    end        
    
    cmpOriSpfType_str = switchh(cmpOriSpfType, {'ori', 'spf', ''}, {'(ORIENTATION ONLY)', '(SPATIAL FREQ. ONLY)', ''});
    cmpOriSpfType_fileext = switchh(cmpOriSpfType, {'ori', 'spf', ''}, {'_ori', '_spf', ''});
        
    if strcmp(cmpType, 'psth') && strcmp(gratingType, 'drifting')
        error('Can''t do psth comparisons with drifting grating data');
    end
    
%     if onlyDoDegreeCmpAtSamePenetration && strcmp(cmpType, 'degree')
%         filename_ext = [filename_ext '_samePen'];
%     end    
    
    if strcmp(cmpType, 'clusters')
        nResamples = 0;
    end
    
    subtractSpont_str = iff(curSubtractSpont, '(SPONT Subtracted)', '');
    bccType_str = iff(bccWP && strcmp(cmpType, 'degree'), sprintf('(BCC distribution : %s)', bccType), '');
    
    fprintf('\n *** Creating pair data file for %s grating stimuli %s %s... *** %s', gratingType, cmpOriSpfType_str, subtractSpont_str, bccType_str);

    % a. Retrieve list of OSPs    
    S1 = load(getFileName('osps', cmpOriSpfType_fileext));
    allCells = S1.allCells;

    
    
% %     nSpf = arrayfun(@(s) length(s.sp), allCells);
% %     allCells = allCells(nSpf > 1);
    % b. Generate list of pairs of cells to be compared.
%         allCells = sortStructArray(allCells, 'Gid');

    cellGroupIds  = [ allCells.Gid    ]';
    cellIds       = [ allCells.cellId ]';
    
    
    locData = [allCells.locData];
    PenIds  = [locData.PenId]';
    CatIds  = [locData.CatId]';
        
    allF1oDCsTypes = [allCells.F1oDC];
    F1oDCs = [allF1oDCsTypes.(F1oDC_field)];    

    assert(issorted(cellGroupIds));
%     nCells = length(cellIds);
        
    
    [uGroupIds, groupIdList, nUnitsEachSite] = uniqueList(cellGroupIds);
    
    indOfFirstCellInGroups = cellfun(@(x) x(1), groupIdList);
    sitesWithMU = cellIds(indOfFirstCellInGroups) == 0;    
    nSites = length(uGroupIds);
    nCellsEachSite = nUnitsEachSite(:) - sitesWithMU(:);
    indExpandToAllUnits = arrayfun(@(nU, site_id) ones(1, nU)*site_id, nUnitsEachSite, 1:nSites, 'un', 0);
    indExpandToAllUnits = [indExpandToAllUnits{:}];
    
        
    muIdxs   = find(cellIds == 0);  nMUs = length(muIdxs);
    cellIdxs = find(cellIds ~=0);   nCells = length(cellIdxs);
    fprintf('\nUsing data from %d units (%d cells and %d multi-units)\n', length(cellIds), nCells, nMUs);        
%         nUnits = nMUs + nCells;

    oris      = {allCells(indOfFirstCellInGroups).ori}'; % arrayfun(@(s) s.ori,  allCells(indOfFirstCellInGroups), 'un', 0);
    sps_pix   = {allCells(indOfFirstCellInGroups).sp}';  %  arrayfun(@(s) s.sp,   allCells(indOfFirstCellInGroups), 'un', 0);
    phs       = {allCells(indOfFirstCellInGroups).ph}';
    groupPenIds = PenIds( indOfFirstCellInGroups );
    groupCatIds = CatIds( indOfFirstCellInGroups );
    
    assert(isequal(PenIds, groupPenIds(indExpandToAllUnits)))
    assert(isequal(CatIds, groupCatIds(indExpandToAllUnits)))
%         grpAnimal = arrayfun(@(s) s.animal, allCells(indOfFirstCellInGroups));

%     degPerPixs = {allCells(indOfFirstCellInGroups).degPerPix};
%     sps_deg = cellfun(@times, sps_pix, num2cell(degPerPixs), 'un', false);
    sps_deg = sps_pix;

    if showOriSpPhCombinations
        grpStrs = cellfun(@(a,b,c) ['[' vec2str(a) '] X [' vec2str(b) '] X [' vec2str(c) ']'], oris, sps_deg, phs, 'un', false);     
        disp('Ori/Sp/Phase combinations:');
        varBreakdown(grpStrs);    
    end

    rows_arg = {};
    switch cmpType
        case 'phase',
            nOriSpPhs = [cellfun(@length, oris), cellfun(@length, sps_deg), cellfun(@length, phs)];
    %         spfTypes = cell2mat(sps_pix(nOriSpPhs(:,2) == 10) );
            if strcmp(gratingType, 'flashed') && ~mixDiffFlashedGratingSpfSequences
                
                comparableStimuliSets = [cell2mat( sps_pix), nOriSpPhs(:, 3) ];
                rows_arg = {'rows'};
            elseif strcmp(gratingType, 'drifting') && ~mixDiffDriftingGratingSpfSequences
                comparableStimuliSets = cellfun(@(x) num2str(x', '% .1f'), sps_pix, 'un', 0);
                rows_arg = {};
            else % are willing to compare phase tuning curves when sptial frequencies were sampled differently
                comparableStimuliSets = nOriSpPhs;
                rows_arg = {'rows'};
            end            

        case 'degree', % ori & spf degree tuning just uses interpolated tuning curves / statistics of tuning curves, so any ori/ori or spf/spf pair is fine
            % ori/dir tuning curves are of length, 30 (N = 3), 36 (N = 142), or 72 (N = 267)
            % spf tuning curves are of length 10 (N = 56) or 12 (N = 124)

            %             fullOriSet = (nOriSpPhs(:,1) >= 30);            
            % compare only spf-spf (1 vs 1),  and ori-ori (0 vs 0)
            
            if bccWP
                % compare only at the same penetration
                comparableStimuliSets = groupPenIds;
                
            elseif bccWA
                % compare only within the same animal
                comparableStimuliSets = groupCatIds;
                
            elseif bccBP_WA
                comparableStimuliSets = {groupCatIds, groupPenIds};
                
            else
                if permuteCellsToSitesWithSameNCells                
                    comparableStimuliSets = nCellsEachSite;  
                else                    
                    comparableStimuliSets = ones(1, nSites);  
                end
            end
            
        
        case 'clusters',            
            % for within-cell comparisons - we take clusters that belong to the same cell.
            % (we will collect these below)
            
            % for between-cell comparisons - we take clusters that are definitely from different
            % cells - by limiting to only pairs from different sites.            
            comparableStimuliSets = ones(length(uGroupIds), 1);
            
            
        case 'psth'
            PSTHs         = [ allCells(:).PSTH];
            frameLengths_ms = cat(1, PSTHs.frameLength_ms);

            comparableStimuliSets = frameLengths_ms(indOfFirstCellInGroups);        
        
        
    end
    

    
    
    if any( strcmp(cmpType, {'degree', 'phase'}) )
        
        cellinfo = struct('nUnitsEachSite', nUnitsEachSite, 'sitesWithMU', sitesWithMU);
        if ~iscell(comparableStimuliSets) %
            [~,subRanges] = uniqueList(comparableStimuliSets, rows_arg{:});
            [Wcc_idxs, Wcm_idxs, Bcc_idxs, Bcm_idxs, Bmm_idxs] = generateListOfCellPairs(cellinfo, 1, subRanges, 'originalOrder');
            cellinfo.subRanges = subRanges;
        else % can't do randomization control in this case
            for i = 1:2
                [~,subRanges{i}] = uniqueList(comparableStimuliSets{i}, rows_arg{:}); %#ok<AGROW>
                [Wcc_idxs_C{i}, Wcm_idxs_C{i}, Bcc_idxs_C{i}, Bcm_idxs_C{i}] = generateListOfCellPairs(cellinfo, 1, subRanges{i}, 'originalOrder'); %#ok<AGROW>
            end
            Wcc_idxs = Wcc_idxs_C{1}; assert(isequal(Wcc_idxs_C{1}, Wcc_idxs_C{2}));
            Wcm_idxs = Wcm_idxs_C{1}; assert(isequal(Wcm_idxs_C{1}, Wcm_idxs_C{2}));
            Bcc_idxs = setdiff(Bcc_idxs_C{1}, Bcc_idxs_C{2});
            Bcm_idxs = setdiff(Bcm_idxs_C{1}, Bcm_idxs_C{2});
            cellinfo.subRanges = subRanges;
        end
        
        nBcc_total = length(Bcc_idxs);
        if ~addCellMUPairs
            Wcm_idxs = [];
            Bcm_idxs = [];
        end

        S = struct('Wcc_idxs', Wcc_idxs, 'Wcm_idxs', Wcm_idxs, 'Bcc_idxs', Bcc_idxs, 'Bcm_idxs', Bcm_idxs, 'Bmm_idxs', Bmm_idxs);
        pairLabels = {'Within-site cell-cell pairs', 'Within-site cell-multi-unit pairs', ...
                      'Between-site cell-cell pairs', 'Between-site cell-multi-unit pairs', 'Between-site multiunit pairs'};
        grpType = 'cells';

        
    elseif strcmp(cmpType, 'clusters') % 2. Generate within-cell cluster-cluster pairs (same site)
        % We want:
        %   Wcc_idxs   (within-cell cluster-cluster pairs)
        %   Bcc_idxs   (between-cell (same site) cluster-cluster pairs)
        %   BScc_idxs  (Between-Site cluster-cluster pairs) - just computed above                
                
        BScc_idxs = [];%Bcc_idxs; % Between-cell (different site)
        f = 100;
        clusterIds       = [ allCells.clustId ]';
        effGroupIds = cellGroupIds*f + cellIds; % effective "groups" (are really cells)        
        
        % 3. Generate between-cell cluster-cluster pairs (same site)
        [uGroupCellIds, groupCellIdList, nClustersEachCell] = uniqueList(effGroupIds);                                
        clustGroupIds = floor(uGroupCellIds/f);            
        [~, subRanges_withinSite] = uniqueList(clustGroupIds);        
        
        [Wcc_idxs, ~, Bcc_idxs] = generateListOfCellPairs(nClustersEachCell, 1, subRanges_withinSite, 'originalOrder');
        dblCheck = 1;
        if dblCheck
            % BCds should be different Gids
            if ~isempty(BScc_idxs)
                BScc_pairs = ind2subV([nCells, nCells], BScc_idxs);        
                BScc_Gids = cellGroupIds(BScc_pairs);
                assert(all(BScc_Gids(:,1) ~= BScc_Gids(:,2)));
            end
            if any(Bcc_idxs == 49296)
                3;
            end
                
            % Bcc should be have same Gids, different cellIds 
            Bcc_pairs = ind2subV([nCells, nCells], Bcc_idxs);
            Bcc_Gids = cellGroupIds(Bcc_pairs);
            Bcc_cellIds = cellIds(Bcc_pairs);
            assert(  all(Bcc_Gids(:,1)    == Bcc_Gids(:,2)) );
            assert(  all(Bcc_cellIds(:,1) ~= Bcc_cellIds(:,2)) );
            
            % Wccs should have same Gids, same cellIds, different clusterIds
            Wcc_pairs = ind2subV([nCells, nCells], Wcc_idxs);        
            Wcc_Gids = cellGroupIds(Wcc_pairs);
            Wcc_cellIds = cellIds(Wcc_pairs);
            Wcc_clusterIds = clusterIds(Wcc_pairs);
            assert(  all(Wcc_Gids(:,1)       == Wcc_Gids(:,2)) );
            assert(  all(Wcc_cellIds(:,1)    == Wcc_cellIds(:,2)) );
            assert(  all(Wcc_clusterIds(:,1) ~= Wcc_clusterIds(:,2)) );                        
        end        
        
        S = struct('Wcc_idxs', Wcc_idxs,  'Wcm_idxs', [],  'Bcc_idxs', Bcc_idxs,  'Bcm_idxs', [],  'BScc_idxs', BScc_idxs);
        pairLabels = {'Within-cell cluster-cluster pairs', '', 'Between-cell (same site) cluster-cluster pairs', '', 'Between-cell (different site) cluster-cluster pairs'};
        grpType = 'clusters';
        
    end
        
        
     %   Wcc_idxs    (within-cell cluster-cluster pairs)
        %   Bcc_idxs  (between-cell (same site) cluster-cluster pairs)
        %   BScc_idxs
    pairNs = structfun(@length, S);

    if randomizeCellPairOrder
        pairTypes = fieldnames(S);
        for i = 1:length(pairTypes)            
            S.(pairTypes{i}) = S.(pairTypes{i})(randperm(pairNs(i)));
        end        
    end
        
    if ~isempty(maxNPairs)
        pairTypes = fieldnames(S);
        for i = 1:length(pairTypes)
            S.(pairTypes{i}) = reduceNPairsTo(S.(pairTypes{i}), maxNPairs);
        end
        if ~isempty(maxBcmPairs)
            S.Bcm_idxs = reduceNPairsTo(S.Bcm_idxs, maxBcmPairs);
        end
    end
    pairNs = structfun(@length, S);
    
    
    idx_ne = pairNs > 0;
    lbls = cellfun(@(n, s) sprintf('   %7d %s \n', n, s), num2cell(pairNs(idx_ne)'), pairLabels(idx_ne), 'un', 0);
    fprintf('From %d %s, found : \n%s', nCells, grpType, [lbls{:}]);
       
    if visualizeCellPairs
        %%
        nUnits = length(allCells);
        A = zeros(nUnits);
        for i = 1:nUnits
            A(i,i) = 1;
        end    
        A(1,1) = -1;

%     allPairs = sort([WccPairs; WcmPairs; BccPairs; BcmPairs]);
%     idxMtx = zeros(nUnits, nUnits);
%     idxMtx(sort(allPairs)) = 1:length(allPairs);    
%     figure(44);
%     imagesc(idxMtx);
%     3;


        tk.Wcc = .3;
        tk.Wcm = -.3;
        tk.Bcc = .8;
        tk.Bcm = -.8;
        tk.Bmm = -.5;
        
        tk_labels = fieldnames(tk);
        tk_vals = cell2mat( struct2cell(tk) );
        
        [tk_vals, new_ord] = sort(tk_vals);
        tk_labels = tk_labels(new_ord);

        A( Wcc_idxs ) =  .3;
        A( Wcm_idxs ) = -.3;
        A( Bcc_idxs ) =  .8;
        A( Bcm_idxs ) = -.8;    
        A( Bmm_idxs ) = -.5;    

        figure(1);
        imagesc(A); axis equal tight;
        title('All Pairs');
        h_col = colorbar;
        set(h_col, 'ytick', tk_vals, 'yticklabel', tk_labels);
    end        
    
    
    if nResamples > 0
        sc_str = iff(preserveNSimpleComplexAtASite, '(Preserving # simple/complex at a site)', '');
        ab_str = iff(preserveNAlignedAtASite,       '(Preserving # aligned/unaligned cells at a site)', '');                
        
        fprintf('Generating resampling list (N = %d) of within-site distributions %s%s: \n   |', nResamples, sc_str, ab_str);        
        cellinfo.randseed = randSeed;
        cellinfo.F1oDC_field = F1oDC_field;
        if preserveNSimpleComplexAtASite
            SimpleComplexIndices = {find(F1oDCs < 1), find(F1oDCs >= 1), find(isnan(F1oDCs))};
            cellinfo.subpopulations = SimpleComplexIndices;
        elseif preserveNAlignedAtASite
            %%
            statField = iff(subtractSpont, 'oriStats_ss', 'oriStats_si');
            useAllSpkMU = 0;
            dOriField = iff(useAllSpkMU, 'Ddir_pref_allSpkMU', 'Ddir_pref_smlSpkMU');            
            
            dOriMU = nestedFields(allCells, 'tuningStats', statField, dOriField, 1);
            isUnaligned = ibetween(dOriMU, 45, 135);

            UnalignedAlignedIndices = {find(isUnaligned == 1), find(isUnaligned == 0 | isnan(isUnaligned))};
            cellinfo.subpopulations = UnalignedAlignedIndices;
        end
        
%         nUnits = length(allCells);
%         allPairs = sort([Wcc_idxs; Wcm_idxs; Bcc_idxs; Bcm_idxs]);
%         idxMtx = zeros(nUnits, nUnits);
%         idxMtx(sort(allPairs)) = 1:length(allPairs);    
%         3;
        
        
        
         if ~iscell(comparableStimuliSets)
             [Wrcc_idxs, Wrcm_idxs] = generateListOfCellPairs(cellinfo, nResamples, subRanges);
         else
             %%
             % sample from between-site distribution.
             nWcc = length(Wcc_idxs);
             nWcm = length(Wcm_idxs);
             [Wrcc_idxs, Wrcm_idxs] = deal( cell(1, nResamples) );
             progressBar('init-', nResamples, 30);
             for i = 1:nResamples
                 Wrcc_idxs{i} = randsample(Bcc_idxs, nWcc);
                 if nWcm > 0
                    Wrcm_idxs{i} = randsample(Bcm_idxs, nWcm); 
                 end
                 progressBar(i);
             end
             
             %%
             
             
         end
        
        
        if ~addCellMUPairs
            Wrcm_idxs = cell(1,nResamples);
        end               
        if strcmp(cmpType, 'degree')
            %%
            nUniqueIdxs = length( uniqueInts(Wrcc_idxs) );
    %         if any(strcmp(Wcc_idxs) missing_str
            nMissing = (length(Bcc_idxs) + length(Wcc_idxs)) - nUniqueIdxs;
            fprintf('\nThese contain %d unique pairs (missing %d) \n', nUniqueIdxs, nMissing);
        end
        
        chkOnlyCompareComparableStimuli = 1 && ~iscell(comparableStimuliSets);
        if chkOnlyCompareComparableStimuli
            %%
%             nUnits = length(allCells);
            cell_cmpVal = comparableStimuliSets(indExpandToAllUnits);
            cmpStim_M = bsxfun(@eq, cell_cmpVal(:), cell_cmpVal(:)');
            
%             [cell_cmpVal_x, cell_cmpVal_y] = meshgrid(cell_cmpVal);
            
            for i = 1:nResamples
%                 [xi, yi] = ind2sub([nUnits, nUnits], Wrcc_idxs{i});
                Wrcc_cmp_match  = ( cmpStim_M(Wrcc_idxs{i}) );
                assert(all(Wrcc_cmp_match));
            end    
            
        end
        
        chkPreservationOfNSimpleNComplex = 1;
        if chkPreservationOfNSimpleNComplex && preserveNSimpleComplexAtASite
            
            
            
            simple_idx = (F1oDCs>=1);
            nSimpMtx = bsxfun(@plus, simple_idx(:), simple_idx(:)');
                        
            [uSC, sc_count] = uniqueCount( nSimpMtx(Wcc_idxs) );
            for i = 1:nResamples
                [uSC_i, sc_count_i] = uniqueCount( nSimpMtx(Wrcc_idxs{i}) );
                assert(isequal(sc_count, sc_count_i));         
                assert(~isequal(Wrcc_idxs{i}, Wcc_idxs))
                
                if 0
                    %%
                    nUnits = length(allCells);
                    A = zeros(nUnits);
                    for i = 1:nUnits
                        A(i,i) = 1;
                    end
                    A(1,1) = -1;
                                        
                    A( Wcc_idxs ) =  .5;
                    A( Bcc_idxs ) =  1;
                    
                    figure(1);
                    imagesc(A); axis equal tight;
                    
%                     SimpleComplexIndices                   
                    
                    
                    
                end
                
            end            
            
        end
        
    else
        [Wrcc_idxs, Wrcm_idxs] = deal({});
    end
    S.Wrcc_idxs = Wrcc_idxs; 
    S.Wrcm_idxs = Wrcm_idxs; 
    S.nBcc_total = nBcc_total;
    S.Wrcc_details = cellinfo;  %#ok<STRNU>
    assert(isfield(cellinfo, 'subRanges'));
        
    [pathname, filename] = getFileName('pairs', cmpOriSpfType_fileext);
    save([pathname, filename], '-struct', 'S', '-v6');
    fprintf(['\nSaved data for ' gratingType ' pairs stimuli to the file ' filename '\n']);
    
    toc;

end

% function u = uniqueCellInts(C)
%     c_max = max(cellfun(@max, C));
%     tf = false(1,c_max);
%     for i = 1:length(C)
%         tf(C{i}) = 1;
%     end
%     u = find(tf);
% end



%                 phCmp = cellfun(@(V) V == phs, phSets);                
%                 if any (phCmp),  
%                     phType = find(phCmp, 1);
%                 else
%                     phSets = [phSets; phs]; 
%                     phType = length(phSets); 
%                 end
%                 if (any(nPhases == 4) && isequal(phs, phsSet4)) || (any(nPhases == 8) && isequal(phs, phsSet8))


% % get sample ori/sp/ph from each group
%     sampleOris = cell(1,length(uniqueGroupIds));
%     sampleSps  = cell(1,length(uniqueGroupIds));
%     samplePhs  = cell(1,length(uniqueGroupIds));
%     for Gid_i = 1:length(uniqueGroupIds)
%         OSPi = allCells(inds{Gid_i}(1));
%         sampleOris{Gid_i} = OSPi.oris;
%         sampleSps{Gid_i}  = OSPi.sps;
%         samplePhs{Gid_i}  = OSPi.phs;
%     end


function pairIdxs = reduceNPairsTo(pairIdxs, n)
    rand('state', 0); %#ok<RAND>
    if length(pairIdxs) > n
        idx = randperm(length(pairIdxs));
        pairIdxs = pairIdxs(idx(1:n));
    end        
end


    

    % make 'contrib' variable
%     nMaxCellsPerSite = max(nCellsEachSite);        
%     nPairsForNCells = [0, 0, arrayfun(@(n) nchoosek(n, 2), 2:nMaxCellsPerSite)]; % at 2 extra 0s, for 0 and 1 cell.
%     nPairsEachSite =  [nPairsForNCells(nCellsEachSite+1)]; % add 1 to avoid an index of 0.
%     contrib = arrayfun(@(nCells, nPairs) ones(1,nPairs)* (2/nCells), nCellsEachSite, nPairsEachSite, 'un', false);
%     contrib = [contrib{:}]';


%         indOfFirstCellInGroups = find(diff([0 cellGroupIds]) ~= 0);
%         nUnitsEachSite = diff([indOfFirstCellInGroups, (nCells+1)]);
%         sitesWithMU = cellIds(indOfFirstCellInGroups) == 0;                

%{
            if strcmp(gratingType, 'flashed')
                oriOK = find( nestedFields(allCells, 'stats', 'tuningStats', 'oriStats_si', 'cellOK') );                
                spfOK = find( nestedFields(allCells, 'stats', 'tuningStats', 'spfStats_si', 'cellOK') );                                                
                
            elseif strcmp(gratingType, 'drifting')
                nSpfs = cellfun(@length, sps_deg);
                isSpfBatch = nSpfs > 1; 
                isOriBatch = nSpfs == 1; 
                rows_arg = {};

                if ~onlyDoDegreeCmpAtSamePenetration
                    % compare only spf-spf (1 vs 1),  and ori-ori (0 vs 0)
                    comparableStimuliSets = isSpfBatch;  
                else
                    % compare only spf-spf,  and ori-ori   *at the same penetration*
                    spfBatchPenIds = isSpfBatch .* groupPenIds;
                    oriBatchPenIds = 1000*isOriBatch .* groupPenIds;                
                    assert( isempty( setdiff( intersect(spfBatchPenIds, oriBatchPenIds), 0) ) ); % don't overlap individually
                    assert( isempty( intersect(intersect(spfBatchPenIds, oriBatchPenIds), spfBatchPenIds+oriBatchPenIds) )); % adding doesn't overlap

                    comparableStimuliSets = spfBatchPenIds + oriBatchPenIds;

                end
            end
%}


%     
%     if strcmp(cmpType, 'degree') 
%         switch gratingType, 
%             case 'flashed',
%                 is_ori = find( nestedFields(allCells, 'stats', 'tuningStats', 'oriStats_si', 'cellOK') );                
%                 is_spf = find( nestedFields(allCells, 'stats', 'tuningStats', 'spfStats_si', 'cellOK') );
%             case 'drifting',                 
%                 is_ori =  find( strncmp({allCells.stimType}, 'Grating:Orientation', 12));
%                 is_spf =  find( strncmp({allCells.stimType}, 'Grating:Spatial Frequency', 12));
%         end
%                 
%         switch cmpOriSpfType
%             case 'ori', select_idx = is_ori;                
%             case 'spf', select_idx = is_spf;
%         end
%         
%         allCells = allCells(select_idx);        
%     end    


%{

nPerLoop = 20;
             nLoops = nResamples / nPerLoop;
             [Bcc_idxs_C, Bcm_idxs_C] = deal(cell(1, nLoops));
             for loop_i = 1:nLoops
                 for i = 1:2
                     [~, ~, Bcc_idxs_range_C{i}, Bcm_idxs_range_C{i}] = generateListOfCellPairs(cellinfo, nPerLoop, subRanges{i}); %#ok<AGROW>
                 end
%                  Wcc_idxs_C{loop_i} = Wcc_idxs_range_C{1};
%                  Wcm_idxs_C{loop_i} = Wcm_idxs_range_C{1};
                 for j = 1:nPerLoop
                     Bcc_idxs_C{loop_i}{j} = setdiff(Bcc_idxs_range_C{1}{j}, Bcc_idxs_range_C{2}{j});
                     Bcm_idxs_C{loop_i}{j} = setdiff(Bcm_idxs_range_C{1}{j}, Bcm_idxs_range_C{2}{j});
                 end
             end
             %%
             Wrcc_idxs = [Bcc_idxs_C{:}];
             Wrcm_idxs = [Bcm_idxs_C{:}];
%}