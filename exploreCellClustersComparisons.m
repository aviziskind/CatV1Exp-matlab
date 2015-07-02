function exploreCellClustersComparisons

%     gratingType = curGratingType('');  % FLASHED_GRATING = 1; DRIFTING_GRATING = 2;
    
    curCmpType('clusters');
    curGroupingType('clustersPruned');
    curGratingType('drifting');    
    curPairTypes('Wcc', 'Bcc');    
        
    ospDatafile = getFileName('osps');    
    pairDatafile= getFileName('pairs');
    cmpDatafile = getFileName('comparisons');

    
%     fprintf('Loading ... '); tic;    
    S1 = load(ospDatafile);        
    allCells = S1.allCells;
    clear S1;
    nUnits = length(allCells);    
    
    S2 = load(pairDatafile);    
    [Wcc_pairIdxs, Wcm_pairIdxs, Bcc_pairIdxs, Bcm_pairIdxs, Wrcc_pairIdxs, Wrcm_pairIdxs] = ...
        deal(S2.Wcc_idxs, S2.Wcm_idxs, S2.Bcc_idxs, S2.Bcm_idxs, S2.Wrcc_idxs, S2.Wrcm_idxs);   
    clear S2;
    
    S3 = load(cmpDatafile); 
    [pairData, S, pairTypes, measures, locations] = deal(S3.pairData, S3.allStatsC, S3.pairTypes, S3.measures, S3.locations);
    clear S3;
%     flds_pd = fieldnames(pairData(1));


    pairGids = cat(1, pairData.Gids);
    pairCellIds = cat(1, pairData.cellIds);
    pairClustIds = cat(1, pairData.clustIds);

    
    Wcc_pairIdxs = sort(Wcc_pairIdxs);
    allPairIdxs = {Wcc_pairIdxs, Wrcc_pairIdxs, Bcc_pairIdxs, Wcm_pairIdxs, Wrcm_pairIdxs, Bcm_pairIdxs};
    allPairIdxs_list = {Wcc_pairIdxs, unique([Wrcc_pairIdxs{:}]), Bcc_pairIdxs, Wcm_pairIdxs, unique([Wrcm_pairIdxs{:}]), Bcm_pairIdxs};
    allPairTypes = {'Wcc', 'Wrcc', 'Bcc',  'Wcm', 'Wrcm', 'Bcm'};
    pairTypes = pairTypes(ord(cellfun(@(s) find(strcmp(s, allPairTypes)), pairTypes)));
    
    minNSpikes = min(cat(2, pairData.nSpikes), [], 2);
    
    3;
            
    pairTypesAvailable = cellfun(@(pr) any ( strcmp(pairTypes, pr)), allPairTypes);
    pairIdxs = cellfun( @(pr) allPairIdxs{find(strcmp(allPairTypes, pr))},  pairTypes, 'un', false);

    pairIdxList = unique(cat(1, allPairIdxs_list{pairTypesAvailable}));
    nPairs = length(pairIdxList);
    
    idxMtx = zeros(nUnits, nUnits, 'uint32');
    idxMtx(pairIdxList) = 1:length(pairIdxList);
    
    4;    
    
    Wcc_idxs = double( sort(idxMtx(Wcc_pairIdxs)) );
    Bcc_idxs = double( sort(idxMtx(Bcc_pairIdxs)) );
    
    assert( all (pairGids(:,1) == pairGids(:,2) ) );
    assert( all (pairCellIds(Wcc_idxs,1) == pairCellIds(Wcc_idxs,2) ) );
    assert( all (pairCellIds(Bcc_idxs,1) ~= pairCellIds(Bcc_idxs,2) ) );

        
    pairPrunePct = [pairData.pairPrunePct];
    pairRefrPct = [pairData.pairRefrPct];
    pairMinIsis = [pairData.pairMinIsis];
    pairRefrPeriod_ms = [pairData.pairRefrPeriod_ms];

    doPCAonCCG = 0;
    if doPCAonCCG
        3;
        allCCGs = cat(1, pairData.pairCCG16ms);
        normFactor = 1./normV(allCCGs, 2); normFactor(isinf(normFactor)) = 0;
        allCCGs_norm = bsxfun(@times, allCCGs, normFactor);
        [CCG_PCA_coeff, CCG_PCA_comps] = doPCA(allCCGs_norm, 16);

        binE = linspace(-16, 16, 80+1); binC = binEdge2cent(binE);
        tic;
        [CCG_GLF_coeff, CCG_GLF_comps] = GLF(allCCGs_norm', 4, 15);
        toc;
        3;
        figure(5); plot(binC, CCG_PCA_comps, '.-'); xlim([-16, 16]); 

        legend(legendarray('PCA ', 1:4));
        title('first 4 PCA components')

        figure(1); clf; 
        plot(CCG_PCA_coeff(1,Bcc_idxs), CCG_PCA_coeff(2,Bcc_idxs), 'r.', 'markersize', 1); hold on;
        plot(CCG_PCA_coeff(1,Wcc_idxs), CCG_PCA_coeff(2,Wcc_idxs), 'b.', 'markersize', 1);
        axis tight; xlabel('PCA 1'); ylabel('PCA 2');

        figure(2); clf; 
        plot(CCG_PCA_coeff(3,Bcc_idxs), CCG_PCA_coeff(4,Bcc_idxs), 'r.', 'markersize', 1); hold on;
        plot(CCG_PCA_coeff(3,Wcc_idxs), CCG_PCA_coeff(4,Wcc_idxs), 'b.', 'markersize', 1);
        axis tight; xlabel('PCA 3'); ylabel('PCA 4');
        figure(6); plot(binC, CCG_GLF_comps, '.-'); xlim([-16, 16]); legend(legendarray('GLF ', 1:4));
        title('first 4 GLF components')

        3;

        figure(3); clf; 
        plot(CCG_GLF_coeff(1,Bcc_idxs), CCG_GLF_coeff(2,Bcc_idxs), 'r.', 'markersize', 1); hold on;
        plot(CCG_GLF_coeff(1,Wcc_idxs), CCG_GLF_coeff(2,Wcc_idxs), 'b.', 'markersize', 1);
        axis tight; xlabel('GLF 1'); ylabel('GLF 2');

        figure(4); clf; 
        plot(CCG_GLF_coeff(3,Bcc_idxs), CCG_GLF_coeff(4,Bcc_idxs), 'r.', 'markersize', 1); hold on;
        plot(CCG_GLF_coeff(3,Wcc_idxs), CCG_GLF_coeff(4,Wcc_idxs), 'b.', 'markersize', 1);
        axis tight; xlabel('GLF 3'); ylabel('GLF 4');
    end
    
%     allFeatureDists = [pairData.featureDists];
%     selectIdx = find( (pairMinIsis < .5) & ([pairData.fullWaveform_cc] < .3) & (pairRefrPeriod_ms > 1.1));
%     selectIdx = find( minNSpikes > 1000 );
    
%     selectIdx = find( [pairData.fullWaveform_cc] > .1 );
    selectIdx = [];
    
%     allClustIds = cat(1, pairData.clustIds);
%     allClustIds = sort(allClustIds, 2, 'ascend');

%     tf_28 = (allClustIds(:,1) == 2) & (allClustIds(:,2) == 8);
%     tf_82 = (allClustIds(:,1) == 8) & (allClustIds(:,2) == 2);
%     
%     selectIdx = find( tf_28 | tf_82 );
    
    
    if ~isempty(selectIdx)
        Bcc_idxs_select = intersectq(selectIdx, Bcc_idxs);
        Wcc_idxs_select = intersectq(selectIdx, Wcc_idxs);
    else
        Bcc_idxs_select = Bcc_idxs;
        Wcc_idxs_select = Wcc_idxs;
    end
    
%     allPCA = log10([allFeatureDists.PCA_overlap]);
%     allGLF = log10([allFeatureDists.GLF_overlap]);
%     displayNHist({allPCA, allGLF})
%     displayNHist({allPCA(Wcc_idxs), allPCA(Bcc_idxs)})
%     displayNHist({allGLF(Wcc_idxs), allGLF(Bcc_idxs)})
    
    3;
%     idxL
    
    
%     XXfldName = 'pairRefrPct';   XXname = '# refractory spikes';
%     YYfldName = 'pairPrunePct';  YYname = '# pruned spikes';

%     XXfldName = 'negAmps_dist';  XXname = 'negAmps_dist';
%     YYfldName = 'fullWaveform_cc';  YYname = 'Waveform CC';

    XXfldName = 'fullWaveform_cc';
    YYfldName = 'negAmps_cc';

%     XXfldName = 'pairCcgCenterRatio4';
%     YYfldName = 'pairCcgLrRatio16';
    
    XXname = XXfldName;
    YYname = YYfldName;


    XX = [pairData.(XXfldName)];
    YY = [pairData.(YYfldName)];  

    showBCCs = 0;
    showWCCs = 1;
    
    figId = 1;
    figure(figId); clf; hold on;
    if showBCCs
        plot(XX(Bcc_idxs_select), YY(Bcc_idxs_select), 'ro', 'markersize', 2, 'markerfacecolor', 'none'); 
    end
    if showWCCs
        plot(XX(Wcc_idxs_select), YY(Wcc_idxs_select), 'b+', 'markersize', 2, 'markerfacecolor', 'none'); 
    end
    
%     set(gca, 'xscale', 'log', 'yscale', 'log')
    h_tit = title('');
    set(figId, 'windowButtonDownFcn', {@selectPointsInFigure, @updateTitleWithId});
    xlabel(XXname); ylabel(YYname)
    3;
    
   
    
    
    doPruneVsRefr = 0;
    if doPruneVsRefr    
        
        idxuse = find(pairPrunePct > 0 & pairRefrPct > 0);
        prune_samp = pairPrunePct(idxuse);
        refr_samp = pairRefrPct(idxuse);
        
        cmd = 'cum'; mkr = '.'; Nbin = 100;
        lims_prune = log10(lims(prune_samp, 0));  prune_binE = linspace(lims_prune(1), lims_prune(2), Nbin);
        lims_refr = log10(lims(refr_samp, 0));    refr_binE = linspace(lims_refr(1), lims_refr(2), Nbin);
        N_prune_Bcc = histcnt(log10(pairPrunePct(Bcc_idxs)), prune_binE);  N_prune_Bcc = ff(N_prune_Bcc, cmd);        
        N_prune_Wcc = histcnt(log10(pairPrunePct(Wcc_idxs)), prune_binE);  N_prune_Wcc = ff(N_prune_Wcc, cmd);

        N_refr_Bcc = histcnt(log10(pairRefrPct(Bcc_idxs)), refr_binE); N_refr_Bcc = ff(N_refr_Bcc, cmd);
        N_refr_Wcc = histcnt(log10(pairRefrPct(Wcc_idxs)), refr_binE); N_refr_Wcc = ff(N_refr_Wcc, cmd);

        figure(2);
        prune_binC = binEdge2cent(prune_binE);
        d_prune = max(abs(N_prune_Bcc-N_prune_Wcc));
        plot(prune_binC, N_prune_Bcc, ['r' mkr '-'], prune_binC, N_prune_Wcc, ['b' mkr '-']); 
        title(sprintf('prune (dmax = %.2f)', d_prune));


        figure(3);
        refr_binC = binEdge2cent(refr_binE);
        d_refr = max(abs(N_refr_Bcc-N_refr_Wcc));
        plot(refr_binC, N_refr_Bcc, ['r' mkr '-'], refr_binC, N_refr_Wcc, ['b' mkr '-']); title('refr');
        title(sprintf('refr (dmax = %.2f)', d_refr));


        3;
    %     N_prune_Wcc = histcnt(log10(pairRefrPct(Wcc_idxs)), linspace(lims_prune(1), lims_prune(2), 30));
        3;
    end    
    
    function updateTitleWithId(glob_id, grp_id, loc_id, pt_x, pt_y) %#ok<INUSL>
        
        isBCC = showBCCs && (grp_id==1);
        isWCC = (~showBCCs && (grp_id==1)) || (showBCCs && (grp_id==2));
        assert(xor(isBCC, isWCC));
        
        if isBCC
            pd = pairData(Bcc_idxs_select(loc_id));                     
            assert(pd.cellIds(1) ~= pd.cellIds(2));
        elseif isWCC
            pd = pairData(Wcc_idxs_select(loc_id));            
            assert(pd.cellIds(1) == pd.cellIds(2));
        end
        assert( (pd.(XXfldName) == pt_x) && (pd.(YYfldName) == pt_y) );
        
        refrS = iff(pt_x == round(pt_x), '%d', '%.2f');
        pruneS = iff(pt_y == round(pt_y), '%d', '%.2f');
        
        Gid = pd.Gids(1);
        cellIds = pd.cellIds;
        clustIds = pd.clustIds;
        
        if isWCC
            s = sprintf(['Gid = %d, (cellId = %d), clustIds %d & %d. [%s = ' refrS ', %s = ' pruneS ' ]'], Gid, cellIds(1), clustIds, XXname, pt_x, YYname, pt_y);
        elseif isBCC
            s = sprintf(['Gid = %d, clustIds %d & %d [%s = ' refrS ', %s = ' pruneS ' ]'], Gid, clustIds, XXname, pt_x, YYname, pt_y);
        end      
        compareWvfms(Gid, clustIds(1), clustIds(2))
        
        fprintf([s '\n']);
        set(h_tit, 'string', s);
        3;

    end

    return;
    
    
    
%     allPairIdxs = {Wcc_pairIdxs, Wrcc_pairIdxs, Bcc_pairIdxs, Wcm_pairIdxs, Wrcm_pairIdxs, Bcm_pairIdxs};
%     allPairIdxs_list = {Wcc_pairIdxs, unique([Wrcc_pairIdxs{:}]), Bcc_pairIdxs, Wcm_pairIdxs, unique([Wrcm_pairIdxs{:}]), Bcm_pairIdxs};
%     allPairTypes = {'Wcc', 'Wrcc', 'Bcc',  'Wcm', 'Wrcm', 'Bcm'};
%     pairTypes = pairTypes(ord(cellfun(@(s) find(strcmp(s, allPairTypes)), pairTypes)));
    
    % 
    Wcc_pairs = ind2subV([nUnits, nUnits], Wcc_pairIdxs);
    Bcc_pairs = ind2subV([nUnits, nUnits], Bcc_pairIdxs);
    
    % Get Ori_si, Ori_ss, and Spf_si stats
    oriCells_tf = strncmp({allCells.stimType}, 'Grating:Orientation', 19);
    oriCells_idx = find(oriCells_tf);
    spfCells_idx = find(~oriCells_tf);

    allOriCells = allCells( oriCells_tf);
    allSpfCells = allCells(~oriCells_tf);    
    
    allOriStats_si = structNestedFields(allOriCells, {'stats', 'tuningStats', 'oriStats_si'});
    allOriStats_ss = structNestedFields(allOriCells, {'stats', 'tuningStats', 'oriStats_ss'});
    allSpfStats_si = structNestedFields(allSpfCells, {'stats', 'tuningStats', 'spfStats_si'});
    
    % Renumber Wcc & Bcc pairs to index the subsets of only orientation batch cells, or only spatial-freq
    % batch cells
    [Wcc_pairs_ori, Wcc_oo_pairIdxs] = renumberSubsetOfCellPairs(Wcc_pairs, Wcc_pairIdxs, oriCells_idx, idxMtx);
    [Wcc_pairs_spf, Wcc_ss_pairIdxs] = renumberSubsetOfCellPairs(Wcc_pairs, Wcc_pairIdxs, spfCells_idx, idxMtx);
              
    [Bcc_pairs_ori, Bcc_oo_pairIdxs] = renumberSubsetOfCellPairs(Bcc_pairs, Bcc_pairIdxs, oriCells_idx, idxMtx);
    [Bcc_pairs_spf, Bcc_ss_pairIdxs] = renumberSubsetOfCellPairs(Bcc_pairs, Bcc_pairIdxs, spfCells_idx, idxMtx);
    
    allOriGids = [allOriCells.Gid];
    allOriCellIds = [allOriCells.cellId];
    allSpfGids = [allSpfCells.Gid];
    allSpfCellIds = [allSpfCells.cellId];

    % Double check that the renumbering is all ok.
    checkRenumbering = 1;
    if checkRenumbering
        assert( all(allOriGids(Wcc_pairs_ori(:,1)) == allOriGids(Wcc_pairs_ori(:,2))));
        assert( all(allSpfGids(Wcc_pairs_spf(:,1)) == allSpfGids(Wcc_pairs_spf(:,2))));    

        assert( all(allOriGids(Bcc_pairs_ori(:,1)) ~= allOriGids(Bcc_pairs_ori(:,2))));
        assert( all(allSpfGids(Bcc_pairs_spf(:,1)) ~= allSpfGids(Bcc_pairs_spf(:,2))));    
        
        
%         pairGids    = cat(1, pairData.Gids);
%         pairCellIds = cat(1, pairData.cellIds);
%         
%         Wcc_oo_pairIdxs
%         
%         pairGids(
        
    end    
    
    
    
    
    3;
    
    


end


function y = ff(x, cmd)
    switch cmd
        case 'cum', y = cumsum(x); y = y / y(end);
        case 'norm', y = x/sum(x);
    end
end
            

