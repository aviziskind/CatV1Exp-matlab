function printDegreeOfTuningComparisons
    doFigs = [1:11];
%     doFigs = [14];

    doPlots = 1             && 1;
    doPrintStats = 1        && 1;
    doPrintSigTests = 1     && 1;
    dispPlotStats = 0;

%     differentiateSimpleComplexCellsInStats = 1;
%     onlyCompareSpfAtSamePenetration = curSamePen;
    filename_ext = '';
    

    doOri = 1;
    doSpf = 1;    
    
    saveControlDistribs = 0;
    addProbSubplots = 1;
    addFlashedOrDriftingToFigTitles = 1;
    curCmpType('degree');
    
    addLegend_drifting = 1;
    addLegend_flashed = 0;
    
    title_fsize = 12;
    
    gratingType = curGratingType('');
%     if doPrintSigTests || doPrintStats
        curPairTypes('Wcc', 'Wrcc', 'Bcc');            
%     elseif doPrintStats
%         curPairTypes('Wcc', 'Bcc');            
%     end
%     curPairTypes('Wcc', 'Wrcc', 'Bcc');            
    
    statsDatafile = getFileName('controls', filename_ext);
        
%     fprintf('Loading ... '); tic;    
    % load ORI data

%     function [allUnits, pairData_os, S_os, pairTypes, measures_os, pairIdxs_os, pairIdxList_os, idxMtx_os, ...
%     ]    = loadDegreeOfTuningData(osType)
%         
% 
% 
%         
%     end
%     
%     
%     [pairData_ori, S_ori, pairTypes, measures_ori, pairIdxs_ori, pairIdxList_ori, idxMtx_ori,...
%         = loadDegreeOfTuningData('ori')
%  
    
    % load ORI data
    S1_ori = load(getFileName('osps', '_ori'));        
    allOriUnits = S1_ori.allCells;
    clear S1;
    nUnits_ori = length(allOriUnits);    
    
    S2_ori = load(getFileName('pairs', '_ori'));            
    S3_ori = load(getFileName('comparisons', '_ori')); 
    [pairData_ori, S_ori, pairTypes, measures_ori] = deal(S3_ori.pairData, S3_ori.allStatsC, S3_ori.pairTypes, S3_ori.measures);    
    [pairIdxs_ori, pairIdxList_ori, idxMtx_ori] = useSubsetOfPairIdxs(S2_ori, pairTypes, nUnits_ori);
    clear S2;
    clear S3;
        
        
    
    % load SPF data
    S1_spf = load(getFileName('osps', '_spf'));        
    allSpfUnits = S1_spf.allCells;
    clear S1;
    nUnits_spf = length(allSpfUnits);    
    
    S2_spf = load(getFileName('pairs', '_spf'));            
    S3_spf = load(getFileName('comparisons', '_spf')); 
    [pairData_spf, S_spf, pairTypes, measures_spf] = deal(S3_spf.pairData, S3_spf.allStatsC, S3_spf.pairTypes, S3_spf.measures);    
    [pairIdxs_spf, pairIdxList_spf, idxMtx_spf] = useSubsetOfPairIdxs(S2_spf, pairTypes, nUnits_spf);
    clear S2;
    clear S3;

    measures = uniqueInOrder([measures_ori, measures_spf]);
    
    
    limitToNPermutes = 100;
    
    % ori reindexing
    Wcc_oo_pairIdxs_M = pairIdxs_ori{ find(strcmp(pairTypes, 'Wcc'), 1) };
    Wcc_oo_pairIdxs = idxMtx_ori( Wcc_oo_pairIdxs_M );
    Wcc_pairs_ori = ind2subV([nUnits_ori, nUnits_ori], Wcc_oo_pairIdxs_M);
    if any(strcmp(pairTypes, 'Wrcc'))
        if ~isempty(limitToNPermutes)
            i_wrcc = find(strcmp(pairTypes, 'Wrcc'), 1);
            pairIdxs_ori{ i_wrcc } = pairIdxs_ori { i_wrcc }(1:limitToNPermutes);
        end

        Wrcc_oo_pairIdxs_M = pairIdxs_ori{ find(strcmp(pairTypes, 'Wrcc'), 1) };        
        Wrcc_oo_pairIdxs =  cellfun(@(idxs_M) idxMtx_ori( idxs_M ), Wrcc_oo_pairIdxs_M, 'un', 0);                
        Wrcc_pairs_ori = cellfun(@(idxs) ind2subV([nUnits_ori, nUnits_ori], idxs), Wrcc_oo_pairIdxs_M, 'un', 0);                
    end
    Bcc_oo_pairIdxs_M  = pairIdxs_ori{ find(strcmp(pairTypes, 'Bcc'), 1) };
    Bcc_oo_pairIdxs = idxMtx_ori( Bcc_oo_pairIdxs_M );
    Bcc_pairs_ori = ind2subV([nUnits_ori, nUnits_ori], Bcc_oo_pairIdxs_M);

        
    % spf reindexing
    Wcc_ss_pairIdxs_M = pairIdxs_spf{ find(strcmp(pairTypes, 'Wcc'), 1) };
    Wcc_ss_pairIdxs = idxMtx_spf( Wcc_ss_pairIdxs_M );
    Wcc_pairs_spf = ind2subV([nUnits_spf, nUnits_spf], Wcc_ss_pairIdxs_M);

    if any(strcmp(pairTypes, 'Wrcc'))
        if ~isempty(limitToNPermutes)
            i_wrcc = find(strcmp(pairTypes, 'Wrcc'), 1);
            pairIdxs_spf{ i_wrcc } = pairIdxs_spf{ i_wrcc }(1:limitToNPermutes);
        end        
        
        Wrcc_ss_pairIdxs_M = pairIdxs_spf{ find(strcmp(pairTypes, 'Wrcc'), 1) };        
        Wrcc_ss_pairIdxs =  cellfun(@(idxs_M) idxMtx_spf( idxs_M ), Wrcc_ss_pairIdxs_M, 'un', 0);                
        Wrcc_pairs_spf = cellfun(@(idxs) ind2subV([nUnits_spf, nUnits_spf], idxs), Wrcc_ss_pairIdxs_M, 'un', 0);                
    end
    Bcc_ss_pairIdxs_M  = pairIdxs_spf{ find(strcmp(pairTypes, 'Bcc'), 1) };
    Bcc_ss_pairIdxs = idxMtx_spf( Bcc_ss_pairIdxs_M );
    Bcc_pairs_spf = ind2subV([nUnits_spf, nUnits_spf], Bcc_ss_pairIdxs_M);
    
    
    % Choose cells (cellID > 0) that are OK
    oriIsCell = [allOriUnits.cellId] > 0;
    spfIsCell = [allSpfUnits.cellId] > 0;
    if strcmp(gratingType, 'flashed')                
        oriCells_use = oriIsCell & nestedFields(allOriUnits, 'stats', 'tuningStats', 'oriStats_si', 'cellOK');
        spfCells_use = spfIsCell & nestedFields(allSpfUnits, 'stats', 'tuningStats', 'spfStats_si', 'cellOK');                        
    elseif strcmp(gratingType, 'drifting')
        oriCells_use = oriIsCell & strncmp({allCells.stimType}, 'Grating:Orientation', 19);
        spfCells_use = spfIsCell & strncmp({allCells.stimType}, 'Grating:Spatial Freq', 19);
    end
        
    if doOri
        oriCells_idx = find(oriCells_use);
        allOriCells = allOriUnits(oriCells_idx);
                
        allOriStats_si = nestedFields(allOriUnits, {'stats', 'tuningStats', 'oriStats_si'});
        allOriSpkFeatures = [allOriUnits.spkFeatures];
        
        doSpontSubtractedMeasures = any(cellfun(@(s) ~isempty(strfind(s, '_ss')), measures)) && ...
            isfield(allOriCells(1).stats.tuningStats, 'oriStats_ss');

        if doSpontSubtractedMeasures
            allOriStats_ss = nestedFields(allOriCells, {'stats', 'tuningStats', 'oriStats_ss'});
        end            
        assert(all([allOriStats_si.cellOK]));        
    else
        doSpontSubtractedMeasures = 0;
    end
    
    if doSpf
        spfCells_idx = find(spfCells_use);
        allSpfCells = allSpfUnits( spfCells_idx );    
        allSpfStats_si = nestedFields(allSpfUnits, {'stats', 'tuningStats', 'spfStats_si'});
    
        assert(all([allSpfStats_si.cellOK]));
        
    end
    
%     if doOri && doSpf
%         ori_spf_cells_idx = find(oriCells_use & spfCells_use);
%         oriSpfCells = allCells(ori_spf_cells_idx);
%         
%     end
        
    
        
    if doOri                
        allOriGids = [allOriUnits.Gid];
        ori_cellF1oDC = [allOriStats_si.F1oDC];
        
        allOriUnitStats_si = nestedFields(allOriUnits, {'stats', 'tuningStats', 'oriStats_si'});
        ori_unitF1oDC = [allOriUnitStats_si.F1oDC];
        
        ori_pairF1oDC_Wcc = ori_unitF1oDC(Wcc_pairs_ori);       
        ori_pairF1oDC_Bcc = ori_unitF1oDC(Bcc_pairs_ori);       
    else
        [Wcc_oo_pairIdxs, Wrcc_oo_pairIdxs, Bcc_oo_pairIdxs, ori_cellF1oDC, ori_pairF1oDC_Wcc] = deal([]);
    end

    if doSpf
        
        allSpfGids = [allSpfUnits.Gid];
        spf_cellF1oDC = [allSpfStats_si.F1oDC];
        
        allSpfUnitStats_si = nestedFields(allSpfUnits, {'stats', 'tuningStats', 'spfStats_si'});
        spf_unitF1oDC = [allSpfUnitStats_si.F1oDC];
        
        spf_pairF1oDC_Wcc = spf_unitF1oDC(Wcc_pairs_spf);       
        spf_pairF1oDC_Bcc = spf_unitF1oDC(Bcc_pairs_spf);     
        
        
%         allSpfGids = [allSpfCells.Gid];
%         spf_cellF1oDC = [allSpfStats_si.F1oDC];
%         spf_pairF1oDC_Wcc = spf_cellF1oDC(Wcc_pairs_spf);
%         spf_pairF1oDC_Bcc = spf_cellF1oDC(Bcc_pairs_spf);
%         
%         spf_nSimp_Wcc = sum(spf_pairF1oDC_Wcc > 1, 2);
%         spf_nSimp_Bcc = sum(spf_pairF1oDC_Bcc > 1, 2);
    else
        [Wcc_ss_pairIdxs, Wrcc_ss_pairIdxs, Bcc_ss_pairIdxs, spf_cellF1oDC, spf_pairF1oDC_Wcc] = deal([]);
    end
    
        
    % get ori outlier idxs
    if doOri
    
        ori_is_outlier = [allOriStats_si.Dori_pref_allSpkMU] > 45;
        ori_is_outlier_idx = find(ori_is_outlier);
    %     oriNorm_idx    = oriCells_idx(~ori_is_outlier);
    %     oriOutlier_idx = oriCells_idx(ori_is_outlier);    

        pair_is_outlier = binarySearch(ori_is_outlier_idx, Wcc_pairs_ori, [], 'exact');
        idx_0outliers =  ~any(pair_is_outlier, 2) ;
        idx_1outlier  =   xor(pair_is_outlier(:,1), pair_is_outlier(:,2)) ;            
        idx_2outliers =   all(pair_is_outlier, 2) ;


        Wcc_oo_pairIdxs_norm      = Wcc_oo_pairIdxs(idx_0outliers);  % Wcc_pairs_ori_norm      = Wcc_pairs_ori(idx_0outliers,:); 
        Wcc_oo_pairIdxs_1outlier  = Wcc_oo_pairIdxs(idx_1outlier);   % Wcc_pairs_ori_1outlier  = Wcc_pairs_ori(idx_1outlier,:);  
        Wcc_oo_pairIdxs_2outliers = Wcc_oo_pairIdxs(idx_2outliers);  % Wcc_pairs_ori_2outliers = Wcc_pairs_ori(idx_2outliers,:);      
    end    
        
        
    
%     allOriCellIds = [allOriCells.cellId];
%     allSpfCellIds = [allSpfCells.cellId];

    % Double check that the renumbering is all ok.
    checkRenumbering = 1;
    if checkRenumbering
        if doOri
            assert( all(allOriGids(Wcc_pairs_ori(:,1)) == allOriGids(Wcc_pairs_ori(:,2))));
            assert( all(allOriGids(Bcc_pairs_ori(:,1)) ~= allOriGids(Bcc_pairs_ori(:,2))));
        end        
        if doSpf
            assert( all(allSpfGids(Wcc_pairs_spf(:,1)) == allSpfGids(Wcc_pairs_spf(:,2))));    
            assert( all(allSpfGids(Bcc_pairs_spf(:,1)) ~= allSpfGids(Bcc_pairs_spf(:,2))));                    
        end
    end    
            
        
    if doPrintStats
        fprintf('********************* TABLE 2: STATISTICS FOR PAIRWISE DIFFERENCES **********************\n');
        fprintf('    Parameter       |    Mean     |      Std    |    Median   |     P25     |     P75      | N\n');        
        ori_pairTypeIdxs = {Wcc_oo_pairIdxs, Bcc_oo_pairIdxs};        
        spf_pairTypeIdxs = {Wcc_ss_pairIdxs, Bcc_ss_pairIdxs};        
    
%         measures_all = {'Dw_ori_glob_si', 'Dw_ori_glob_ss', 'Dw_ori_loc_si'  'Dw_ori_loc_ss', 'D_ori_pref' ...
%                         'D_dsi_si', 'D_dsi_ss', 'D_dir_pref', ...
%                         'Dw_spf'    'D_spf_pref' };        
        measures_all = measures;

        idx_Dori = find(strcmp(measures_all, 'D_ori_pref'), 1);        
        measures_all = [measures_all(1:idx_Dori), {'D_ori_pref*', 'D_ori_pref_MU', 'D_ori_pref_MU*'}, measures_all(idx_Dori+1:end)];
                
        if strcmp(gratingType, 'drifting')        
            idx_Ddir = find(strcmp(measures_all, 'D_dir_pref'), 1);
            measures_all = [measures_all(1:idx_Ddir), {'D_dir_pref_MU'}, measures_all(idx_Ddir+1:end)];
        end
        
        for i = 1:length(measures_all)            
            
            measure_si = measures_all{i};            
            measure_si_noStar = strrep(measure_si, '*', '');
            isOriMeasure = any(strcmp(measure_si_noStar, measures_ori));                                        
            
            if strcmp(measure_si(end-2:end), '_ss')
                continue;
            end            
            if strcmp(measure_si(end-2:end), '_si')
                measure_ss = strrep(measure_si, '_si', '_ss');                
            else
                measure_ss = '';
            end
            measure_full = strrep(measure_si, '_si', '');
            measure = measure_full;
            
            mu_measure = ~isempty(strfind(measure, 'MU'));
            measure = strrep(measure, '_MU', '');
            
            
            measures_ofSameType = iff(isOriMeasure,measures_ori, measures_spf);            
            
            idx_si = find(strcmp(measures_ofSameType, strrep(measure_si, '*', '') ));
            idx_ss = find(strcmp(measures_ofSameType, strrep(measure_ss, '*', '') ));
            spont_idxs = unique([idx_si, idx_ss]);
            pairTypes_str = {'Wcc', 'Bcc'};            
            
            
            removeOutliers = ~isempty(strfind(measure, '*'));            
            
            if isOriMeasure && ~doOri
                continue;
            end
            if ~isOriMeasure && ~doSpf
                continue;
            end
            
            if isOriMeasure
                pairTypeIdxs = ori_pairTypeIdxs;
                S = S_ori;
                pairData = pairData_ori;
                if removeOutliers
                    pairTypeIdxs{1} = Wcc_oo_pairIdxs_norm;
                end
                cellF1oDCs = ori_cellF1oDC;
                pairF1oDCs = ori_pairF1oDC_Wcc;
            else
                S = S_spf;
                pairData = pairData_spf;
                pairTypeIdxs = spf_pairTypeIdxs;
                cellF1oDCs = spf_cellF1oDC;
                pairF1oDCs = spf_pairF1oDC_Wcc;
            end
            
            
            isOriDirPref = ~isempty(strfind(measure, 'ori_pref')) || ~isempty(strfind(measure, 'dir_pref'));            
            
            if isOriDirPref
                pairTypes_str = {'Wcc'}; % don't bother with 'Bcc'
            end
            
            if mu_measure
                pairTypes_str = {'Wcm'};
                spont_idxs = nan;
            end    
            if removeOutliers
                measure = strrep(measure, '*', '');
            end
            nSpont = length(spont_idxs);
            
            switch measure
                case 'Dw_ori_glob', fprintf('***** ORIENTATION ******\n');
                case 'D_dsi',       fprintf('***** DIRECTION  ******\n');
                case 'Dw_spf',      
                    fprintf('***** SPATIAL FREQUENCY ******\n');
                    2;
                case 'dR_spont_abs',
                    fprintf('***** RESPONSES AT NULL ORIENTATION ******\n');
            end
                    
            
            for pt_i = 1:length(pairTypes_str)
                [vals_mean, vals_std, vals_median, vals_P25, vals_P75, Npr, Ncl] = deal(cell(1, nSpont));
                                    
                
                for spont_i = 1:nSpont
                    if ~mu_measure % Cell-Cell measures
                        vals = S{ spont_idxs(spont_i) }.val(  pairTypeIdxs{pt_i}   );                                                
                    else
                        % Cell-MU measures:
                        switch measure 
                            case 'D_ori_pref', vals = [allOriStats_si.Dori_pref_allSpkMU];
                            case 'D_dir_pref', vals = [allOriStats_si.Ddir_pref_allSpkMU];
                            otherwise, error('!');
                        end
                        if removeOutliers
                           vals = vals(~ori_is_outlier); 
                        end
                    end
                    
                    if (pt_i == 1) && ~mu_measure % only for within-site, cell-cell
                        idx_nonnans = ~isnan(vals);
                        gids = pairData.Gids(pairTypeIdxs{pt_i}(idx_nonnans),:);
                        cids = pairData.cellIds(pairTypeIdxs{pt_i}(idx_nonnans),:);                        
                        allGC = unique([gids(:), cids(:)], 'rows');
                        Ncl{spont_i} = size(allGC,1);
                    end
                    
                    if mu_measure % cell-multiunit                        
                        Ncl{spont_i} = length( vals);
                    end

                    [vals_mean{spont_i}, vals_std{spont_i}, vals_median{spont_i}, ...
                        vals_P25{spont_i}, vals_P75{spont_i}, Npr{spont_i}] = getMeanStdPctile(vals);  
                end                        
                
                w = num2str( switchh( nSpont, [1 2], [11, 5]));
                if ~isempty(strfind(measure_si, 'ori')) || ~isempty(strfind(measure_si, 'dir'))
                    fmt = ['%' w '.1f'];
                elseif ~isempty(strfind(measure_si, 'dsi')) || ~isempty(strfind(measure_si, 'spf'))
                    fmt = ['%' w '.2f'];
                end
                
                if nSpont == 1
                    pairTmp_str = ['| ' fmt ' '];
                    num_tmp = '%d';
                else
                    pairTmp_str = ['| ' fmt '/' fmt ' '];
                    num_tmp = '%d/%d';
                end
                nPairs_str = sprintf([num_tmp ' Pr; '], Npr{:});
                
                if pt_i == 1
                    nCells_str = [sprintf(num_tmp, Ncl{:}) ' Cl; '];
                else
                    nCells_str = '';
                end
                   
                if mu_measure
                    nMU_str = sprintf('%d MU', length(unique(allOriGids)) );
                    nPairs_str = '';
                else
                    nMU_str = '';
                end
                nCells_str = [nCells_str nMU_str];  %#ok<AGROW>
                
                if removeOutliers
                    backspace(1); % remove the <enter> to join with the previous paragraph
                end
                
                str_template = ['%14s (' pairTypes_str{pt_i} ')' repmat( pairTmp_str, 1, 5) ' | ' nPairs_str '%s \n'];        

                fprintf(str_template, ...
                    measure_full, vals_mean{:}, vals_std{:}, vals_median{:}, vals_P25{:}, vals_P75{:}, nCells_str);
            end        
            fprintf('\n');        
        end
        fprintf('*******************************************************************************************\n\n\n');
        3;
        
        
    end
    
    if doPrintSigTests                        

        Npermutes = max(length(Wrcc_oo_pairIdxs), length(Wrcc_ss_pairIdxs));
        fprintf('\n\n ************************** SIGNIFICANCE TESTS (Npermute = %d) ********************************* \n', Npermutes);
        fprintf('  Parameter    | Median Ratio| Median Prob |  Mean Ratio |  Mean Prob  |   KS_stat   |   KS prob   |      KS prob(indep)    |\n') 
        
        ori_pairTypeIdxs = {Wcc_oo_pairIdxs, Wrcc_oo_pairIdxs, Bcc_oo_pairIdxs};
        spf_pairTypeIdxs = {Wcc_ss_pairIdxs, Wrcc_ss_pairIdxs, Bcc_ss_pairIdxs};        
    
        measures_all = measures;        
        idx_Dori = find(strcmp(measures_all, 'D_ori_pref'), 1);        
        measures_all = [measures_all(1:idx_Dori), {'D_ori_pref*'}, measures_all(idx_Dori+1:end)];
        
%         measures_all = {'D_F1pair'};
        
        for i = 1:length(measures_all)
                        
            
            measure_si = measures_all{i};                                    
            if strcmp(measure_si(end-2:end), '_ss')
                continue;
            end            
            if strcmp(measure_si(end-2:end), '_si')
                measure_ss = strrep(measure_si, '_si', '_ss');                
            else
                measure_ss = '';
            end
            measure  = strrep(measure_si, '_si', '');
            
            if strcmp(measure, 'D_F1pair')
                continue
                3;
            end
            
            isOriMeasure = isempty(strfind(measure, 'spf')) && isempty(strfind(measure, 'F1oDC'));
            
            if isOriMeasure && ~doOri
                continue;
            end
            if ~isOriMeasure && ~doSpf
                continue;
            end
%             isOriMeasure = isempty(strfind(measure, 'spf'));
            removeOutliers = ~isempty(strfind(measure, '*'));            
            

            measures_ofSameType = iff(isOriMeasure,measures_ori, measures_spf);                        
            idx_si = find(strcmp(measures_ofSameType, strrep(measure_si, '*', '') ));
            idx_ss = find(strcmp(measures_ofSameType, strrep(measure_ss, '*', '') ));
            spont_idxs = unique([idx_si, idx_ss]);
                        
            nSpont = length(spont_idxs);
%             pairTypes_str = {'Wcc', 'Bcc'};    
            pairTypesHere = {'Wcc', 'Bcc', 'Wrcc'};            
            Wcc_idx = find(strcmp(pairTypes, 'Wcc'), 1);
            Bcc_idx = find(strcmp(pairTypes, 'Bcc'), 1);
            Wrcc_idx = find(strcmp(pairTypes, 'Wrcc'), 1);
            
            isOriMeasure = isempty(strfind(measure, 'spf')) && isempty(strfind(measure, 'F1'));
            if isOriMeasure
                pairTypeIdxs = ori_pairTypeIdxs;
                if removeOutliers
                    pairTypeIdxs{1} = Wcc_oo_pairIdxs_norm;
                end
                S = S_ori;
            else
                pairTypeIdxs = spf_pairTypeIdxs;
                S = S_spf;
            end
            
            [vals_mean, vals_median, vals_KS, N] = deal(cell(length(pairTypes), nSpont));
            [medianRatio, medianProb, meanRatio, meanProb, ...
                ksStat, ksProb, ksProb2, vals_Bcc] = deal( cell(1, nSpont) );
            
            %%
            if 0 %&& strcmp(measures_all, 'D_F1pair')
                Frac_ss = @(x) nnz(x == 2)/length(x);
                Frac_sc = @(x) nnz(x == 1)/length(x);
                Frac_cc = @(x) nnz(x == 0)/length(x);
                for pt_i = 1:length(pairTypes)
                    allVals = S{ spont_idxs(1) }.val;                                        
                    pairIdxs = pairTypeIdxs{pt_i};
                    if ~iscell(pairIdxs)
%                         vals = allVals ( pairIdxs );
                        frac_ss{pt_i} = Frac_ss(allVals ( pairIdxs ));
                        frac_sc{pt_i} = Frac_sc(allVals ( pairIdxs ));
                        frac_cc{pt_i} = Frac_cc(allVals ( pairIdxs ));
                    else
                        frac_ss{pt_i} = cellfun(@(idxs) Frac_ss( allVals ( idxs ) ), pairIdxs );
                        frac_sc{pt_i} = cellfun(@(idxs) Frac_sc( allVals ( idxs ) ), pairIdxs );
                        frac_cc{pt_i} = cellfun(@(idxs) Frac_cc( allVals ( idxs ) ), pairIdxs );
                    end
                    
                end
                
                
                
            end                
              %%  

            for spont_i = 1:nSpont
                for pt_i = 1:length(pairTypes)                                                        
                    
                    allVals = S{ spont_idxs(spont_i) }.val;
                    ctrl_dist = getCtrlDist(allVals);
                    
                    pairIdxs = pairTypeIdxs{pt_i};
                    if ~iscell(pairIdxs)
                        vals = allVals ( pairIdxs );                        
                        [vals_mean{pt_i, spont_i}, vals_median{pt_i, spont_i}, vals_KS{pt_i, spont_i}, N{pt_i, spont_i}] = getMeanMedianKS(vals, ctrl_dist);  
                    else
                        
                        [vals_mean{pt_i, spont_i}, vals_median{pt_i, spont_i}, vals_KS{pt_i, spont_i}, N{pt_i, spont_i}] = ...
                        cellfun(@(idxs) getMeanMedianKS( allVals ( idxs ), ctrl_dist), pairIdxs ) ;                                              
                    end
                                        
                end
                
                %%% Median Data 
                Med_wcc = vals_median{Wcc_idx, spont_i};
                Med_bcc = vals_median{Bcc_idx, spont_i};
                Med_permute = vals_median{Wrcc_idx, spont_i};
                medianRatio{spont_i} = Med_bcc/Med_wcc;
                medianProb{spont_i}  = getRandomizedProb(Med_wcc, Med_permute, 'left');
                                                                                                               
                %%% Mean Data 
                Mean_wcc     = vals_mean{Wcc_idx, spont_i};
                Mean_bcc     = vals_mean{Bcc_idx, spont_i};
                Mean_permute = vals_mean{Wrcc_idx, spont_i};                
                
                meanRatio{spont_i} = Mean_bcc/Mean_wcc;
                meanProb{spont_i}  = getRandomizedProb(Mean_wcc, Mean_permute, 'left');
                
                
                %%% KS Data
                KS_wcc = vals_KS{Wcc_idx, spont_i};
                KS_permute = vals_KS{Wrcc_idx, spont_i};                                

                ksStat{spont_i} = KS_wcc;                
                ksProb{spont_i} = getRandomizedProb(KS_wcc, KS_permute, 'right');
                
                
                vals_Wcc          = nonnans( S{ spont_idxs(spont_i) }.val(  pairTypeIdxs{Wcc_idx}  )); 
                vals_Bcc{spont_i} = nonnans( S{ spont_idxs(spont_i) }.val(  pairTypeIdxs{Bcc_idx}  )); 
                [h, ksProb2{spont_i}] = kstest2(vals_Wcc, vals_Bcc{spont_i}, .05, 'larger');
%                 mwwProb{spont_i} = ranksum(vals_Wcc, vals_Bcc{spont_i});    
                
                showDistribs = 0;
                if showDistribs
                                        
                    L_permute = lims(Med_permute);
                    L = lims([vals_median{:}], .05);
                    L(1) = roundToNearest(L(1), .1, 'down');
                    L(2) = roundToNearest(L(2), .1, 'up');
                    
                    [n, xout] = hist(Med_permute, 35 );
                    distM = mean(Med_permute); distS= std(Med_permute);
                    G = @(b, x) b*gaussian(x, distM, distS);
                    A = nlinfit(xout, n, G, max(n));
                    [g_x,g_y] = fplot(@(x) G(A, x), L_permute, 'r');
                                        
                    erfc1 = @(x) erfc(x/sqrt(2)); % integral of gaussian with variance 1 (instead of 1/2)                            
                    nStdDevFromM = abs(Med_wcc - distM)/distS;
                    pval_est = erfc1(nStdDevFromM);

                    figure(76); clf;
                    subplotGap(1,1, 1, [], [.01 0 .15]);
                    hbar = bar(xout, n, 1);
                    set(hbar, 'faceColor', [0 .92 0]); hold on;
                    xlim(L);

                    harr = text(Med_wcc, 0, '\uparrow');
                    set(harr, 'verticalAlignment', 'top', 'horizontalAlignment', 'center', 'fontSize', 25, 'color', 'b', 'fontWeight', 'bold');
                    hline = drawVerticalLine(Med_wcc, 'color', 'b');
                    set(hline, 'lineWidth', 2);

                    h_gauss = plot(g_x, g_y, 'r');
                    set(h_gauss, 'linewidth', 2);
                    
                    
%                     harr = text(Med_wcc, 0, '\downarrow');
%                     set(harr, 'verticalAlignment', 'baseline', 'horizontalAlignment', 'center', 'fontSize', 25, 'color', 'b');
                    
                    3;
                end

                
            end        

            
%             w = num2str( switchh( nSpont, [1 2], [11, 5]));
            if nSpont == 1
                fmt_val = ['%11.2f'];                                
                fmt_pval = ['%11.4f'];
                fmt_pval_w = ['%21.3g'];
            elseif nSpont == 2
                fmt_pval = ['%5.4f'];                
                fmt_pval_w = ['%10.3g'];
            end
            
            if nSpont == 1                
                value_str  = ['| ' fmt_val ' '];
                pval_str   = ['| ' fmt_pval  ' '];                
                pval_str_w = ['| ' fmt_pval_w ' '];                                                
            else
                pval_str   = ['| ' fmt_pval   '/' fmt_pval   ' '];                
                pval_str_w = ['| ' fmt_pval_w '/' fmt_pval_w ' '];                
            end
            value_pval_str = [value_str, pval_str];
            
                                    
%             str_template = ['%14s ' repmat( pval_str, 1, 4)  repmat( pval_str_w, 1, 2)  ' | \n'];        
%             str_template = ['%14s ' repmat( pval_str, 1, 6), repmat( pval_str_w, 1, 1) ' | \n'];        
            str_template = ['%14s ' repmat( value_pval_str, 1, 3) repmat( pval_str_w, 1, 1) ' | \n'];        

            
            fprintf(str_template, ...
                measure, medianRatio{:}, medianProb{:}, meanRatio{:}, meanProb{:}, ksStat{:}, ksProb{:}, ksProb2{:} );
    
            3;
            if saveControlDistribs && ~removeOutliers
                if nSpont == 1
                    S_save.(measure_si) = struct('vals_Bcc',  vals_Bcc{1} , ...
                                              'vals_Wrcc_median', vals_median{Wrcc_idx, 1}, ...
                                              'vals_Wrcc_mean', vals_mean{Wrcc_idx, 1});
                else
                    S_save.(measure_si) = struct('vals_Bcc',  vals_Bcc{1} , ...
                                                 'vals_Wrcc_median', vals_median{Wrcc_idx, 1}, ...
                                                 'vals_Wrcc_mean', vals_mean{Wrcc_idx, 1});
                    S_save.(measure_ss) = struct('vals_Bcc',  vals_Bcc{2} , ...
                                                 'vals_Wrcc_median', vals_median{Wrcc_idx, 2}, ...
                                                 'vals_Wrcc_mean', vals_mean{Wrcc_idx, 2});
                end
            end
            
        end        
        fprintf('\n\n');
        
        if saveControlDistribs
            save(statsDatafile, '-struct', 'S_save');
        end        
        
        
    end
    
    if doPlots
        addLegends = (strcmp(gratingType, 'flashed') &&  addLegend_flashed) || ...
                     (strcmp(gratingType, 'drifting') &&  addLegend_drifting);
        if addFlashedOrDriftingToFigTitles
            fg_str = iff(strcmp(gratingType, 'flashed'), '(Flashed Gratings)', '(Drifting Gratings)');
            add_fg = @(s) {s, fg_str};
        else
            add_fg = @(s) s;
        end
        
        if any(doFigs == 1) && doOri
            % FIGURE 1 --> scatter plots of Global/local ORI width
            w_ori_global_si1 = [allOriStats_si(Wcc_pairs_ori(:,1)).w_ori_global];
            w_ori_global_si2 = [allOriStats_si(Wcc_pairs_ori(:,2)).w_ori_global];
            w_ori_local_si1 = [allOriStats_si(Wcc_pairs_ori(:,1)).w_ori_local];
            w_ori_local_si2 = [allOriStats_si(Wcc_pairs_ori(:,2)).w_ori_local];

            [w_ori_global_si1, w_ori_global_si2] = xyyx(w_ori_global_si1, w_ori_global_si2);
            [w_ori_local_si1,  w_ori_local_si2]  = xyyx(w_ori_local_si1, w_ori_local_si2);

            if doSpontSubtractedMeasures
                w_ori_global_ss1 = [allOriStats_ss(Wcc_pairs_ori(:,1)).w_ori_global];
                w_ori_global_ss2 = [allOriStats_ss(Wcc_pairs_ori(:,2)).w_ori_global];
                w_ori_local_ss1 = [allOriStats_ss(Wcc_pairs_ori(:,1)).w_ori_local];
                w_ori_local_ss2 = [allOriStats_ss(Wcc_pairs_ori(:,2)).w_ori_local];            

                [w_ori_global_ss1, w_ori_global_ss2] = xyyx(w_ori_global_ss1, w_ori_global_ss2);
                [w_ori_local_ss1,  w_ori_local_ss2]  = xyyx(w_ori_local_ss1, w_ori_local_ss2);
                
                subM = 2;
                spont_incl_str = ' Spont Included';                    
            else
                subM = 1;
                spont_incl_str = '';    
            end
            ori_w_types = {'GLOBAL Orientation width', 'LOCAL Orientation width'};
%             ori_w_types = {'w_{Global}^{ORI}', 'w_{Local}^{ORI}'};
            xylims = [0 61];
            figure(1); clf; % plots of orientation width 
            h_ax1(1) = subplotGap(subM,2, 1);  plot(w_ori_global_si1, w_ori_global_si2, '+', 'markersize', 2); ht1(1) = title(add_fg(sprintf('%s%s', ori_w_types{1}, spont_incl_str)));   xlabel('Cell 1'); ylabel('Cell 2'); axis([xylims xylims]); axis square; 
            h_ax1(2) = subplotGap(subM,2, 2);  plot(w_ori_local_si1,  w_ori_local_si2,  '+', 'markersize', 2); ht1(2) = title(add_fg(sprintf('%s%s', ori_w_types{2}, spont_incl_str)));   xlabel('Cell 1'); ylabel('Cell 2'); axis([xylims xylims]); axis square; 
            if doSpontSubtractedMeasures
                h_ax1(3) = subplotGap(subM,2, 3);  plot(w_ori_global_ss1, w_ori_global_ss2, '+', 'markersize', 2); ht1(3) = title(add_fg(sprintf('%s Spont Subtracted', ori_w_types{1}))); xlabel('Cell 1'); ylabel('Cell 2'); axis([xylims xylims]); axis square; 
                h_ax1(4) = subplotGap(subM,2, 4);  plot(w_ori_local_ss1,  w_ori_local_ss2,  '+', 'markersize', 2); ht1(4) = title(add_fg(sprintf('%s Spont Subtracted', ori_w_types{2}))); xlabel('Cell 1'); ylabel('Cell 2'); axis([xylims xylims]); axis square; 
            end
            set(ht1, 'fontsize', title_fsize);
            3;
        end    




        if any(doFigs == 2) && doOri
        % FIGURE 2 --> cumulative differences in Global/local ORI width

            ori_glob_si_idx = find(strcmp(measures_ori, 'Dw_ori_glob_si'));
            ori_loc_si_idx = find(strcmp(measures_ori, 'Dw_ori_loc_si'));
            ori_glob_ss_idx = find(strcmp(measures_ori, 'Dw_ori_glob_ss'));
            ori_loc_ss_idx = find(strcmp(measures_ori, 'Dw_ori_loc_ss'));
        
            plots = {'Dw_ori_glob_si', 'Dw_ori_loc_si'; 'Dw_ori_glob_ss', 'Dw_ori_loc_ss'};
            fig2_ms_idxs = [ori_glob_si_idx, ori_loc_si_idx; ori_glob_ss_idx, ori_loc_ss_idx];
%             ori_w_types = {'w_{Global}^{ORI}', 'w_{Local}^{ORI}'};
            ori_w_types = {'GLOBAL Orientation width', 'LOCAL Orientation width'};
            spont_types = {' Spont Included', ' Spont Subtracted'};
            ori_w_x_labels = {'Differences in w_{Global}^{ORI}, degrees', 'Differences in w_{Local}^{ORI}, degrees'};
            
            if ~doSpontSubtractedMeasures
                plots = plots(1,:);
                spont_types = {''};
            end        
            nSpont = size(plots,1);
            
            nPlots = numel(plots);
            figure(2); clf;
            binIdCloseTo1 = 0;
            maxDiffTh = 0.98;
            [h_ax2, h_ax2_inset] = deal( zeros(1,nPlots) );
            [Wcc_binVals, Bcc_binVals] = deal( cell(1, nPlots) );

            nAcross = 2;
            

            for plot_i = 1:nPlots
                [spont_idx,ori_w_type_idx] = ind2sub(size(plots), plot_i);
                
                ms_idx = fig2_ms_idxs(plot_i);
                h_ax2(plot_i) = subplotGap(subM,nAcross,plot_i, [], [0 0 0], [0 0, 0]);
    %             binEdges = S_ori{ms_idx}.binEdges;
    %             nBins = length(binEdges)-1;
                binEdges = [0:2.5:90];
                binCents = binEdge2cent(binEdges);
    
                [Wcc_cumFracBins, Wcc_binVals{plot_i}] = cumFracBinsFromVals(S_ori{ms_idx}.val(Wcc_oo_pairIdxs), binEdges);
                [Bcc_cumFracBins, Bcc_binVals{plot_i}] = cumFracBinsFromVals(S_ori{ms_idx}.val(Bcc_oo_pairIdxs), binEdges);            

                binIdCloseTo1 = max([binIdCloseTo1, find(Wcc_cumFracBins >= maxDiffTh, 1), find(Bcc_cumFracBins >= maxDiffTh, 1)] );
    %             Wcc_cumFracBins2 = cumFracBinsFromBinIds(S_ori{ms_idx}.bins(Wcc_oo_pairIdxs), nBins);
    %             Bcc_cumFracBins2 = cumFracBinsFromBinIds(S_ori{ms_idx}.bins(Bcc_oo_pairIdxs), nBins);

                plot(binEdges, [0 Wcc_cumFracBins], 'bo-'); hold on;
                plot(binEdges, [0 Bcc_cumFracBins], 'rs:');
                xlabel(ori_w_x_labels{ori_w_type_idx});
                ylabel('Cumulative fraction of pairs');

                ylim([0 1]);
                title(add_fg(sprintf('%s%s', ori_w_types{ori_w_type_idx}, spont_types{spont_idx})), 'fontsize', title_fsize);
                if (plot_i == 2) && addLegends
                    h_leg2 = legend({'Within-site diffs', 'Between-site diffs'}, 'location', 'NW', 'fontsize', 8);
                end            

            end    

            if addProbSubplots
                drawnow;
                for plot_i = 1:nPlots
                    pos_inset = insetAxesPosition(get(h_ax2(plot_i), 'position'), [.4, .02, .6, .6]);
                    h_ax2_inset(plot_i) = axes('outerposition', pos_inset, 'ytick', [], 'fontsize', 8, 'nextplot', 'add');
                    plot(binCents, Wcc_binVals{plot_i}/sum(Wcc_binVals{plot_i}), 'b.-');
                    plot(binCents, Bcc_binVals{plot_i}/sum(Bcc_binVals{plot_i}), 'r.-');                
                    axis(h_ax2_inset(plot_i), 'tight'); 
                    
                    ylims = get(h_ax2_inset(plot_i), 'ylim');
                    set(h_ax2_inset(plot_i), 'ylim', [0, ylims(2)*1.05]);                                        
                    
                    if (plot_i == 2) && addLegends
                        ax_inset_pos = get(h_ax2_inset(plot_i), 'position');                        
                        leg_pos = get(h_leg2, 'position');
                        LB = [ax_inset_pos(1), ax_inset_pos(2)+ax_inset_pos(4)+.02];
                        set(h_leg2, 'position', [LB, leg_pos(3:4)]);                        
                    end
                end
                
                
            end

            xlim2 = 45; %roundToNearest( binEdges(binIdCloseTo1+1), 10, 'up');
            set(h_ax2, 'xlim', [0 xlim2]);
            if addProbSubplots
                set(h_ax2_inset, 'xlim', [0 xlim2]);
            end
    3;

        end
        3;
        if any(doFigs == 3) && doOri
           % 3a. scatterplot plot preferred directions
            ori_pref1 = [allOriStats_si(Wcc_pairs_ori(:,1)).ori_pref_deg];
            ori_pref2 = [allOriStats_si(Wcc_pairs_ori(:,2)).ori_pref_deg];
            oriMax = 180;
            diffGt90_12 = (ori_pref1 - ori_pref2)> oriMax/2;
            diffGt90_21 = (ori_pref2 - ori_pref1)> oriMax/2;
            ori_pref1(diffGt90_12) = ori_pref1(diffGt90_12)-oriMax;
            ori_pref1(diffGt90_21) = ori_pref1(diffGt90_21)+oriMax;

    %         ori_pref2(diffGt90_21) = ori_pref2(diffGt90_21)-180;

            [ori_pref1, ori_pref2] = xyyx(ori_pref1, ori_pref2);        

            figure(3); clf; % plots of orientation width             
            subplotGap(1,2,1);
            plot(ori_pref1, ori_pref2, '+', 'markersize', 2); 
            title(add_fg('Preferred Orientation'), 'fontsize', title_fsize);
            axis(oriMax*[-1/2, 3/2, -1/2, 3/2]);
            tks = oriMax*[-1/2:1/4:3/2];
            set(gca, 'xtick', tks, 'ytick', tks);        
            xlabel('Cell 1'); ylabel('Cell 2');
            axis square;

            hax3b = subplotGap(1,2,2);
            pref_ori_idx = find(strcmp(measures_ori, 'D_ori_pref'));
            3;
            binEdges = S_ori{pref_ori_idx}.binEdges;
            binCents = binEdge2cent(binEdges);
            nBins = length(binEdges)-1;


            [Wcc_cumFracBins, Wcc_binVals] = cumFracBinsFromBinIds(S_ori{pref_ori_idx}.bins(Wcc_oo_pairIdxs), nBins);
    %         if showBccPairs
            [Bcc_cumFracBins, Bcc_binVals] = cumFracBinsFromBinIds(S_ori{pref_ori_idx}.bins(Bcc_oo_pairIdxs), nBins);       
    %         end

            plot(binEdges, [0 Wcc_cumFracBins], 'bo-'); hold on;
    %         if showBccPairs
                plot(binEdges, [0 Bcc_cumFracBins], 'rs:');
    %         end
            xlim(binEdges([1, end]));
            set(gca, 'xtick', [0:15:90])
            title(add_fg('Preferred Orientation'), 'fontsize', title_fsize);
            xlabel('Difference in preferred orienation, degrees');
            ylabel('Cumulative fraction of pairs');
            if addLegends
                legend({'Within-site diffs', 'Between-site diffs'}, 'location', 'SE', 'fontsize', 8);
            end

            
            if addProbSubplots && 0
                drawnow;
                for plot_i = 1:nPlots
                    pos_inset = insetAxesPosition(get(hax3b, 'position'), [.4, .02, .6, .6]);
                    h_ax3b_inset = axes('outerposition', pos_inset, 'ytick', [], 'fontsize', 8, 'nextplot', 'add');
                    plot(binCents, Wcc_binVals/sum(Wcc_binVals), 'b.-');
                    plot(binCents, Bcc_binVals/sum(Bcc_binVals), 'r.-');                
                    axis(h_ax3b_inset(plot_i), 'tight'); 
                    
                    ylims = get(h_ax2_inset(plot_i), 'ylim');
                    set(h_ax3b_inset(plot_i), 'ylim', [0, ylims(2)*1.05]);
                end
                
                set(h_ax3b_inset, 'xlim', binEdges([1, end]));
            end
            3;
           % 3b. cumulative difference between preferred directions 

        end

        if any(doFigs == 4) && doOri

            figure(4); clf;        
            subplotGap(1,2,1);
            pref_ori_idx = find(strcmp(measures_ori, 'D_ori_pref'), 1);        
            binE = S_ori{pref_ori_idx}.binEdges;
            binC = binEdge2cent(binE);
            binVals_cc_all = histcnt(S_ori{pref_ori_idx}.val(Wcc_oo_pairIdxs), binE);
            
            binVals_cc_norm = histcnt(S_ori{pref_ori_idx}.val(Wcc_oo_pairIdxs_norm), binE);
            binVals_cc_1outlier = histcnt(S_ori{pref_ori_idx}.val(Wcc_oo_pairIdxs_1outlier), binE);
            binVals_cc_2outliers = histcnt(S_ori{pref_ori_idx}.val(Wcc_oo_pairIdxs_2outliers), binE);                                   
            assert(all (binVals_cc_all == binVals_cc_norm + binVals_cc_1outlier + binVals_cc_2outliers ));
            
            h_bar4a = bar(binC, [binVals_cc_norm, (binVals_cc_1outlier + binVals_cc_2outliers)], 1, 'stacked');
            set(h_bar4a(1), 'facecolor', 'b');
            set(h_bar4a(2), 'facecolor', 'g');
            
            xlim(binE([1, end]));
            set(gca, 'xtick', [0:15:90]);
            title(add_fg('Preferred Orientation : Pairwise differences'), 'fontsize', title_fsize-1);
            xlabel('Difference in Preferred Orientation, degrees');
            ylabel('Number of cell-cell pairs');
            legend({'Typical cells', 'Outlier cells'})
            if dispPlotStats
                th_deg = 60;
                n_tot = nnz(S_ori{pref_ori_idx}.val(Wcc_oo_pairIdxs) > th_deg);
                n_norm = nnz(S_ori{pref_ori_idx}.val(Wcc_oo_pairIdxs_norm) > th_deg);
                n_out = n_tot-n_norm;
                pct_out = n_out/n_tot*100;
                fprintf('Outlier cells contribute to %d out of %d (%.1f%%) of the pairs above %d degrees: \n', ...
                    n_out, n_tot, pct_out, th_deg);               
            end

            subplotGap(1,2,2);        
            binVals_cm_norm     = histcnt([allOriStats_si(~ori_is_outlier).Dori_pref_allSpkMU], binE);
            binVals_cm_outliers = histcnt([allOriStats_si( ori_is_outlier).Dori_pref_allSpkMU], binE);
            
            h_bar4b = bar(binC, [binVals_cm_norm(:), binVals_cm_outliers(:)], 1, 'stacked');
            set(h_bar4b(1), 'facecolor', 'b');
            set(h_bar4b(2), 'facecolor', 'g');                        
            
            xlim(binE([1, end]));
            set(gca, 'xtick', [0:15:90]);
            title(add_fg('Preferred Orientation : Differences from Multi-units'), 'fontsize', title_fsize-1);
            xlabel('Difference in Preferred Orientation, degrees');
            ylabel('Number of cell-multiunit pairs');            
    %         
    %         binVals_cm2 = histcnt([allOriStats_si.Dori_pref_smlSpkMU], binE);
    %         subplot(1,3,3);        
    %         bar(binC, binVals_cm2, 1);
    %         xlim(binE([1, end]));
    %         title('Preferred Orientation : Differences from Multi-units');        
    %         xlabel('Difference in Preferred Orientation, degrees');
    %         ylabel('Number of Cell pairs');
    %         
            3;



        end

        if any(doFigs == 5) && doOri
            % comparison with multi-units
            gratingStr = [titleCase(gratingType) ' gratings'];
            Dori_pref_MU = [allOriStats_si(oriCells_idx).Dori_pref_allSpkMU];  
            spkAmps = [allOriSpkFeatures(oriCells_idx).spikeAmp];                
            w_ori_global = [allOriStats_si(oriCells_idx).w_ori_global];
            w_ori_local = [allOriStats_si(oriCells_idx).w_ori_local];
            F1oDCs = [allOriStats_si(oriCells_idx).F1oDC];

            idx_mark = Dori_pref_MU > 45;
            figure(5); clf;
            subplotGap(1,3,1); h1 = plot(Dori_pref_MU(~idx_mark), spkAmps(~idx_mark), 'b+', ...
                                      Dori_pref_MU( idx_mark), spkAmps( idx_mark), 'r+');              
            xlabel(' '); %xlabel('Difference from multiunits, degrees'); 
            ylabel('spike amplitude, mV');        
            xlim([0 90]); set(gca, 'xtick', [0:15:90]);
            pvalU_spkAmp = ranksum(spkAmps(~idx_mark), spkAmps( idx_mark) );
            [h, pvalT_spkAmp] = ttest2(spkAmps(~idx_mark), spkAmps( idx_mark) );
            title({' ', sprintf('p_U = %.2g, p_T = %.2g', pvalU_spkAmp, pvalT_spkAmp)});
            med_amps_marked   = nanmedian( spkAmps( idx_mark) );
            med_amps_unmarked = nanmedian( spkAmps( ~idx_mark) );
            fprintf('Median amp of of outliers: %.2f. Median amp of non-outliers: %.2f\n', ...
                med_amps_marked, med_amps_unmarked);
            
            subplotGap(1,3,2); h2 = plot(Dori_pref_MU(~idx_mark),w_ori_global(~idx_mark), 'b+', ...
                                      Dori_pref_MU( idx_mark),w_ori_global( idx_mark), 'r+');             
            xlabel('Difference from multiunits, degrees'); 
            ylabel('w_{Global}^{ORI}');
            xlim([0 90]); set(gca, 'xtick', [0:15:90]);                                  
            pvalU_oriGlobal = ranksum(w_ori_global(~idx_mark), w_ori_global( idx_mark) );
            [h, pvalT_oriGlobal] = ttest2(w_ori_global(~idx_mark), w_ori_global( idx_mark) );
            title({gratingStr, sprintf('p_U = %.2g, p_T = %.2g', pvalU_oriGlobal, pvalT_oriGlobal)});
            med_w_ori_glob_marked   = nanmedian( w_ori_global( idx_mark) );
            med_w_ori_glob_unmarked = nanmedian( w_ori_global( ~idx_mark) );
            fprintf('Median w_ori_global of of outliers: %.2f. Median w_ori_global of non-outliers: %.2f\n',...
                med_w_ori_glob_marked, med_w_ori_glob_unmarked);
            
%             w_ori_local = F1oDCs;
            subplotGap(1,3,3); h3 = plot(Dori_pref_MU(~idx_mark),w_ori_local(~idx_mark), 'b+', ...
                                         Dori_pref_MU( idx_mark),w_ori_local( idx_mark), 'r+');         
            xlabel(' '); %xlabel('Difference from multiunits, degrees'); 
            ylabel('w_{Local}^{ORI}');
            xlim([0 90]); set(gca, 'xtick', [0:15:90]);                                  
            pvalU_oriLocal = ranksum( nonnans( w_ori_local(~idx_mark)), nonnans( w_ori_local( idx_mark) ));
            [h, pvalT_oriLocal] = ttest2(w_ori_local(~idx_mark), w_ori_local( idx_mark) );
            title({' ', sprintf('p_U = %.2g, p_T = %.2g', pvalU_oriLocal, pvalT_oriLocal)});            
            med_w_ori_loc_marked   = nanmedian( w_ori_local( idx_mark) );
            med_w_ori_loc_unmarked = nanmedian( w_ori_local( ~idx_mark) );
            fprintf('Median amp of of outliers: %.2f. Median amp of non-outliers: %.2f\n', ...
                med_w_ori_loc_marked, med_w_ori_loc_unmarked);
            3;
            
%             subplotGap(1,4,); h3 = plot(Dori_pref_MU(~idx_mark),w_ori_local(~idx_mark), 'b+', ...
%                                       Dori_pref_MU( idx_mark),w_ori_local( idx_mark), 'r+');         
%             xlabel(' '); %xlabel('Difference from multiunits, degrees'); 
%             ylabel('w_{Local}^{ORI}');
%             xlim([0 90]); set(gca, 'xtick', [0:15:90]);                                  
%             pvalU_oriLocal = ranksum( nonnans( w_ori_local(~idx_mark)), nonnans( w_ori_local( idx_mark) ));
%             [h, pvalT_oriLocal] = ttest2(w_ori_local(~idx_mark), w_ori_local( idx_mark) );
%             title({' ', sprintf('p_U = %.2g, p_T = %.2g', pvalU_oriLocal, pvalT_oriLocal)});            
%             med_w_ori_loc_marked   = nanmedian( w_ori_local( idx_mark) );
%             med_w_ori_loc_unmarked = nanmedian( w_ori_local( ~idx_mark) );
%             fprintf('Median amp of of outliers: %.2f. Median amp of non-outliers: %.2f\n', ...
%                 med_w_ori_loc_marked, med_w_ori_loc_unmarked);
            
            set([h1, h2, h3], 'markersize', 2);
            3;
        end    

        if any(doFigs == 6) && strcmp(gratingType, 'drifting') && doOri
            % scatterplot of DSI's
            dsi_si1 = [allOriStats_si(Wcc_pairs_ori(:,1)).DSI_global];
            dsi_si2 = [allOriStats_si(Wcc_pairs_ori(:,2)).DSI_global];
            [dsi_si1, dsi_si2] = xyyx(dsi_si1, dsi_si2);

            if doSpontSubtractedMeasures
                dsi_ss1 = [allOriStats_ss(Wcc_pairs_ori(:,1)).DSI];
                dsi_ss2 = [allOriStats_ss(Wcc_pairs_ori(:,2)).DSI];
                [dsi_ss1, dsi_ss2] = xyyx(dsi_ss1, dsi_ss2);
                subN = 2;
                spont_incl_str = ' Spont Included';
            else
                spont_incl_str = '';
                subN = 1;
            end
            
%             dsi_str = 'DSI';
            dsi_str = 'Direction Selectivity Index';
            
            figure(6); clf; 
            
%             h_ax6(1) = subplotGap(1,subN, 1);  
            h_ax6(1) = subplotGap(1,2, 1);  
            plot(dsi_si1, dsi_si2, '+', 'markersize', 2); ht6(1) = title(sprintf('%s%s', dsi_str, spont_incl_str));
            xlabel('Cell 1'); ylabel('Cell 2'); axis([0 1 0 1]); axis square;
%             if doSpontSubtractedMeasures
%                 h_ax6(2) = subplotGap(1,subN, 2);  plot(dsi_ss1, dsi_ss2, '+', 'markersize', 2); ht6(2) = title(sprintf('%s Spont Subtracted', dsi_str)); xlabel('Cell 1'); ylabel('Cell 2'); axis([0 1 0 1]); axis square;
%             end
            set(ht6, 'fontsize', title_fsize);
            
        end

        if any(doFigs == 7) && strcmp(gratingType, 'drifting') && doSpf

            DSI_measures = {'D_dsi_si', 'D_dsi_ss'};
            dsi_si_idx = find(strcmp(measures_ori, 'D_dsi_glob_si'), 1);
            dsi_ss_idx = find(strcmp(measures_ori, 'D_dsi_glob_ss'), 1);
            dsi_ms_idxs = [dsi_si_idx, dsi_ss_idx];

            if ~doSpontSubtractedMeasures
                DSI_measures = DSI_measures(1);
                spont_incl_str = '';
            else
                spont_incl_str = ' Spont Included';
            end
            nPlots = length(DSI_measures);
            [h_ax7, h_ax7_inset] = deal( zeros(1,nPlots) );
            [Wcc_binVals, Bcc_binVals] = deal(cell(1,nPlots));
            
%             figure(7); clf;
            for i = 1:nPlots
                
            
                binEdges = [0:.05:1]; %S_ori{dsi_si_idx}.binEdges;
                binCents = binEdge2cent(binEdges);
        %         nBins = length(binEdges)-1;

                [Wcc_cumFracBins, Wcc_binVals{i}] = cumFracBinsFromVals(S_ori{dsi_ms_idxs(i)}.val(Wcc_oo_pairIdxs), binEdges);
                [Bcc_cumFracBins, Bcc_binVals{i}] = cumFracBinsFromVals(S_ori{dsi_ms_idxs(i)}.val(Bcc_oo_pairIdxs), binEdges);
           
            
%                 h_ax7(i) = subplotGap(1,nPlots,i);
                figure(6);
                h_ax7 = subplotGap(1,2, 2);  
                            
                plot(binEdges, [0 Wcc_cumFracBins], 'o-'); hold on;
                plot(binEdges, [0 Bcc_cumFracBins], 'rs:'); 
                xlim(binEdges([1, end]));
                xlabel('Difference in DSI');
                ylabel('Cumulative fraction of pairs');
                ht7 = title(sprintf('DSI%s', spont_incl_str));
                set(ht7, 'fontsize', title_fsize)

                if (i == nPlots) && addLegends
                    h_leg7 = legend({'Within-site diffs', 'Between-site diffs'}, 'location', 'NW', 'fontsize', 8);
                end
                
                
                
            end
            
            if addProbSubplots
                drawnow;
                for plot_i = 1:nPlots
                    pos_inset = insetAxesPosition(get(h_ax7(plot_i), 'position'), [.4, .02, .6, .6]);
                    h_ax7_inset(plot_i) = axes('outerposition', pos_inset, 'ytick', [], 'fontsize', 8, 'nextplot', 'add');
                    plot(binCents, Wcc_binVals{plot_i}/sum(Wcc_binVals{plot_i}), 'b.-');
                    plot(binCents, Bcc_binVals{plot_i}/sum(Bcc_binVals{plot_i}), 'r.-');                
                    axis(h_ax7_inset(plot_i), 'tight'); 
                    xlim([0 1]);
                    ylims = get(h_ax7_inset(plot_i), 'ylim');
                    set(h_ax7_inset(plot_i), 'ylim', [0, ylims(2)*1.05]);
                    
                    if (plot_i == nPlots) && addLegends
                        ax_inset_pos = get(h_ax7_inset(plot_i), 'position');                        
                        leg_pos = get(h_leg7, 'position');
                        LB = [ax_inset_pos(1), ax_inset_pos(2)+ax_inset_pos(4)+.02];
                        set(h_leg7, 'position', [LB, leg_pos(3:4)]);                        
                    end

                    
                end
            end
            
            
        end



        if any(doFigs == 8) && strcmp(gratingType, 'drifting') && doSpf
           % 3a. scatterplot plot preferred directions
            dirMax = 360;
            dir_pref1 = [allOriStats_si(Wcc_pairs_ori(:,1)).dir_pref_deg];
            dir_pref2 = [allOriStats_si(Wcc_pairs_ori(:,2)).dir_pref_deg];

            diffGt90_12 = (dir_pref1 - dir_pref2)> dirMax/2;
            diffGt90_21 = (dir_pref2 - dir_pref1)> dirMax/2;
            dir_pref1(diffGt90_12) = dir_pref1(diffGt90_12)-dirMax;
            dir_pref1(diffGt90_21) = dir_pref1(diffGt90_21)+dirMax;

    %         dir_pref2(diffGt90_21) = dir_pref2(diffGt90_21)-180;

            [dir_pref1, dir_pref2] = xyyx(dir_pref1, dir_pref2);        

            figure(8); clf; % plots of orientation width 
            subplotGap(1,2,1);
            plot(dir_pref1, dir_pref2, '+', 'markersize', 2); 
            title('Preferred Direction', 'fontsize', title_fsize);
            axis(dirMax*[-1/2, 3/2, -1/2, 3/2]);
            tks = dirMax*[-1/2:1/4:3/2];
            set(gca, 'xtick', tks, 'ytick', tks);
            xlabel('Cell 1'); ylabel('Cell 2');
            axis square;
           3;
           
            subplotGap(1,2,2);
            pref_ori_idx = find(strcmp(measures_ori, 'D_ori_pref'));
            3;
            binEdges = S_ori{pref_ori_idx}.binEdges;
    
            allBinIds = S_ori{pref_ori_idx}.bins(Wcc_oo_pairIdxs);
            cumFracBins = cumFracBinsFromBinIds(allBinIds, length(binEdges)-1);
                    
            Wcc_cumFracBins = cumFracBinsFromVals(S_ori{pref_ori_idx}.val(Wcc_oo_pairIdxs), binEdges);
            Bcc_cumFracBins = cumFracBinsFromVals(S_ori{pref_ori_idx}.val(Bcc_oo_pairIdxs), binEdges);            
            
            plot(binEdges, [0 Wcc_cumFracBins], 'bo-'); hold on;
            plot(binEdges, [0 Bcc_cumFracBins], 'rs:');
            xlim(binEdges([1, end]));
            title('Distribution of Differences', 'fontsize', title_fsize);
            ylabel('Cumulative fraction of pairs');
            xlabel('Difference in preferred direction, degrees');

            title('Preferred Direction');
            if addLegends
                legend({'Within-site diffs', 'Between-site diffs'}, 'location', 'best', 'fontsize', 8);
            end            

            3;
           % 3b. cumulative difference between preferred directions 

        end


        if any(doFigs == 9) && strcmp(gratingType, 'drifting') && doSpf

            pref_dir_idx = find(strcmp(measures_ori, 'D_dir_pref'));

            binEdges = S_ori{pref_dir_idx}.binEdges;
            binCents = binEdge2cent(binEdges);
            
            binVals_cc_all = histcnt(S_ori{pref_dir_idx}.val(Wcc_oo_pairIdxs), binEdges);                                   
            
            binVals_cc_norm = histcnt(S_ori{pref_dir_idx}.val(Wcc_oo_pairIdxs_norm), binEdges);
            binVals_cc_1outlier = histcnt(S_ori{pref_dir_idx}.val(Wcc_oo_pairIdxs_1outlier), binEdges);
            binVals_cc_2outliers = histcnt(S_ori{pref_dir_idx}.val(Wcc_oo_pairIdxs_2outliers), binEdges);                                   
            assert(all (binVals_cc_all == binVals_cc_norm + binVals_cc_1outlier + binVals_cc_2outliers ));
            
            figure(9); clf;
            subplotGap(1,2,1);
            h_bar9a = bar(binCents, [binVals_cc_norm, (binVals_cc_1outlier + binVals_cc_2outliers)], 1, 'stacked');
            set(h_bar9a(1), 'facecolor', 'b');
            set(h_bar9a(2), 'facecolor', 'g');
                       
            set(gca, 'xtick', [0:30:180]);
            xlim(binEdges([1, end]));
            title('Preferred Direction : Pairwise differences');        
            xlabel('Difference in Preferred Direction, degrees');
            ylabel('Number of cell-cell pairs');
            legend({'Typical cells', 'Outlier cells'})

            if dispPlotStats
                range_deg = 90+[-20, 30];
                n_tot = nnz(ibetween( S_ori{pref_dir_idx}.val(Wcc_oo_pairIdxs), range_deg));
                n_norm = nnz(ibetween( S_ori{pref_dir_idx}.val(Wcc_oo_pairIdxs_norm), range_deg));
                n_out = n_tot-n_norm;
                pct_out = n_out/n_tot*100;
                fprintf('Outlier cells contribute to %d out of %d (%.1f%%) of the pairs in the range %d-%d: \n', ...
                    n_out, n_tot, pct_out, range_deg);               
            end
                                                            
            subplotGap(1,2,2);        
            dori_mu_all = [allOriStats_si.Ddir_pref_allSpkMU];
            dori_mu_norm = [allOriStats_si(~ori_is_outlier).Ddir_pref_allSpkMU];
            dori_mu_out = [allOriStats_si(ori_is_outlier).Ddir_pref_allSpkMU];
            binVals_cm_norm     = histcnt(dori_mu_norm, binEdges);
            binVals_cm_outliers = histcnt(dori_mu_out, binEdges);
            
            h_bar9b = bar(binCents, [binVals_cm_norm(:), binVals_cm_outliers(:)], 1, 'stacked');
            set(h_bar9b(1), 'facecolor', 'b');
            set(h_bar9b(2), 'facecolor', 'g');
            
            xlim(binEdges([1, end]));
            set(gca, 'xtick', [0:30:180]);
            title('Preferred Direction : Differences from Multi-units');        
            xlabel('Difference in Preferred Orientation, degrees');
            ylabel('Number of cell-multiunit pairs');

            if dispPlotStats
                
                dori_ori_all = S_ori{pref_dir_idx}.val(Wcc_oo_pairIdxs);
                dori_ori_norm = S_ori{pref_dir_idx}.val(Wcc_oo_pairIdxs_norm);
                assert(all(ibetween(S_ori{pref_dir_idx}.val(Wcc_oo_pairIdxs), 0, 180)));
                n_tot = length(dori_mu_all);
                n_norm = length(dori_mu_norm);
                n_out = length(dori_mu_out);
                pct_i = @(n) n/n_tot*100;
                pct_e = @(n) n/n_norm*100;
                n_within_30     = nnz( dori_mu_all < 30) ;
                n_within_30_opp = nnz( dori_mu_all > 150) ;
                n_within_30_45       = nnz( ibetween(dori_mu_all, 30, 45) );
                n_within_135_150       = nnz( ibetween(dori_mu_all, 135, 150) );
                n_aligned = nnz( dori_mu_norm < 45);
                n_anti_aligned = nnz( dori_mu_norm > 135);
                p_aligned = n_aligned/n_norm;
                p_anti_aligned = n_anti_aligned/n_norm;
                p_gt90_exp = 2*p_aligned*p_anti_aligned*100;
                
                
                pctPair_gt90_obs_all = nnz( dori_ori_all > 90)/length(dori_ori_all) * 100;
                pctPair_gt90_obs_norm = nnz( dori_ori_norm > 90)/length(dori_ori_norm) * 100;
                
                n_gt90_obs_norm = nnz( dori_mu_norm > 90);
                ci = bayesConfInt(n_aligned, n_norm, .95);
                fprintf('Cells within 30 of site preferred: %.1f%% (%.1f%% with outliers excluded.\n', pct_i(n_within_30), pct_e(n_within_30) );
                fprintf('Cells within 30 of site opposite: %.1f%% (%.1f%% with outliers excluded.\n', pct_i(n_within_30_opp), pct_e(n_within_30_opp) );
                fprintf('Cells within 30-45 of site preferred: %.1f%% (%.1f%% with outliers excluded.\n', pct_i(n_within_30_45), pct_e(n_within_30_45) );
                fprintf('Cells within 135-150 of site opposite: %.1f%% (%.1f%% with outliers excluded.\n', pct_i(n_within_135_150), pct_e(n_within_135_150) );
                fprintf('Outlier cells (between 45 and 135): %d (%.2f%%)\n', n_out, pct_i(n_out))
                fprintf('Cells aligned: %d. (%.2f%%). Anti-aligned : %d (%.2f%%)\n', n_aligned, pct_e(n_aligned), n_anti_aligned, pct_e(n_anti_aligned) );
                fprintf('N pairs greater than 90 apart (expected): 2 x %.3f x %.3f = %.2f%% \n', p_aligned, p_anti_aligned, p_gt90_exp );
                fprintf('N pairs greater than 90 apart (observed): (%.1f%%). (outliers excluded) (%.1f%%) (outliers included)\n', pctPair_gt90_obs_norm, pctPair_gt90_obs_all );
                fprintf('Bayes confidence interval: %.3f in same direction, %.3f (%d-%d%% in opp direction)\n', ci, round((1-ci([2, 1]))*100));
                3;
            end

        end



        if any(doFigs == 10) && doSpf
            % 10a. Scatterplot of spatial frequency tuning widths
            w_spf_si1 = [allSpfStats_si(Wcc_pairs_spf(:,1)).w_spf];
            w_spf_si2 = [allSpfStats_si(Wcc_pairs_spf(:,2)).w_spf];

            [w_spf_si1, w_spf_si2] = xyyx(w_spf_si1, w_spf_si2);

            figure(10); clf; % plots of orientation width 
            h1 = subplotGap(1,2, 1); 
            plot(w_spf_si1, w_spf_si2, '+', 'markersize', 2); 
            ht10a = title(add_fg('Spatial frequency tuning width (Octaves)'));
            set(ht10a, 'fontsize', title_fsize);
            
            xlabel('Cell 1'); ylabel('Cell 2');
            axis square;
            3;


            h_ax10 = subplotGap(1,2, 2); 
            maxDiffTh = .99;
            w_spf_idx = find(strcmp(measures_spf, 'Dw_spf'));
            binEdges = S_spf{w_spf_idx}.binEdges;
            binEdges = binEdges(1):0.25:binEdges(end);
            binCents = binEdge2cent(binEdges);

            [Wcc_cumFracBins, Wcc_binVals] = cumFracBinsFromVals(S_spf{w_spf_idx}.val(Wcc_ss_pairIdxs), binEdges);
            [Bcc_cumFracBins, Bcc_binVals] = cumFracBinsFromVals(S_spf{w_spf_idx}.val(Bcc_ss_pairIdxs), binEdges);
            binIdCloseTo1 = max(find(Wcc_cumFracBins >= maxDiffTh, 1), find(Bcc_cumFracBins >= maxDiffTh, 1) );
    %         Wcc_cumFracBins = cumFracBinsFromBinIds(S_spf{w_spf_idx}.bins(Wcc_ss_pairIdxs), nBins);
    %         Bcc_cumFracBins = cumFracBinsFromBinIds(S_spf{w_spf_idx}.bins(Bcc_ss_pairIdxs), nBins);                    

            plot(binEdges, [0 Wcc_cumFracBins], 'o-'); hold on
            plot(binEdges, [0 Bcc_cumFracBins], 'rs:'); 
    %         xlim(binEdges([1, end]));
            xlims = [0, binEdges(binIdCloseTo1+1)];
            xlim(xlims);
%             title('Spatial Frequency Tuning Width');
            ht10 = title(add_fg('Distribution of differences'));
            set(ht10, 'fontsize', title_fsize)
            ylabel('Cumulative fraction of pairs');
            xlabel('Difference in spatial frequency tuning width, octaves');
            if addLegends
                h_leg10 = legend({'Within-site diffs', 'Between-site diffs'}, 'location', 'SE', 'fontsize', 8);
            end            
            
            if addProbSubplots
                drawnow;
                
                pos_inset = insetAxesPosition(get(h_ax10, 'position'), [.4, .02, .6, .6]);
                h_ax10_inset = axes('outerposition', pos_inset, 'ytick', [], 'fontsize', 8, 'nextplot', 'add');
                plot(binCents, Wcc_binVals/sum(Wcc_binVals), 'b.-');
                plot(binCents, Bcc_binVals/sum(Bcc_binVals), 'r.-');                
                axis(h_ax10_inset, 'tight');     
                xlim(xlims);

                ylims = get(h_ax10_inset, 'ylim');
                set(h_ax10_inset, 'ylim', [0, ylims(2)*1.05]);

                if addLegends
                    ax_inset_pos = get(h_ax10_inset, 'position');
                    leg_pos = get(h_leg10, 'position');
                    LB = [ax_inset_pos(1), ax_inset_pos(2)+ax_inset_pos(4)+.02];
                    set(h_leg10, 'position', [LB, leg_pos(3:4)]);
                end
                
            end
            
            
        end


        if any(doFigs == 11) && doSpf
            % 10a. Scatterplot of preferred spatial frequency 
            f_opt_si1 = [allSpfStats_si(Wcc_pairs_spf(:,1)).f_opt];
            f_opt_si2 = [allSpfStats_si(Wcc_pairs_spf(:,2)).f_opt];

            [f_opt_si1, f_opt_si2] = xyyx(f_opt_si1, f_opt_si2);

            figure(11); clf; % plots of orientation width 
            h_ax11a = subplotGap(1,2, 1); 
            h_plot11 = loglog(f_opt_si1, f_opt_si2, '+', 'markersize', 2); 
            title(add_fg('Preferred Spatial frequency (cyc/deg)'), 'fontsize', title_fsize);
            xlabel('Cell 1'); ylabel('Cell 2');
            
            f_range =  [0.08, 3];
            all_f_opt = [allSpfStats_si.f_opt];
            all_f_opt = all_f_opt(  ibetween(all_f_opt, f_range) ); % exclude the one outlier for flashed gratings (2504:3), with f_opt = 0.06 c/deg
            
            switch gratingType, 
                case 'drifting', xticks_mark = [.1, .2, .3, .4, .6, 1, 2, 3];         ext_frac = 0.08; n_extra = 0;
                case 'flashed',  xticks_mark = [.1, .2, .3, .4, .6, .8, 1, 1.2, 1.6]; ext_frac = 0.04; n_extra = 6;
            end
            
            absLims = lims(all_f_opt, ext_frac, [], 'log');
            d_log = floor(log10(absLims(1))) : ceil(log10(absLims(2)));                        
                
            all_ticksC = {}; 
            for i = 1:length(d_log)-1
                new_ticks = [10^d_log(i) : 10^d_log(i) : 10^d_log(i+1) + n_extra*10^d_log(i) ];
                all_ticksC{i} = [new_ticks];                
            end
                        
            all_ticks = unique([all_ticksC{:}, xticks_mark]);
%             i_range = find(ibetween(all_ticks, absLims'));
%             all_ticks = all_ticks( [ max(i_range(1)-1, 1) : min(i_range(end)+1, length(all_ticks))  ]  );
            
            all_ticks = all_ticks( ibetween(all_ticks, absLims') );            
            
            idx_noShow = binarySearch(xticks_mark, all_ticks, [], 'exact') == 0;
            
            ticks = all_ticks;
            tick_labels = arrayfun(@(tk) sprintf('%.2g', tk), all_ticks, 'un', 0);
            tick_labels(idx_noShow) = {''};
            set(gca, 'xtick', ticks, 'ytick', ticks, ...
                'xticklabel', tick_labels, 'yticklabel', tick_labels, 'fontsize', 8);
            set(gca, 'xlim', absLims, 'ylim', absLims);
            axis square;
                 %%       
            
            h_ax11b = subplotGap(1,2, 2); 
            maxDiffTh = 0.99;        
            spf_idx = find(strcmp(measures_spf, 'D_spf_pref'), 1);
            binEdges = S_spf{spf_idx}.binEdges;
            binEdges = binEdges(1):0.25:binEdges(end);
            binCents = binEdge2cent(binEdges);
    %         nBins = length(binEdges)-1;

    %         Wcc_cumFracBins = cumFracBinsFromBinIds(S_spf{spf_idx}.bins(Wcc_ss_pairIdxs), nBins);
    %         Bcc_cumFracBins = cumFracBinsFromBinIds(S_spf{spf_idx}.bins(Bcc_ss_pairIdxs), nBins);

            [Wcc_cumFracBins, Wcc_binVals] = cumFracBinsFromVals(S_spf{spf_idx}.val(Wcc_ss_pairIdxs), binEdges);
            [Bcc_cumFracBins, Bcc_binVals] = cumFracBinsFromVals(S_spf{spf_idx}.val(Bcc_ss_pairIdxs), binEdges);        
            binIdCloseTo1 = max(find(Wcc_cumFracBins >= maxDiffTh, 1), find(Bcc_cumFracBins >= maxDiffTh, 1) );

            plot(binEdges, [0 Wcc_cumFracBins], 'o-'); hold on;
            plot(binEdges, [0 Bcc_cumFracBins], 'rs:');
            xlims = [0, binEdges(binIdCloseTo1+1)];
            xlim(xlims);
    %         xlim(binEdges([1, end]));
            title(add_fg('Preferred Spatial frequency '), 'fontsize', title_fsize);
            xlabel('Difference in preferred spatial frequency, octaves');
            if addLegends
                h_leg11 = legend({'Within-site diffs', 'Between-site diffs'}, 'location', 'SE', 'fontsize', 8);
            end
            ylabel('Cumulative fraction of pairs');
            3;
            %%
            if addProbSubplots
                drawnow;                
                pos_inset = insetAxesPosition(get(h_ax11b, 'position'), [.4, .02, .6, .6]);
                h_ax11_inset = axes('outerposition', pos_inset, 'ytick', [], 'fontsize', 8, 'nextplot', 'add');
                plot(binCents, Wcc_binVals/sum(Wcc_binVals), 'b.-');
                plot(binCents, Bcc_binVals/sum(Bcc_binVals), 'r.-');                
                axis(h_ax11_inset, 'tight');     
                xlim(xlims);

                ylims = get(h_ax11_inset, 'ylim');
                set(h_ax11_inset, 'ylim', [0, ylims(2)*1.05]);
                
                if addLegends
                    ax_inset_pos = get(h_ax11_inset, 'position');
                    leg_pos = get(h_leg11, 'position');
                    LB = [ax_inset_pos(1), ax_inset_pos(2)+ax_inset_pos(4)+.02];
                    set(h_leg11, 'position', [LB, leg_pos(3:4)]);
                end
                
            end   
            
            if dispPlotStats
                f_opt = [allSpfStats_si.f_opt];
                fprintf('Range of spatial frequencies is [%.2f, %.2f] c/deg\n', min(f_opt), max(f_opt));
            end
            
            
        end

                
        % supplementary figure;
        if any(doFigs == 12) && strcmp(gratingType, 'drifting') && doOri
            %%
            % FIGURE S5 --> scatter plots of DSI vs. Global/local ORI width
            w_ori_global_si = [allOriStats_si.w_ori_global]; w_ori_global_si = w_ori_global_si(:);
            w_ori_local_si  = [allOriStats_si.w_ori_local];  w_ori_local_si = w_ori_local_si(:);
            dsi_si          = [allOriStats_si.DSI_global];   dsi_si = dsi_si(:);

            [cc_dsi_wglob, p_dsi_wglob] = corr(dsi_si(:), w_ori_global_si(:), 'rows', 'complete');
            [cc_dsi_wloc, p_dsi_wloc] = corr(dsi_si(:), w_ori_local_si(:), 'rows', 'complete');                        
            
            if doSpontSubtractedMeasures
                w_ori_global_ss = [allOriStats_ss.w_ori_global]; w_ori_global_ss = w_ori_global_ss(:);
                w_ori_local_ss  = [allOriStats_ss.w_ori_local];  w_ori_local_ss = w_ori_local_ss(:);
                dsi_ss          = [allOriStats_ss.DSI];    dsi_ss = dsi_ss(:);

                [cc_dsi_wglob_ss, p_dsi_wglob_ss] = corr(dsi_ss(:), w_ori_global_ss(:), 'rows', 'complete');
                [cc_dsi_wloc_ss, p_dsi_wloc_ss] = corr(dsi_ss(:), w_ori_local_ss(:), 'rows', 'complete');                
                
                Nrows = 2;
                spont_incl_str = ' Spont Included';
            else
                Nrows = 1;
                spont_incl_str = '';
            end
                        
            figure(12); clf; % plots of orientation width 
            h_ax12(1) = subplotGap(Nrows,2, 1);      plot(dsi_si, w_ori_global_si, '+', 'markersize', 2); 
            title(sprintf('%s. r = %.2f', spont_incl_str, cc_dsi_wglob));  xlabel('DSI'); ylabel('w_{Global}^{ORI}');
            h_ax12(2) = subplotGap(Nrows,2, 2);      plot(dsi_si, w_ori_local_si,  '+', 'markersize', 2); 
            title(sprintf('%s. r = %.2f', spont_incl_str, cc_dsi_wloc));   xlabel('DSI'); ylabel('w_{Local}^{ORI}');
            
            if doSpontSubtractedMeasures
                h_ax12(3) = subplotGap(Nrows,2, 3);  plot(dsi_ss, w_ori_global_ss, '+', 'markersize', 2);
                title(sprintf('Spontaneous Subtracted. r = %.2f', cc_dsi_wglob_ss)); xlabel('DSI'); ylabel('w_{Global}^{ORI}');            
                h_ax12(4) = subplotGap(Nrows,2, 4);  plot(dsi_ss, w_ori_local_ss,  '+', 'markersize', 2); 
                title(sprintf('Spontaneous Subtracted. r = %.2f', cc_dsi_wloc_ss));  xlabel('DSI'); ylabel('w_{Local}^{ORI}');
            end

            
            fprintf('CC between DSI and global ori width: %.3f (p = %.2g)\n', cc_dsi_wglob, p_dsi_wglob);
            fprintf('CC between DSI and local ori width: %.3f (p = %.2g)\n', cc_dsi_wloc, p_dsi_wloc);
            
            idx = ~ori_is_outlier;
            [cc_dsi_wglob_norm, p_dsi_wglob_norm] = corr(dsi_si(idx), w_ori_global_si(idx), 'rows', 'complete');
            [cc_dsi_wloc_norm, p_dsi_wloc_norm] = corr(dsi_si(idx), w_ori_local_si(idx), 'rows', 'complete');                        
            
            fprintf('CC between DSI and global ori width (Outliers removed): %.3f (p = %.2g)\n', cc_dsi_wglob_norm, p_dsi_wglob_norm);
            fprintf('CC between DSI and local ori width (Outliers removed): %.3f (p = %.2g)\n', cc_dsi_wloc_norm, p_dsi_wloc_norm);
            
            % with outliers removed            
            
            if doSpontSubtractedMeasures
                fprintf('CC between DSI and global ori width (spont subtracted): %.3f (p = %.2g)\n', cc_dsi_wglob_ss, p_dsi_wglob_ss);
                fprintf('CC between DSI and local ori width (spont subtracted): %.3f (p = %.2g)\n', cc_dsi_wloc_ss, p_dsi_wloc_ss);
                                
                [cc_dsi_wglob_ss_norm, p_dsi_wglob_ss_norm] = corr(dsi_ss(idx), w_ori_global_ss(idx), 'rows', 'complete');
                [cc_dsi_wloc_ss_norm, p_dsi_wloc_ss_norm] = corr(dsi_ss(idx), w_ori_local_ss(idx), 'rows', 'complete');                        

                fprintf('CC between DSI and global ori width (spont subtracted; Outliers removed): %.3f (p = %.2g)\n', cc_dsi_wglob_ss_norm, p_dsi_wglob_ss_norm);
                fprintf('CC between DSI and local ori width (spont subtracted; Outliers removed): %.3f (p = %.2g)\n', cc_dsi_wloc_ss_norm, p_dsi_wloc_ss_norm);
                
            end
            
            
            
            3;
        end
        
        % verify that no outliers for spatial frequency
        if any(doFigs == 13) 
            % 1. plot of diff in pref spatial frequency vs max spatial freq tuning width of the pair
            % 2. plot of diff in pref spatial frequency vs min spike amplitude of the pair

            spf_idx = find(strcmp(measures_spf, 'D_spf_pref'), 1);
            dspf_pref = S_spf{spf_idx}.val(Wcc_ss_pairIdxs);

            w_spf1 = [allSpfStats_si(Wcc_pairs_spf(:,1)).w_spf];
            w_spf2 = [allSpfStats_si(Wcc_pairs_spf(:,2)).w_spf];
            max_w_spf = max(w_spf1, w_spf2)';
            
            allSpfSpkFeatures = [allSpfCells.spkFeatures];
            spkAmp1 = [allSpfSpkFeatures(Wcc_pairs_spf(:,1)).spikeAmp];
            spkAmp2 = [allSpfSpkFeatures(Wcc_pairs_spf(:,2)).spikeAmp];
            min_spkAmp = min(spkAmp1, spkAmp2)';
    
            figure(13); clf;
            subplotGap(1,2,1);
            plot(dspf_pref, max_w_spf, '.');
            [r1, p1] = corr(dspf_pref, max_w_spf, 'type', 'spearman', 'rows', 'complete');
            xlabel('Difference in preferred spatial frequency');
            ylabel('Maximum spatial frequency tuning width of the pair');
            title(sprintf('p_s = %.2g', p1));
            
            subplotGap(1,2,2);
            plot(dspf_pref, min_spkAmp, '.');
            xlabel('Difference in preferred spatial frequency');
            ylabel('Minimum spike amplitude of the pair');
            [r2, p2] = corr(dspf_pref, min_spkAmp, 'type', 'spearman', 'rows', 'complete');
            title(sprintf('p_s = %.2g', p2));
            3;
            
        end
        
        if any(doFigs == 14)  % tuning width vs preferred.
            
           % ori width vs preferred ori
            
           % ori width vs spf width ?
            if 0 && strcmp(gratingType, 'flashed');
                %%
                os_oriStats = nestedFields(oriSpfCells, 'stats', 'tuningStats', 'oriStats_si');
                os_spfStats = nestedFields(oriSpfCells, 'stats', 'tuningStats', 'spfStats_si');
                
                os_ori_glob_w = [os_oriStats.w_ori_global];
                os_ori_loc_w = [os_oriStats.w_ori_global];
                os_spf_w = [os_spfStats.w_spf];
                
                plot(os_ori_glob_w, os_spf_w)                
            end                
                
           
           % dori width vs dspf width. ?
           
           
           subplot(1,2,1); 
           3;
           
           pref_ori_idx = find(strcmp(measures_ori, 'D_ori_pref'), 1);           
           ori_glob_si_idx = find(strcmp(measures_ori, 'Dw_ori_glob_si'), 1);
           ori_loc_si_idx = find(strcmp(measures_ori, 'Dw_ori_loc_si'), 1);

           dOri_pref = S_spf{pref_ori_idx}.val(Wcc_oo_pairIdxs);
           dOri_w_glob = S_spf{ori_glob_si_idx}.val(Wcc_oo_pairIdxs);
           dOri_w_loc = S_spf{ori_loc_si_idx}.val(Wcc_oo_pairIdxs);
           
           %%
           nsimp_oo = sum(ori_pairF1oDC_Wcc > 1, 2);
           ori_idx_ss_pair = nsimp_oo== 2;
           ori_idx_cc_pair = nsimp_oo == 0;
           ori_idx_sc_pair = nsimp_oo == 1;

           corrType = 'pearson';
           %%
            figure(14); clf;
            subplotGap(1,2,1);  hold on; box on;
            plot(dOri_pref(ori_idx_sc_pair), dOri_w_glob(ori_idx_sc_pair), 'mo', 'markersize', 2);
            plot(dOri_pref(ori_idx_ss_pair), dOri_w_glob(ori_idx_ss_pair), 'bo', 'markersize', 2);
            plot(dOri_pref(ori_idx_cc_pair), dOri_w_glob(ori_idx_cc_pair), 'ro', 'markersize', 2);
            xlabel('Difference in preferred orientation');
            ylabel('Difference in Global Ori Width');
            ori_glob_str = getSimpleComplexPairCorr(ori_pairF1oDC_Wcc, dOri_pref, dOri_w_glob, corrType); % cc_ori_glob_all, cc_ori_glob_ss, cc_ori_glob_sc, cc_ori_glob_cc, str
            title([{sprintf('\\bf %s Gratings : \\Delta Pref Ori vs \\DeltaOri Global Width\\rm', titleCase(gratingType))}, ori_glob_str]);

            
            subplotGap(1,2,2);  hold on; box on;
            plot(dOri_pref(ori_idx_sc_pair), dOri_w_loc(ori_idx_sc_pair), 'mo', 'markersize', 2);
            plot(dOri_pref(ori_idx_ss_pair), dOri_w_loc(ori_idx_ss_pair), 'bo', 'markersize', 2);
            plot(dOri_pref(ori_idx_cc_pair), dOri_w_loc(ori_idx_cc_pair), 'ro', 'markersize', 2);
            xlabel('Difference in preferred orientation');
            ylabel('Difference in Local Ori Width');
            ori_loc_str = getSimpleComplexPairCorr(ori_pairF1oDC_Wcc, dOri_pref, dOri_w_loc, corrType); % cc_ori_glob_all, cc_ori_glob_ss, cc_ori_glob_sc, cc_ori_glob_cc, str
            title([{sprintf('\\bf %s Gratings : \\Delta Pref Ori vs \\DeltaOri Local Width\\rm', titleCase(gratingType))}, ori_loc_str]);

            
            if strcmp(gratingType, 'drifting')
                %%
                pref_dir_idx = find(strcmp(measures_ori, 'D_dir_pref'), 1);
                dsi_si_idx = find(strcmp(measures_ori, 'D_dsi_glob_si'), 1);
                
                dDir_pref = S_spf{pref_dir_idx}.val(Wcc_oo_pairIdxs);
                dDSI = S_spf{dsi_si_idx}.val(Wcc_oo_pairIdxs);                               
                dDir_pref = dOri_pref;
                
                figure(15); clf;
                hold on; box on;
                plot(dDir_pref(ori_idx_sc_pair), dDSI(ori_idx_sc_pair), 'mo', 'markersize', 2);
                plot(dDir_pref(ori_idx_ss_pair), dDSI(ori_idx_ss_pair), 'bo', 'markersize', 2);
                plot(dDir_pref(ori_idx_cc_pair), dDSI(ori_idx_cc_pair), 'ro', 'markersize', 2);
                xlabel('Difference in preferred Direction');
                ylabel('Difference in DSI');
                dsi_str = getSimpleComplexPairCorr(ori_pairF1oDC_Wcc, dDir_pref, dDSI, corrType); % cc_spf_glob_all, cc_spf_glob_ss, cc_spf_glob_sc, cc_spf_glob_cc, str
                title([{'\\bf\\DeltaPref Direction vs \\Delta DSI\\rm'}, dsi_str]);
                title([{sprintf('\\bf %s Gratings : \\DeltaPref Direction vs \\Delta DSI\\rm', titleCase(gratingType))}, dsi_str]);
                
%                 xlim([-1, 181]); set(gca, 'xtick', [0:45:180])
                xlim([-1, 91]); set(gca, 'xtick', [0:30:90])
                
                3;
            end
            
            

            
            %%
            3;
            
            nsimp_ss = sum(spf_pairF1oDC_Wcc > 1, 2);
           spf_idx_ss_pair = nsimp_ss == 2;
           spf_idx_cc_pair = nsimp_ss == 0;
           spf_idx_sc_pair = nsimp_ss == 1;
            
            
           pref_spf_idx = find(strcmp(measures_ori, 'D_spf_pref'), 1);           
           spf_w_idx = find(strcmp(measures_ori, 'Dw_spf'), 1);

           dSpf_pref = S_spf{pref_spf_idx}.val(Wcc_ss_pairIdxs);
           dSpf_w = S_spf{spf_w_idx}.val(Wcc_ss_pairIdxs);

            figure(16); clf;
            hold on; box on;
            plot(dSpf_pref(spf_idx_sc_pair), dSpf_w(spf_idx_sc_pair), 'mo', 'markersize', 2);
            plot(dSpf_pref(spf_idx_ss_pair), dSpf_w(spf_idx_ss_pair), 'bo', 'markersize', 2);
            plot(dSpf_pref(spf_idx_cc_pair), dSpf_w(spf_idx_cc_pair), 'ro', 'markersize', 2);
            xlabel('Difference in preferred Spatial Frequency');
            ylabel('Difference in Spatial Frequency Width');
            spf_glob_str = getSimpleComplexPairCorr(spf_pairF1oDC_Wcc, dSpf_pref, dSpf_w, corrType); % cc_spf_glob_all, cc_spf_glob_ss, cc_spf_glob_sc, cc_spf_glob_cc, str
            title([{'\\bf\\DeltaPref Spf vs \\Delta Spf Width\\rm'}, spf_glob_str]);
            title([{sprintf('\\bf %s Gratings : \\DeltaPref Spf vs \\Delta Spf Width\\rm', titleCase(gratingType))}, spf_glob_str]);            
            
            
            
        end
        
    end
    
    if dispPlotStats
        
        % dependency on simple/complex        
        idx_simple_spf = spf_cellF1oDC >= 1;
        idx_complex_spf = spf_cellF1oDC < 1;
        nSimpleTot = nnz(idx_simple_spf);
        nComplexTot = nnz(idx_complex_spf);
        fprintf('For %s gratings: %d simple cells and %d complex cells (total of %d)\n', gratingType, nSimpleTot, nComplexTot, length(spf_cellF1oDC))
        
        spf_nSimp_Wcc = sum(spf_pairF1oDC_Wcc > 1, 2);
        spf_nSimp_Bcc = sum(spf_pairF1oDC_Bcc > 1, 2);        

        
        %% Preferred spatial frequency
        % A. Does S/C affect preferred spatial frequency?
        compareSimpleComplexCellStats([allSpfStats_si.f_opt], idx_simple_spf, idx_complex_spf, 'Preferred Spatial Frequency', 51);

        % B. Does SS/SC/CC affect differences in preferred spatial frequency?
        dspf_pref_idx = find(strcmp(measures_spf, 'D_spf_pref'), 1);
        allD_SpfPref_Wcc = S_spf{dspf_pref_idx}.val(Wcc_ss_pairIdxs);
        allD_SpfPref_Bcc = S_spf{dspf_pref_idx}.val(Bcc_ss_pairIdxs);
        
        compareSimpleComplexPairingStats(allD_SpfPref_Wcc, allD_SpfPref_Bcc, spf_nSimp_Wcc, spf_nSimp_Bcc, 'Preferred Spatial Frequency', 52);
        
        
        %% Spatial frequency tuning width
        % A. Does S/C affect single cell tuning?
        compareSimpleComplexCellStats([allSpfStats_si.w_spf], idx_simple_spf, idx_complex_spf, 'Spatial Frequency Tuning Width', 53);
                
        % B. Does SS/SC/CC affect differences in spatial frequency tuning width?
        dw_spf_idx = find(strcmp(measures_spf, 'Dw_spf'), 1);
        allDw_Spf_Wcc = S_spf{dw_spf_idx}.val(Wcc_ss_pairIdxs);
        allDw_Spf_Bcc = S_spf{dw_spf_idx}.val(Bcc_ss_pairIdxs);
        compareSimpleComplexPairingStats(allDw_Spf_Wcc, allDw_Spf_Bcc, spf_nSimp_Wcc, spf_nSimp_Bcc, 'Spatial Frequency Tuning Width', 54);
        
        
        if 1 || strcmp(gratingType, 'flashed')  % can test whether F1/DC affects orientation tuning width
            %% Load F1/DCs%             
            idx_simple_ori = ori_cellF1oDC >= 1;
            idx_complex_ori = ori_cellF1oDC < 1;
            
            ori_nSimp_Wcc = sum(ori_pairF1oDC_Wcc > 1, 2);
            ori_nSimp_Bcc = sum(ori_pairF1oDC_Bcc > 1, 2);
            
            nSimpleTot = nnz(idx_simple_ori);
            nComplexTot = nnz(idx_complex_ori);
            fprintf('For %s gratings: %d simple cells and %d complex cells (total of %d)\n', gratingType, nSimpleTot, nComplexTot, length(ori_cellF1oDC));
                                        

            
            %% Preferred Orientation 
            % Single cell Preferred Orientation
            compareSimpleComplexCellStats([allOriStats_si.ori_pref_deg], idx_simple_ori, idx_complex_ori, 'Preferred Orientation', 55);
            
            % Differences in preferred Orientation
            dori_pref_idx = find(strcmp(measures_ori, 'D_ori_pref'), 1);
            allD_OriPref_Wcc = S_ori{dori_pref_idx}.val(Wcc_oo_pairIdxs);
            allD_OriPref_Bcc = S_ori{dori_pref_idx}.val(Bcc_oo_pairIdxs);
                                   
            
            compareSimpleComplexPairingStats(allD_OriPref_Wcc, allD_OriPref_Bcc, ori_nSimp_Wcc, ori_nSimp_Bcc, 'Preferred Orientation', 56);

            %% Orientation Tuning Width
            % Single cell global/local orientation tuning width?
            compareSimpleComplexCellStats([allOriStats_si.w_ori_global], idx_simple_ori, idx_complex_ori, 'Global Orientation Width', 57);                        
            
            dori_w_glob_idx = find(strcmp(measures_ori, 'Dw_ori_glob_si'), 1);
            allD_OriW_glob_Wcc = S_ori{dori_w_glob_idx}.val(Wcc_oo_pairIdxs);
            allD_OriW_glob_Bcc = S_ori{dori_w_glob_idx}.val(Bcc_oo_pairIdxs);
            compareSimpleComplexPairingStats(allD_OriW_glob_Wcc, allD_OriW_glob_Bcc, ori_nSimp_Wcc, ori_nSimp_Bcc, 'Global Orientation width', 58);

            
            compareSimpleComplexCellStats([allOriStats_si.w_ori_local],  idx_simple_ori, idx_complex_ori, 'Local Orientation Width', 59);
            dori_w_loc_idx = find(strcmp(measures_ori, 'Dw_ori_loc_si'), 1);
            allD_OriW_loc_Wcc = S_ori{dori_w_loc_idx}.val(Wcc_oo_pairIdxs);
            allD_OriW_loc_Bcc = S_ori{dori_w_loc_idx}.val(Bcc_oo_pairIdxs);
            compareSimpleComplexPairingStats(allD_OriW_loc_Wcc, allD_OriW_loc_Bcc, ori_nSimp_Wcc, ori_nSimp_Bcc, 'Local Orientation width', 60);                                    

            if strcmp(gratingType, 'drifting')
                %% DSI
                % Single cell DSI
                compareSimpleComplexCellStats([allOriStats_si.DSI_global], idx_simple_ori, idx_complex_ori, 'Direction Selectivity Index', 61);

                d_dsi_idx = find(strcmp(measures_ori, 'D_dsi_glob_si'), 1);
                allD_DSI_Wcc = S_ori{d_dsi_idx}.val(Wcc_oo_pairIdxs);
                allD_DSI_Bcc = S_ori{d_dsi_idx}.val(Bcc_oo_pairIdxs);
                compareSimpleComplexPairingStats(allD_DSI_Wcc, allD_DSI_Bcc, ori_nSimp_Wcc, ori_nSimp_Bcc, 'Direction Selectivity Index', 62);                
            end
            3;
            
        end
        
                
        
        %% probability of simple/complex pairing
        spf_pairTypeIdxs = {Wcc_ss_pairIdxs, Wrcc_ss_pairIdxs, Bcc_ss_pairIdxs};        
        allGC = [allSpfCells.Gid]*10000 + [allSpfCells.cellId];
        idx_used = binarySearch(allGC, unique(allGC(Wcc_pairs_spf(:))), [], 0);
        
        ncells_used = length(idx_used);
        nsimp_used = nnz( spf_cellF1oDC(idx_used) >=1 );
        ncomp_used = nnz( spf_cellF1oDC(idx_used) < 1 );
        
        p_simple = nsimp_used/ncells_used;
        p_complex = 1-p_simple;
        
        frac_ss_exp = p_simple^2;
        frac_cc_exp = p_complex^2;
        frac_sc_exp = 2*p_simple*p_complex;
        
        Frac_ss = @(x) nnz(x == 2)/length(x);
        Frac_sc = @(x) nnz(x == 1)/length(x);
        Frac_cc = @(x) nnz(x == 0)/length(x);
        s_idx = find(strcmp(measures_spf, 'D_F1pair'), 1);
        for pt_i = 1:length(pairTypes)
            allVals = S_spf{ s_idx }.val;
            pairIdxs = spf_pairTypeIdxs{pt_i};
            if ~iscell(pairIdxs)
                %                         vals = allVals ( pairIdxs );
                frac_ss{pt_i} = Frac_ss(allVals ( pairIdxs ));
                frac_sc{pt_i} = Frac_sc(allVals ( pairIdxs ));
                frac_cc{pt_i} = Frac_cc(allVals ( pairIdxs ));
            else
                frac_ss{pt_i} = cellfun(@(idxs) Frac_ss( allVals ( idxs ) ), pairIdxs ) ;
                frac_sc{pt_i} = cellfun(@(idxs) Frac_sc( allVals ( idxs ) ), pairIdxs ) ;
                frac_cc{pt_i} = cellfun(@(idxs) Frac_cc( allVals ( idxs ) ), pairIdxs ) ;
            end
            
        end
        
        
        %             gids = pairdata.gids(pairtypeidxs{pt_i}(idx_nonnans),:);
        %             cids = pairdata.cellids(pairtypeidxs{pt_i}(idx_nonnans),:);
        %             allGC = unique([gids(:), cids(:)], 'rows');
        %             ncl{spont_i} = size(allGC,1);
        %%
        nwcc = length(spf_pairTypeIdxs{1});
        n_sc_obs = frac_sc{1}*nwcc;
        n_sc_exp = frac_sc_exp*nwcc;
        
        p_sc_pair = 2*p_simple*p_complex;
        p_sc_obs = binocdf(n_sc_obs, nwcc, p_sc_pair);
                
        fprintf('of a total of %d cells used, %d were simple, %d were complex\n', ncells_used, nsimp_used, ncomp_used)
        fprintf('expected: %.1f %% ss pairs, %.1f %% sc pairs, %.1f cc pairs\n', frac_ss_exp*100, frac_sc_exp*100, frac_cc_exp*100);
        fprintf('observed: %.1f %% (%d) ss pairs, %.1f %% (%d) sc pairs, %.1f (%d) cc pairs\n', frac_ss{1}*100, frac_ss{1}*nwcc, frac_sc{1}*100, frac_sc{1}*nwcc, frac_cc{1}*100, frac_cc{1}*nwcc);
        fprintf('probability of observing %d sc pairs (or fewer) instead of (expected) %.1f if random (with p_sc = %.3f): p = %.3g\n', n_sc_obs, n_sc_exp, p_sc_pair, p_sc_obs);
        
        3;
        
        
        
    end
    %     S, pairTypes, measures
% 
%     'Dw_ori_loc_si'
%     'Dw_ori_loc_ss'
%     'Dw_ori_glob_si'
%     'Dw_ori_glob_ss'
    
%             'Dw_ori_glob_si'        'D_dsi_si'  
%     'D_dsi_ss'    'D_ori_pref'    'D_dir_pref'    'Dw_spf', 'D_spf_pref'




end


function [XY, YX] = xyyx(x, y)
    XY = [x(:); y(:)];
    YX = [y(:); x(:)];
end

function [pairs_renumbered, pairsIdxs_selected] = renumberSubsetOfCellPairs(curPairs, curPairIdxs, selectedCells_idx, idxMtx)
    
    if iscell(curPairs)
%         [pairs_renumbered, pairsIdxs_selected] = cellfun(@(curP, curPidxs) ...
%             renumberSubsetOfCellPairs(curP, curPidxs, selectedCells_idx, idxMtx), curPairs, curPairIdxs, 'un', 0);
        nPairs = length(curPairIdxs);
        pairs_renumbered = cell(1, nPairs);
        pairsIdxs_selected = cell(1, nPairs);
        fprintf('Renumbering cell pairs ... ');
        progressBar('init-', nPairs);
        for p_i = 1:nPairs            
            [pairs_renumbered{p_i}, pairsIdxs_selected{p_i}] = ...
            renumberSubsetOfCellPairs(curPairs{p_i}, curPairIdxs{p_i}, selectedCells_idx, idxMtx);
            progressBar(p_i);
        end
        progressBar('done');
        return;
    end

    pairs_renumbered = binarySearch(selectedCells_idx, curPairs, [], 'exact');
    idxPairs_inSubset = find( all(pairs_renumbered, 2) > 0);
    pairs_renumbered = pairs_renumbered(idxPairs_inSubset, :);  % remove entries with '0'
    pairsIdxs_selected = idxMtx(curPairIdxs(idxPairs_inSubset));    
    
end


function [cumFracBins, binVals] = cumFracBinsFromVals(allVals, binEdges)
    
    binVals = histcnt(allVals, binEdges)';
    cumFracBins = cumsum(binVals);
    cumFracBins = cumFracBins / cumFracBins(end);    

end

function [cumFracBins, binVals] = cumFracBinsFromBinIds(allBinIds, nBins)

    [binIds, binCounts] = uniqueCount(allBinIds);  %         [binIds, binIdxs] = uniqueList(allBinIds);
    idx_keep = binIds>0;
    binVals = zeros(1, nBins);
    binVals(binIds(idx_keep)) = binCounts(idx_keep);
    cumsumBins = cumsum(binVals);
    cumFracBins = cumsumBins/cumsumBins(end);

end

function s = fix_exp_str(s)
    s = strrep(s, 'e-0', 'e-');
    s = strrep(s, 'e-0', 'e-');
end

function [x_mean, x_std, x_median, x_p25, x_p75, n] = getMeanStdPctile(x)
     x_mean = nanmean(x);   
     x_std = nanstd(x);
     x_median = nanmedian(x);
     x_p25 = prctile(x,25);
     x_p75 = prctile(x,75);
     n = nnz(~isnan(x));    
end

function [x_mean, x_median, ks_stat, n] = getMeanMedianKS(x, ctrl_dist)
     x_mean = nanmean(x);   
     x_median = nanmedian(x);
     n = nnz(~isnan(x));    
     if (nargin > 1) && ~isempty(ctrl_dist)
        ks_stat = getKSstat(x, ctrl_dist);
     else
        ks_stat = 0;
     end
end

function ctrl_dist = getCtrlDist(X)
    binEdges = [-inf ; unique(X(:)); inf];    
    binCounts  =  histc (X(:), binEdges, 1);
    sumCounts  =  cumsum(binCounts)./sum(binCounts);
    
    CDF  =  sumCounts(1:end-1);    
    ctrl_dist = struct('binEdges', binEdges, 'CDF', CDF);
end
   

function KSstatistic = getKSstat(x, ctrl_dist, tail)

    % binEdges    =  [-inf ; sort([x1;x2]) ; inf];

    binCounts  =  histc (x , ctrl_dist.binEdges, 1);
    sumCounts  =  cumsum(binCounts)./sum(binCounts);
    sampleCDF  =  sumCounts(1:end-1);

    if nargin < 3
        tail = 0;
    end

    switch tail
       case  0      %  2-sided test: T = max|F1(x) - F2(x)|.
          deltaCDF  =  abs(sampleCDF - ctrl_dist.CDF);

       case -1      %  1-sided test: T = max[F2(x) - F1(x)].
          deltaCDF  =  ctrl_dist.CDF - sampleCDF;

       case  1      %  1-sided test: T = max[F1(x) - F2(x)].
          deltaCDF  =  sampleCDF - ctrl_dist.CDF;
    end

    KSstatistic   =  max(deltaCDF);
end

function pval = getRandomizedProb(val_wcc, val_permute, tail)
                    
    nPermutes = length(val_permute);
   
    useCorrectProbs = 1;
    if useCorrectProbs
        probFunc = @(L,N) iff(L>0, (L+1)/(N+2),  0);
    else
        probFunc = @(L,N) L/N;
    end

    if nargin < 3
        tail = 'both';
    end
    
    if ischar(tail)
        tail = switchh(tail, {'left', 'right', 'both'}, [-1, 1, 0]);
    end
    
    switch tail
        case -1,  %e.g.  mean/medians/ks-stats
            L = nnz( val_permute <= val_wcc); % left-tailed test
        case 1,
            L = nnz( val_permute >= val_wcc); % right-tailed test
        case 0,   %e.g.  cc
            M = mean(val_permute);            
            L = nnz( abs(val_permute-M) >= abs(val_wcc-M) ); % two-sided test            
            
    end
    
    pval = probFunc(L, nPermutes);
                    
end



function pos_inset = insetAxesPosition(pos, relDist)
    [L, B, W, H] = dealV(relDist);
%     [.5, .5, .05, .05]

    pos_inset = [pos(1) + pos(3)*L;
                 pos(2) + pos(4)*B;
                 pos(3)*W;
                 pos(4)*H];
end

function compareSimpleComplexCellStats(p, idx_simple, idx_complex, pname, fig_id)
%%
    p_simple = p(idx_simple);
    p_complex = p(idx_complex);
    mean_simple = mean(p_simple);   mean_complex = mean(p_complex);
    median_simple = median(p_simple);   median_complex = median(p_complex);
    [~, pval_ks] = kstest2(p_simple, p_complex); % ks test:
    [~, pval_t] = ttest2(p_simple, p_complex);
    pval_U = ranksum(p_simple, p_complex);
    fprintf(' ** %s  (N simple = %d. N complex = %d. N tot = %d): \n', pname, nnz(idx_simple), nnz(idx_complex), nnz(p));
    fprintf('    Simple Cells : Mean : %.3g.  Std: %.3g.  Median %.3g\n', mean_simple, std(p_simple), median_simple )
    fprintf('    Complex Cells: Mean : %.3g.  Std: %.3g.  Median %.3g\n', mean_complex, std(p_complex), median_complex )
    fprintf('    T-test : p = %.3g, U-test: p = %.3g.  KS test p = %.3g \n', pval_t, pval_U, pval_ks);    
    if nargin >= 5
        %%
        figure(fig_id); clf;
        h = hist2({p_simple, p_complex}, 20, 'norm', 'line');
        set(h(2), 'color', 'r');
        set(h, 'linewidth', 2)
        gratingStr = curGratingType('');
        title({sprintf('%s (%s gratings)', pname, gratingStr), ...
               sprintf('Md: Simple = %.3g, Complex = %.3g', median_simple, median_complex ), ...
               sprintf('p_U = %.3g. p_t = %.3g', pval_U, pval_t)}, 'fontsize', 10); 
        legend('Simple', 'Complex');
        ylims = ylim;
%         if pval_U < .05
        line(median_simple*[1, 1],  [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'b', 'linestyle', '--');
        line(median_complex*[1, 1], [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'r', 'linestyle', '--');

        line(mean_simple*[1, 1],  [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'b', 'linestyle', ':');
        line(mean_complex*[1, 1], [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'r', 'linestyle', ':');
            
        
    end
    3;
    
end

function compareSimpleComplexPairingStats(dval_wcc, dval_bcc, nsimp_wcc, nsimp_bcc, dname, fig_id)
% spf_nSimp_Bcc, allD_SpfPref_Bcc
% spf_pairF1oDC_Wcc, allD_SpfPref_Wcc, spf_pairF1oDC_Bcc, allD_SpfPref_Bcc

    %%
    idx_ss_pair_wcc = nsimp_wcc == 2;
    idx_cc_pair_wcc = nsimp_wcc == 0;
    idx_sc_pair_wcc = nsimp_wcc == 1;

    idx_ss_pair_bcc = nsimp_bcc == 2;
    idx_cc_pair_bcc = nsimp_bcc == 0;
    idx_sc_pair_bcc = nsimp_bcc == 1;
        
    dval_ss_wcc = dval_wcc(idx_ss_pair_wcc);
    dval_sc_wcc = dval_wcc(idx_sc_pair_wcc);
    dval_cc_wcc = dval_wcc(idx_cc_pair_wcc);

    dval_ss_bcc = dval_bcc(idx_ss_pair_bcc);
    dval_sc_bcc = dval_bcc(idx_sc_pair_bcc);
    dval_cc_bcc = dval_bcc(idx_cc_pair_bcc);
    
    % standard analysis
    pu_ss_sc = ranksum(dval_ss_wcc, dval_sc_wcc);
    pu_ss_cc = ranksum(dval_ss_wcc, dval_cc_wcc);
    pu_sc_cc = ranksum(dval_cc_wcc, dval_sc_wcc);

    [~, pt_ss_sc] = ttest2(dval_ss_wcc, dval_sc_wcc);
    [~, pt_ss_cc] = ttest2(dval_ss_wcc, dval_cc_wcc);
    [~, pt_sc_cc] = ttest2(dval_cc_wcc, dval_sc_wcc);

    [~, pks_ss_sc] = kstest2(dval_ss_wcc, dval_sc_wcc);
    [~, pks_ss_cc] = kstest2(dval_ss_wcc, dval_cc_wcc);
    [~, pks_sc_cc] = kstest2(dval_cc_wcc, dval_sc_wcc);

    fprintf('Differences in %s \n', dname);
    fprintf('  SS : median %.3g, mean %.3g, rms = %.3g\n', median(dval_ss_wcc), mean(dval_ss_wcc), rms(dval_ss_wcc));
    fprintf('  SC : median %.3g, mean %.3g, rms = %.3g\n', median(dval_sc_wcc), mean(dval_sc_wcc), rms(dval_sc_wcc));
    fprintf('  CC : median %.3g, mean %.3g, rms = %.3g\n', median(dval_cc_wcc), mean(dval_cc_wcc), rms(dval_cc_wcc));
    fprintf('  SS/SC : U-test: p = %6.2g.  T-test p = %6.2g. KS-test: p = %6.g\n', pu_ss_sc, pt_ss_sc, pks_ss_sc);
    fprintf('  SC/CC : U-test: p = %6.2g.  T-test p = %6.2g. KS-test: p = %6.g\n', pu_sc_cc, pt_sc_cc, pks_sc_cc);
    fprintf('  SS/CC : U-test: p = %6.2g.  T-test p = %6.2g. KS-test: p = %6.g\n\n', pu_ss_cc, pt_ss_cc, pks_ss_cc);
    
    % clustering index approach.
    median_bw_ratio_ss = median(dval_ss_bcc)/median(dval_ss_wcc);
    median_bw_ratio_sc = median(dval_sc_bcc)/median(dval_sc_wcc);
    median_bw_ratio_cc = median(dval_cc_bcc)/median(dval_cc_wcc);
    median_bw_ratio = median(dval_bcc)/median(dval_wcc);

    mean_bw_ratio_ss = mean(dval_ss_bcc)/mean(dval_ss_wcc);
    mean_bw_ratio_sc = mean(dval_sc_bcc)/mean(dval_sc_wcc);
    mean_bw_ratio_cc = mean(dval_cc_bcc)/mean(dval_cc_wcc);
    mean_bw_ratio = mean(dval_bcc)/mean(dval_wcc);
    
    medianRatio_str = sprintf('  Median ratios: All pairs: %.3f.  SS: %.3f   SC: %.3f   CC: %.3f', median_bw_ratio, median_bw_ratio_ss, median_bw_ratio_sc, median_bw_ratio_cc);
    meanRatio_str =  sprintf('  Mean ratios:   All pairs: %.3f.  SS: %.3f   SC: %.3f   CC: %.3f', mean_bw_ratio, mean_bw_ratio_ss, mean_bw_ratio_sc, mean_bw_ratio_cc);
    fprintf(' * Clustering indices\n');
    fprintf('%s\n', medianRatio_str);
    fprintf('%s\n', meanRatio_str);
    
    
    3;
    
    if nargin >= 4
        %%
        figure(fig_id); clf;
        h = hist2({dval_ss_wcc, dval_sc_wcc, dval_cc_wcc}, 15, 'line', 'norm');
        set(h(1), 'color', 'b');
        set(h(2), 'color', 'm');
        set(h(3), 'color', 'r');
        title(dname); 
        legend('SS', 'SC', 'CC');        
        
        set(h, 'linewidth', 2)
        gratingStr = curGratingType('');
        title({sprintf('%s (%s gratings)', dname, gratingStr), ...
               sprintf('Md: SS= %.3g, SC = %.3g, CC = %.3g', median(dval_ss_wcc), median(dval_sc_wcc), median(dval_cc_wcc) ), ...
               sprintf('p_U: SS/SC=%.3g,   SC/CC = %.3g,   SS/CC = %.3g', pu_ss_sc, pu_sc_cc, pu_ss_cc), ...
               medianRatio_str}, 'fontsize', 10) 
        ylims = ylim;
%         if pval_U < .05
        line(median(dval_ss_wcc)*[1, 1],  [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'b', 'linestyle', '--');
        line(median(dval_sc_wcc)*[1, 1],  [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'm', 'linestyle', '--');
        line(median(dval_cc_wcc)*[1, 1],  [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'r', 'linestyle', '--');

        line(mean(dval_ss_wcc)*[1, 1],  [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'b', 'linestyle', ':');
        line(mean(dval_sc_wcc)*[1, 1],  [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'm', 'linestyle', ':');
        line(mean(dval_cc_wcc)*[1, 1],  [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'r', 'linestyle', ':');
        
        
    end    
    3;
    
end


function [pairIdxs, pairIdxList, idxMtx] = useSubsetOfPairIdxs(allPT, pairTypes, nUnits)
    
    if isstruct(allPT)
        [Wcc_pairIdxs, Wcm_pairIdxs, Bcc_pairIdxs, Bcm_pairIdxs, Wrcc_pairIdxs, Wrcm_pairIdxs] = ...
        deal(allPT.Wcc_idxs, allPT.Wcm_idxs, allPT.Bcc_idxs, allPT.Bcm_idxs, allPT.Wrcc_idxs, allPT.Wrcm_idxs);   
        allPairIdxs = {Wcc_pairIdxs, Wrcc_pairIdxs, Bcc_pairIdxs, Wcm_pairIdxs, Wrcm_pairIdxs, Bcm_pairIdxs};
    elseif iscell(allPT)
        assert(length(allPT) == 6);
    end
    allPairTypes = {'Wcc', 'Wrcc', 'Bcc',  'Wcm', 'Wrcm', 'Bcm'};
    
%     skipWrccUnique = exist('allInBccFlag', 'var') && ~isempty('allInBccFlag'), 
    
%     allPairIdxs_list = cellfun(@catIfCell, allPairIdxs, 'un', 0);
    
%     allPairIdxs_list = cellfun(@catIfCell, allPairIdxs, 'un', 0);
    
    pairTypes = pairTypes(  ord(cellfun(@(s) find(strcmp(s, allPairTypes)), pairTypes)) );
    
    pairTypesAvailable = cellfun(@(pr) any ( strcmp(pairTypes, pr)), allPairTypes);
    
    pairIdxs = cellfun( @(pr) allPairIdxs{ strcmp(allPairTypes, pr) },  pairTypes, 'un', false);
        
%     pairIdxList = unique(cat(1, allPairIdxs_list{pairTypesAvailable}));
    pairIdxList = uniqueInts( allPairIdxs (pairTypesAvailable) );        

    idxMtx = zeros(nUnits, nUnits, 'uint32');
    idxMtx(pairIdxList) = 1:length(pairIdxList);

end

function [str, v_all, v_ss, v_sc, v_cc] = getSimpleComplexPairCorr(pairF1oDCs, x1, x2, corrType)
    nsimp = sum(pairF1oDCs > 1, 2);
    idx_ss_pair = nsimp == 2;
    idx_cc_pair = nsimp == 0;
    idx_sc_pair = nsimp == 1;
    cc_str = switchh(corrType, {'pearson', 'spearman'}, {'cc', '\\rho'});
    
    x1 = x1(:);
    x2 = x2(:);
    [v_all.cc, v_all.p] = corr(x1, x2, 'type', 'spearman');
    [v_ss.cc, v_ss.p] = corr(x1(idx_ss_pair), x2(idx_ss_pair), 'type', corrType);
    [v_sc.cc, v_sc.p] = corr(x1(idx_sc_pair), x2(idx_sc_pair), 'type', corrType);
    [v_cc.cc, v_cc.p] = corr(x1(idx_cc_pair), x2(idx_cc_pair), 'type', corrType);
    str1 = sprintf('All (%d) : %s = %.2f. p = %.2g', nnz(x1), cc_str, v_all.cc, v_all.p);
    str2 = sprintf('S/S (%d) : %s = %.2f. p = %.2g', nnz(idx_ss_pair), cc_str, v_ss.cc, v_ss.p);
    str3 = sprintf('C/C (%d) : %s = %.2f. p = %.2g', nnz(idx_cc_pair), cc_str, v_cc.cc, v_cc.p);
    str4 = sprintf('S/C (%d) : %s = %.2f. p = %.2g', nnz(idx_sc_pair), cc_str, v_sc.cc, v_sc.p);
    str = {str1, str2, str3, str4};

end


% function u = uniqueCellInts(C)
%     C = C(~cellfun(@isempty, C));
%     if isempty(C)
%         u = [];
%         return;
%     end
%     
%     idx_subcell = find(cellfun(@iscell, C));
%     if ~isempty(idx_subcell)
%         C(idx_subcell) = cellfun(@uniqueCellInts, C(idx_subcell));
%     end
%     
%     c_max = max(cellfun(@max, C));
%     tf = false(1,c_max);
%     for i = 1:length(C)
%         tf(C{i}) = 1;
%     end
%     u = find(tf(:));
% end

    

%                 erfc1 = @(x) erfc(x/sqrt(2)); % integral of gaussian with variance 1 (instead of 1/2)                            
%                 nStdDevFromM = abs(Med_wcc - permuteDistribM)/permuteDistribS;
%                 pval_GaussEst = erfc1(nStdDevFromM);
% 
%                 mwwProb{spont_i} = pval_GaussEst;    



%{
    this doesn't save much time (maybe a second or two?)
                        [v_mean, v_median, v_KS, v_N] = deal(zeros(1,Npermutes));
                        for perm_i = 1:Npermutes
                            [v_mean(perm_i), v_median(perm_i), v_KS(perm_i), v_N(perm_i)] = ...
                            getMeanMedianKS( allVals ( pairIdxs{perm_i} ), ctrl_dist) ;  
                        end                    
                        [vals_mean{pt_i, spont_i}, vals_median{pt_i, spont_i}, vals_KS{pt_i, spont_i}, N{pt_i, spont_i}] = ...
                            deal(v_mean, v_median, v_KS, v_N);
%}



%{
%%
N = 100000;
sigma = 2;
mu = 9;

A = sigma*randn(1, N)+mu;

idx = randi(N, 2, 100000);
prs = A(idx);
ds = diff(prs, [], 1);


exp_var = 2*sigma^2
act_var = var(ds)

s1 = prs(1,:);
s2 = prs(2,:);
s1sqr = sqrt( mean(s1.^2) )
s2sqr = sqrt( mean(s2.^2) )
s1s2 = 2*sqrt( mean(2.*s1.*s2) )

s1sqr+s2sqr-s1s2

sqrt( mean( s1.^2 + s2.^2 - 2*s1.*s2 ) )
%%
s1b = s1-mu;
s2b = s1-mu;

std(s1b)
mean(s2b)


% mean( s1.^2 + s2.^2 - 2*s1.*s2 ) - (mean(s1.^2) + mean(s2.^2) - mean(2*s1.*s2) )

% m = mean(ds)
% s = std(ds)
% rms = sqrt(mean(ds.^2))

%}


%{    
    cells_tf = [allCells.cellId] > 0;
    if strcmp(gratingType, 'flashed')        
        allOri_ok = nestedFields(allCells, 'stats', 'tuningStats', 'oriStats_si', 'cellOK');
        allSpf_ok = nestedFields(allCells, 'stats', 'tuningStats', 'spfStats_si', 'cellOK');
        oriCells_tf = cells_tf & allOri_ok;
        spfCells_tf = cells_tf & allSpf_ok;           
                
    elseif strcmp(gratingType, 'drifting')
        oriCells_tf = strncmp({allCells.stimType}, 'Grating:Orientation', 19) & cells_tf;
        spfCells_tf = strncmp({allCells.stimType}, 'Grating:Spatial Freq', 19) & cells_tf;        
    end
    
    doOri = nnz(oriCells_tf) > 0;
    doSpf = nnz(spfCells_tf) > 0;
    %}


%{
% Renumber Wcc & Bcc pairs to index the subsets of only orientation batch cells, or only spatial-freq
    % batch cells
    if doOri
%         [Wcc_pairs_ori, Wcc_oo_pairIdxs] = renumberSubsetOfCellPairs(Wcc_pairs, Wcc_pairIdxs, oriCells_idx, idxMtx);    
%         if any(strcmp(pairTypes, 'Wrcc'))        
%             [Wrcc_pairs_ori, Wrcc_oo_pairIdxs] = renumberSubsetOfCellPairs(Wrcc_pairs, Wrcc_pairIdxs, oriCells_idx, idxMtx);
%             assert( all(cellfun(@length, Wrcc_oo_pairIdxs) == length(Wcc_oo_pairIdxs) ) )
%         end
%         [Bcc_pairs_ori, Bcc_oo_pairIdxs] = renumberSubsetOfCellPairs(Bcc_pairs, Bcc_pairIdxs, oriCells_idx, idxMtx);
                
        allOriGids = [allOriCells.Gid];
        ori_cellF1oDC = [allOriStats_si.F1oDC];
        ori_pairF1oDC_Wcc = ori_cellF1oDC(Wcc_pairs_ori);       
        ori_pairF1oDC_Bcc = ori_cellF1oDC(Bcc_pairs_ori);       
    else
        [Wcc_oo_pairIdxs, Wrcc_oo_pairIdxs, Bcc_oo_pairIdxs, ori_cellF1oDC, ori_pairF1oDC_Wcc] = deal([]);
    end

    if doSpf
%         [Wcc_pairs_spf, Wcc_ss_pairIdxs] = renumberSubsetOfCellPairs(Wcc_pairs, Wcc_pairIdxs, spfCells_idx, idxMtx);       
%         if any(strcmp(pairTypes, 'Wrcc'))        
%             [Wrcc_pairs_spf, Wrcc_ss_pairIdxs] = renumberSubsetOfCellPairs(Wrcc_pairs, Wrcc_pairIdxs, spfCells_idx, idxMtx);    
%             assert( all(cellfun(@length, Wrcc_ss_pairIdxs) == length(Wcc_ss_pairIdxs) ) )
%         end
%         [Bcc_pairs_spf, Bcc_ss_pairIdxs] = renumberSubsetOfCellPairs(Bcc_pairs, Bcc_pairIdxs, spfCells_idx, idxMtx);
        
        allSpfGids = [allSpfCells.Gid];
        spf_cellF1oDC = [allSpfStats_si.F1oDC];
        spf_pairF1oDC_Wcc = spf_cellF1oDC(Wcc_pairs_spf);
        spf_pairF1oDC_Bcc = spf_cellF1oDC(Bcc_pairs_spf);
        
        spf_nSimp_Wcc = sum(spf_pairF1oDC_Wcc > 1, 2);
        spf_nSimp_Bcc = sum(spf_pairF1oDC_Bcc > 1, 2);
    else
        [Wcc_ss_pairIdxs, Wrcc_ss_pairIdxs, Bcc_ss_pairIdxs, spf_cellF1oDC, spf_pairF1oDC_Wcc] = deal([]);
    end
    
%}