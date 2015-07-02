function printDegreeOfTuningComparisons(suppressOutput_flag)
    doFigs = [1:13, 106, 107];
%%
%     doFigs = [1:13];

    doPlots = 1             && 1;
        doMainFigurePlots = 0;
        doDepthAndDistancePlots = 1;
        doBetweenAnimalTests = 0;
        doSpikeSortingCheckFigure = 0;
        doOutlierStatsFigure = 0;
        doDSIvsOriWidthFigure = 0;
        doOriVsSpikeAmplitudeFigure = 0;
        doCheckSpatialFrequencyOutliersFig = 0;
        doPrefVsWidthFigure = 0;
        doScatterVsWidthPlots = 0;
    doPairStats = 0        && 1;
    doSigTests = 0     && 1;
%     doMiscPairData = 1;
%     doMiscPairData = 1;
    dispPlotStats = 0;
    dispMiscStats = 0;
    doSimpleComplexStats = 0;
    
    plotInColor = 0;
    doBootstrapsOfClusterIndices = 1;
        nBoots = 1000;
        bootMedianRatios = 1;
        bootMeanRatios = 1;

    subtractSpont = curSubtractSpont;
    preserveSC = curPreserveSimpleComplex;
    preserveAB = curPreserveAligned;
    bccType = curBccType;
        opt.bccType = bccType;
        opt.minNSitesPerPen = 3;
        opt.minNSitesPerAnimal = 3;
        opt.excludeSameGroupPairs = 1;
    
    addSubplotLetters = 1;
    printMUproperties = 1;

    spikeAmpFactor = 1/4;
    
    modMedianRatioForW_spf = 1;
    
%     differentiateSimpleComplexCellsInStats = 1;
    filename_ext = '';

    F1oDC_field = 'F1oDC_maxR_avP_sm';

    doOri = 1;
    doSpf = 1; %~subtractSpont;

%     multiUnitSpikes = 'all spikes';
    multiUnitSpikes = 'small spikes';
    
    
    opt.applyStdErrorThresholds = true;
    opt.maxOriStdErr_deg = 5;
    opt.maxDSIStdErr = 0.1;
    opt.maxSpfStdErr = 0.5;
    
    showWorking = nargin == 0 || isempty(suppressOutput_flag);
    if ~showWorking
        [doPlots, doPairStats, doSigTests, dispPlotStats, dispMiscStats] = deal(0);
        doPairStats = any(suppressOutput_flag == 1);
        doSigTests = any(suppressOutput_flag == 2);
        doSimpleComplexStats = any(suppressOutput_flag == 3);
        doDepthAndDistancePlots = any(suppressOutput_flag == 4);
    end
    opt.showWorking = showWorking;
    comments_C = {};
%     X = rand(50,1);
%     CD = ks_getCtrlDist(X);
%     y = X(randi(length(X), 10, 1));
%     ks = getKSstat(y(:), CD);
%     3;
%     return;
    
    labelSpontIfIncluded = 1;
    labelSpontIfSubtracted = 0;

    
    saveControlDistribs = 0;
    addProbSubplots = 1;
    addFlashedOrDriftingToFigTitles = 1;
    curCmpType('degree');
    
    addLegend_drifting = 1;
    addLegend_flashed = 0;
    
    
    gratingType = curGratingType('');
%     if doSigTests || doPairStats
        curPairTypes('Wcc', 'Wrcc', 'Bcc');            
%     elseif doPairStats
%         curPairTypes('Wcc', 'Bcc');            
%     end
%     curPairTypes('Wcc', 'Wrcc', 'Bcc');            
    
%     statsDatafile = getFileName('controls', filename_ext);
        
    
%     opt.limitToNPermutes = 500;
    opt.limitToNPermutes = [];

    spont_s = iff(subtractSpont, 'Spont Subtracted', 'Spont Included');
    
    switch multiUnitSpikes
        case 'small spikes', dOriMU_field = 'Dori_pref_smlSpkMU';
                             dDirMU_field = 'Ddir_pref_smlSpkMU';
        case 'all spikes',   dOriMU_field = 'Dori_pref_allSpkMU';
                             dDirMU_field = 'Ddir_pref_allSpkMU';
    end
    
    
    
    useCriteria = 1;
    criteria = struct;
    if useCriteria
        %%
        minID = 10;
                      
        criteria.minID = minID;
        
%         criteria.MID_fit_cc.min_ovlp = 0.3;        
    else
        %%
        criteria = struct;
    end
    
    opt.criteria = criteria;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    miscPairStats = struct;
    if doOri
        [allOriUnits, pairData_ori, S_ori, pairTypes, measures_ori, pairIdxs_ori, pairIdxList_ori, idxMtx_ori, ...
            Wcc_oo_pairIdxs, Wcc_oo_pairIdxs_M, Wcc_pairs_ori, ...
            Wrcc_oo_pairIdxs, Wrcc_oo_pairIdxs_M, Wrcc_pairs_ori, ...
            Bcc_oo_pairIdxs, Bcc_oo_pairIdxs_M, Bcc_pairs_ori, ...
            oriCells_use, allOriUnitStats, Wrcc_params_ori, ...
            idx_cellsAtEachSite_ori, idx_pairsAtEachSite_ori, miscPairStats] ...
                = loadDegreeOfTuningData('ori', opt, miscPairStats);            

            allOriSpkFeatures = [allOriUnits.spkFeatures];
            ori_unitF1oDC = [allOriUnits.(F1oDC_field)];
    %         ori_pairF1oDC_Wcc = ori_unitF1oDC(Wcc_pairs_ori);       
    %         ori_pairF1oDC_Bcc = ori_unitF1oDC(Bcc_pairs_ori);                           
    %         ori_pairF1oDC_Wrcc = cellfun(@(idxs) ori_unitF1oDC(idxs), Wrcc_pairs_ori, 'un', 0);
            allOriCellStats = allOriUnitStats(oriCells_use);
            allOriCells = allOriUnits(oriCells_use);
            allOriCellGids = [allOriUnits(oriCells_use).Gid];
            allOriCellIDs = [allOriSpkFeatures(oriCells_use).IsolationDistance];
            allOriErrors = [allOriCellStats.error_jack];
            
            nOriUnits = length( allOriUnits );
            
            randInfo_ori = Wrcc_params_ori;
            randInfo_ori.nUnits = nOriUnits;
            randInfo_ori.idxMtx = idxMtx_ori;
    else
        [Wcc_oo_pairIdxs, Wcc_oo_pairIdxs_M, Wcc_pairs_ori, Bcc_oo_pairIdxs, Bcc_oo_pairIdxs_M, Bcc_pairs_ori] = deal([]);        
        [Wrcc_oo_pairIdxs, Wrcc_oo_pairIdxs_M, Wrcc_pairs_ori, measures_ori] = deal({});
        
        
    end
    
   
    
    
    %%
    
    if doSpf
        [allSpfUnits, pairData_spf, S_spf, pairTypes, measures_spf, pairIdxs_spf, pairIdxList_spf, idxMtx_spf, ...
            Wcc_ss_pairIdxs, Wcc_ss_pairIdxs_M, Wcc_pairs_spf, ...
            Wrcc_ss_pairIdxs, Wrcc_ss_pairIdxs_M, Wrcc_pairs_spf, ...
            Bcc_ss_pairIdxs, Bcc_ss_pairIdxs_M, Bcc_pairs_spf, ...
            spfCells_use, allSpfUnitStats, Wrcc_params_spf, ...
            idx_cellsAtEachSite_spf, idx_pairsAtEachSite_spf, miscPairStats] ...
                = loadDegreeOfTuningData('spf', opt, miscPairStats);        

            spf_unitF1oDC = [allSpfUnits.(F1oDC_field)];
    %         spf_pairF1oDC_Wcc = spf_unitF1oDC(Wcc_pairs_spf);       
    %         spf_pairF1oDC_Bcc = spf_unitF1oDC(Bcc_pairs_spf);
    %         spf_pairF1oDC_Wrcc = cellfun(@(idxs) spf_unitF1oDC(idxs), Wrcc_pairs_spf, 'un', 0);
            allSpfCellStats = allSpfUnitStats(spfCells_use);
            allSpfCells = allSpfUnits(spfCells_use);
            allSpfErrors = [allSpfCellStats.error_jack];
            
            nSpfUnits = length( allSpfUnits );
            
            
            randInfo_spf = Wrcc_params_spf;
            randInfo_spf.nUnits = nSpfUnits;
            randInfo_spf.idxMtx = idxMtx_spf;
    else
        [Wcc_ss_pairIdxs, Wcc_ss_pairIdxs_M, Wcc_pairs_spf, Bcc_ss_pairIdxs, Bcc_ss_pairIdxs_M, Bcc_pairs_spf] = deal([]);        
        [Wrcc_ss_pairIdxs, Wrcc_ss_pairIdxs_M, Wrcc_pairs_spf, measures_spf] = deal({});        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    measures = uniqueInOrder([measures_ori, measures_spf]);    

    % get ori outlier idxs
    if doOri        
        
        cell_dOriMU = [allOriCellStats.(dOriMU_field)];
        unit_dOriMU = [allOriUnitStats.(dOriMU_field)];
        
        ori_cell_is_outlier = cell_dOriMU >  45;
%         ori_cell_is_typical = cell_dOriMU <= 45;
        
        ori_unit_is_outlier = unit_dOriMU > 45;
        ori_unit_is_outlier_idx = find(ori_unit_is_outlier);
        
        ori_unit_is_typical = unit_dOriMU <= 45;
        ori_unit_is_typical_idx = find(ori_unit_is_typical);
        
    %     oriNorm_idx    = oriCells_idx(~ori_unit_is_outlier);
    %     oriOutlier_idx = oriCells_idx(ori_unit_is_outlier);    

        %%
        bothGoodCells = ( all(~isnan(unit_dOriMU(Wcc_pairs_ori)),2) );
    
        pair_is_outlier = binarySearch(ori_unit_is_outlier_idx, Wcc_pairs_ori, [], 'exact');
        pair_is_typical = binarySearch(ori_unit_is_typical_idx, Wcc_pairs_ori, [], 'exact');
        %%
        idx_noOutliers =  ~any(pair_is_outlier, 2) & bothGoodCells; % ~any(pair_is_outlier, 2);
        idx_1outlier  =   xor(pair_is_outlier(:,1), pair_is_outlier(:,2)) & bothGoodCells ;            
        idx_2outliers =   all(pair_is_outlier, 2) & bothGoodCells ;
        idx_withOutliers = ~idx_noOutliers & bothGoodCells;
        assert(nnz(idx_noOutliers) + nnz(idx_1outlier) + nnz(idx_2outliers) + nnz(~bothGoodCells) == length(bothGoodCells))
%%
        Wcc_oo_pairIdxs_norm      = Wcc_oo_pairIdxs(idx_noOutliers);  % Wcc_pairs_ori_norm      = Wcc_pairs_ori(idx_0outliers,:); 
        Wcc_oo_pairIdxs_1outlier  = Wcc_oo_pairIdxs(idx_1outlier);   % Wcc_pairs_ori_1outlier  = Wcc_pairs_ori(idx_1outlier,:);  
        Wcc_oo_pairIdxs_2outliers = Wcc_oo_pairIdxs(idx_2outliers);  % Wcc_pairs_ori_2outliers = Wcc_pairs_ori(idx_2outliers,:);      
        
        nOutliers = nnz(ori_cell_is_outlier);
        nCellsTotal = length(cell_dOriMU);
        nGoodOriCells = nnz(~isnan(cell_dOriMU));
        miscPairStats.outlierStats = struct('nOutliers', nOutliers, 'nGoodOriCells', nGoodOriCells);
        3;
                
        if showWorking
            fprintf('\n%s gratings, %s : %d / %d (%.2f%%) cells are outliers (d_ori from MU > 45)\n\n', gratingType, spont_s, ...
                nOutliers, nCellsTotal, nOutliers/nCellsTotal*100)
            3;
        end
       
        
%         nResamplesTotal = 10000;
        
      


    end    
            
    Npermutes = max(length(Wrcc_oo_pairIdxs), length(Wrcc_ss_pairIdxs));
    
    
    
    
    
    if doPairStats
        
        pairDiffs_matFile = getFileName('pairDiffs');
        pairDiffs_S.cmpType = 'degree';
        pairDiffs_S.gratingType = gratingType;
        pairDiffs_S.subtractSpont = subtractSpont;
        pairDiffs_S.bccType = bccType;
        pairDiffs_S.preserveSC = preserveSC;
        pairDiffs_S.preserveAB = preserveAB;

        if showWorking
            fprintf('********************* TABLE 2: STATISTICS FOR PAIRWISE DIFFERENCES **********************\n');
            fprintf('    Parameter       |    Mean     |      Std    |    Median   |     P25     |     P75      |  N\n');        
        end
        ori_pairTypeIdxs = {Wcc_oo_pairIdxs, Bcc_oo_pairIdxs};
        spf_pairTypeIdxs = {Wcc_ss_pairIdxs, Bcc_ss_pairIdxs};        
    
%         measures_all = {'Dw_ori_glob', 'Dw_ori_loc'  'D_ori_pref' ...
%                         'D_dsi', 'D_dir_pref', ...
%                         'Dw_spf'    'D_spf_pref' };        
        measures_all = measures;
                        
        if printMUproperties
%             is_MU_measure = cellfun(@(s) ~isempty(strfind(s, '_MU')), measures_all);
%             measures_all = measures_all(~is_MU_measure);
        
            % add D_ori_pref_MU and D_dir_pref_MU
            %%
            idx_dOri = find(strcmp(measures_all, 'D_ori_pref'));
            measures_all = [measures_all(1:idx_dOri), 'D_ori_pref_MU', measures_all(idx_dOri+1:end)];
            
            if strcmp(gratingType, 'drifting')
                idx_dDir = find(strcmp(measures_all, 'D_dir_pref'));
                measures_all = [measures_all(1:idx_dDir), 'D_dir_pref_MU', measures_all(idx_dDir+1:end)];                
            end
            
            
        end
        allMeasureNames = cell(1, length(measures_all)*length(pairTypes)); all_ms_idx = 1;
        
        for i = 1:length(measures_all)            
            %%
            measure = measures_all{i};         
            
            isMUmeasure = ~isempty(strfind(measure, 'MU'));
            isOriMeasure = any(strcmp(measure, measures_ori)) || isMUmeasure;                                        
                                                            
            measures_ofSameType = iff(isOriMeasure, measures_ori, measures_spf);            
            
            idx_S = find(strcmp(measures_ofSameType, measure));                        
            
            if ~isMUmeasure
                pairTypes_str = {'Wcc', 'Bcc'};            
            else
                pairTypes_str = {'Wcm'};            
            end
                        
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
%                 cellF1oDCs = ori_cellF1oDC;
%                 pairF1oDCs = ori_pairF1oDC_Wcc;
            else
                S = S_spf;
                pairData = pairData_spf;
                pairTypeIdxs = spf_pairTypeIdxs;
%                 cellF1oDCs = spf_cellF1oDC;
%                 pairF1oDCs = spf_pairF1oDC_Wcc;
            end
            
            
            isOriDirPref = ~isempty(strfind(measure, 'ori_pref')) || ~isempty(strfind(measure, 'dir_pref'));            
            
%             if isOriDirPref
% %                 pairTypes_str = {'Wcc'}; % don't bother with 'Bcc'
%             end
            
            
%             mu_measure = 0;
%             if mu_measure
% %                 pairTypes_str = {'Wcm'};                
%             end                
            
            if showWorking
                switch measure
                    case 'Dw_ori_glob', fprintf('***** ORIENTATION ******\n');
                    case 'D_dsi',       fprintf('***** DIRECTION  ******\n');
                    case 'Dw_spf',      
                        fprintf('***** SPATIAL FREQUENCY ******\n');
                        2;
                    case 'dR_spont_abs',
                        fprintf('***** RESPONSES AT NULL ORIENTATION ******\n');
                end
            end
            
            for pt_i = 1:length(pairTypes_str)
                                                    
                if ~isMUmeasure % Cell-Cell measures
                    vals = S{ idx_S }.val(  pairTypeIdxs{pt_i}   );                                                 %#ok<FNDSB>
                    idx_use = true(size(vals));
                     
                else
                    %%
                    % Cell-MU measures:
                    switch measure 
                        case 'D_ori_pref_MU', 
                            vals = [allOriCellStats.(dOriMU_field)];
                            
                            if opt.applyStdErrorThresholds
                                errs = nestedFields(allOriCellStats, 'error_jack', 'ori_pref');
                                idx_use = (errs <= opt.maxOriStdErr_deg);
                            else
                                idx_use = true(size(vals));
                            end
                           
                        case 'D_dir_pref_MU', 
                            vals = [allOriCellStats.(dDirMU_field)];
                            
                            if opt.applyStdErrorThresholds
                                errs = nestedFields(allOriCellStats, 'error_jack', 'ori_pref', 1);
                                idx_use = (errs <= opt.maxOriStdErr_deg);
                            else
                                idx_use = true(size(vals));
                            end
                        otherwise, error('!');
                    end
                    
                    
                    
                    
                end
                  

                if (pt_i == 1) && ~isMUmeasure % only for within-site, cell-cell
                    idx_nonnans = ~isnan(vals);
                    gids = pairData.Gids(pairTypeIdxs{pt_i}(idx_nonnans & idx_use),:);
                    cids = pairData.cellIds(pairTypeIdxs{pt_i}(idx_nonnans & idx_use),:);                        
                    allGC = unique([gids(:), cids(:)], 'rows');
                    Ncl = size(allGC,1);
                end
                vals = vals(idx_use); 

                if isMUmeasure % cell-multiunit                        
                    Ncl = length( vals);
                end

                [vals_mean, vals_std, vals_median, ...
                    vals_P25, vals_P75, Npr] = getMeanStdPctile(vals);  
                
                
                err_str = '';
                if strcmp(curDegreeOEmode, 'oe_diff') % also calculate median ratios for variability 
                        
                    ms_fld = switchh(measure, {'Dw_ori_glob', 'Dw_ori_loc', 'D_ori_pref', 'D_dsi_glob', 'Dw_spf', 'D_spf_pref'}, ...
                                              {'w_ori_global', 'w_ori_local', 'ori_pref', 'DSI_global', 'w_spf', 'f_opt', ''});
                    if ~isempty(ms_fld)
                        if isOriMeasure
                            allErrMs = [allOriErrors.(ms_fld)];
                        else
                            allErrMs = [allSpfErrors.(ms_fld)];
                        end
                        medianErr = nanmedian(allErrMs);
                        meanErr = nanmean(allErrMs);
                        stdErr = nanstd(allErrMs);
                        pctls = prctile(allErrMs, [25,75]);
%                         err_str = sprintf(' [err: med = %.2f; %.2f � %.2f]', medianErr, meanErr, stdErr);
                        err_str = sprintf(' [err: med = %.2f; %.2f - %.2f]', medianErr, pctls);
                        
                    else
                        err_str = '';
                    end
                    
%                     err_str = sprintf(' [med_err = %.2f]', medianErr);
                end
                
                if strcmp(measure, 'cc_ori')
                    3;
                end
                
                w = num2str( 11 );
                if (~isempty(strfind(measure, 'ori')) || ~isempty(strfind(measure, 'dir'))) && isempty(strfind(measure, 'cc_'))
                    fmt = ['%' w '.1f�'];
                elseif ~isempty(strfind(measure, 'dsi')) || ~isempty(strfind(measure, 'spf')) || ~isempty(strfind(measure, 'cc_'))
                    fmt = ['%' w '.2f'];
                end
                
                pairTmp_str = ['| ' fmt ' '];
                num_tmp = '%d';
                nPairs_str = sprintf([num_tmp ' Pr; '], Npr);
                
                if pt_i == 1
                    nCells_str = [sprintf(num_tmp, Ncl) ' Cl'];
                else
                    nCells_str = '';
                end
                   
                if isMUmeasure
                    nMU_str = sprintf('%d MU', length(unique(allOriCellGids(idx_use) )) );
                    nMU_str_disp = ['; ' nMU_str];
                    nPairs_str = '';
                    N_fields = {'N1', nMU_str, 'N2', nCells_str};
                else
                    nMU_str = '';
                    nMU_str_disp = '';
                    if strcmp(pairTypes_str{pt_i}, 'Wcc')
                        N_fields = {'N1', nPairs_str, 'N2', nCells_str};
                    elseif strcmp(pairTypes_str{pt_i}, 'Bcc')
                        N_fields = {'N1', nPairs_str, 'N2', ''};
                    else
                        error('!');
                    end
                    
                end
                nCells_str = [nCells_str nMU_str_disp];  %#ok<AGROW>
                                
                str_template = ['%14s (' pairTypes_str{pt_i} ')' repmat( pairTmp_str, 1, 5) ' | ' nPairs_str '%s ' err_str ' \n'];        

                if showWorking
                    fprintf(str_template, ...
                        measure, vals_mean, vals_std, vals_median, vals_P25, vals_P75, nCells_str);
                end
                
                measure_fld = [measure '_' pairTypes_str{pt_i}];
                pairDiffs_S.(measure_fld) = struct('mean', vals_mean, 'std', vals_std, 'median', vals_median, 'P25', vals_P25, 'P75', vals_P75, ...
                    N_fields{:} ); 
                allMeasureNames{all_ms_idx} = measure_fld; all_ms_idx = all_ms_idx+1;
                
            end        
            if showWorking
                fprintf('\n');
            end               
        end
        %%
        pairingMeasures = {'D_aligned_pair_Wcc', 'D_aligned_pair_Bcc',   'D_F1pair_ori_Wcc', 'D_F1pair_ori_Bcc',   'D_F1pair_spf_Wcc', 'D_F1pair_spf_Bcc'};
        allMeasureNames_save = allMeasureNames( ~cellfun(@isempty, allMeasureNames) );
        allMeasureNames_save = allMeasureNames_save( ~strCcmp(allMeasureNames_save, pairingMeasures)  );
        
        pairDiffs_S.allMeasureNames_orig = measures_all;
        pairDiffs_S.allMeasureNames = allMeasureNames_save;
        pairDiffs_S.measures_ori = measures_ori;
        pairDiffs_S.measures_spf = measures_spf;
        pairDiffs_S.columns = fieldnames(pairDiffs_S.(measure_fld));
%         pairDiffs_S.miscStats = miscStats;
        save(pairDiffs_matFile, '-struct', 'pairDiffs_S');             

        if showWorking
            fprintf('*******************************************************************************************\n\n\n');
        end
        3;
        
        
    end
    
%     Wrcc_oo_pairIdxs_cat = [Wrcc_oo_pairIdxs{:}];
%     Wrcc_ss_pairIdxs_cat = [Wrcc_ss_pairIdxs{:}];
    
    gauss1stddev_pct = erf(1/sqrt(2)) * 100; % integral(@(x) gaussian(x, 0, 1), -1, 1) * 100; % 

    if doSigTests                        
        doKStest = 1;
        
        

        printResultsForSC_AAPairing = 1;
        
        [num_ss, num_sc, num_cc,   num_aa, frac_ab, frac_bb, frac_ua, frac_ub, frac_uu] = deal(cell(1, 3));

        sigTestComments_C = {};
        pairStats_matFile = getFileName('pairStats');
        pairStats_S.cmpType = 'degree';
        pairStats_S.gratingType = gratingType;
        pairStats_S.subtractSpont = subtractSpont;
        pairStats_S.bccType = bccType;
        pairStats_S.preserveSC = preserveSC;                
        pairDiffs_S.preserveAB = preserveAB;
        
        saveSigStats = ~preserveAB;
        miscSigStats = struct;
        if showWorking
            fprintf('\n\n ************************** SIGNIFICANCE TESTS (Npermute = %d) ********************************* \n', Npermutes);
            fprintf('   Parameter     | Median Ratio| Median Prob | MedianR(lo) | MedianR(hi) |  Mean Ratio |  Mean Prob  |   KS_stat   |   KS prob   |      KS prob(indep)    |\n') 
        end
        
        ori_pairTypeIdxs = {Wcc_oo_pairIdxs, Wrcc_oo_pairIdxs, Bcc_oo_pairIdxs};
        spf_pairTypeIdxs = {Wcc_ss_pairIdxs, Wrcc_ss_pairIdxs, Bcc_ss_pairIdxs};        
    
        measures_all = measures;
%         measures_all(strcmp(measures_all, 'D_F1pair')) = [];
        
        pairingMeasures = {'D_aligned_pair', 'D_F1pair_ori', 'D_F1pair_spf'};
        spont_null_measures_idx = find(strncmp(measures_all, 'dR90', 4) | strncmp(measures_all, 'dR_spont', 6));
        spont_null_measures = measures_all(spont_null_measures_idx); % {'dR_spont_abs'    'dR90_total_abs'    'dR90_stim_abs', 'dR_spont_rel', 'dR90_total_rel'    'dR90_stim_rel'};
%         measures_all = pairingMeasures;
        if preserveAB
            measures_all = {'D_aligned_pair'};
        else
%             measures_all(strcmp(measures_all, 'D_aligned_pair')) = [];
        end
        
        idxsToDo = 1:length(measures_all);
%         idxsToDo = spont_null_measures_idx;

        for i = idxsToDo
            %%
            measure_i = measures_all{i};                                    
                                    
%             isOriMeasure = isempty(strfind(measure_i, 'spf')) && isempty(stind(measure_i, 'F1')) || ~isempty(strfind(measure_i, 'aligned'));
            isOriMeasure = any(strcmp(measure_i, measures_ori)) ;
            isSpfMeasure = any(strcmp(measure_i, measures_spf));
            
            isPairingMeasure = any(strcmp(measure_i, pairingMeasures));
            
            if isOriMeasure && ~doOri
                continue;
            end
            if isSpfMeasure && ~doSpf
                continue;
            end
            
%             isOriMeasure = isempty(strfind(measure, 'spf'));
            %%
            measures_ofSameType = iff(isOriMeasure,measures_ori, measures_spf);                        
            idx_S = find(strcmp(measures_ofSameType, measure_i ));                        
            if isempty(idx_S)
                continue;
            end
                                    
%             pairTypes_str = {'Wcc', 'Bcc'};    
            pairTypesHere = {'Wcc', 'Bcc', 'Wrcc'};            
            Wcc_idx = find(strcmp(pairTypes, 'Wcc'), 1);
            Bcc_idx = find(strcmp(pairTypes, 'Bcc'), 1);
            Wrcc_idx = find(strcmp(pairTypes, 'Wrcc'), 1);
            
            
            if isOriMeasure
                pairTypeIdxs = ori_pairTypeIdxs;
                S = S_ori;
            else
                pairTypeIdxs = spf_pairTypeIdxs;
                S = S_spf;
            end
            
            [vals_mean, vals_median, vals_KS, N] = deal(cell(length(pairTypes), 1));
            
            
            %%
            
              %%  
            if isPairingMeasure
                
                if strncmp(measure_i, 'D_F1pair', 8)
                    ori_spf_type = measure_i(end-2:end);
                    C_code = 0; S_code = 1;
                    CC_code = C_code+C_code;
                    SC_code = S_code+C_code;
                    SS_code = S_code+S_code;
                    
                    Num_ss = @(x) nnz(x == SS_code);
                    Num_sc = @(x) nnz(x == SC_code);
                    Num_cc = @(x) nnz(x == CC_code);
%                     N_ss = @(x) nnz(x == 2);
%                     N_sc = @(x) nnz(x == 1);
%                     N_cc = @(x) nnz(x == 0);
                    for pt_i = 1:length(pairTypes)
                        allVals = S{ idx_S }.val;
                        pairIdxs = pairTypeIdxs{pt_i};
                        if ~iscell(pairIdxs)
                            num_ss{pt_i} = Num_ss(allVals ( pairIdxs ));
                            num_sc{pt_i} = Num_sc(allVals ( pairIdxs ));
                            num_cc{pt_i} = Num_cc(allVals ( pairIdxs ));
                        else
                            num_ss{pt_i} = cellfun(@(idxs) Num_ss( allVals ( idxs ) ), pairIdxs );
                            num_sc{pt_i} = cellfun(@(idxs) Num_sc( allVals ( idxs ) ), pairIdxs );
                            num_cc{pt_i} = cellfun(@(idxs) Num_cc( allVals ( idxs ) ), pairIdxs );
                        end
                    end                 
                    %%
%                     L = nnz(num_sc{Wrcc_idx} <= num_sc{Wcc_idx}); N = length(num_sc{Wrcc_idx});
%                     pval = (L+1)/(N+2);
                    
%                     pval = getRandomizedProb(num_sc{Wcc_idx}, num_sc{Wrcc_idx}, 'right');
                    
                    
                    nPairsWcc = nnz(~isnan( allVals( pairTypeIdxs{Wcc_idx} ) ) );
                    nPairsBcc = nnz(~isnan( allVals( pairTypeIdxs{Bcc_idx} ) ) );
                    nPairsWrcc = cellfun( @(idxs)   nnz(~isnan( allVals( idxs ) ) ), pairTypeIdxs{Wrcc_idx});

                    frac_sc_observed = num_sc{Wcc_idx} / nPairsWcc;
                    frac_sc_rand     = num_sc{Wrcc_idx} ./ nPairsWrcc;
                    %%
                    L = nnz(frac_sc_observed >= frac_sc_rand); N = length(num_sc{Wrcc_idx});
                    [pval, two_sided_sgn] = getRandomizedProb(frac_sc_observed, frac_sc_rand, 'left');
%%
%                     nWcc = nnz(~isnan( allVals(pairTypeIdxs{Wcc_idx}) ) );
%                     nBcc = length(pairTypeIdxs{Bcc_idx});
%                     nnW = @(pct) round(pct*nWcc);
%                     nnB = @(pct) round(pct*nBcc);
                    sc_str{1} = sprintf('In Pw  (total %d pairs), had %d (%.2f%%) SS pairs, %d (%.2f%%) SC pairs and %d (%.2f%%) CC pairs', ...
                        nPairsWcc, num_ss{Wcc_idx}, num_ss{Wcc_idx}/nPairsWcc*100,   num_sc{Wcc_idx}, num_sc{Wcc_idx}/nPairsWcc*100, num_cc{Wcc_idx}, num_cc{Wcc_idx}/nPairsWcc*100 ); 
                    sc_str{2} = sprintf('In Pb  (total %d pairs), had %d (%.2f%%) SS pairs, %d (%.2f%%) SC pairs and %d (%.2f%%) CC pairs', ...
                        nPairsBcc, num_ss{Bcc_idx}, num_ss{Bcc_idx}/nPairsBcc*100,   num_sc{Bcc_idx}, num_sc{Bcc_idx}/nPairsBcc*100, num_cc{Bcc_idx}, num_cc{Bcc_idx}/nPairsBcc*100 ); 
                    sc_str{3} = sprintf('In PW-r (av %.1f pairs), had average of %.1f (%.3f%%) SS pairs, %.1f (%.3f%%) SC pairs and %.1f (%.3f%%) CC pairs', ...
                        mean(nPairsWrcc),  mean(num_ss{Wrcc_idx}), mean(num_ss{Wrcc_idx}./nPairsWrcc*100),  mean(num_sc{Wrcc_idx}), mean(num_sc{Wrcc_idx}./nPairsWrcc*100), mean(num_cc{Wrcc_idx}), mean(num_cc{Wrcc_idx}./nPairsWrcc*100) ); 
                    sc_str{4} = sprintf('nMixed in Wrcc <= Wcc: %d/%d. p = %.5f', L, N, pval );
                    sc_descrip = sprintf('***%s for %s gratings', measure_i, gratingType);
                    if showWorking
                        fprintf('\n%s\n', sc_descrip);
                        cellfun(@(str) fprintf('%s\n', str), sc_str);                                              
                    end
                    sigTestComments_C = [sigTestComments_C, ' ', sc_descrip, sc_str]; 
                    
                    simpComplexShuffStats = struct('numSSpairs', num_ss{Wcc_idx}, 'numSCpairs', num_sc{Wcc_idx}, 'numCCpairs', num_cc{Wcc_idx}, ...
                                                   'nSC_rand_mean', mean(num_sc{Wrcc_idx}./nPairsWrcc), 'L', L, 'simpComp_p', pval);
                    miscSigStats.(sprintf('simpComplexShuffStats_%s', ori_spf_type)) = simpComplexShuffStats;
                    
                elseif strcmp(measure_i, 'D_aligned_pair')
                    %%
%                     Frac_aa = @(x) nnz(x == 2)/length(x);
%                     Frac_ab = @(x) nnz(x == 11)/length(x);
%                     Frac_bb = @(x) nnz(x == 20)/length(x);
%                     
%                     Frac_ua = @(x) nnz(x == 101)/length(x);   
%                     Frac_ub = @(x) nnz(x == 110)/length(x);   
%                     Frac_uu = @(x) nnz(x == 200)/length(x);   
                    a_code = 1; b_code = 10; u_code = 100;
                    aa_code = a_code + a_code;
                    ab_code = a_code + b_code;
                    bb_code = b_code + b_code;
                    
                    ua_code = u_code + a_code;
                    ub_code = u_code + b_code;
                    uu_code = u_code + u_code;

                    Num_aa = @(x) nnz(x == aa_code);
                    Num_ab = @(x) nnz(x == ab_code);
                    Num_bb = @(x) nnz(x == bb_code);
                    
                    Num_ua = @(x) nnz(x == ua_code);   
                    Num_ub = @(x) nnz(x == ub_code);   
                    Num_uu = @(x) nnz(x == uu_code);   
                    
                    for pt_i = 1:length(pairTypes)
                        allVals = S{ idx_S }.val;
                        pairIdxs = pairTypeIdxs{pt_i};
                        if ~iscell(pairIdxs)
                            %                         vals = allVals ( pairIdxs );
%                             n_aa{pt_i} = N_aa(allVals ( pairIdxs ));
%                             n_ab{pt_i} = N_ab(allVals ( pairIdxs ));
%                             n_bb{pt_i} = N_bb(allVals ( pairIdxs ));
%                             nPairs = nnz(isnan(pairIdxs);
                            num_aa{pt_i} = Num_aa(allVals ( pairIdxs ));
                            num_ab{pt_i} = Num_ab(allVals ( pairIdxs ));
                            num_bb{pt_i} = Num_bb(allVals ( pairIdxs ));
                            num_ua{pt_i} = Num_ua(allVals ( pairIdxs ));
                            num_ub{pt_i} = Num_ub(allVals ( pairIdxs ));
                            num_uu{pt_i} = Num_uu(allVals ( pairIdxs ));
                        else
%                             n_aa{pt_i} = cellfun(@(idxs) Num_aa( allVals ( idxs ) ), pairIdxs );
%                             n_ab{pt_i} = cellfun(@(idxs) Num_ab( allVals ( idxs ) ), pairIdxs );
%                             n_bb{pt_i} = cellfun(@(idxs) Num_bb( allVals ( idxs ) ), pairIdxs );
                            num_aa{pt_i} = cellfun(@(idxs) Num_aa( allVals ( idxs ) ), pairIdxs );
                            num_ab{pt_i} = cellfun(@(idxs) Num_ab( allVals ( idxs ) ), pairIdxs );
                            num_bb{pt_i} = cellfun(@(idxs) Num_bb( allVals ( idxs ) ), pairIdxs );
                            num_ua{pt_i} = cellfun(@(idxs) Num_ua( allVals ( idxs ) ), pairIdxs );
                            num_ub{pt_i} = cellfun(@(idxs) Num_ub( allVals ( idxs ) ), pairIdxs );
                            num_uu{pt_i} = cellfun(@(idxs) Num_uu( allVals ( idxs ) ), pairIdxs );
                        end
                    end                   
                    %%                    
                    nPairsWcc = nnz(~isnan( allVals( pairTypeIdxs{Wcc_idx} ) ) );
                    nPairsBcc = nnz(~isnan( allVals( pairTypeIdxs{Bcc_idx} ) ) );
                    nPairsWrcc = cellfun( @(idxs)   nnz(~isnan( allVals( idxs ) ) ), pairTypeIdxs{Wrcc_idx});
                    
                    frac_ab_observed = num_ab{Wcc_idx} / nPairsWcc;
                    frac_ab_rand     = num_ab{Wrcc_idx} ./ nPairsWrcc;

%                     L = nnz(frac_ab_observed >= frac_ab_rand); N = length(num_ab{Wrcc_idx});
                    [pval, two_sided_sgn, L, N] = getRandomizedProb(frac_ab_observed, frac_ab_rand);
%                     [pval, two_sided_sgn] = getRandomizedProb(num_ab{Wcc_idx}, num_ab{Wrcc_idx});
                    
                    
                    aa_comment_str = {};
                    for jj = 1:3
                        %%
                        switch jj
                            case 1, 
                                name = 'within-site'; nPr_str = sprintf(' (total %d pairs)', nPairsWcc);
                                idx = Wcc_idx;
%                                 nn_str = nnW;
                                doAv = false;
                                n_fld = '%d';
                                nPairs = nPairsWcc;
                            case 2,
                                name = 'between-site'; nPr_str = sprintf(' (total %d pairs)', nPairsBcc);
                                idx = Bcc_idx;
%                                 nn_str = nnB;
                                doAv = false;
                                n_fld = '%d';
                                nPairs = nPairsBcc;
                            case 3
                                name = 'randomized within-site'; nPr_str =  sprintf(' (each with %d pairs)', nPairsWcc);
                                idx = Wrcc_idx;
%                                 nn_str = nnWr;
                                doAv = true;
                                n_fld = '%.1f';
                                nPairs = nPairsWcc;
                        end
                        aa_str = sprintf([n_fld ' (%.2f%%) AA pairs'], mean(num_aa{idx}), mean( num_aa{idx}/nPairs*100) );
                        ab_str = sprintf([n_fld ' (%.2f%%) AB pairs'], mean(num_ab{idx}), mean( num_ab{idx}/nPairs*100) );
                        bb_str = sprintf([n_fld ' (%.2f%%) BB pairs'], mean(num_bb{idx}), mean( num_bb{idx}/nPairs*100) );
                        ua_str = sprintf([n_fld ' (%.2f%%) UA pairs'], mean(num_ua{idx}), mean( num_ua{idx}/nPairs*100) );
                        ub_str = sprintf([n_fld ' (%.2f%%) UB pairs'], mean(num_ub{idx}), mean( num_ub{idx}/nPairs*100) );
                        uu_str = sprintf([n_fld ' (%.2f%%) UU pairs'], mean(num_uu{idx}), mean( num_uu{idx}/nPairs*100) );
                        
                        aa_comment_str{jj} = sprintf('In %s distribution %s, had \n %25s,  %25s,  %25s \n %25s,  %25s,  %25s', name, nPr_str, aa_str, ab_str, bb_str, ua_str, ub_str, uu_str);
                    end                                        
%                     aa_str{2} = sprintf('In between-site distribution (total %d pairs), had \n %d (%.2f%%) AA pairs, %d (%.2f%%) AB pairs and %d (%.2f%%) BB pairs', ...
%                         nBcc, nnB(frac_aa{Bcc_idx}), frac_aa{Bcc_idx}*100,  nnB(frac_ab{Bcc_idx}), frac_ab{Bcc_idx}*100, nnB(frac_bb{Bcc_idx}), frac_bb{Bcc_idx}*100 );
%                     aa_str{3} = sprintf('In randomized within-site distribution, had \n %.1f (%.3f%%) AA pairs, %.1f (%.3f%%) AB pairs and %.1f (%.3f%%) BB pairs', ...
%                         mean(nnW(frac_aa{Wrcc_idx})), mean(frac_aa{Wrcc_idx}*100), mean(nnW(frac_ab{Wrcc_idx})), mean(frac_ab{Wrcc_idx}*100), mean(nnW(frac_bb{Wrcc_idx})), mean(frac_bb{Wrcc_idx}*100) );

                    aa_comment_str{end+1} = sprintf('nMixed in Wrcc >= Wcc: %d/%d. p = %.4f', L, N, (L+1)/(N+2) );                    
                    aa_descrip = sprintf('***%s for %s gratings', measure_i, gratingType);
                    if showWorking
                        %%
                        fprintf('\n%s\n', aa_descrip);
                        cellfun(@(str) fprintf('%s\n', str), aa_comment_str);                        
                    end                    
                    %%
                    nAA = num_aa{Wcc_idx};
                    nAB = num_ab{Wcc_idx};
                    nBB = num_bb{Wcc_idx};
                    nWithU = num_ua{Wcc_idx} + num_ub{Wcc_idx} + num_uu{Wcc_idx};
                    assert(nAA + nAB + nBB + nWithU == nPairsWcc);

                    %%
%                     nWithU = nWcc - (nAA + nAB + nBB);      
                    alignedAntiAlignedStats = struct('numAApairs_Wcc', nAA, ...
                                                     'numABpairs_Wcc', nAB, ...
                                                     'numBBpairs_Wcc', nBB, ...
                                                     'numWithU_Wcc', nWithU, ...  
                                                     'fracABpairs_Wrcc', mean( frac_ab_rand ), ...
                                                     'clustering_pval', pval);
                    miscSigStats.alignedAntiAlignedStats = alignedAntiAlignedStats;
%                     
%                     miscSigStats.simpComplexShuffStats = simpComplexShuffStats;
                    
                    %%
                    sigTestComments_C = [sigTestComments_C, ' ', aa_descrip, aa_comment_str]; %#ok<AGROW>
                    3;
                    
                end
                3;
                continue;
                
            end
                        
            allVals = S{ idx_S }.val;
            if doKStest && any(~isnan(allVals(:)))
                ctrl_dist = ks_getCtrlDist(allVals);
            else
                ctrl_dist = [];
            end
            for pt_i = 1:length(pairTypes)                                                                        

                pairIdxs = pairTypeIdxs{pt_i};
                if ~iscell(pairIdxs)                    
                    [vals_mean{pt_i}, vals_median{pt_i}, vals_KS{pt_i}, N{pt_i}] = ...
                        getMeanMedianKS( tovector(allVals ( pairIdxs,:,: )), ctrl_dist);  
                else                    
                    [vals_mean{pt_i}, vals_median{pt_i}, vals_KS{pt_i}, N{pt_i}] = ...
                        getMeanMedianKS( allVals, ctrl_dist, pairIdxs);                                              
                end

            end
            
            vals_wcc = allVals( pairTypeIdxs {1} );
           
            
            
            %%% Median Data 
            Med_wcc = vals_median{Wcc_idx};
            Med_bcc = vals_median{Bcc_idx};
            Med_permute = vals_median{Wrcc_idx};
            medianRatio = Med_bcc/Med_wcc;
            medianRatio_rand = Med_bcc./Med_permute;
            medianProb  = getRandomizedProb(medianRatio, medianRatio_rand, 'right');

            
            %%% Mean Data 
            Mean_wcc     = vals_mean{Wcc_idx};
            Mean_bcc     = vals_mean{Bcc_idx};
            Mean_permute = vals_mean{Wrcc_idx};                

            meanRatio = Mean_bcc/Mean_wcc;
            meanRatio_rand = Mean_bcc./Mean_permute;
            meanProb  = getRandomizedProb(meanRatio, meanRatio_rand, 'right');


            %%% KS Data
            KS_wcc = vals_KS{Wcc_idx};
            KS_permute = vals_KS{Wrcc_idx};                                

            ksStat = KS_wcc;                
            ksProb = getRandomizedProb(KS_wcc, KS_permute, 'right');


            vals_Wcc = nonnans( S{ idx_S  }.val(  pairTypeIdxs{Wcc_idx},:,:  )); 
            vals_Bcc = nonnans( S{ idx_S }.val(  pairTypeIdxs{Bcc_idx},:,: )); 
            haveValues = ~isempty(vals_Wcc) && ~isempty(vals_Bcc);
            if haveValues
                [h, ksProb2] = kstest2(vals_Wcc, vals_Bcc, .05, 'larger');
            else
                ksProb2 = nan;
            end 
            
            
             %%
             
            bootCI_medianRatio_fields = {};
            bootCI_meanRatio_fields = {};
            
            if doBootstrapsOfClusterIndices
%                 fprintf('Bootstrapping cluster indices ...'); tic;
%                 dist_Wcc = allVals( pairTypeIdxs { Wcc_idx } );
%                 dist_Bcc = allVals( pairTypeIdxs { Bcc_idx } );
                
                %%
                if bootMedianRatios && 1
                    if haveValues
%                         nBoots = 10000;
                        medians_wcc_boot = bootstrp(nBoots, @median, vals_Wcc);
                        medians_bcc_boot = bootstrp(nBoots, @median, vals_Bcc);
                        medianRatios_boot = medians_bcc_boot ./ medians_wcc_boot;
                        medianRatio_ci_68 = getBootCI(medianRatios_boot, gauss1stddev_pct);
                        
                        show = 0;
                        
                        if show
                            %%
                            allNBootTest = [10, 50, 100, 500, 1000, 5000, 10000];
                            for nbi = 1:length(allNBootTest)
                                medianRatio_ci_68_v_n{nbi} = getBootCI(medianRatios_boot(1:allNBootTest(nbi)), gauss1stddev_pct);
                            end
                            figure(59); 
                            ci = cat(1, medianRatio_ci_68_v_n{:});
                            plot(allNBootTest, ci, 'o-')
                            set(gca, 'xscale', 'log')
                        end
                        3;
%                         medianRatio_ci_5_95 = getBootCI(medianRatios_boot, 90);
                        bootCI_medianRatio_fields = { 'medianRatio_lo', medianRatio_ci_68(1), 'medianRatio_hi', medianRatio_ci_68(2) };
                    else
                        bootCI_medianRatio_fields = { 'medianRatio_lo', nan, 'medianRatio_hi', nan };
                    end
                end
                
                
                if bootMeanRatios && 1
                    if haveValues
                        means_wcc_boot = bootstrp(nBoots, @mean, vals_Wcc);
                        means_bcc_boot = bootstrp(nBoots, @mean, vals_Bcc);
                        meanRatios_boot = means_bcc_boot ./ means_wcc_boot;
                        
                        meanRatio_ci = getBootCI(meanRatios_boot, gauss1stddev_pct);
                        bootCI_meanRatio_fields = { 'meanRatio_lo', meanRatio_ci(1), 'meanRatio_hi', meanRatio_ci(2) };
                    else
                        bootCI_meanRatio_fields = { 'meanRatio_lo', nan, 'meanRatio_hi', nan };
                    end
                end
                
%                 fprintf('done. '); toc;
                %%                
                
            end
            %%
            
            
            
            if strcmp(curDegreeOEmode, 'oe_diff') % also calculate median ratios for variability 
                3;
                ms_fld = switchh(measure_i, {'Dw_ori_glob', 'Dw_ori_loc', 'D_ori_pref', 'D_DSI', 'Dw_spf', 'D_spf_pref'}, ...
                                            {'w_ori_global', 'w_ori_local', 'ori_pref', 'DSI_global', 'w_spf', 'f_opt', ''});
                if ~isempty(ms_fld)
                    if isOriMeasure
                        allErrMs = [allOriErrors.(ms_fld)];
                    else
                        allErrMs = [allSpfErrors.(ms_fld)];
                    end
                    medianErr = nanmedian(allErrMs);
                    medianRatio_var = Med_bcc / medianErr;
                else
                    medianRatio_var = nan;
                end

                medianRatio_var_flds = {'medianRatio_var', medianRatio_var};
                var_str = sprintf(' [var = %.2f]', medianRatio_var);
            else
                medianRatio_var_flds = {};
                var_str = '';
            end
            
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


            fmt_val = ['%11.2f'];                                
            fmt_pval = ['%11.4f'];
            fmt_pval_w = ['%21.3g'];
            
            value_str  = ['| ' fmt_val ' '];
            pval_str   = ['| ' fmt_pval  ' '];                
            pval_str_w = ['| ' fmt_pval_w ' '];                                                

            value_pval_str = [value_str, pval_str];
                                                
            str_template = ['%16s '      value_pval_str,          value_str value_str,   value_pval_str,          value_pval_str     repmat( pval_str_w, 1, 1) ' | '];        
            str_print = sprintf(str_template, ...
                           measure_i,    medianRatio, medianProb, medianRatio_ci_68,     meanRatio, meanProb,     ksStat, ksProb,    ksProb2 );
            
            
            if showWorking
                
                fprintf('%s%s\n', str_print, var_str)
                
            end

            
            pairStats_S.(measure_i) = struct('medianRatio', medianRatio, 'medianProb', medianProb, 'meanRatio', meanRatio, 'meanProb', meanProb, ...
                'ksStat', ksStat, 'ksProb', ksProb, bootCI_medianRatio_fields{:}, bootCI_meanRatio_fields{:}, ...
                'Wcc_mean', Mean_wcc, 'Bcc_mean', Mean_bcc, 'Wrcc_mean', single(Mean_permute), ... 
                'Wcc_median', Med_wcc, 'Bcc_median', Med_bcc, 'Wrcc_median', single(Med_permute), medianRatio_var_flds{:} ); 
            
            3;
            if saveControlDistribs 
                S_save.(measure_i) = struct('vals_Bcc',  vals_Bcc{1} , ...
                                             'vals_Wrcc_median', vals_median{Wrcc_idx, 1}, ...
                                             'vals_Wrcc_mean', vals_mean{Wrcc_idx, 1});
            end
            
        end        
        if showWorking
            fprintf('\n\n');
        end
        
        if saveControlDistribs
            save(statsDatafile, '-struct', 'S_save');
        end          
        
        if saveSigStats
            allMeasureNames_save = measures_all( ~strCcmp(measures_all, pairingMeasures) );
            pairStats_S.allMeasureNames = allMeasureNames_save;
            pairStats_S.measures_ori = measures_ori;
            pairStats_S.measures_spf = measures_spf;
            pairStats_S.columns = fieldnames(pairStats_S.(measures_all{1}));
            pairStats_S.comments = sigTestComments_C;
        else
            pairStats_S.allMeasureNames = {};
            pairStats_S.columns = {};
        end
        pairStats_S.miscStats = miscSigStats;
        save(pairStats_matFile, '-struct', 'pairStats_S');             
3;
        
        
    end
    
    
    
    
    

  
    ori_glob_idx = find(strcmp(measures_ori, 'Dw_ori_glob'));
    ori_loc_idx = find(strcmp(measures_ori, 'Dw_ori_loc'));
    ori_pref_idx = find(strcmp(measures_ori, 'D_ori_pref'));

    spf_w_idx = find(strcmp(measures_spf, 'Dw_spf'), 1);
    spf_pref_idx = find(strcmp(measures_spf, 'D_spf_pref'), 1);

    
    dOriW_global_cc = S_ori{ori_glob_idx}.val;
    dOriW_global_bcc = dOriW_global_cc(Bcc_oo_pairIdxs);
    dOriW_global_wcc = dOriW_global_cc(Wcc_oo_pairIdxs);
    
    dOriW_local_cc = S_ori{ori_loc_idx}.val;
    dOriW_local_bcc = dOriW_local_cc(Bcc_oo_pairIdxs);
    dOriW_local_wcc = dOriW_local_cc(Wcc_oo_pairIdxs);

    dOriPref_cc = S_ori{ori_pref_idx}.val;
    dOriPref_wcc = dOriPref_cc(Wcc_oo_pairIdxs);
    dOriPref_bcc = dOriPref_cc(Bcc_oo_pairIdxs);
    dOriPref_wcc_norm = dOriPref_cc(Wcc_oo_pairIdxs_norm);

    if strcmp(gratingType, 'drifting')
        dsi_idx = find(strcmp(measures_ori, 'D_dsi_glob'), 1);
        dDSI_cc = S_ori{dsi_idx}.val;
        dDSI_wcc = dDSI_cc(Wcc_oo_pairIdxs);
        dDSI_bcc = dDSI_cc(Bcc_oo_pairIdxs);
    end

    
    dSpf_w_cc = S_spf{spf_w_idx}.val;
    dSpf_w_wcc = dSpf_w_cc(Wcc_ss_pairIdxs);
    dSpf_w_bcc = dSpf_w_cc(Bcc_ss_pairIdxs);

    dSpfPref_cc = S_spf{spf_pref_idx}.val;
    dSpfPref_wcc = dSpfPref_cc(Wcc_ss_pairIdxs);
    dSpfPref_bcc = dSpfPref_cc(Bcc_ss_pairIdxs);
        
    
    if doPlots || doPairStats
        
        
        if plotInColor
            scat_col = 'b';
            w_col = 'b';  b_col = 'r';
            bar_norm = 'b'; bar_out = 'g';
            lineWidth_main = 1;
            lineWidth_inset = 1;
        else
            %%
            scat_col = 'k';
            w_col = 'k'; b_col = [.4 .4 .4];
            lineWidth_main = 2;
            lineWidth_inset = 2;
            bar_norm = [.4, .4, .4]; bar_out = [.85, .85, .85];
        end
        
        
        
%             mk = '+'; sz = 2;
        mk = 'o'; sz = 3;
        
        printPlotsToFilesForPaper = 1;
            printSmallSize = 0;
        
            
        cumMkSize = 5;
        distMkSize = 2;
        dashedLine = ':';
        
        ori_scatter_tick_spacing = 90;
        ori_cum_tick_spacing = 15;
        
        if printPlotsToFilesForPaper
            
            set(0,'DefaultFigureWindowStyle','normal') 
%              set(0,'DefaultFigureWindowStyle','docked') 
            fig_LB = [150, 280];
            fig_offset = [70, -10];
            printWith2RowsForFlashedDrifting = 1;
            if printSmallSize
                title_fsize = 8;
                xlabel_fsize = 6;
                ylabel_fsize = 6;
                legend_fsize = 5;
                fig_height_row = 210;
                fig_width_row = 420;

                lineWidth_main = 1;
                cumMkSize = 3;
                dashedLine = ':';
                
                ori_scatter_tick_spacing = 90;
            else
                title_fsize = 11; % 9
                xlabel_fsize = 12; %9;
                ylabel_fsize = 12; %9;
                legend_fsize = 7.5;
                fig_height_row = 350;
                fig_width_row = 730;
                
            end
        else
            set(0,'DefaultFigureWindowStyle','docked') 
            printWith2RowsForFlashedDrifting = 0;
            
            title_fsize = 12;
            xlabel_fsize = 10;
            ylabel_fsize = 10;
            legend_fsize = 8;

        end

        saveFilesNow = 1 && strcmp(bccType, 'full') && ~curPreserveSimpleComplex && ~curPreserveAligned;
        
        subSpcM = [0.01 0.02, 0.01];
        subSpcN = [0.02 0.00, 0.02];
        
%         within_between_legend_str = {'Within-site diffs', 'Between-site diffs'};
        within_between_legend_str = {'Within-site', 'Between-site'};
        
        printWith2RowsForFlashedDrifting = 1;
        plotFolder = 'C:\Users\Avi\Documents\MATLAB\columbia\miller\CatV1Exp\Figures\DegreePaper\';
%         fig_format = '-deps';
        fig_format = 'pdf';
        
        subplotLetter_offset = [-0.02, -0.0]; % offset (right, up)
%         subplotLetter_up_gap = 0.02;   
        subplotLetFont = 'Helvetica';
%         subplotLetFont = 'Arial';
%         subplotLetFont = 'Georgia';

%         subplotLetterOffset = [.04, subplotLetter_up_gap];
        
        addLegends = (strcmp(gratingType, 'flashed') &&  addLegend_flashed) || ...
                     (strcmp(gratingType, 'drifting') &&  addLegend_drifting);
        if addFlashedOrDriftingToFigTitles
            fg_str = iff(strcmp(gratingType, 'flashed'), '(Flashed gratings)', '(Drifting gratings)');
            add_fg = @(s) {s, fg_str};
        else
            add_fg = @(s) s;
        end
        grating_offset = switchh(gratingType, {'drifting', 'flashed'}, [0, 2]);
        
        if printWith2RowsForFlashedDrifting
            subRows = 2;
            row_idx = switchh(gratingType, {'drifting', 'flashed'}, [1, 2]);
            if row_idx == 1
                clf_ifFirstRow = @() clf;
            else
                clf_ifFirstRow = @() [];
            end
        else
            subRows = 1;
            row_idx = 1;
        end
        
        %%
        
        oriWidthScatter_figId = 3;
        oriWidthCumProb_figId = 4;
        oriPrefScatter_figId = 5;
        oriPrefProbDist_figId = 6;
        DSI_figId = 7;
        dirPrefScatter_figId = 8;
        dirPrefProbDist_figId = 9;
        spfWidth_figId = 10;
        spfPref_figId = 11;
        %%
        
        if doPlots && doMainFigurePlots && any(doFigs == 1) && doOri
            %%
            % FIGURE --> scatter plots of Global/local ORI width
            w_ori_global1 = [allOriUnitStats(Wcc_pairs_ori(:,1)).w_ori_global];
            w_ori_global2 = [allOriUnitStats(Wcc_pairs_ori(:,2)).w_ori_global];
            w_ori_local1 = [allOriUnitStats(Wcc_pairs_ori(:,1)).w_ori_local];
            w_ori_local2 = [allOriUnitStats(Wcc_pairs_ori(:,2)).w_ori_local];

            [w_ori_global1, w_ori_global2] = xyyx(w_ori_global1, w_ori_global2);
            [w_ori_local1,  w_ori_local2]  = xyyx(w_ori_local1, w_ori_local2);

            if subtractSpont
                spont_str = iff(labelSpontIfSubtracted, ' (Spont Subtr.)', '');    
            else
                spont_str = iff(labelSpontIfIncluded,   ' (Spont Incl.)', '');                    
            end
            ori_w_types = {'Global orientation width', 'Local orientation width'};
%             ori_w_types = {'w_{ORI}^{Global}', 'w_{ORI}^{Local}'};
            xylims = [0 61];
            figure(oriWidthScatter_figId); clf_ifFirstRow(); % plots of orientation width 
            if printPlotsToFilesForPaper
                set(gcf, 'position', [fig_LB + (fig_offset * oriWidthScatter_figId), [fig_width_row, fig_height_row*2]], 'color', 'w')
            end                 
            
            h_ax1(1) = subplotGap(subRows,2, row_idx, 1, subSpcM, subSpcN);  plot(w_ori_global1, w_ori_global2, [scat_col mk], 'markersize', sz); 
            ht1(1) = title(add_fg(sprintf('%s%s', ori_w_types{1}, spont_str)));   
            xlabel('Cell 1', 'fontsize', xlabel_fsize); ylabel('Cell 2', 'fontsize', ylabel_fsize); axis([xylims xylims]); axis square; 
            h_ax1(2) = subplotGap(subRows,2, row_idx, 2, subSpcM, subSpcN);  plot(w_ori_local1,  w_ori_local2,  [scat_col mk], 'markersize', sz); 
            ht1(2) = title(add_fg(sprintf('%s%s', ori_w_types{2}, spont_str)));   
            xlabel('Cell 1', 'fontsize', xlabel_fsize); ylabel('Cell 2', 'fontsize', ylabel_fsize); axis([xylims xylims]); axis square; 
            
            set(ht1, 'fontsize', title_fsize);
            if addSubplotLetters
                h_let(1) = addSubplotLetter(subRows, 2, row_idx, 1, subSpcM, subSpcN, char('A'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont );
                h_let(2) = addSubplotLetter(subRows, 2, row_idx, 2, subSpcM, subSpcN, char('B'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont );            
            end
            3;
            
            3;
        end    
        
        
%         dOriW_local_bcc = S_ori{ori_loc_idx}.val(Bcc_oo_pairIdxs);
%         dOriW_local_wcc = S_ori{ori_loc_idx}.val(Wcc_oo_pairIdxs);
% 
%         dOriW_global_bcc = S_ori{ori_glob_idx}.val(Bcc_oo_pairIdxs);
%         dOriW_global_wcc = S_ori{ori_glob_idx}.val(Wcc_oo_pairIdxs);
        
        
        if doPlots && doMainFigurePlots && any(doFigs == 2) && doOri
        % FIGURE 2 --> cumulative differences in Global/local ORI width
            %%
        
            plots = {'Dw_ori_glob', 'Dw_ori_loc'};
            fig2_ms_idxs = [ori_glob_idx, ori_loc_idx];
%             ori_w_types = {'w_{ORI}^{Global}', 'w_{ORI}^{Local}'};
            ori_w_types = {'Global orientation width', 'Local orientation width'};
            
            ori_w_x_labels = {'Diff. in w_{ORI}^{Global}, degrees', 'Diff. in w_{ORI}^{Local}, degrees'};
            
            if subtractSpont
                spont_str = iff(labelSpontIfSubtracted, ' (Spont Subtracted)', '');    
            else
                spont_str = iff(labelSpontIfIncluded,   ' (Spont Included)', '');                    
            end
            Wcc_vals = {dOriW_global_wcc, dOriW_local_wcc};                
            Bcc_vals = {dOriW_global_bcc, dOriW_local_bcc};
                        
            nPlots = 2;
            figure(oriWidthCumProb_figId); clf_ifFirstRow();
            if printPlotsToFilesForPaper
                set(gcf, 'position', [fig_LB + (fig_offset * oriWidthCumProb_figId), [fig_width_row, fig_height_row*2]], 'color', 'w')
            end                 

            binIdCloseTo1 = 0;
            maxDiffTh = 0.98;
            [h_ax2, h_ax2_inset] = deal( zeros(1,nPlots) );
            [Wcc_binVals, Bcc_binVals] = deal( cell(1, nPlots) );            
            

            
            for plot_i = 1:nPlots
                %%
                ori_w_type_idx = plot_i;
                subM = 1;
                ms_idx = fig2_ms_idxs(plot_i);
                h_ax2(plot_i) = subplotGap(subRows,2, row_idx, plot_i, subSpcM, subSpcN);
                if addSubplotLetters
                    h_let(plot_i) = addSubplotLetter(subRows,2, row_idx, plot_i, subSpcM, subSpcN, char('A'+grating_offset+plot_i-1), subplotLetter_offset, 'fontname', subplotLetFont);
                end

    %             binEdges = S_ori{ms_idx}.binEdges;
    %             nBins = length(binEdges)-1;
                binEdges = [0:2.5:90];
                binCents = binEdge2cent(binEdges);
    
                [Wcc_cumFracBins, Wcc_binVals{plot_i}] = cumFracBinsFromVals(Wcc_vals{ori_w_type_idx}, binEdges);
                [Bcc_cumFracBins, Bcc_binVals{plot_i}] = cumFracBinsFromVals(Bcc_vals{ori_w_type_idx}, binEdges);            

                binIdCloseTo1 = max([binIdCloseTo1, find(Wcc_cumFracBins >= maxDiffTh, 1), find(Bcc_cumFracBins >= maxDiffTh, 1)] );
    %             Wcc_cumFracBins2 = cumFracBinsFromBinIds(S_ori{ms_idx}.bins(Wcc_oo_pairIdxs), nBins);
    %             Bcc_cumFracBins2 = cumFracBinsFromBinIds(S_ori{ms_idx}.bins(Bcc_oo_pairIdxs), nBins);
        
                plot(binEdges, [0 Wcc_cumFracBins], 'o-',              'color', w_col, 'lineWidth', lineWidth_main, 'markersize', cumMkSize); hold on;
                plot(binEdges, [0 Bcc_cumFracBins], ['s' dashedLine] , 'color', b_col, 'lineWidth', lineWidth_main, 'markersize', cumMkSize);
                xlabel(ori_w_x_labels{ori_w_type_idx}, 'fontsize', xlabel_fsize);
                ylabel('Cumulative fraction of pairs', 'fontsize', ylabel_fsize);

                ylim([0 1]);
                title(add_fg(sprintf('%s%s', ori_w_types{ori_w_type_idx}, spont_str)), 'fontsize', title_fsize);
                if (plot_i == 2) && addLegends
                    h_leg2 = legend(within_between_legend_str, 'location', 'NW', 'fontsize', legend_fsize);
                end            

            end    

            if addProbSubplots
                drawnow;
                for plot_i = 1:nPlots
                    pos_inset = insetAxesPosition(get(h_ax2(plot_i), 'position'), [.4, .02, .6, .6]);
                    h_ax2_inset(plot_i) = axes('outerposition', pos_inset, 'ytick', [], 'fontsize', 8, 'nextplot', 'add');
                    plot(binCents, Wcc_binVals{plot_i}/sum(Wcc_binVals{plot_i}), '-', 'color', w_col, 'markersize', distMkSize, 'linewidth', lineWidth_inset);
                    plot(binCents, Bcc_binVals{plot_i}/sum(Bcc_binVals{plot_i}), ':', 'color', b_col, 'markersize', distMkSize, 'linewidth', lineWidth_inset);
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

        
        if doPlots && doMainFigurePlots && any(doFigs == 3) && doOri
            %%
           % 3a. scatterplot plot preferred directions
            ori_pref1 = [allOriUnitStats(Wcc_pairs_ori(:,1)).ori_pref_deg];
            ori_pref2 = [allOriUnitStats(Wcc_pairs_ori(:,2)).ori_pref_deg];
            oriMax = 180;
            diffGt90_12 = (ori_pref1 - ori_pref2)> oriMax/2;
            diffGt90_21 = (ori_pref2 - ori_pref1)> oriMax/2;
            ori_pref1(diffGt90_12) = ori_pref1(diffGt90_12)-oriMax;
            ori_pref1(diffGt90_21) = ori_pref1(diffGt90_21)+oriMax;

    %         ori_pref2(diffGt90_21) = ori_pref2(diffGt90_21)-180;

            [ori_pref1, ori_pref2] = xyyx(ori_pref1, ori_pref2);        
            
            figure(oriPrefScatter_figId); clf_ifFirstRow(); % plots of orientation width             
            if printPlotsToFilesForPaper
                set(gcf, 'position', [fig_LB + (fig_offset * oriPrefScatter_figId), [fig_width_row, fig_height_row*2]], 'color', 'w')
            end                             
            subplotGap(subRows,2, row_idx, 1, subSpcM, subSpcN);            

            plot(ori_pref1, ori_pref2, mk, 'markersize', sz, 'color', scat_col); 
            title(add_fg('Preferred orientation'), 'fontsize', title_fsize);
%             axis(oriMax*[-1/2, 3/2, -1/2, 3/2]);
            ori_ticks = -90 : ori_scatter_tick_spacing : 270;
%             tks = oriMax*[-1/2:1/4:3/2];
            set(gca, 'xtick', ori_ticks, 'ytick', ori_ticks, 'xlim', lims(ori_ticks), 'ylim', lims(ori_ticks));        
            xlabel('Cell 1', 'fontsize', xlabel_fsize);
            ylabel('Cell 2', 'fontsize', ylabel_fsize);
            axis square;

            hax3b = subplotGap(subRows,2, row_idx ,2, subSpcM, subSpcN);
            
            3;
            binEdges = S_ori{ori_pref_idx}.binEdges;
            binCents = binEdge2cent(binEdges);
            nBins = length(binEdges)-1;


            [Wcc_cumFracBins, Wcc_binVals] = cumFracBinsFromBinIds(S_ori{ori_pref_idx}.bins(Wcc_oo_pairIdxs), nBins);
    %         if showBccPairs
            [Bcc_cumFracBins, Bcc_binVals] = cumFracBinsFromBinIds(S_ori{ori_pref_idx}.bins(Bcc_oo_pairIdxs), nBins);       
    %         end

            plot(binEdges, [0 Wcc_cumFracBins], 'o-',             'color', w_col, 'lineWidth', lineWidth_main, 'markersize', cumMkSize); hold on;
    %         if showBccPairs
            plot(binEdges, [0 Bcc_cumFracBins], ['s' dashedLine], 'color', b_col, 'lineWidth', lineWidth_main, 'markersize', cumMkSize);
    %         end
            xlim(binEdges([1, end]));
            set(gca, 'xtick', [0:15:90])
            title(add_fg('Preferred orientation'), 'fontsize', title_fsize);
            xlabel('Diff. in pref. orientation, degrees', 'fontsize', xlabel_fsize);
            ylabel('Cumulative fraction of pairs', 'fontsize', ylabel_fsize);
            if addLegends
                legend(within_between_legend_str, 'location', 'SE', 'fontsize', legend_fsize);
            end

            if addSubplotLetters
                addSubplotLetter(subRows,2, row_idx, 1, subSpcM, subSpcN, char('A'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont);
                addSubplotLetter(subRows,2, row_idx, 2, subSpcM, subSpcN, char('B'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont);            
            end
            %%
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

        %%
      
        
        % figure 4 and diff ori stats:        
        th_deg = 60;
        nPair_tot_gtTh = nnz(dOriPref_wcc > th_deg);         
        nPair_norm_gtTh = nnz(dOriPref_wcc_norm > th_deg);
        nPair_out_gtTh = nPair_tot_gtTh-nPair_norm_gtTh;
        pctPair_out_gtTh = nPair_out_gtTh/nPair_tot_gtTh*100;
                
        miscPairStats.outlierStats.nPairsGt60Outliers = nPair_out_gtTh;       
        miscPairStats.outlierStats.nPairsGt60Total = nPair_tot_gtTh;

        
        %%
        if showWorking
            fprintf('\n\n*** Outliers in Diff Ori Distribution (%s gratings, %s):\nOutlier cells contribute to %d out of %d (%.1f%%) of the pairs above %d degrees: \n', ...
                gratingType, spont_s, nPair_out_gtTh, nPair_tot_gtTh, pctPair_out_gtTh, th_deg);     
        end
        
        if doPlots && doMainFigurePlots && (any(doFigs == 4) && doOri) 
%             pref_dir_idx = find(strcmp(measures_ori, 'D_dir_pref'));
            %%
            figure(oriPrefProbDist_figId); clf_ifFirstRow();        
            if printPlotsToFilesForPaper
                set(gcf, 'position', [fig_LB + (fig_offset * oriPrefProbDist_figId), [fig_width_row, fig_height_row*2]], 'color', 'w')
            end                             

            subplotGap(subRows,2, row_idx, 1, subSpcM, subSpcN);
            ori_pref_idx = find(strcmp(measures_ori, 'D_ori_pref'), 1);      
            if addSubplotLetters
                addSubplotLetter(subRows,2, row_idx, 1, subSpcM, subSpcN, char('A'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont);
            end
            
            
            binE = S_ori{ori_pref_idx}.binEdges;
            binC = binEdge2cent(binE);
            dOriPref_cc = S_ori{ori_pref_idx}.val;
            dOriPref_wcc = dOriPref_cc(Wcc_oo_pairIdxs);
            binVals_cc_all = histcnt(dOriPref_cc(Wcc_oo_pairIdxs), binE);

            binVals_cc_norm = histcnt(dOriPref_cc(Wcc_oo_pairIdxs_norm), binE);
            binVals_cc_1outlier = histcnt(dOriPref_cc(Wcc_oo_pairIdxs_1outlier), binE);
            binVals_cc_2outliers = histcnt(dOriPref_cc(Wcc_oo_pairIdxs_2outliers), binE);  binVals_cc_2outliers = binVals_cc_2outliers(:);                             
%             assert(all (binVals_cc_all == binVals_cc_norm + binVals_cc_1outlier + binVals_cc_2outliers ));
                        
            
            
            
            
            h_bar4a = bar(binC, [binVals_cc_norm, (binVals_cc_1outlier + binVals_cc_2outliers)], 1, 'stacked');
            set(h_bar4a(1), 'facecolor', bar_norm);
            set(h_bar4a(2), 'facecolor', bar_out);
            
            xlim(binE([1, end]));
            set(gca, 'xtick', [0:15:90]);
            title(add_fg('Pref. orientation: pairwise differences'), 'fontsize', title_fsize);
            xlabel('Diff. in pref. orientation, degrees', 'fontsize', xlabel_fsize);
            ylabel('Number of cell-cell pairs', 'fontsize', ylabel_fsize);
            legend({'Typical cells', 'Outlier cells'}, 'fontsize', legend_fsize);
            if dispPlotStats
                %%
                          
            end

%             medians = nan(1, 6);
%             diams = nan(1,6);
%             diams(2) = 12.7; diams(4) = 12.7; diams(5) = 25;
            %%
%                 allCellLocData = [allOriUnits(oriCells_use).locData];
%                 allCellElecType = {allCellLocData.ElectrodeType};
%                 [uElecType, elec_idxs] = uniqueList(allCellElecType);
%                 
%                 elec_idx = 1;
%                 idxs_use = elec_idxs{elec_idx};
% 
%             figure(55); clf;
            subplotGap(subRows,2, row_idx,2, subSpcM, subSpcN);                    
            dori_mu_all = [allOriCellStats.(dOriMU_field)]; 
            dori_mu_norm = dori_mu_all(~ori_cell_is_outlier);
            dori_mu_out  = dori_mu_all(ori_cell_is_outlier);
            binVals_cm_norm     = histcnt(dori_mu_norm, binE);
            binVals_cm_outliers = histcnt(dori_mu_out, binE);            
            
            h_bar4b = bar(binC, [binVals_cm_norm(:), binVals_cm_outliers(:)], 1, 'stacked');
            set(h_bar4b(1), 'facecolor', bar_norm);
            set(h_bar4b(2), 'facecolor', bar_out);                        
            
            xlim(binE([1, end]));
            set(gca, 'xtick', [0:15:90]);
            title(add_fg('Pref. orientation: diff. from multi-units'), 'fontsize', title_fsize);
            xlabel('Diff. in pref. orientation, degrees', 'fontsize', xlabel_fsize);
            ylabel('Number of cell-multiunit pairs', 'fontsize', ylabel_fsize);            
            if addSubplotLetters
                addSubplotLetter(subRows,2, row_idx, 2, subSpcM, subSpcN, char('B'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont);
            end
            
            
%             showMedian = 0;
            
                
                
                
%                 allUnitLocData = [allOriUnits.locData];
%                 allUnitElecData = {allUnitLocData.ElectrodeType};
%                 allPairElecType = allUnitElecData(Wcc_pairs_ori(:,1));
%                 allUnitLocData = [allOriUnits.locData];
                
%                 med_c_mu = median(dori_mu_all(idxs_use));
%                 drawVerticalLine(med_c_mu, 'linestyle', ':', 'color', 'r')
%                 title(sprintf('elecType = %s (N = %d). median = %.2f', uElecType{elec_idx}, length(idxs_use), med_c_mu));
                
                
            %%
    %         
    %         binVals_cm2 = histcnt([allOriUnitStats.Dori_pref_smlSpkMU], binE);
    %         subplot(1,3,3);        
    %         bar(binC, binVals_cm2, 1);
    %         xlim(binE([1, end]));
    %         title('Preferred Orientation: Differences from Multi-units');        
    %         xlabel('Diff. in Preferred Orientation, degrees');
    %         ylabel('Number of Cell pairs');
    %         
            3;

            estimateShapeOfMUdistribution = 0;
            if estimateShapeOfMUdistribution
                %%
                binE = 0:5:90;
                binC = binEdge2cent(binE);
                dori_mu_all = [allOriUnitStats(oriCells_use).(dOriMU_field)];
                binVals_cm_all  = histcnt(dori_mu_all, binE);
                binVals_cm_all = binVals_cm_all/(sum(binVals_cm_all)*diff(binE(1:2)));
                
                
                % fit exponential curve
                f_gauss = @(beta, x) (1./beta(1))*exp(-x/beta(1));
                beta0 = 1./[binVals_cm_all(1)];
                b = nlinfit(binC(:), binVals_cm_all(:), f_gauss, beta0);
                
                b2 = expfit(binVals_cm_all);

                
                %
                b_val = 11.7;
                diffFromMU = -b_val*log(1-rand(1,10000) );
                n_samp = histcnt(dori_mu_all, binE);
                n_samp = n_samp/(sum(n_samp)*diff(binE(1:2)));                
                
                figure(102); clf; hold on;
                plot(binC, binVals_cm_all, 'o-');                                
                plot(binC, f_gauss(b, binC), 'r:');
                plot(binC, n_samp, 'gs-');
                b, pearsonR(binVals_cm_all(:), f_gauss(b, binC))^2
                
                
                
            end
            
            showComparisonToAlbusFigure = 0;
            if showComparisonToAlbusFigure
                %%
                figure(504); clf;
                binE = [0:15:90];
                binC = binEdge2cent(binE);
                
                binVals_cc_all = histcnt(dOriPref_cc(Wcc_oo_pairIdxs), binE);                        
                binVals_cc_all_pct = binVals_cc_all / sum(binVals_cc_all) * 100;
                        
                h_albus = bar(binC, binVals_cc_all_pct, 1);
                set(h_albus(1), 'facecolor', 'none', 'linewidth', 2);
            
                xlim(lims(binE, .01));
                set(gca, 'xtick', [0:30:90]);                
                title('Differences in pref. ori');
                ylabel('% of pairs');
                set(gca, 'tickDir', 'out');
                box off;
                3;
                
            end
            3;


        end

        
        if doPlots && doMainFigurePlots && any(doFigs == 5) && strcmp(gratingType, 'drifting') && doOri
            %% (5a) scatterplots of DSI's
            dsi1 = [allOriUnitStats(Wcc_pairs_ori(:,1)).DSI_global];
            dsi2 = [allOriUnitStats(Wcc_pairs_ori(:,2)).DSI_global];
            [dsi1, dsi2] = xyyx(dsi1, dsi2);

            if subtractSpont
                spont_str = iff(labelSpontIfSubtracted, ' (Spont Subtracted)', '');    
            else
                spont_str = iff(labelSpontIfIncluded,   ' (Spont Included)', '');                    
            end            
%             dsi_str = 'DSI';
            dsi_str = 'Direction selectivity index';
            
            figure(DSI_figId); clf;
            if printPlotsToFilesForPaper
                set(gcf, 'position', [fig_LB + (fig_offset * DSI_figId) + [0, fig_height_row], [fig_width_row, fig_height_row]], 'color', 'w')
            end                             

            
%             h_ax5(1) = subplotGap(1,subN, 1);  
            h_ax5(1) = subplotGap(1,2, 1, 1);              
            plot(dsi1, dsi2, mk, 'markersize', sz, 'color', scat_col); ht5(1) = title({[], sprintf('%s%s', dsi_str, spont_str)});
            xlabel('Cell 1', 'fontsize', xlabel_fsize); 
            ylabel('Cell 2', 'fontsize', ylabel_fsize); 
            axis([0 1 0 1]); axis square;
%             if subtractSpont
%                 h_ax5(2) = subplotGap(1,subN, 2);  plot(dsi_ss1, dsi_ss2, mk, 'markersize', sz); ht5(2) = title(sprintf('%s Spont Subtracted', dsi_str)); xlabel('Cell 1'); ylabel('Cell 2'); axis([0 1 0 1]); axis square;
%             end
            set(ht5, 'fontsize', title_fsize);
            
            
            %%% (5b) CDF of DSI
            
%             DSI_measures = {'D_dsi'};
%             dsi_idx = find(strcmp(measures_ori, 'D_dsi_glob'), 1);


            nPlots = 1;
            [h_ax7, h_ax7_inset] = deal( zeros(1,nPlots) );
            [Wcc_binVals, Bcc_binVals] = deal(cell(1,nPlots));
            
            for i = 1:nPlots                
            
                binEdges = [0:.05:1]; %S_ori{dsi_idx}.binEdges;
                binCents = binEdge2cent(binEdges);
        %         nBins = length(binEdges)-1;

                [Wcc_cumFracBins, Wcc_binVals{i}] = cumFracBinsFromVals(dDSI_wcc, binEdges);
                [Bcc_cumFracBins, Bcc_binVals{i}] = cumFracBinsFromVals(dDSI_bcc, binEdges);
           
            
%                 h_ax7(i) = subplotGap(1,nPlots,i);                
                h_ax7 = subplotGap(1,2, 1, 2);            
                            
                plot(binEdges, [0 Wcc_cumFracBins], 'o-',             'color', w_col, 'linewidth', lineWidth_main, 'markersize', cumMkSize); hold on;
                plot(binEdges, [0 Bcc_cumFracBins], ['s' dashedLine], 'color', b_col, 'linewidth', lineWidth_main, 'markersize', cumMkSize); 
                xlim(binEdges([1, end]));
                xlabel('Difference in DSI', 'fontsize', xlabel_fsize);
                ylabel('Cumulative fraction of pairs', 'fontsize', ylabel_fsize);
%                 ht7 = title({[], sprintf('Differences in DSI%s', spont_str)});
                ht7 = title({[], sprintf(' %s', spont_str)});
                set(ht7, 'fontsize', title_fsize)

                if (i == nPlots) && addLegends
                    h_leg7 = legend(within_between_legend_str, 'location', 'NW', 'fontsize', legend_fsize);
                end
                
                
                
            end
            
            if addProbSubplots
                drawnow;
                for plot_i = 1:nPlots
                    pos_inset = insetAxesPosition(get(h_ax7(plot_i), 'position'), [.4, .02, .6, .6]);
                    h_ax7_inset(plot_i) = axes('outerposition', pos_inset, 'ytick', [], 'fontsize', 8, 'nextplot', 'add');
                    plot(binCents, Wcc_binVals{plot_i}/sum(Wcc_binVals{plot_i}), '-', 'color', w_col, 'markersize', distMkSize, 'linewidth', lineWidth_inset);
                    plot(binCents, Bcc_binVals{plot_i}/sum(Bcc_binVals{plot_i}), ':', 'color', b_col, 'markersize', distMkSize, 'linewidth', lineWidth_inset);                
                                     
                    axis(h_ax7_inset(plot_i), 'tight'); 
                    xlim([0 1]);
                    ylims = get(h_ax7_inset(plot_i), 'ylim');
                    set(h_ax7_inset(plot_i), 'ylim', [0, ylims(2)*1.05]);
                    
                    if (plot_i == nPlots) && addLegends
                        ax_inset_pos = get(h_ax7_inset(plot_i), 'position');                        
                        leg_pos = get(h_leg7, 'position');
                        LB = [ax_inset_pos(1), ax_inset_pos(2)+ax_inset_pos(4)+.02];
%                         set(h_leg7, 'position', [LB, leg_pos(3:4)]);                        
                    end

                    
                end
            end
            if addSubplotLetters
                addSubplotLetter(1,2, 1, 1, subSpcM, subSpcN, 'A', subplotLetter_offset, 'fontname', subplotLetFont);
                addSubplotLetter(1,2, 1, 2, subSpcM, subSpcN, 'B', subplotLetter_offset, 'fontname', subplotLetFont);
            end
            
            

            
        end



        if doPlots && doMainFigurePlots && any(doFigs == 6) && strcmp(gratingType, 'drifting') && doSpf
            %% 6a. scatterplot plot preferred directions
            dirMax = 360;
            dir_pref1 = [allOriUnitStats(Wcc_pairs_ori(:,1)).dir_pref_deg];
            dir_pref2 = [allOriUnitStats(Wcc_pairs_ori(:,2)).dir_pref_deg];

            diffGt90_12 = (dir_pref1 - dir_pref2)> dirMax/2;
            diffGt90_21 = (dir_pref2 - dir_pref1)> dirMax/2;
            dir_pref1(diffGt90_12) = dir_pref1(diffGt90_12)-dirMax;
            dir_pref1(diffGt90_21) = dir_pref1(diffGt90_21)+dirMax;

    %         dir_pref2(diffGt90_21) = dir_pref2(diffGt90_21)-180;

            [dir_pref1, dir_pref2] = xyyx(dir_pref1, dir_pref2);        

            figure(dirPrefScatter_figId); clf; % plots of orientation width 
            if printPlotsToFilesForPaper
                set(gcf, 'position', [fig_LB + (fig_offset * dirPrefScatter_figId)  + [0, fig_height_row], [fig_width_row, fig_height_row]], 'color', 'w')
            end                             
            subplotGap(1,2, 1, 1, subSpcM, subSpcN);            
            plot(dir_pref1, dir_pref2, mk, 'markersize', sz, 'color', scat_col); 
            title({[], 'Preferred direction'}, 'fontsize', title_fsize);
%             axis(dirMax*[-1/2, 3/2, -1/2, 3/2]);
%             tks = dirMax*[-1/2:1/4:3/2];        
            dir_ticks = -180 : ori_scatter_tick_spacing*2 : 540;
%             tks = oriMax*[-1/2:1/4:3/2];
            set(gca, 'xtick', dir_ticks, 'ytick', dir_ticks, 'xlim', lims(dir_ticks), 'ylim', lims(dir_ticks));        

%             set(gca, 'xtick', tks, 'ytick', tks);
            xlabel('Cell 1', 'fontsize', xlabel_fsize); 
            ylabel('Cell 2', 'fontsize', ylabel_fsize);
            axis square;
           3;
           
           % 6b. cumulative difference between preferred directions 
            subplotGap(1,2, 1, 2, subSpcM, subSpcN);
            dir_pref_idx = find(strcmp(measures_ori, 'D_dir_pref'));
            
            dDirPref_cc = S_ori{dir_pref_idx}.val;
            dDirPref_wcc = dDirPref_cc(Wcc_oo_pairIdxs);
%             binVals_cc_all = histcnt(dDirPref_cc(Wcc_oo_pairIdxs), binE);

            
            3;
%             binEdges = S_ori{dir_pref_idx}.binEdges;
            binEdges = [0:10:180];
    
%             allBinIds = S_ori{dir_pref_idx}.bins(Wcc_oo_pairIdxs);
%             cumFracBins = cumFracBinsFromBinIds(allBinIds, length(binEdges)-1);
                    

            Wcc_cumFracBins = cumFracBinsFromVals(dDirPref_cc(Wcc_oo_pairIdxs), binEdges);
            Bcc_cumFracBins = cumFracBinsFromVals(dDirPref_cc(Bcc_oo_pairIdxs), binEdges);            
            
            plot(binEdges, [0 Wcc_cumFracBins], 'o-',             'color', w_col, 'linewidth', lineWidth_main, 'markersize', cumMkSize); hold on;
            plot(binEdges, [0 Bcc_cumFracBins], ['s' dashedLine], 'color', b_col, 'linewidth', lineWidth_main, 'markersize', cumMkSize);
            xlim(binEdges([1, end]));
            set(gca, 'xtick', [0:30:180])
            title({[], 'Distribution of differences'}, 'fontsize', title_fsize);
            ylabel('Cumulative fraction of pairs', 'fontsize', ylabel_fsize);
            xlabel('Diff. in pref. direction, degrees', 'fontsize', xlabel_fsize);

%             title({[], 'Differences in preferred direction'}, 'fontsize', title_fsize);
            title({[], ' '}, 'fontsize', title_fsize);
            if addLegends
                legend(within_between_legend_str, 'location', 'SE', 'fontsize', legend_fsize);
            end            
            if addSubplotLetters
                addSubplotLetter(1,2, 1, 1, subSpcM, subSpcN, 'A', subplotLetter_offset, 'fontname', subplotLetFont);
                addSubplotLetter(1,2, 1, 2, subSpcM, subSpcN, 'B', subplotLetter_offset, 'fontname', subplotLetFont);
            end
            
            3;
           
           

        end

        if strcmp(gratingType, 'drifting')
            pref_dir_idx = find(strcmp(measures_ori, 'D_dir_pref'), 1);
            dDirPref_cc = S_ori{pref_dir_idx}.val;
            dDirPref_wcc = dDirPref_cc(Wcc_oo_pairIdxs);
        end
        
        doFig7 = doPlots && doMainFigurePlots && any(doFigs == 7);
        if (doFig7 || doPairStats) && strcmp(gratingType, 'drifting') && doOri                        
            
            % Figure 7 & d(preferred direction) stats
            restrictCMtoCellsWithMoreThanOnePair = 0;
            
            binEdges = S_ori{pref_dir_idx}.binEdges;
            binCents = binEdge2cent(binEdges);

            
            binVals_cc_all = histcnt(dDirPref_cc(Wcc_oo_pairIdxs), binEdges);                                               
            binVals_cc_norm = histcnt(dDirPref_cc(Wcc_oo_pairIdxs_norm), binEdges);
            binVals_cc_1outlier = histcnt(dDirPref_cc(Wcc_oo_pairIdxs_1outlier), binEdges);
            binVals_cc_2outliers = histcnt(dDirPref_cc(Wcc_oo_pairIdxs_2outliers), binEdges);                                   
%             assert(all (binVals_cc_all == binVals_cc_norm + binVals_cc_1outlier + binVals_cc_2outliers )); this may no longer be true, since some pairs are excluded even though 
% cells have well defined preferred orientations, but the ori of the MU is not well-defined.;
                
            dori_mu_all = [allOriCellStats.(dDirMU_field)];
            dori_mu_norm = dori_mu_all(~ori_cell_is_outlier);
            dori_mu_out  = dori_mu_all(ori_cell_is_outlier);
            
                
                
%                     [allOriUnits, pairData_ori, S_ori, pairTypes, measures_ori, pairIdxs_ori, pairIdxList_ori, idxMtx_ori, ...
%             Wcc_oo_pairIdxs, Wcc_oo_pairIdxs_M, Wcc_pairs_ori, ...
%             Wrcc_oo_pairIdxs, Wrcc_oo_pairIdxs_M, Wrcc_pairs_ori, ...
%             Bcc_oo_pairIdxs, Bcc_oo_pairIdxs_M, Bcc_pairs_ori, ...
%             oriCells_use, allOriUnitStats];
        
               

            
            if restrictCMtoCellsWithMoreThanOnePair
                idx_oriCells_use = find(oriCells_use);

                unitIsInWccPair = false(1, length(allOriUnits));
                unitIsInWccPair(unique(Wcc_pairs_ori(:))) = 1;

                glob_idx_cells_usedInPairs = find(oriCells_use & unitIsInWccPair);

                idx_cells_usedInPairs = binarySearch(idx_oriCells_use, glob_idx_cells_usedInPairs);                    
                idx_C_use = idx_cells_usedInPairs;
            else
                idx_C_use = 1:length(allOriCellStats);
            end
            %%
            dori_mu_all = [allOriCellStats(idx_C_use).(dDirMU_field)];
            dori_mu_norm = dori_mu_all(~ori_cell_is_outlier(idx_C_use));
            dori_mu_out  = dori_mu_all(ori_cell_is_outlier(idx_C_use));


            dori_ori_all = nonnans(  dDirPref_cc(Wcc_oo_pairIdxs) );
            dori_ori_norm = nonnans( dDirPref_cc(Wcc_oo_pairIdxs_norm) );
            assert(all(ibetween( nonnans( dDirPref_cc(Wcc_oo_pairIdxs) ), 0, 180)));
            
            n_cells_tot = nnz(~isnan(dori_mu_all));
            n_cells_aligned_MU = nnz( dori_mu_all < 45);  
            n_cells_anti_aligned_MU = nnz( dori_mu_all > 135);                
            n_cells_a_aa = n_cells_aligned_MU + n_cells_anti_aligned_MU; % == length(dori_mu_norm);
            n_cells_unaligned_MU = nnz( ibetween(dori_mu_all, [45, 135])); % == length(dori_mu_out);
            n_cell_pairs_tot = length(dori_ori_all);
            n_cell_pairs_a_aa = length(dori_ori_norm);

            n_a_o_u = [n_cells_aligned_MU, n_cells_anti_aligned_MU, n_cells_unaligned_MU];
            assert(sum(n_a_o_u) == n_cells_tot);

            pct_i = @(n) n/n_cells_tot*100;
            pct_e = @(n) n/n_cells_a_aa*100;
            n_within_30_MU     = nnz( dori_mu_all < 30) ;
            n_within_30_opp_MU = nnz( dori_mu_all > 150) ;
            n_within_30_45_MU       = nnz( ibetween(dori_mu_all, 30, 45) );
            n_within_135_150_MU       = nnz( ibetween(dori_mu_all, 135, 150) );
            p_cell_aligned = n_cells_aligned_MU/n_cells_tot;
            p_cell_anti_aligned = n_cells_anti_aligned_MU/n_cells_tot;
%             p_cc_a_aa_gt90_exp = 2*p_cell_aligned*p_cell_anti_aligned;

%                 nnz(idx_0outliers)
            pctPair_gt90_obs_all = nnz( dori_ori_all > 90)/length(dori_ori_all) * 100;
            pctPair_gt90_obs_norm = nnz( dori_ori_norm > 90)/length(dori_ori_norm) * 100;

            n_cc_a_aa_gt90_obs = nnz( dori_ori_norm > 90);
%                 ci = bayesConfInt(n_cells_aligned_MU, n_cells_a_aa, .95);
%                 assert(sum([n_cells_aligned_MU, n_cells_anti_aligned_MU n_cells_unaligned_MU]) == n_cells_tot);
            ci = bayesConfInt([n_cells_aligned_MU, n_cells_anti_aligned_MU], n_cells_tot, .95);

%                 pct_cc_gt90 = nnz( dDirPref_cc(Wcc_oo_pairIdxs) > 90 ) / length(Wcc_oo_pairIdxs) * 100;

%             assert(nnz(idx_noOutliers) + nnz(idx_1outlier) + nnz(idx_2outliers) == length(idx_noOutliers))

            p_aligned_ci = ci(1,:)/n_cells_tot;
            p_anti_aligned_ci = ci(2,:)/n_cells_tot;
            p_unaligned_ci = ci(3,:)/n_cells_tot;

%             n_cc_a_aa_gt90_exp = n_cell_pairs_a_aa* p_cc_a_aa_gt90_exp;
%             probObserving_n_CC_pairs_gt90 = 1-binocdf(n_cc_a_aa_gt90_obs, n_cell_pairs_a_aa, p_cc_a_aa_gt90_exp );
                
            dDirMUStats = struct('nDCellMU_0_30', n_within_30_MU, 'nDCellMU_30_45', n_within_30_45_MU, ...
                               'nDCellMU_135_150', n_within_135_150_MU, 'nDCellMU_0_30_opp', n_within_30_opp_MU, ...
                               'nAligned', n_cells_aligned_MU, 'nAntiAligned', n_cells_anti_aligned_MU, 'nUnaligned', n_cells_unaligned_MU, ...
                               'pAligned_ci', p_aligned_ci, 'pAntiAligned_ci', p_anti_aligned_ci, 'pUnaligned_ci', p_unaligned_ci);
            miscPairStats.dDirMUStats = dDirMUStats;
                            
            
            
            if dispPlotStats ;
                %%
                str{1} = sprintf('\n\n***Differences in Direction Statistics (%s gratings, %s)', gratingType, spont_s);
                str{2} = sprintf('Total # of cells: %d. Without outliers: %d.  (Number of outliers : %d)', n_cells_tot, n_cells_a_aa, n_cells_unaligned_MU);
                str{3} = sprintf('Cells within 30 of site preferred: %d : %.1f%% (%.1f%% with outliers excluded.)', n_within_30_MU, pct_i(n_within_30_MU), pct_e(n_within_30_MU) );
                str{4} = sprintf('Cells within 30 of site opposite: %d : %.1f%% (%.1f%% with outliers excluded.)', n_within_30_opp_MU, pct_i(n_within_30_opp_MU), pct_e(n_within_30_opp_MU) );
                str{5} = sprintf('Cells within 30-45 of site preferred: %d : %.1f%% (%.1f%% with outliers excluded.)', n_within_30_45_MU, pct_i(n_within_30_45_MU), pct_e(n_within_30_45_MU) );
                str{6} = sprintf('Cells within 135-150 of site opposite: %d : %.1f%% (%.1f%% with outliers excluded.)', n_within_135_150_MU, pct_i(n_within_135_150_MU), pct_e(n_within_135_150_MU) );
                str{7} = sprintf('Outlier cells (between 45 and 135): %d (%.2f%%)', n_cells_unaligned_MU, pct_i(n_cells_unaligned_MU));
                
                str{8} = sprintf('Cells aligned: %d. (%.2f%%). Anti-aligned : %d (%.2f%%)', n_cells_aligned_MU, pct_i(n_cells_aligned_MU), n_cells_anti_aligned_MU, pct_i(n_cells_anti_aligned_MU) );
%                 str{9} = sprintf('N pairs greater than 90 apart (expected): 2 x %.3f x %.3f = %.2f%% ', p_cell_aligned, p_cell_anti_aligned, p_cc_a_aa_gt90_exp*100 );
%                 str{10} = sprintf('N pairs greater than 90 apart (observed): (%.1f%%). (outliers excluded) (%.1f%%) (outliers included)', pctPair_gt90_obs_norm, pctPair_gt90_obs_all );
                str{9} = sprintf('*Aligned & anti-aligned with multiunits:');
                str{10} = sprintf('Empirically: %d are aligned, %d are anti-aligned and %d are unaligned (Tot = %d)',n_cells_aligned_MU, n_cells_anti_aligned_MU, n_cells_unaligned_MU, n_cells_tot);
                str{11} = sprintf('Empirical probabilities: aligned: p = %.1f, anti-aligned, p = %.1f. unaligned, p = %.1f', n_a_o_u/ n_cells_tot);                
                str{12} = sprintf('Confidence intervals: aligned: [%.1f - %.1f%%], anti-aligned: [%.1f-%.1f%%], unaligned [%.1f-%.1f%%]',ci'/n_cells_tot*100);
                %%
%                 fprintf('Bayes confidence interval: [%.3f, %.3f] in same direction, [%d-%d%%] in opp direction)\n', ci, round((1-ci([2, 1]))*100));
                3;
                cellfun(@(str_i) fprintf('%s\n', str_i), str);
            end
            
        
            if  doFig7
                %%

                restrictCMtoCellsWithMoreThanOnePair = 1;            


                figure(dirPrefProbDist_figId); clf;
                subplotGap(1,2, 1, 1, subSpcM, subSpcN);
                if printPlotsToFilesForPaper
                    set(gcf, 'position', [fig_LB + (fig_offset * dirPrefProbDist_figId) + [0, fig_height_row], [fig_width_row, fig_height_row]], 'color', 'w')
                end

                h_bar7a = bar(binCents, [binVals_cc_norm, (binVals_cc_1outlier + binVals_cc_2outliers)], 1, 'stacked');
                set(h_bar7a(1), 'facecolor', bar_norm);
                set(h_bar7a(2), 'facecolor', bar_out);

                set(gca, 'xtick', [0:30:180]);
                xlim(binEdges([1, end]));
                title({[], 'Pref. direction: pairwise differences'}, 'fontsize', title_fsize);        
                xlabel('Diff. in pref. direction, degrees', 'fontsize', xlabel_fsize);
                ylabel('Number of cell-cell pairs', 'fontsize', ylabel_fsize);
                legend({'Typical cells', 'Outlier cells'}, 'fontsize', legend_fsize)
%%
                if dispPlotStats && 0
                    %%
                    range_deg = 90+[-30, 30];
                    n_tot = nnz(ibetween( dDirPref_cc(Wcc_oo_pairIdxs), range_deg));
                    nPair_norm_gtTh = nnz(ibetween( dDirPref_cc(Wcc_oo_pairIdxs_norm), range_deg));
                    n_out = n_tot-nPair_norm_gtTh;
                    pct_out = n_out/n_tot*100;
                    fprintf('Outlier cells contribute to %d out of %d (%.1f%%) of the pairs in the range %d-%d: \n', ...
                        n_out, n_tot, pct_out, range_deg);               
                end
%%
                subplotGap(1,2, 1, 2, subSpcM, subSpcN);
                binVals_cm_norm     = histcnt(dori_mu_norm, binEdges);
                binVals_cm_outliers = histcnt(dori_mu_out, binEdges);

                h_bar7b = bar(binCents, [binVals_cm_norm(:), binVals_cm_outliers(:)], 1, 'stacked');
                set(h_bar7b(1), 'facecolor', bar_norm);
                set(h_bar7b(2), 'facecolor', bar_out);

                xlim(binEdges([1, end]));
                set(gca, 'xtick', [0:30:180]);
                title({[], 'Pref. direction: diff. from multi-units'}, 'fontsize', title_fsize);        
                xlabel('Diff. in pref. direction, degrees', 'fontsize', xlabel_fsize);
                ylabel('Number of cell-multiunit pairs', 'fontsize', ylabel_fsize);
                if addSubplotLetters
                    addSubplotLetter(1,2, 1, 1, subSpcM, subSpcN, 'A', subplotLetter_offset, 'fontname', subplotLetFont);
                    addSubplotLetter(1,2, 1, 2, subSpcM, subSpcN, 'B', subplotLetter_offset, 'fontname', subplotLetFont);
                end
                3;
                

            end
        end

% 
%         nwcc = length(spf_pairTypeIdxs{1});
%         n_sc_obs = frac_sc{1}*nwcc;
%         n_sc_exp = frac_sc_exp*nwcc;
%         
%         p_sc_pair = 2*p_simple*p_complex;
%         p_sc_obs = binocdf(n_sc_obs, nwcc, p_sc_pair);
%                 

        if doPlots && doMainFigurePlots && any(doFigs == 8) && doSpf
            %%
            % 8a. Scatterplot of spatial frequency tuning widths
            w_spf1 = [allSpfUnitStats(Wcc_pairs_spf(:,1)).w_spf];
            w_spf2 = [allSpfUnitStats(Wcc_pairs_spf(:,2)).w_spf];

            [w_spf1, w_spf2] = xyyx(w_spf1, w_spf2);
            idx_use = ~isnan(w_spf1) & ~isnan(w_spf2);
            w_spf1 = w_spf1(idx_use);
            w_spf2 = w_spf2(idx_use);

            figure(spfWidth_figId); clf_ifFirstRow(); % plots of orientation width 
            if printPlotsToFilesForPaper
                set(gcf, 'position', [fig_LB + (fig_offset * spfWidth_figId), [fig_width_row, fig_height_row*2]], 'color', 'w')
            end
            h1 = subplotGap(subRows,2, row_idx,  1, subSpcM, subSpcN);
            plot(w_spf1, w_spf2, mk, 'markersize', sz, 'color', scat_col); 
            ht8a = title(add_fg('Spatial frequency tuning width'));
            set(ht8a, 'fontsize', title_fsize);
            
            xlabel('Cell 1', 'fontsize', xlabel_fsize); 
            ylabel('Cell 2', 'fontsize', ylabel_fsize);
            axis square;
            3;


            h_ax8 = subplotGap(subRows,2, row_idx,  2, subSpcM, subSpcN); 
            maxDiffTh = .99;            
            binEdges = S_spf{spf_w_idx}.binEdges;
            binEdges = binEdges(1):0.25:binEdges(end);
            binCents = binEdge2cent(binEdges);

            [Wcc_cumFracBins, Wcc_binVals] = cumFracBinsFromVals(dSpf_w_wcc, binEdges);
            [Bcc_cumFracBins, Bcc_binVals] = cumFracBinsFromVals(dSpf_w_bcc, binEdges);
            binIdCloseTo1 = max(find(Wcc_cumFracBins >= maxDiffTh, 1), find(Bcc_cumFracBins >= maxDiffTh, 1) );
    %         Wcc_cumFracBins = cumFracBinsFromBinIds(S_spf{spf_w_idx}.bins(Wcc_ss_pairIdxs), nBins);
    %         Bcc_cumFracBins = cumFracBinsFromBinIds(S_spf{spf_w_idx}.bins(Bcc_ss_pairIdxs), nBins);                    

            plot(binEdges, [0 Wcc_cumFracBins], 'o-',             'color', w_col, 'linewidth', lineWidth_main, 'markersize', cumMkSize); hold on
            plot(binEdges, [0 Bcc_cumFracBins], ['s' dashedLine], 'color', b_col, 'linewidth', lineWidth_main, 'markersize', cumMkSize); 
            
            
    %         xlim(binEdges([1, end]));
            xlims = [0, binEdges(binIdCloseTo1+1)];
            xlim(xlims);
%             title('Spatial Frequency Tuning Width');
%             ht8 = title(add_fg('Differences in spat. freq. tuning width'));
            ht8 = title(add_fg('Spatial frequency tuning width'));
            set(ht8, 'fontsize', title_fsize)
            ylabel('Cumulative fraction of pairs', 'fontsize', ylabel_fsize);
            xlabel('Diff. in spatial freq. tuning width', 'fontsize', xlabel_fsize);
            if addLegends
                h_leg8 = legend(within_between_legend_str, 'location', 'SE', 'fontsize', legend_fsize);
            end            
            
            if addProbSubplots
                drawnow;
                
                pos_inset = insetAxesPosition(get(h_ax8, 'position'), [.4, .02, .6, .6]);
                h_ax8_inset = axes('outerposition', pos_inset, 'ytick', [], 'fontsize', 8, 'nextplot', 'add');
                plot(binCents, Wcc_binVals/sum(Wcc_binVals), '-', 'color', w_col, 'markersize', distMkSize, 'linewidth', lineWidth_inset);
                plot(binCents, Bcc_binVals/sum(Bcc_binVals), ':', 'color', b_col, 'markersize', distMkSize, 'linewidth', lineWidth_inset);            
                                
                axis(h_ax8_inset, 'tight');     
                xlim(xlims);

                ylims = get(h_ax8_inset, 'ylim');
                set(h_ax8_inset, 'ylim', [0, ylims(2)*1.05]);

                if addLegends
                    ax_inset_pos = get(h_ax8_inset, 'position');
                    leg_pos = get(h_leg8, 'position');
                    LB = [ax_inset_pos(1), ax_inset_pos(2)+ax_inset_pos(4)+.02];
                    set(h_leg8, 'position', [LB, leg_pos(3:4)]);
                end
                
            end
            if addSubplotLetters
                addSubplotLetter(subRows,2, row_idx, 1, subSpcM, subSpcN, char('A'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont);            
                addSubplotLetter(subRows,2, row_idx, 2, subSpcM, subSpcN, char('B'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont);
            end
            

        end


        f_opt = [allSpfCellStats.f_opt];
        if doPairStats            
            [p25, p75] = dealV( prctile(f_opt, [25, 75]) );
            range_octaves = log2( p75 / p25);
            miscPairStats.spfCellsStats.spf_min = min(f_opt);
            miscPairStats.spfCellsStats.spf_max = max(f_opt);            
            miscPairStats.spfCellsStats.range_octaves = range_octaves;
            miscPairStats.spfCellsStats.nRatiosLT2 = nnz(dSpfPref_wcc < 1);
            
        end
        if showWorking
            fprintf('Range of spatial frequencies is [%.2f, %.2f] c/deg\n', min(f_opt), max(f_opt));
        end
        
        
        % checking whether isolation distance affects any of the measures
        % of clustering
        if doPlots && doSpikeSortingCheckFigure
            %%
            [ccs,pvals_cc, pvals_cmp] = deal(cell(1,2)); %zeros(nDists, nMsCmp));
            for os_i = 1:2
                if os_i == 1 && doOri

                    minIDinPair = pairData_ori.minID(Wcc_oo_pairIdxs);
                    GLF_overlap = pairData_ori.GLF_overlap(Wcc_oo_pairIdxs);
                    PCA_overlap = pairData_ori.PCA_overlap(Wcc_oo_pairIdxs);
                    waveform_ed = pairData_ori.fullWaveform_ed(Wcc_oo_pairIdxs);
                    negAmps_overlap = pairData_ori.negAmps_overlap(Wcc_oo_pairIdxs);

                    measures_cmp_names = {'D_ori_pref', 'Dw_ori_glob', 'Dw_ori_loc'};
                    measures_cmp = {dOriPref_wcc, dOriW_global_wcc, dOriW_local_wcc};

                    if strcmp(gratingType, 'drifting')
                        measures_cmp_names = [measures_cmp_names, 'D_dsi_glob'];
                        measures_cmp = [measures_cmp, dDSI_wcc];
                    end
                    
                elseif os_i == 2 && doSpf
                     minIDinPair = pairData_spf.minID(Wcc_ss_pairIdxs);
                    GLF_overlap = pairData_spf.GLF_overlap(Wcc_ss_pairIdxs);
                    PCA_overlap = pairData_spf.GLF_overlap(Wcc_ss_pairIdxs);
                    waveform_ed = pairData_spf.fullWaveform_ed(Wcc_ss_pairIdxs);
                    negAmps_overlap = pairData_spf.negAmps_overlap(Wcc_ss_pairIdxs);
                    %%
                    
                    measures_cmp_names = {'D_spfW', 'D_spf_pref'};
                    measures_cmp = {dSpf_w_wcc, dSpfPref_wcc};
                end
                                
                dists = {minIDinPair, GLF_overlap, waveform_ed };..., PCA_overlap, negAmps_overlap};
                distNames = {'minID', 'GLFoverlap', 'waveform_ed' }; ..., 'PCAoverlap', 'negAmpsOverlap'};
                nDists = length(dists);

                    
                nMsCmp = length(measures_cmp);
                                
%                 ms_idx_use = 1;
%                 measures_ori_cmp = measures_ori_cmp(1);

                
                
                gratingOffset = iff(strcmp(gratingType, 'flashed'), 100, 0);
                pctile_frac = 1/4;
                show = 1;
                figure(350+ os_i + gratingOffset); clf;
                
                for di = 1:nDists
                    for mi = 1:nMsCmp          
                        
                        
                        [ccs{os_i}(di, mi), pvals_cc{os_i}(di, mi)] = corr(  dists{di}(:), measures_cmp{mi}(:), 'type', 'spearman', 'rows', 'complete', 'tail', 'right');                        
                                   
                        dist_pctles = prctile(dists{di}, [pctile_frac, 1-pctile_frac]*100);
                        idx_lower_frac = find(dists{di} < dist_pctles(1));
                        idx_higher_frac = find(dists{di} > dist_pctles(2));
                        [pvals_cmp{os_i}(di, mi), ~, stats(di,mi)] = ranksum(measures_cmp{mi}(idx_lower_frac), measures_cmp{mi}(idx_higher_frac), 'tail', 'left');

                        
                        if show
                            subplotGap(nMsCmp, nDists,  mi, di); hold on;
                            x = prctile(dists{di}, [pctile_frac/2, (1-pctile_frac/2)]*100);
                            y = [nanmedian(measures_cmp{mi}(idx_lower_frac)), nanmedian(measures_cmp{mi}(idx_higher_frac))];
                                
                            plot(dists{di}(:), measures_cmp{mi}(:), '.'); 
                            plot(x, y, 'ro-', 'linewidth', 3);
                            drawVerticalLine([min(dists{di}), dist_pctles, max(dists{di})], 'linewidth', 2, 'color', 'g')
                            xlabel(distNames{di}, 'fontsize', xlabel_fsize);
                            ylabel(measures_cmp_names{mi}, 'interp', 'none', 'fontsize', ylabel_fsize);
                            s1 = sprintf('cc = %.2f, pcc = %.2g', ccs{os_i}(di, mi), pvals_cc{os_i}(di, mi));
                            s2 = sprintf('d = %.2f, pcmp = %.2g', diff(y), pvals_cmp{os_i}(di, mi));
                            title({s1, s2})
                        end
                    end
                end
                   
                
            end
            
            3;
            
            
            
        end
        
        
        if doPlots && doMainFigurePlots && any(doFigs == 9) && doSpf
            %%
            % 9a. Scatterplot of preferred spatial frequency 
            f_opt1 = [allSpfUnitStats(Wcc_pairs_spf(:,1)).f_opt];
            f_opt2 = [allSpfUnitStats(Wcc_pairs_spf(:,2)).f_opt];

%             [f_opt1, f_opt2] = xyyx(f_opt1, f_opt2);
%             lims_plot = [.08, 5];
%             idx_use = ibetween(f_opt1, lims_plot) & ibetween(f_opt2, lims_plot);
%             [f_opt1, f_opt2] = deal(f_opt1(idx_use), f_opt2(idx_use));

            figure(spfPref_figId); clf_ifFirstRow(); % plots of orientation width 
            if printPlotsToFilesForPaper
                set(gcf, 'position', [fig_LB + (fig_offset * spfPref_figId), [fig_width_row, fig_height_row*2]], 'color', 'w')
            end
            h_ax9a = subplotGap(subRows,2, row_idx,  1, subSpcM, subSpcN); 
            h_plot9 = loglog(f_opt1, f_opt2, mk, 'markersize', sz, 'color', scat_col); 
            title(add_fg('Preferred spatial frequency'), 'fontsize', title_fsize);
            xlabel('Cell 1', 'fontsize', xlabel_fsize); 
            ylabel('Cell 2', 'fontsize', ylabel_fsize);
            
            f_range =  [0.08, 3];
            all_f_opt = [f_opt1];
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
                   
            
            h_ax9b = subplotGap(subRows,2, row_idx,  2, subSpcM, subSpcN); 
            maxDiffTh = 0.99;                    
            spf_pref_idx = find(strcmp(measures_spf, 'D_spf_pref'), 1);
            binEdges = S_spf{spf_pref_idx}.binEdges;
            binEdges = binEdges(1):0.25:binEdges(end);
            binCents = binEdge2cent(binEdges);
    %         nBins = length(binEdges)-1;

    %         Wcc_cumFracBins = cumFracBinsFromBinIds(S_spf{spf_pref_idx}.bins(Wcc_ss_pairIdxs), nBins);
    %         Bcc_cumFracBins = cumFracBinsFromBinIds(S_spf{spf_pref_idx}.bins(Bcc_ss_pairIdxs), nBins);
            
            [Wcc_cumFracBins, Wcc_binVals] = cumFracBinsFromVals(dSpfPref_wcc, binEdges);
            [Bcc_cumFracBins, Bcc_binVals] = cumFracBinsFromVals(dSpfPref_bcc, binEdges);        
            binIdCloseTo1 = max(find(Wcc_cumFracBins >= maxDiffTh, 1), find(Bcc_cumFracBins >= maxDiffTh, 1) );

            plot(binEdges, [0 Wcc_cumFracBins], 'o-',             'color', w_col, 'linewidth', lineWidth_main, 'markersize', cumMkSize); hold on;
            plot(binEdges, [0 Bcc_cumFracBins], ['s' dashedLine], 'color', b_col, 'linewidth', lineWidth_main, 'markersize', cumMkSize);
            xlims = [0, binEdges(binIdCloseTo1+1)];
            xlim(xlims);
    %         xlim(binEdges([1, end]));
%             title(add_fg('Differences in pref. spatial frequency'), 'fontsize', title_fsize);
            title(add_fg('Preferred spatial frequency'), 'fontsize', title_fsize);
            
            xlabel('Diff. in pref. spatial freq., octaves', 'fontsize', xlabel_fsize);
            if addLegends
                h_leg9 = legend(within_between_legend_str, 'location', 'SE', 'fontsize', legend_fsize);
            end
            ylabel('Cumulative fraction of pairs', 'fontsize', ylabel_fsize);
            3;
            
            if addProbSubplots
                drawnow;                
                pos_inset = insetAxesPosition(get(h_ax9b, 'position'), [.4, .02, .6, .6]);
                h_ax9_inset = axes('outerposition', pos_inset, 'ytick', [], 'fontsize', 8, 'nextplot', 'add');
               
                plot(binCents, Wcc_binVals/sum(Wcc_binVals), '-', 'color', w_col, 'markersize', distMkSize, 'linewidth', lineWidth_inset);
                plot(binCents, Bcc_binVals/sum(Bcc_binVals), ':', 'color', b_col, 'markersize', distMkSize, 'linewidth', lineWidth_inset);            
                
                axis(h_ax9_inset, 'tight');     
                xlim(xlims);

                ylims = get(h_ax9_inset, 'ylim');
                set(h_ax9_inset, 'ylim', [0, ylims(2)*1.05]);
                
                if addLegends
                    ax_inset_pos = get(h_ax9_inset, 'position');
                    leg_pos = get(h_leg9, 'position');
                    LB = [ax_inset_pos(1), ax_inset_pos(2)+ax_inset_pos(4)+.02];
                    set(h_leg9, 'position', [LB, leg_pos(3:4)]);
                end
                
            
            end
            if addSubplotLetters
                addSubplotLetter(subRows,2, row_idx, 1, subSpcM, subSpcN, char('A'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont);            
                addSubplotLetter(subRows,2, row_idx, 2, subSpcM, subSpcN, char('B'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont);
            end
            
            if dispPlotStats
                
                pctPairsLessThan1 = nnz(dSpfPref_wcc < 1)/length(dSpfPref_wcc) * 100; % for comparison with other studies in supplement.                
            end
            

            
            
            
        end
        
        
        figureFolder = [CatV1Path 'Figures' filesep 'DegreePaper' filesep];
        fig_format = 'pdf';
        
        if printPlotsToFilesForPaper && doMainFigurePlots &&  saveFilesNow && (row_idx == 2) && (subtractSpont == 1)
            %%
            
            input('Press Return to Export to PDF');
            %%
            
            
            %%
            figure(oriWidthScatter_figId);
            fig_filename_oriWidthScatter = sprintf('%sFigure%d_oriWidthScatter.pdf', figureFolder, oriWidthScatter_figId);
            export_fig(oriWidthScatter_figId, fig_format, fig_filename_oriWidthScatter);
%%
            figure(oriWidthCumProb_figId);
            fig_filename_oriCumWidth = sprintf('%sFigure%d_oriWidthCumProb.pdf', figureFolder, oriWidthCumProb_figId);
            export_fig(oriWidthCumProb_figId, fig_filename_oriCumWidth);

            figure(oriPrefScatter_figId);
            fig_filename_oriPrefScatter = sprintf('%sFigure%d_oriPrefScatter.pdf', figureFolder, oriPrefScatter_figId);
            export_fig(oriPrefScatter_figId, fig_format, fig_filename_oriPrefScatter);
        
            figure(oriPrefProbDist_figId);
            fig_filename_oriPrefProb = sprintf('%sFigure%d_oriPrefProbDist.pdf', figureFolder, oriPrefProbDist_figId);
            export_fig(oriPrefProbDist_figId, fig_format, fig_filename_oriPrefProb);

            figure(DSI_figId);
            fig_filename_DSI = sprintf('%sFigure%d_DSI.pdf', figureFolder, DSI_figId);
            export_fig(DSI_figId, fig_format, fig_filename_DSI);

            figure(dirPrefScatter_figId);
            fig_filename_dirPrefScatter = sprintf('%sFigure%d_dirPrefScatter.pdf', figureFolder, dirPrefScatter_figId);
            export_fig(dirPrefScatter_figId, fig_format, fig_filename_dirPrefScatter);

            figure(dirPrefProbDist_figId);
            fig_filename_dirPrefProb = sprintf('%sFigure%d_dirPrefProbDist.pdf', figureFolder, dirPrefProbDist_figId);
            export_fig(dirPrefProbDist_figId, fig_format, fig_filename_dirPrefProb);

            figure(spfWidth_figId);
            fig_filename_spfWidth = sprintf('%sFigure%d_spfWidth.pdf', figureFolder, spfWidth_figId);
            export_fig(spfWidth_figId, fig_format, fig_filename_spfWidth);
            
            figure(spfPref_figId);
            fig_filename_spfPref = sprintf('%sFigure%d_spfPref.pdf', figureFolder, spfPref_figId);
            export_fig(spfPref_figId, fig_format, fig_filename_spfPref);

        end
        
        
        
        
        
        if doBetweenAnimalTests
            %%
            fig_i = 400 + curGratingType*100;
            allCatIds_used = [2     4     5     9    10    12    16    18    19    33    34    36    43    44    48    52    53    54];
            nCatsTot = length(allCatIds_used);
            
            allCellLocData_ori = [allOriCells.locData];
            cellCatIds_ori_orig = [allCellLocData_ori.CatId];
            cellCatIds_ori = binarySearch(allCatIds_used, cellCatIds_ori_orig, 1, 0); assert(all(cellCatIds_ori > 0));
            F1oDCs_ori = [allOriCells.F1oDC_maxR_avP_sm];
            
            
            allCellLocData_spf = [allSpfCells.locData];
            cellCatIds_spf_orig = [allCellLocData_spf.CatId];
            cellCatIds_spf = binarySearch(allCatIds_used, cellCatIds_spf_orig, 1, 0); assert(all(cellCatIds_spf > 0));
            F1oDCs_spf = [allSpfCells.F1oDC_maxR_avP_sm];
            
            cell_w_ori_global = [allOriCellStats.w_ori_global];
            cell_w_ori_local = [allOriCellStats.w_ori_local];
            cell_ori_pref = [allOriCellStats.ori_pref_deg];
            cell_DSI = [allOriCellStats.DSI_global];

            cell_w_spf = [allSpfCellStats.w_spf];
            cell_spf_pref = [allSpfCellStats.f_opt];

            set(0,'DefaultFigureWindowStyle','docked') 
            for stim_type = 1:2
                %%
                if stim_type == 1
                    all_ms = {cell_w_ori_global, cell_w_ori_local, cell_ori_pref, cell_DSI, F1oDCs_ori};
                    all_ms_names = {'Global Ori Width', 'Local Ori Width', 'Preferred Ori', 'DSI', 'F1/DC (ori)'};
%                     all_ms = {F1oDCs_ori};
%                     all_ms_names = {'F1/DC (ori)'};
                    cellCatIds = cellCatIds_ori;
                    cellGroupIds = [allOriCells.Gid];
                    cellPenIds = [allCellLocData_ori.PenId];
                    F1oDCs = F1oDCs_ori;
                else
                    all_ms = {cell_w_spf, cell_spf_pref, F1oDCs_spf};
                    all_ms_names = {'Spatial freq. tuning width', 'Preferred Spatial freq.', 'F1/DC (spf)'};
%                     all_ms = {F1oDCs_spf};
%                     all_ms_names = {'F1/DC (spf)'};
                    cellCatIds = cellCatIds_spf;
                    cellGroupIds = [allSpfCells.Gid];
                    cellPenIds = [allCellLocData_spf.PenId];
                    F1oDCs = F1oDCs_spf;
                end

                
                
                for ms_i = 1:length(all_ms)
                    
                    if strcmp(gratingType, 'flashed') && strcmp(all_ms_names{ms_i}, 'DSI')
                        fig_i = fig_i+1;
                        continue;
                    end

                    %%
%                     f1odc_ranges = {[0, 2], [1, 2], [0, 1]};
%                     f1odc_ranges = {[0, 2]};
                    f1odc_ranges = {[1, 2], [0, 1]};
                    nRanges = length(f1odc_ranges);

                    for ri = 1:nRanges
          %             f1odc_range = [0, 1];
         %             f1odc_range = [0, 2];
                        f1odc_range = f1odc_ranges{ri};
                        idx_ok = find(ibetween(F1oDCs, f1odc_range));
                        [uCatIds, catIdxs] = uniqueList(cellCatIds);
                        %%
                        for i = 1:length(catIdxs)
                            catIdxs{i} = intersect(catIdxs{i}(:), idx_ok(:));
                        end
                        idx_allCats = cat(1, catIdxs{:});
                        %%
                        min_nCells = 5;
                        idx_use = cellfun(@length, catIdxs) >= min_nCells;
                        uCatIds = uCatIds(idx_use);
                        catIdxs = catIdxs(idx_use);
                        nCats = length(catIdxs);
                        nCellsEachCat = cellfun(@length, catIdxs);
                        for i = 1:nCats
                            nSitesEachCat(i) = length(unique(cellGroupIds(catIdxs{i})));
                            nPensEachCat(i) = length(unique(cellPenIds(catIdxs{i})));
                        end

                        catIdxs_withFull = [catIdxs, idx_allCats];

                    
            %             nCats = length(uCatIds);
                        X = all_ms{ms_i};
                        pvals = nan(nCats, nCats);
                        pvals_withFull = nan(nCats, 1);
                        for i = 1:nCats
                            idxs_i = catIdxs{i};
                            pvals_withFull(i) = ranksum(X(idxs_i), X(idx_allCats));
                            for j = 1:i-1
                                idxs_j = catIdxs{j};
                                pvals(i,j) = ranksum(X(idxs_i), X(idxs_j) );

                            end
                        end
                        pvals = [pvals, pvals_withFull];

                        figure(fig_i); clf;

                        all_markers = marker( floor([1:nCatsTot] / 6)+2, 'os^');

    %                     all_colors = color_s(
                        f1odc_str = switchh(f1odc_range, {[0, 2], [0, 1], [1, 2]}, {'[ALL CELLS]', '[COMPLEX CELLS]', '[SIMPLE CELLS]', sprintf(' [F1/DC in %.1f - %.1f]', f1odc_range)});
                        title_str = sprintf('%s (%s) %s', all_ms_names{ms_i}, gratingType, f1odc_str);

                        binArg = 10;
                        if strncmp(all_ms_names{ms_i}, 'F1', 2)
                            binArg = [0:.2:2];
                        end
                        h_ax(1) = subplotGap(1,3,1,1);
                        vals = cellfun(@(idx) X(idx), catIdxs_withFull, 'un', 0);
                        h1 = hist2(vals, binArg, 'line', 'norm');
                        title({'', 'Distributions'});

                        h_ax(2) = subplotGap(1,3,1,2);
                        h2 = hist2(vals, binArg, 'line', 'norm', 'cum');
                        col_order = mat2cell(get(gca, 'colorOrder'), ones(1,7), 3);
                        col_order = col_order(1:end-1);
    %                     linestyles = ['osx'];
                        for i = 1:length(h1)-1
                            set([h1(i) h2(i)], 'marker', all_markers( uCatIds(i)), 'color', color_s( uCatIds(i), col_order ), 'markersize', 10, 'linewidth', 2)
                        end
                         set([h1(end) h2(end)], 'marker', '*', 'color', 'k', 'markersize', 10, 'linewidth', 2)
                        legend_strs = [arrayfun(@(i) sprintf('Cat %d (%dc, %ds, %dp)', uCatIds(i), nCellsEachCat(i), nSitesEachCat(i), nPensEachCat(i)), 1:nCats, 'un', 0), ...
                            sprintf('All cats (%dc, %ds, %dp)', length(idx_allCats), sum(nSitesEachCat), sum(nPensEachCat)) ];
                        h_leg = legend(legend_strs, 'fontsize', 7, 'location', 'SE');
                        
                        %x = get(h1(1), 'xdata');
                        xlims = get(h_ax(1), 'xlim'); %lims(x, .1);
                        xlims_ext = lims(xlims, .4);
%                         if length(binArg) > 1
%                             xlims_ext = lims(binArg, .1);
%                         end
                        xlim([xlims(1), xlims_ext(2)]);


                        title({title_str, 'Cumulative Distributions'}, 'interp', 'none');

                        h_ax(3) = subplotGap(1,3,1,3);
                        imagesc(-log10(pvals));
                        axis square;
                        colorbar;
                        caxis([0, 4]);


                        xlabel('Cat j'); h_ylab = ylabel('Cat i'); title([title_str ':  Comparisons'], 'interp', 'none');
                        set(gca, 'fontSize', 7)
                        set(gca, 'xtick', 1:nCats+1, 'xticklabel', [num2cell(uCatIds), 'All'], 'ytick', 1:nCats+1, 'yticklabel', uCatIds)
                        title({'', 'Pairwise comparisons (-log10(p))'}, 'fontsize', 10);
%                         p_ylab_ds = get(h_ylab, 'position'); p_ylab = ds2nfu(h_ax(2), p_ylab_ds(1), p_ylab_ds(2));
                        p_ylab = get(h_ax(3), 'outerposition');
                        px_ylab = p_ylab(1);
                        drawnow;
                        
                        p_leg = get(h_leg, 'position');
                        new_p_leg = p_leg; new_p_leg(1) = px_ylab - p_leg(3)+.008;
                                           new_p_leg(2) = .01;         
                        set(h_leg, 'position', new_p_leg);
                        
                        3;
                        fig_i = fig_i + 1;
                    end % for f1/dc range
                    3;
                end % for measure
            
            end % for ori/spf
            3;
       
            %%
            
%             allOriCellStats = allOriUnitStats(oriCells_use);
%     %         allOriCells = allOriUnits(oriCells_use);
%             allOriCellGids = [allOriUnits(oriCells_use).Gid];
%             allOriCellIDs = [allOriSpkFeatures(oriCells_use).IsolationDistance];
%             allOriErrors = [allOriCellStats.error_jack];
%  
            
            
            
        end
        % 
        %%        
        % supplementary figures;
        doOutlierStatsFigure = doPlots && doOutlierStatsFigure && doOri;
        clear h5
        if doOutlierStatsFigure || doPairStats
            %%
            printData = showWorking;
            % comparison with multi-units
            spkFeatures = allOriSpkFeatures(oriCells_use);
            spkWH = [spkFeatures.spkWidthHeight];
            
            gratingStr = [titleCase(gratingType) ' gratings'];
            Dori_pref_MU = [allOriCellStats.(dOriMU_field)];  
            
            w_ori_global = [allOriCellStats.w_ori_global];
            w_ori_local = [allOriCellStats.w_ori_local];
            F1oDCs = ori_unitF1oDC(oriCells_use);
            ID = [spkFeatures.IsolationDistance];
            spkAmps = [spkFeatures.spikeAmp] * spikeAmpFactor;
            FWHM = [spkWH.FWHM];
            PtP_height = [spkWH.PtP_height];
            PtP_width = [spkWH.PtP_width];
            error_jack = [allOriCellStats.error_jack];
            w_ori_global_jack = [error_jack.w_ori_global];
            

            Y_vals       = {spkAmps, w_ori_global, w_ori_local, FWHM, PtP_height, PtP_width, F1oDCs, ID, w_ori_global_jack};            
            Y_labels     = {'spike amplitude, \muV', 'w_{ORI}^{Global}', 'w_{ORI}^{Local}', 'Full-Width at Half-Max (FWHM)', 'Peak-to-Peak height', 'Peak-to-Peak width', 'F1/DC', 'Isolation Distance', 'w_ori_global_jack'};
            Y_name_short = {'Spike amplitude', 'Global Ori width', 'Local Ori Width', 'FWHM', 'PtP height', 'PtP width','F1/DC', 'Isolation Distance', 'w_ori_global_jack'};
            Y_name_fld   = {'SpikeAmps', 'Ori_W_Global', 'Ori_W_Local', 'FWHM', 'PtP_height', 'PtP_width','F1_DC', 'IsolationDistance', 'w_ori_global_jack'};
            
%             idx_use = [1 2 3 4 5 6 7 8];
            idx_plot = [2 3 6];
%             Y_vals_plot = Y_vals(idx_use); Y_labels = Y_labels(idx_use); Y_name_short = Y_name_short(idx_use);
                        
            nVariables = length(Y_vals);
            nPlots = length(idx_plot);
            separateWideOri = 0;
                            
            idx_outlier = Dori_pref_MU > 45;     
            idx_norm = ~idx_outlier;
%             sym1 = 'b+'; sym2 = 'ks'; sym3 = 'r^';
            sym1 = 'bo'; sym2 = 'k^';
            %%
            ranksum_nonnans = @(x,y) ranksum( nonnans(x), nonnans(y) );
            
            if printData
                fprintf('\n\n*** Outlier Statistics for %s Gratings (%s)\n', titleCase(gratingType), spont_s);
            end
            
            if doOutlierStatsFigure
                fig_offset = iff(strcmp(gratingType, 'flashed'), -3, 0);
                figure(106+ fig_offset); clf;
            end
            
            outlierSigStats = struct('ColumnNames', {{'Outl', 'Norm'}});
            for var_i = 1:nVariables
%                 subplotGap(2,nPlots/2,sub_i); 
                Y_i = Y_vals{var_i};
                sub_i = find(idx_plot == var_i);
                                
                pvalU = ranksum_nonnans(Y_i(idx_norm), Y_i( idx_outlier) );                
%                 [h, pvalT_i12] = ttest2(Y_i(idx_norm), Y_i(idx_outlier) );
                
                title_str_now = sprintf('p_U = %.2g', pvalU);
                
                gratingType_str = sprintf('(%s gratings)', titleCase(gratingType) );
                
                
                med_amps_norm = nanmedian( Y_i( idx_norm) );
                med_amps_out = nanmedian( Y_i( idx_outlier) );
                                     
%                 outlierSigStats.(Y_name_fld{var_i}) = [struct('outlier_median', med_amps_out, 'norm_median', med_amps_norm, 'pval_U', pvalU);
                outlierSigStats.(Y_name_fld{var_i})        = [med_amps_out, med_amps_norm];
                outlierSigStats.([Y_name_fld{var_i} '_p']) = [pvalU];                    

                if printData
                    fprintf('%s (Median) : Outliers: %.2f.  Non-outliers: %.2f...  pU(1): pU = %.2g \n', ...
                            Y_name_short{var_i}, med_amps_out,  med_amps_norm, pvalU);                                    
                end
                
                if doOutlierStatsFigure && ~isempty(sub_i)
                    
                    subplotGap(1,nPlots,sub_i);

                        
    %                 h5{sub_i} = plot(Dori_pref_MU(idx_norm), Y_i(idx_norm), sym1, ...
    %                                  Dori_pref_MU(idx_mark1), Y_i(idx_mark1), sym2, ...
    %                                  Dori_pref_MU(idx_mark2), Y_i(idx_mark2), sym3);  %#ok<AGROW>
                        h5{sub_i} = plot(Dori_pref_MU(idx_norm), Y_i(idx_norm), sym1, ...
                                         Dori_pref_MU(idx_outlier), Y_i(idx_outlier), sym2 ); 

                    if sub_i == 1 % floor((nPlots+1)/2)
                        xlabel('Difference from multiunits, degrees');
                        title_str0 = gratingStr;
                    else
                        xlabel(' ');
                        title_str0 = ' ';
                    end
                    ylabel(Y_labels{var_i});

                    xlim([0 90]); set(gca, 'xtick', [0:15:90]);
                    title({Y_name_short{var_i}, gratingType_str, title_str_now});

                    if strcmp(Y_name_short{var_i}, 'PtP width')
                        ylim([0.1, 0.7]); 
                    end
                    
                    if sub_i == nPlots
                        %  legend({'Other cells', 'Outliers #1', 'Outliers #2'}, 'location', 'SE')
                        legend({'Typical cells', 'Outliers'}, 'location', 'NE', 'fontsize', legend_fsize)
                    end
                    if strcmp(Y_labels{sub_i}, 'w_{ORI}^{Global}')
                        % line([45 45], [0 50], 'linestyle', ':', 'color', 'k');
                        % line([45 90], [35 35], 'linestyle', ':', 'color', 'k');
                    end
                    if addSubplotLetters
                        addSubplotLetter(1, nPlots, 1, sub_i, subSpcM, subSpcN, char('A'+(grating_offset/2)*nPlots+sub_i-1), [], 'fontname', subplotLetFont);
                    end
                end
                
               
                
            end
            miscPairStats.outlierSigStats = outlierSigStats;
            3;
            if doOutlierStatsFigure
                set([h5{:}], 'markersize',5);
            end
            3;
        end    
%%
        doFig107 = doPlots && doDSIvsOriWidthFigure;
        if  (doFig107 || doPairStats) && (strcmp(gratingType, 'drifting') && doOri) 
            
            w_ori_global = [allOriCellStats.w_ori_global]; w_ori_global = w_ori_global(:);
            w_ori_local  = [allOriCellStats.w_ori_local];  w_ori_local = w_ori_local(:);
            dsi          = [allOriCellStats.DSI_global];   dsi = dsi(:);

            [cc_dsi_wglob, p_dsi_wglob] = corr(dsi(:), w_ori_global(:), 'rows', 'complete');
            [cc_dsi_wloc, p_dsi_wloc] = corr(dsi(:), w_ori_local(:), 'rows', 'complete');                        

            corr_oriW_DSI = struct('dsi_w_ori_glob_cc', cc_dsi_wglob, ...
                                   'dsi_w_ori_glob_p', p_dsi_wglob, ...
                                   'dsi_w_ori_loc_cc', cc_dsi_wloc, ...
                                   'dsi_w_ori_loc_p' , p_dsi_wloc);

            miscPairStats.corr_oriW_DSI = corr_oriW_DSI;
            if doFig107
                %%
                % FIGURE S5 --> scatter plots of DSI vs. Global/local ORI width

                Nrows = 1;
                if subtractSpont
                    spont_str = '(Spontaneous Subtracted)';
                else
                    spont_str = '(Spontaneous Included)';
                end

                figure(107); clf; % plots of orientation width 
                h_ax101(1) = subplotGap(Nrows,2, 1);      plot(dsi, w_ori_global, mk, 'markersize', sz, 'color', 'k'); 
                t_str1a = 'Global Ori Width  vs  DSI';
                t_str1c = sprintf('r = %.2f (p = %s)', cc_dsi_wglob, removeExtraZerosFromExpStr(sprintf('%.1g', p_dsi_wglob)));
                title({t_str1a, spont_str, t_str1c});  xlabel('DSI'); ylabel('w_{ORI}^{Global}');
                h_ax101(2) = subplotGap(Nrows,2, 2);      plot(dsi, w_ori_local,  mk, 'markersize', sz, 'color', 'k'); 
                t_str2a = 'Local Ori Width  vs  DSI';
                t_str2c = sprintf('r = %.2f (p = %s)', cc_dsi_wloc, removeExtraZerosFromExpStr(sprintf('%.1g', p_dsi_wloc)));
                title({t_str2a, spont_str, t_str2c});   xlabel('DSI'); ylabel('w_{ORI}^{Local}');
                if addSubplotLetters
                    addSubplotLetter(1, 2, 1, 1, subSpcM, subSpcN, char('A'+(1-subtractSpont)*2), subplotLetter_offset, 'fontname', subplotLetFont);
                    addSubplotLetter(1, 2, 1, 2, subSpcM, subSpcN, char('B'+(1-subtractSpont)*2), subplotLetter_offset, 'fontname', subplotLetFont);
                end

    %             if subtractSpont
    %                 h_ax101(3) = subplotGap(Nrows,2, 3);  plot(dsi_ss, w_ori_global_ss, mk, 'markersize', sz);
    %                 title(sprintf('Spontaneous Subtracted. r = %.2f', cc_dsi_wglob_ss)); xlabel('DSI'); ylabel('w_{ORI}^{Global}');            
    %                 h_ax101(4) = subplotGap(Nrows,2, 4);  plot(dsi_ss, w_ori_local_ss,  mk, 'markersize', sz); 
    %                 title(sprintf('Spontaneous Subtracted. r = %.2f', cc_dsi_wloc_ss));  xlabel('DSI'); ylabel('w_{ORI}^{Local}');
    %             end

                fprintf('\n\n*** DSI vs Ori Width (%s)\n', spont_s);
                fprintf('CC between DSI and global ori width: %.3f (p = %.2g)\n', cc_dsi_wglob, p_dsi_wglob);
                fprintf('CC between DSI and local ori width: %.3f (p = %.2g)\n', cc_dsi_wloc, p_dsi_wloc);

    %             idx = ~ori_is_outlier;
    %             [cc_dsi_wglob_norm, p_dsi_wglob_norm] = corr(dsi(idx), w_ori_global(idx), 'rows', 'complete');
    %             [cc_dsi_wloc_norm, p_dsi_wloc_norm] = corr(dsi(idx), w_ori_local(idx), 'rows', 'complete');                        
    %             
    %             fprintf('CC between DSI and global ori width (Outliers removed): %.3f (p = %.2g)\n', cc_dsi_wglob_norm, p_dsi_wglob_norm);
    %             fprintf('CC between DSI and local ori width (Outliers removed): %.3f (p = %.2g)\n', cc_dsi_wloc_norm, p_dsi_wloc_norm);

                % with outliers removed            

    %             if subtractSpont
    %                 fprintf('CC between DSI and global ori width (spont subtracted): %.3f (p = %.2g)\n', cc_dsi_wglob_ss, p_dsi_wglob_ss);
    %                 fprintf('CC between DSI and local ori width (spont subtracted): %.3f (p = %.2g)\n', cc_dsi_wloc_ss, p_dsi_wloc_ss);
    %                                 
    %                 [cc_dsi_wglob_ss_norm, p_dsi_wglob_ss_norm] = corr(dsi_ss(idx), w_ori_global_ss(idx), 'rows', 'complete');
    %                 [cc_dsi_wloc_ss_norm, p_dsi_wloc_ss_norm] = corr(dsi_ss(idx), w_ori_local_ss(idx), 'rows', 'complete');                        
    % 
    %                 fprintf('CC between DSI and global ori width (spont subtracted; Outliers removed): %.3f (p = %.2g)\n', cc_dsi_wglob_ss_norm, p_dsi_wglob_ss_norm);
    %                 fprintf('CC between DSI and local ori width (spont subtracted; Outliers removed): %.3f (p = %.2g)\n', cc_dsi_wloc_ss_norm, p_dsi_wloc_ss_norm);
    %                 
    %             end



                3;
            end
        end
        
        % Supp figure 8 -- dOriPref vs spikeAmplitude
        doFig108 = doPlots && doOriVsSpikeAmplitudeFigure; 
        
        if  (doFig108 && doPairStats) 
            %%
            oriUnitSpikeAmps = [allOriSpkFeatures.spikeAmp] * spikeAmpFactor;            
            
            oriPairSpkAmps = oriUnitSpikeAmps(Wcc_pairs_ori);
            dOriSpkAmps = abs(diff(oriPairSpkAmps, [], 2));
            minOriSpkAmps = max(oriPairSpkAmps, [], 2);
            oriNegAmpsDist = pairData_ori.negAmps_dist(Wcc_oo_pairIdxs);
            
            
            allSpfSpkFeatures = [allSpfUnits.spkFeatures];
            spfUnitSpikeAmps = [allSpfSpkFeatures.spikeAmp] * spikeAmpFactor;
            
            spfPairSpkAmps = spfUnitSpikeAmps(Wcc_pairs_spf);
            dSpfSpkAmps = abs(diff(spfPairSpkAmps, [], 2));
            minSpfSpkAmps = max(spfPairSpkAmps, [], 2);
            spfNegAmpsDist = pairData_spf.negAmps_dist(Wcc_ss_pairIdxs);
            
            %%
            xs = {dOriPref_wcc, dOriW_global_wcc, dOriW_local_wcc, dSpfPref_wcc, dSpf_w_wcc};
            isOri = [1, 1, 1, 0, 0];
%             xlab = {'Ori Pref', 'OriWidth(glob)', 'OriWidth(local)', 'Spf Pref', 'SpfWidth'};
            xlab = {'Preferred Orientation', 'Global Orientation Width', 'Local Orientation Width', 'Preferred Spat Freq', 'Spat Freq Width'};
            if strcmp(gratingType, 'drifting')
                xs = [xs, dDSI_wcc];
                isOri = [isOri, 1];
                xlab = [xlab, 'DSI'];
            end
            %%            
            ori_separateNormAndOutliers = 0;
            separateLowerPctile = 0; pctile_value = 10;
            useAbsAmplitudes = 1;
            testAboveBelowPctile = 0;
            markSignificantTextAsRed = 1;  sig_th = 0.05;
%             xlabels = cellfun(@(s) sprintf('diff (%s [%s])', s, gratingType(1)), xlab, 'un', 0 );            
            xlabels = cellfun(@(s) {sprintf('Difference in %s', s), sprintf('(%s gratings)', titleCase(gratingType))}, xlab, 'un', 0 );            
            
            if useAbsAmplitudes
                minOriSpkAmps_use = abs(minOriSpkAmps);
                minSpfSpkAmps_use = abs(minSpfSpkAmps);
            else
                minOriSpkAmps_use = minOriSpkAmps;
                minSpfSpkAmps_use = minSpfSpkAmps;
            end
            ys_ori = {minOriSpkAmps_use, dOriSpkAmps, oriNegAmpsDist};
            ys_spf = {minSpfSpkAmps_use, dSpfSpkAmps, spfNegAmpsDist};
            
            if useAbsAmplitudes
                pctile_value = 100-pctile_value;
            end
            ori_lower_pct = prctile(minOriSpkAmps_use, pctile_value);
            spf_lower_pct = prctile(minSpfSpkAmps_use, pctile_value);
            if useAbsAmplitudes
                idx_bothLargeAmp_ori = minOriSpkAmps_use > ori_lower_pct;
                idx_bothLargeAmp_spf = minSpfSpkAmps_use > spf_lower_pct;
            else
                idx_bothLargeAmp_ori = minOriSpkAmps_use < ori_lower_pct;
                idx_bothLargeAmp_spf = minSpfSpkAmps_use < spf_lower_pct;                
            end
            
            ylabels = {'max (=min(abs)) SpkAmps', 'Difference in spike amplitudes', '4D cluster centers distance'};
            if useAbsAmplitudes
                ylabels{1} = 'min Spike Amplitude of pair';
            end
            
            y_use_idx = [1];
%             y_use_idx = [1];
            ys_ori = ys_ori(y_use_idx);
            ylabels = ylabels(y_use_idx);
            nY = length(y_use_idx);
            
            
              fig_offset = iff(strcmp(gratingType, 'drifting'), 0, 100);
            for xi = 1:length(xs)
                X = xs{xi};
                figure(108+xi-1 + fig_offset);
                clf;                
                for yi = 1:nY
                    if isOri(xi)
                        Y = ys_ori{yi};
                        idx_bothLarge = idx_bothLargeAmp_ori;
                        pct_low = ori_lower_pct;
                    else
                        Y = ys_spf{yi};
                        idx_bothLarge = idx_bothLargeAmp_spf;
                        pct_low = spf_lower_pct;
                    end   
%                     X(idx_bothLarge) = X(idx_bothLarge)*5;
                    subplot(1,nY,yi); hold on; box on;   
                    if isOri(xi) && ori_separateNormAndOutliers
                        plot(X(idx_noOutliers), Y(idx_noOutliers), 'bo', 'markersize', 3); 
                        plot(X(idx_withOutliers), Y(idx_withOutliers), 'ro', 'markersize', 3);                           
                    else
                        plot(X, Y, 'bo', 'markersize', 3);                            
                    end                    

                    xlabel(xlabels{xi}); 
                    ylabel(ylabels{yi});
                    
                    [cc, cc_p] = corr(X, Y, 'type', 'pearson');
                    [rho, rho_p] = corr(X, Y, 'type', 'spearman');
                    
                    if markSignificantTextAsRed
                        cc_txt_color_str = sprintf('\\color{%s}', iff(cc_p < sig_th, 'red', 'black'));
                        rho_txt_color_str = sprintf('\\color{%s}', iff(rho_p < sig_th, 'red', 'black'));
                    else
                        cc_txt_color_str = '';
                        rho_txt_color_str = '';
                    end
                    
                    str_cc_all = sprintf('%scc = %.2f, p = %.2g', cc_txt_color_str, cc, cc_p);
                    str_rho_all = sprintf('%s\\rho = %.2f, p = %.2g', rho_txt_color_str, rho, rho_p);

                    [cc_nbl, p_nbl] = corr(X(~idx_bothLarge), Y(~idx_bothLarge), 'type', 'pearson');
                    [rho_nbl, rho_p_nbl] = corr(X(~idx_bothLarge), Y(~idx_bothLarge), 'type', 'spearman');
                    
                    str_cc_nbl = sprintf('\\color[rgb]{0 .7 0}cc = %.2f, p = %.2g', cc_nbl, p_nbl);
                    str_rho_nbl = sprintf('\\rho = %.2f, p = %.2g', rho_nbl, rho_p_nbl);                    
                    title_str_C = {str_cc_all, str_rho_all};
                    if separateLowerPctile
                        title_str_C = [title_str_C, str_cc_nbl, str_rho_nbl]; %#ok<*AGROW>
                    end
                    
                    if isOri(xi) && ori_separateNormAndOutliers                        
                        
                        [cc_typ, p_typ] = corr(X(idx_noOutliers), Y(idx_noOutliers), 'type', 'pearson');
                        [rho_typ, rho_p_typ] = corr(X(idx_noOutliers), Y(idx_noOutliers), 'type', 'spearman');
                        str_cc_typ = sprintf('\\color{blue}cc = %.2f, p = %.2g', cc_typ, p_typ);
                        str_rho_typ = sprintf('\\rho = %.2f, p = %.2g', rho_typ, rho_p_typ);
                        
                        [cc_out, p_out] = corr(X(idx_withOutliers), Y(idx_withOutliers), 'type', 'pearson');
                        [rho_out, rho_p_out] = corr(X(idx_withOutliers), Y(idx_withOutliers), 'type', 'spearman');
                        str_cc_out = sprintf('\\color{red}cc = %.2f, p = %.2g', cc_out, p_out);
                        str_rho_out = sprintf('\\rho = %.2f, p = %.2g', rho_out, rho_p_out);
                        
%                         title_str_C = [title_str_C, str_cc_typ, str_rho_typ, str_cc_out, str_rho_out];                        
                    end                    
                                  
                    
                    if separateLowerPctile && yi == 1
                        drawHorizontalLine(pct_low, 'color', [0 .7 0], 'linestyle', '-', 'linewidth', 2);                       
                    end
                    if testAboveBelowPctile
                         med_above = median(X(~idx_bothLarge));
                         med_below = median(X(idx_bothLarge));
                         pval_U = ranksum(X(idx_bothLarge), X(~idx_bothLarge));
                         ylims = ylim;                         
                         
                         line(med_above * [1,1], [pct_low ylims(2)], 'color', 'k', 'linewidth', 2);
                         line(med_below * [1,1], [pct_low ylims(1)], 'color', 'k', 'linewidth', 2);
                         u_test_str = sprintf('\\color{black}Above:%.2f. Below: %.2f. p = %.2g', med_above, med_below, pval_U);
                         title_str_C = [title_str_C, u_test_str];
                    end
                                        
                    if xi == 1
                        set(gca, 'xtick', [0:30:90]); xlim([0 90]);
                    end
                    
                    
%                     title(title_str_C, 'fontsize', 11);      
%                     title(title_str_C, 'fontsize', 11);      
                    xlims = xlim; ylims = ylim;
                    h_t = text( xlims(1) + diff(xlims)*.98, ylims(1) + diff(ylims)*.98, title_str_C, 'fontsize', 9);
                    set(h_t, 'horizontalAlignment', 'right', 'verticalAlignment', 'top');
                    
                    
                    
                end
            end
                    
            3;
                  
%             
%             figure(108);                
%             subplot(1,3,1);
% 
%             subplot(1,3,2);
%             plot(dOriPref_wcc, dOriSpkAmps, 'o', 'markersize', 3);
%             [cc_dOri_dAmps, p_dOri_dAmps] = corr(dOriPref_wcc, dOriSpkAmps, 'type', 'spearman');
%             xlabel(sprintf('diff in Pref. Orientation (%s)', gratingType)); ylabel();
%             title(sprintf('\\rho = %.2f, p = %.2g', cc_dOri_dAmps, p_dOri_dAmps));
%             set(gca, 'xtick', [0:30:90]); xlim([0 90]);
%             
%             subplot(1,3,2);
%             plot(dOriPref_wcc, oriNegAmpsDist, 'o', 'markersize', 3);
%             [cc_dOri_dAmps, p_dOri_dAmps] = corr(dOriPref_wcc, oriNegAmpsDist, 'type', 'spearman');
%             xlabel(sprintf('diff in Pref. Orientation (%s)', gratingType)); ylabel();
%             title(sprintf('\\rho = %.2f, p = %.2g', cc_dOri_dAmps, p_dOri_dAmps));
%             set(gca, 'xtick', [0:30:90]); xlim([0 90]);
%             
%             figure(109);                       
%             subplot(1,2,1);
%             plot(dSpfPref_wcc, minSpfSpkAmps, 'o', 'markersize', 3);
%             [cc_dSpf_minAmps, p_dSpf_minAmps] = corr(dSpfPref_wcc, minSpfSpkAmps, 'type', 'spearman');
%             xlabel(sprintf('diff in Pref. Spat. Frequency (%s)', gratingType)); ylabel('max (=min(abs)) SpkAmps');
%             title(sprintf('\\rho = %.2f, p = %.2g', cc_dSpf_minAmps, p_dSpf_minAmps));
% 
%             subplot(1,2,2);
%             plot(dSpfPref_wcc, spfNegAmpsDist, 'o', 'markersize', 3);
%             [cc_dSpf_dAmps, p_dSpf_dAmps] = corr(dSpfPref_wcc, spfNegAmpsDist, 'type', 'spearman');
%             xlabel(sprintf('diff in Pref. Spat. Frequency (%s)', gratingType)); ylabel('diff (SpkAmps)');
%             title(sprintf('\\rho = %.2f, p = %.2g', cc_dSpf_dAmps, p_dSpf_dAmps));            
%             
%             
            3;
        end
        % verify that no outliers for spatial frequency
        if doPlots && doCheckSpatialFrequencyOutliersFig
            % 1. plot of diff in pref spatial frequency vs max spatial freq tuning width of the pair
            % 2. plot of diff in pref spatial frequency vs min spike amplitude of the pair
            %%
%             dSpfPref_wcc = S_spf{spf_pref_idx}.val(Wcc_ss_pairIdxs);

            w_spf1 = [allSpfUnitStats(Wcc_pairs_spf(:,1)).w_spf];
            w_spf2 = [allSpfUnitStats(Wcc_pairs_spf(:,2)).w_spf];
            max_w_spf = max(w_spf1, w_spf2)';
            
            allSpfSpkFeatures = [allSpfUnits.spkFeatures];
            spkAmp1 = [allSpfSpkFeatures(Wcc_pairs_spf(:,1)).spikeAmp] * spikeAmpFactor;
            spkAmp2 = [allSpfSpkFeatures(Wcc_pairs_spf(:,2)).spikeAmp] * spikeAmpFactor;
            min_spkAmp = min(spkAmp1, spkAmp2)';
    
            figure(201); clf;
            subplotGap(1,2,1);            
            plot(dSpfPref_wcc, max_w_spf, '.');
            [r1, p1] = corr(dSpfPref_wcc, max_w_spf, 'type', 'spearman', 'rows', 'complete');
            xlabel('Difference in preferred spatial frequency');
            ylabel('Maximum spatial frequency tuning width of the pair');
            title(sprintf('p_s = %.2g', p1));
            if addSubplotLetters
                addSubplotLetter(1, 2, 1, 1, subSpcM, subSpcN, char('A'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont);
            end
            
            subplotGap(1,2,2);
            plot(dSpfPref_wcc, min_spkAmp, '.');
            xlabel('Difference in preferred spatial frequency');
            ylabel('Minimum spike amplitude of the pair');
            [r2, p2] = corr(dSpfPref_wcc, min_spkAmp, 'type', 'spearman', 'rows', 'complete');
            title(sprintf('p_s = %.2g', p2));
            if addSubplotLetters
                addSubplotLetter(1, 2, 1, 2, char('B'+grating_offset), subplotLetter_offset, 'fontname', subplotLetFont);
            end
            3;
            
        end

        3;
        doFig14 = doPlots && doPrefVsWidthFigure;
        if doFig14 || doPairStats % tuning width vs preferred.

            % diff in preferred ori vs. difference in ori widths
            [dOriPref_dOriW_global_cc, dOriPref_dOriW_global_p] = corr(dOriPref_wcc, dOriW_global_wcc, 'rows', 'complete');
            [dOriPref_dOriW_local_cc,  dOriPref_dOriW_local_p]  = corr(dOriPref_wcc, dOriW_local_wcc, 'rows', 'complete');
            %%
            [bdOriPref_dOriW_global_cc, bdOriPref_dOriW_global_p] = corr(dOriPref_bcc, dOriW_global_bcc, 'rows', 'complete', 'type', 'pearson');
            [bdOriPref_dOriW_local_cc,  bdOriPref_dOriW_local_p]  = corr(dOriPref_bcc, dOriW_local_bcc, 'rows', 'complete', 'type', 'pearson');
              %%          
            [dSpfPref_dSpfW_cc,  dSpfPref_dSpfW_p]  = corr(dSpfPref_wcc, dSpf_w_wcc, 'rows', 'complete');            
            %%
            [bdSpfPref_dSpfW_cc,  bdSpfPref_dSpfW_p]  = corr(dSpfPref_bcc, dSpf_w_bcc, 'rows', 'complete', 'type', 'pearson');            
            %         dOriVdOriW_stats = struct('dOri_w_ori_global_cc', [cc_dOri_w_ori_global, p_dOri_w_ori_global], ...
            %                                   'dOri_w_ori_local',  [cc_dOri_w_ori_local,  p_dOri_w_ori_local] );
              %%          

            % diff in preferred spf vs. difference in spf width
            

            corr_dPref_dW_stats = struct(...
                'dOriPref_dOriW_global_cc', dOriPref_dOriW_global_cc,...
                'dOriPref_dOriW_global_p', dOriPref_dOriW_global_p,...
                'dOriPref_dOriW_local_cc', dOriPref_dOriW_local_cc,...
                'dOriPref_dOriW_local_p', dOriPref_dOriW_local_p,...        
                'dSpfPref_dSpfW_cc',      dSpfPref_dSpfW_cc, ...
                'dSpfPref_dSpfW_p',      dSpfPref_dSpfW_p);
            
            miscPairStats.dPrefVdW_stats = corr_dPref_dW_stats;
        
            if doFig14
               % ori width vs preferred ori

               % ori width vs spf width ?
                if  strcmp(gratingType, 'flashed') && 0
                    %%
                    oriField = iff(~subtractSpont, 'oriStats_si', 'oriStats_ss');
                    spfField = iff(~subtractSpont, 'spfStats_si', 'spfStats_ss');
                    os_oriStats = nestedFields(oriSpfCells, 'stats', 'tuningStats', oriField);
                    os_spfStats = nestedFields(oriSpfCells, 'stats', 'tuningStats', spfField);

                    os_ori_glob_w = [os_oriStats.w_ori_global];
                    os_ori_loc_w = [os_oriStats.w_ori_global];
                    os_spf_w = [os_spfStats.w_spf];

                    plot(os_ori_glob_w, os_spf_w)    ;            
                end                


               % dori width vs dspf width. ?              

%                dOriPref_wcc = dOriPref_cc(Wcc_oo_pairIdxs);
%                dOriW_glob_wcc = S_ori{ori_glob_idx}.val(Wcc_oo_pairIdxs);
%                dOriW_loc_wcc = S_ori{ori_loc_idx}.val(Wcc_oo_pairIdxs);

               %%
               ori_pairF1oDC_Wcc = ori_unitF1oDC(Wcc_pairs_ori);       
               nsimp_oo = sum(ori_pairF1oDC_Wcc > 1, 2);
               ori_idx_ss_pair = nsimp_oo== 2;
               ori_idx_cc_pair = nsimp_oo == 0;
               ori_idx_sc_pair = nsimp_oo == 1;

               corrType = 'pearson';
               %%
                figure(211); clf;
                subplotGap(1,2,1);  hold on; box on;
                plot(dOriPref_wcc(ori_idx_sc_pair), dOriW_glob_wcc(ori_idx_sc_pair), 'mo', 'markersize', sz);
                plot(dOriPref_wcc(ori_idx_ss_pair), dOriW_glob_wcc(ori_idx_ss_pair), 'bo', 'markersize', sz);
                plot(dOriPref_wcc(ori_idx_cc_pair), dOriW_glob_wcc(ori_idx_cc_pair), 'ro', 'markersize', sz);
                xlabel('Difference in preferred orientation');
                ylabel('Difference in Global Ori Width');
                ori_glob_str = getSimpleComplexPairCorr(ori_pairF1oDC_Wcc, dOriPref_wcc, dOriW_glob_wcc, corrType); % cc_ori_glob_all, cc_ori_glob_ss, cc_ori_glob_sc, cc_ori_glob_cc, str
                title([{sprintf('\\bf %s Gratings : \\Delta Pref Ori vs \\DeltaOri Global Width\\rm', titleCase(gratingType))}, ori_glob_str]);

                fprintf('\n\n*** dPrefOri vs dWOri_global (%s gratings, %s)\n', gratingType, spont_s);
                cellfun(@(str) fprintf('%s\n', str), ori_glob_str);

                subplotGap(1,2,2);  hold on; box on;
                plot(dOriPref_wcc(ori_idx_sc_pair), dOriW_loc_wcc(ori_idx_sc_pair), 'mo', 'markersize', sz);
                plot(dOriPref_wcc(ori_idx_ss_pair), dOriW_loc_wcc(ori_idx_ss_pair), 'bo', 'markersize', sz);
                plot(dOriPref_wcc(ori_idx_cc_pair), dOriW_loc_wcc(ori_idx_cc_pair), 'ro', 'markersize', sz);
                xlabel('Difference in preferred orientation');
                ylabel('Difference in Local Ori Width');
                ori_loc_str = getSimpleComplexPairCorr(ori_pairF1oDC_Wcc, dOriPref_wcc, dOriW_loc_wcc, corrType); % cc_ori_glob_all, cc_ori_glob_ss, cc_ori_glob_sc, cc_ori_glob_cc, str
                title([{sprintf('\\bf %s Gratings : \\Delta Pref Ori vs \\DeltaOri Local Width\\rm', titleCase(gratingType))}, ori_loc_str]);

                fprintf('\n\n*** dPrefOri vs dWOri_local (%s gratings, %s)\n', gratingType, spont_s);
                cellfun(@(str) fprintf('%s\n', str), ori_loc_str);

                if strcmp(gratingType, 'drifting')
                    %%
%                     pref_dir_idx = find(strcmp(measures_ori, 'D_dir_pref'), 1);
%                     dsi_idx = find(strcmp(measures_ori, 'D_dsi_glob'), 1);

%                     dDirPref_wcc = dDirPref_cc(Wcc_oo_pairIdxs);
%                     dDSI = S_ori{dsi_idx}.val(Wcc_oo_pairIdxs);                               
    %                 dDirPref_wcc = dOriPref_wcc;
%%
                    figure(212); clf;
                    hold on; box on;
                    plot(dDirPref_wcc(ori_idx_sc_pair), dDSI_wcc(ori_idx_sc_pair), 'mo', 'markersize', sz);
                    plot(dDirPref_wcc(ori_idx_ss_pair), dDSI_wcc(ori_idx_ss_pair), 'bo', 'markersize', sz);
                    plot(dDirPref_wcc(ori_idx_cc_pair), dDSI_wcc(ori_idx_cc_pair), 'ro', 'markersize', sz);
                    xlabel('Difference in preferred Direction');
                    ylabel('Difference in DSI');
                    dsi_str = getSimpleComplexPairCorr(ori_pairF1oDC_Wcc, dDirPref_wcc, dDSI_wcc, corrType); % cc_spf_glob_all, cc_spf_glob_ss, cc_spf_glob_sc, cc_spf_glob_cc, str
                    title([{'\\bf\\DeltaPref Direction vs \\Delta DSI\\rm'}, dsi_str]);
                    title([{sprintf('\\bf %s Gratings : \\DeltaPref Direction vs \\Delta DSI\\rm', titleCase(gratingType))}, dsi_str]);

    %                 xlim([-1, 181]); set(gca, 'xtick', [0:45:180])
                    xlim([-1, 181]); set(gca, 'xtick', [0:30:90])

                    fprintf('\n\n*** dPrefDir vs dDSI (%s gratings, %s)\n', gratingType, spont_s);
                    cellfun(@(str) fprintf('%s\n', str), dsi_str);
                    3;
                end




                %%
                3;
                spf_pairF1oDC_Wcc = spf_unitF1oDC(Wcc_pairs_spf);     

                nsimp_ss = sum(spf_pairF1oDC_Wcc > 1, 2);
               spf_idx_ss_pair = nsimp_ss == 2;
               spf_idx_cc_pair = nsimp_ss == 0;
               spf_idx_sc_pair = nsimp_ss == 1;




%                dSpfPref_wcc = S_spf{spf_pref_idx}.val(Wcc_ss_pairIdxs);
%                dSpf_w_wcc = S_spf{spf_w_idx}.val(Wcc_ss_pairIdxs);

                figure(213); clf;
                hold on; box on;
                plot(dSpfPref_wcc(spf_idx_sc_pair), dSpf_w_wcc(spf_idx_sc_pair), 'mo', 'markersize', sz);
                plot(dSpfPref_wcc(spf_idx_ss_pair), dSpf_w_wcc(spf_idx_ss_pair), 'bo', 'markersize', sz);
                plot(dSpfPref_wcc(spf_idx_cc_pair), dSpf_w_wcc(spf_idx_cc_pair), 'ro', 'markersize', sz);
                xlabel('Difference in pref. spatial frequency');
                ylabel('Difference in spatial frequency width');
                spf_glob_str = getSimpleComplexPairCorr(spf_pairF1oDC_Wcc, dSpfPref_wcc, dSpf_w_wcc, corrType); % cc_spf_glob_all, cc_spf_glob_ss, cc_spf_glob_sc, cc_spf_glob_cc, str
                title([{'\\bf\\DeltaPref Spf vs \\Delta Spf Width\\rm'}, spf_glob_str]);
                title([{sprintf('\\bf %s Gratings : \\DeltaPref Spf vs \\Delta Spf Width\\rm', titleCase(gratingType))}, spf_glob_str]);            

                fprintf('\n\n*** dPrefSpf vs dWSpf (%s gratings, %s)\n', gratingType, spont_s);
                cellfun(@(str) fprintf('%s\n', str), spf_glob_str);

            end
        end
        3;
        
        if doScatterVsWidthPlots || doPairStats 
            doScatterVsMeanWidthPlots = 1;
            doScatterVsStdWidthPlots = 0;
            minNCells_perSite = 3;
            require_min_pref = 1;
            require_min_width = 1;
            ori_useCircVar = 1;

            %% Orientation
            all_measures = {'orientation', 'direction', 'spatial_freq'};
%             all_measures = {'direction'};
            corrType = 'pearson';
%             corrType = 'spearman';

            set(0,'DefaultFigureWindowStyle','docked') 
            for measure_i = 1:length(all_measures)
                measure = all_measures{measure_i};

                switch measure
                    case 'orientation',
                        all_pref = [allOriUnitStats.ori_pref_deg];
                        all_widths_C = {[allOriUnitStats.w_ori_global], [allOriUnitStats.w_ori_local]};
                        idx_cellsAtEachSite = idx_cellsAtEachSite_ori;
                        if ori_useCircVar
                            circ_std = @(x) sqrt( circVar(ones(size(x)), deg2rad(x)) * 2)*180/(2*pi) ;
                        else
                            circ_std = @(x) calcOriGlobalWidth( ones(size(x)), x);
                        end
                        pref_scatter_func = circ_std;

                        pref_label = 'Preferred Ori';
                        pref_scatter_label = 'Preferred Ori Scatter (circ std)';
                        width_labels = {'Global Ori Tuning Width', 'Local Ori Tuning Width'};
                        
                        pref_label_short = 'ori';
                        pref_scatter_label_short = 'prefOriScatter';
                        width_labels_short = {'wOriGlobal', 'wOriLocal'};

                    case 'direction',
                        if strcmp(gratingType, 'flashed')
                            continue;
                        end

                        all_pref = [allOriUnitStats.Ddir_pref_smlSpkMU];
                        all_widths_C = {[allOriUnitStats.DSI_global]};

                        idx_cellsAtEachSite = idx_cellsAtEachSite_ori;

                        frac_aligned_func = @(x) nnz(x < 45) / nnz(x < 45 | x > 135);
                        aligned_prob = 0.676;
                        binom_prob_func = @(x) (binocdf( nnz(x < 45),  nnz(x < 45 | x > 135), aligned_prob));
                        
%                         pref_scatter_func = frac_aligned_func;
                        pref_scatter_func = binom_prob_func;
                        

                        pref_label = 'Preferred Dir';
                        pref_scatter_label = 'Fraction aligned to site';
                        width_labels = {'DSI'};

                        pref_label_short = 'dir';
                        pref_scatter_label_short = 'fracAligned';
                        width_labels_short = {'DSI'};

                    case 'spatial_freq',

                        all_pref = [allSpfUnitStats.f_opt];
                        all_widths_C = { [allSpfUnitStats.w_spf] };

                        idx_cellsAtEachSite = idx_cellsAtEachSite_spf;
                        pref_scatter_func = @(x) nanstd( log10(x));

                        pref_label = 'Spf Freq';
                        pref_scatter_label = 'Preferred Spf Freq Scatter (std)';
                        width_labels = {'Spf Freq Tuning Width', };
                        
                        
                        pref_label_short = 'spf';
                        pref_scatter_label_short = 'spfPref';
                        width_labels_short = {'spfWidth'};
                end
                nWidths = length(all_widths_C);

                nSites = length(idx_cellsAtEachSite);


                [pref_scatter] = deal(nan(1,nSites));
                [width_std_C, width_mean_C] = deal( repmat({nan(1,nSites)}, 1, nWidths )  );


                for site_i = 1:nSites
                    idx_thisSite = idx_cellsAtEachSite{site_i};

                    pref_thisSite = nonnans([all_pref(idx_thisSite)]);
                    pref_thisSite_ok = nonnans(pref_thisSite);
                    if (length(pref_thisSite) >= minNCells_perSite)  &&  ...            
                        (length(pref_thisSite_ok) >= minNCells_perSite  || ~require_min_pref)
        %                 circ_std = @(x) sqrt( circVar(ones(size(x)), deg2rad(x)) );

                        pref_scatter(site_i) = pref_scatter_func(pref_thisSite_ok );
                    end

                    for wi = 1:nWidths
                    %             orientations_atEachSite = cellfun(@(idx) [allOri_pref(idx)], idx_cellsAtEachSite_ori_use, 'un', 0);
                        all_widths_i_thisSite = all_widths_C{wi}(idx_thisSite);
                        all_widths_i_thisSite_ok = nonnans( all_widths_C{wi}(idx_thisSite) );

                        if (length(all_widths_i_thisSite) >= minNCells_perSite)  &&  ...            
                            (length(all_widths_i_thisSite_ok) >= minNCells_perSite  || ~require_min_width)

                            width_mean_C{wi}(site_i) = mean(all_widths_i_thisSite_ok);
                            width_std_C{wi}(site_i)  = std(all_widths_i_thisSite_ok);                    

                        end
                    end


                end

                for wi = 1:nWidths
                    n(wi) = nnz(all(~isnan([width_mean_C{wi}(:), pref_scatter(:)]), 2));
                    [cc_prefScatter_meanWidth(wi), p_prefScatter_meanWidth(wi)] = corr(width_mean_C{wi}(:), pref_scatter(:), 'rows', 'complete', 'type', corrType);
                    [cc_prefScatter_stdWidth(wi), p_prefScatter_stdWidth(wi)] = corr(width_std_C{wi}(:), pref_scatter(:), 'rows', 'complete', 'type', corrType);
                    

                    miscPairStats.scatterVsMeanWidthStats.(sprintf('%s_%s_vs_mean%s_cc', pref_label_short, pref_scatter_label_short, width_labels_short{wi})) = cc_prefScatter_meanWidth(wi);
                    miscPairStats.scatterVsMeanWidthStats.(sprintf('%s_%s_vs_mean%s_cc_p', pref_label_short, pref_scatter_label_short, width_labels_short{wi})) = p_prefScatter_meanWidth(wi);
                    miscPairStats.scatterVsStdWidthStats.(sprintf('%s_%s_vs_std%s_cc', pref_label_short, pref_scatter_label_short, width_labels_short{wi})) = cc_prefScatter_stdWidth(wi);
                    miscPairStats.scatterVsStdWidthStats.(sprintf('%s_%s_vs_std%s_cc_p', pref_label_short, pref_scatter_label_short, width_labels_short{wi})) = p_prefScatter_stdWidth(wi);                    
                end


                if doScatterVsWidthPlots
                    if doScatterVsMeanWidthPlots
                        figure(100*(1+strcmp(gratingType, 'flashed')) + measure_i*10); clf;
                        for wi = 1:nWidths

                            subplot(1,nWidths,wi); plot(width_mean_C{wi}, pref_scatter, 'o');
                            scatter_type = '';
                            if strcmp(measure, 'orientation')
                                if ori_useCircVar 
                                    scatter_type = '';
                                else
                                    scatter_type = '(Distrib Std Dev)';
                                end
                            end
                            ylabel([pref_scatter_label scatter_type]); 
                            xlabel(sprintf('%s (mean)', width_labels{wi}));
                            grat_str = [titleCase(gratingType) ' Gratings'];
                            cc_str = sprintf('cc = %.2f. p = %.2g (N = %d)', cc_prefScatter_meanWidth(wi), p_prefScatter_meanWidth(wi), n(wi));
                            title({[grat_str], cc_str})

                        end
                    end
                    
                    if doScatterVsStdWidthPlots

                        figure(100*(1+strcmp(gratingType, 'flashed')) + measure_i*10 + 1); clf;
                        for wi = 1:nWidths

                            subplot(1,nWidths,wi); plot(width_std_C{wi}, pref_scatter, 'o');
                            ylabel([pref_scatter_label scatter_type]); 
                            xlabel(sprintf('%s (std)', width_labels{wi}));
                            grat_str = [titleCase(gratingType) ' Gratings'];
                            cc_str = sprintf('cc = %.2f. p = %.2g (N = %d)', cc_prefScatter_stdWidth(wi), p_prefScatter_stdWidth(wi), n(wi));
                            title({grat_str , cc_str})

                        end
                    end
                end
                
                %%


            end
            
        end

    end
    3;
    
    3;
%     
%             % figure 4 and diff ori stats:
%         dOriPref_cc = S_ori{ori_pref_idx}.val;        
%         dOriPref_wcc = dOriPref_cc(Wcc_oo_pairIdxs);
%         dOriPref_wcc_norm = dOriPref_cc(Wcc_oo_pairIdxs_norm);
%         
%                 [allOriUnitStats(Wcc_pairs_ori(:,1)).w_ori_global];
%             w_ori_global2 = [allOriUnitStats(Wcc_pairs_ori(:,2)).w_ori_global];
%             w_ori_local1 = [allOriUnitStats(Wcc_pairs_ori(:,1)).w_ori_local];
%             w_ori_local2 = [allOriUnitStats(Wcc_pairs_ori(:,2)).w_ori_local];
%         
%         dOriW_global_wcc = S_ori{ori_glob_idx}.val(Wcc_oo_pairIdxs);
%         dOriW_local_wcc = S_ori{ori_loc_idx}.val(Wcc_oo_pairIdxs);
% 
%         dOriW_global_bcc = S_ori{ori_glob_idx}.val(Bcc_oo_pairIdxs);
%         dOriW_local_bcc = S_ori{ori_loc_idx}.val(Bcc_oo_pairIdxs);
% 
        
        
    
    
%     idx_cellsAtEachSite_ori, idx_pairsAtEachSite_ori, 
%     idx_cellsAtEachSite_spf, idx_pairsAtEachSite_spf, 
    if dispMiscStats && (subtractSpont == 1)
        %%
        cellsDistrib_matFile = getFileName('cellDistribs');
        cellsDistrib_S.cmpType = 'degree';
        cellsDistrib_S.gratingType = gratingType;
        cellsDistrib_S.subtractSpont = subtractSpont;
        cellsDistrib_S.bccType = bccType;
        cellsDistrib_S.preserveSC = preserveSC;
        cellsDistrib_S.preserveAB = preserveAB;
       3;
       % cells / site; sites/penetration, penetrations/animal.
       stimTypes = {'ori', 'spf'};
       for si = 1:length(stimTypes)
           switch stimTypes{si}
               case 'ori',
                   allLocData = [allOriCells.locData];
                   cellSiteIds = [allOriCells.Gid];
               case 'spf',
                   allLocData = [allSpfCells.locData];
                   cellSiteIds = [allSpfCells.Gid];
           end
           cellPenIds = [allLocData.PenId];
           cellCatIds = [allLocData.CatId];
           %%
           [uSiteIds, cellSiteIdxs, nCellsPerSite] = uniqueList(cellSiteIds);
           [uPenIds, cellPenIdxs, nCellsPerPen] = uniqueList(cellPenIds);
           [uCatIds, cellCatIdxs, nCellsPerCat] = uniqueList(cellCatIds);

           idx_firstCellInSite = cellfun(@(x) x(1), cellSiteIdxs);          
           sitePenIds = cellPenIds(idx_firstCellInSite);
           siteCatIds = cellCatIds(idx_firstCellInSite);
                      
           [uPenIds2, sitePenIdxs, nSitesPerPen] = uniqueList(sitePenIds);
           [uCatIds2, siteCatIdxs, nSitesPerCat] = uniqueList(siteCatIds);
           
           idx_firstCellInPen = cellfun(@(x) x(1), cellPenIdxs);          
           penCatIds = cellCatIds(idx_firstCellInPen);
           [uCatIds3, penCatIdxs, nPensPerCat] = uniqueList(penCatIds);
           
           assert(isequal(uPenIds, uPenIds2));
           assert(isequal(uCatIds, uCatIds2));
           assert(isequal(uCatIds, uCatIds3));

           prctiles_vals = [5, 50, 95];
%            cellsPerSite_medstats = prctile(nCellsPerSite, prctiles_vals);
%            cellsPerPen_medstats = prctile(nCellsPerPen, prctiles_vals);
%            cellsPerCat_medstats = prctile(nCellsPerCat, prctiles_vals);
%            sitesPerPen_medstats = prctile(nSitesPerPen, prctiles_vals);
%            sitesPerCat_medstats = prctile(nSitesPerCat, prctiles_vals);           
%            pensPerCat_medstats = prctile(nPensPerCat, prctiles_vals);
%       
           allMeasureNames_S = {'cellsPerSite', cellsPerSite, ...
                              'cellsPerPen', cellsPerPen, ...
                              'cellsPerCat', cellsPerCat, ...
                              'sitesPerPen', sitesPerPen, ...
                              'sitesPerCat', sitesPerCat, ...
                              'pensPerCat', pensPerCat};
           allMeasureNames = fieldnames(allMeasureNames_S);
                          
           for mi = 1:length(allMeasureNames);
               ms_name = allMeasureNames{mi};
                x = allMeasureNames_S.(ms_name);
                [vals_p5, vals_med, vals_p95] = dealV(prctile(x, [5, 50, 95]));
                vals_mean = mean(x);
                vals_std = std(x);
                cellsDistrib_S.(ms) = struct('mean', vals_mean, 'std', vals_std, 'median', vals_median, 'P5', vals_p5, 'P95', vals_p95 ); 
                idx = strfind(ms_name, 'Per');
                nm1 = titleCase(ms_name(1:idx-1));
                nm2 = ms_name(idx+3:end);
                fprintf('N %s / %s : %.1f : %.1f : %.1f [%.1f +- %.1f]\n', nm1, nm2, vals_p5, vals_med, vals_p95, vals_mean, vals_std);
           
             
                miscDistribStats.(['nCells_' stimType]) = length(cellSiteIds);
                miscDistribStats.(['nSites_' stimType]) = length(uSiteIds);
                miscDistribStats.(['nPens_' stimType]) = length(uPenIds);
                miscDistribStats.(['nCats_' stimType]) = length(uCatIds);
                
                     
%                 [cellsPerSite_medstats, cellsPerSite_meanstats] = getPrctilesMeanStd(nCellsPerSite, prctiles_vals);
%            [cellsPerPen_medstats,  cellsPerPen_meanstats] = getPrctilesMeanStd(nCellsPerPen, prctiles_vals);
%            [cellsPerCat_medstats,  cellsPerCat_meanstats] = getPrctilesMeanStd(nCellsPerCat, prctiles_vals);
%            [sitesPerPen_medstats,  sitesPerPen_meanstats] = getPrctilesMeanStd(nSitesPerPen, prctiles_vals);
%            [sitesPerCat_medstats,  sitesPerCat_meanstats] = getPrctilesMeanStd(nSitesPerCat, prctiles_vals);           
%            [pensPerCat_medstats,   pensPerCat_meanstats] = getPrctilesMeanStd(nPensPerCat, prctiles_vals);
           
%                fprintf('\n\n--------\n%s measures (%s gratings)\n', stimTypes{si}, curGratingType(''));
%                fprintf('N Cells: %d. NSites: %d. NPens: %d. NCats: %d\n', length(cellSiteIds), length(uSiteIds), length(uPenIds), length(uCatIds));
%                fprintf('N cells / site: %.1f : %.1f : %.1f [%.1f +- %.1f]\n', cellsPerSite_medstats, cellsPerSite_meanstats);
%                fprintf('N cells / pen:  %.1f : %.1f : %.1f [%.1f +- %.1f]\n', cellsPerPen_medstats, cellsPerPen_meanstats);
%                fprintf('N cells / cat: %.1f : %.1f : %.1f  [%.1f +- %.1f]\n', cellsPerCat_medstats, cellsPerCat_meanstats);
% 
%                fprintf('N sites / pen: %.1f : %.1f : %.1f  [%.1f +- %.1f]\n', sitesPerPen_medstats, sitesPerPen_meanstats);
%                fprintf('N sites / cat: %.1f : %.1f : %.1f  [%.1f +- %.1f]\n', sitesPerCat_medstats, sitesPerCat_meanstats);
% 
%                fprintf('N pens / cat: %.1f : %.1f : %.1f   [%.1f +- %.1f] \n', pensPerCat_medstats, pensPerCat_meanstats);


            end
            
                   
        cellsDistrib_S.allMeasureNames = allMeasureNames_save;
        cellsDistrib_S.measures_ori = measures_ori;
        cellsDistrib_S.measures_spf = measures_spf;
        cellsDistrib_S.columns = fieldnames(cellsDistrib_S.(measure_fld));
        cellsDistrib_S.miscStats = miscDistribStats;
%         save(cellsDistrib_matFile, '-struct', 'cellsDistrib_S');      
           
           
       end
       
       3;
       
       
       
%%
%{
        %% Direction
        minNCells_perSite = 3;
        require_min_dir_pref = 1;
        require_min_dir_width = 1;

        allDir_pref_dMU = [allOriUnitStats.Ddir_pref_smlSpkMU];
        allDSI = [allOriUnitStats.DSI_global];

        nSites = length(idx_cellsAtEachSite_ori);
%         nCellsPerSite_ori = cellfun(@length, idx_cellsAtEachSite_ori);
%         idx_sitesUse_ori = nCellsPerSite_ori >= minNCells_perSite;
%         idx_cellsAtEachSite_ori_use = idx_cellsAtEachSite_ori(idx_sitesUse_ori);
        
            dir_pref_dMU_atEachSite = cellfun(@(idx) allDir_pref_dMU(idx), idx_cellsAtEachSite_ori, 'un', 0);
%             oriWidths_global_atEachSite = cellfun(@(idx) [allOri_width_global(idx)], idx_cellsAtEachSite_ori_use, 'un', 0);
%             oriWidths_local_atEachSite = cellfun(@(idx) [allOri_width_local(idx)], idx_cellsAtEachSite_ori_use, 'un', 0);

        [dir_fracAligned, mean_DSI] = deal(nan(1,nSites));
        for site_i = 1:nSites
            idx_thisSite = idx_cellsAtEachSite_ori{site_i};
            
            dir_pref_dMU_thisSite = nonnans([allDir_pref_dMU(idx_thisSite)]);
            dir_pref_dMU_thisSite_ok = nonnans(dir_pref_dMU_thisSite);
            if (length(dir_pref_dMU_thisSite) >= minNCells_perSite)  &&  ...            
                (length(dir_pref_dMU_thisSite_ok) >= minNCells_perSite  || ~require_min_ori_pref)
        %         circ_std = @(x) sqrt( circVar(ones(size(x)), deg2rad(x)) );
%                 circ_std = @(x) calcOriGlobalWidth( ones(size(x)), x);
%                 circ_std = @(x) std(x);
                func = @(x) nnz(x < 45) / nnz(x < 45 | x > 135);
                
                dir_fracAligned(site_i) = func(dir_pref_dMU_thisSite_ok );
            end
            
            
            %             orientations_atEachSite = cellfun(@(idx) [allOri_pref(idx)], idx_cellsAtEachSite_dir_use, 'un', 0);
            DSIs_thisSite = allDSI(idx_thisSite);
            DSIs_thisSite_ok = nonnans(DSIs_thisSite);
            if (length(DSIs_thisSite) >= minNCells_perSite)  &&  ...            
                (length(DSIs_thisSite_ok) >= minNCells_perSite  || ~require_min_dir_width)

                mean_DSI(site_i) = mean(DSIs_thisSite_ok);
            end
                        
        end
        
        n_global = nnz(all(~isnan([dir_fracAligned(:), mean_DSI(:)]), 2));
        [cc_oriStd_meanOriWidth_global, p_oriStd_meanOriWidth_global] = corr(dir_fracAligned(:), mean_DSI(:), 'rows', 'complete');
        
        figure(71); clf;
        plot(mean_DSI, dir_fracAligned, 'o');
        ylabel('Proportion of aligned cells'); xlabel(sprintf('mean DSI '));
        grat_str = [titleCase(gratingType) ' Gratings'];
        cc_str = sprintf('cc = %.2f. p = %.2g (N = %d)', cc_oriStd_meanOriWidth_global, p_oriStd_meanOriWidth_global, n_global);
        title({grat_str , cc_str})
        3;
        

       
       
        %% Spatial frequency

        
        minNCells_perSite = 3;
        require_min_spf_pref = 1;
        require_min_spf_width = 1;

        allSpf_pref = [allSpfUnitStats.f_opt];
        allSpf_width = [allSpfUnitStats.w_spf];

        nSites = length(idx_cellsAtEachSite_spf);

        [spfPref_scatter, mean_spfWidth] = deal(nan(1,nSites));
        for site_i = 1:nSites
            idx_thisSite = idx_cellsAtEachSite_spf{site_i};
            
            spfPref_thisSite = nonnans([allSpf_pref(idx_thisSite)]);
            spfPref_thisSite_ok = nonnans(spfPref_thisSite);
            if (length(spfPref_thisSite) >= minNCells_perSite)  &&  ...            
                (length(spfPref_thisSite_ok) >= minNCells_perSite  || ~require_min_spf_pref)
                std_func = @(x) nanstd( log10(x)); 

                spfPref_scatter(site_i) = std_func(spfPref_thisSite_ok );
            end
            
            spfWidths_thisSite = allSpf_width(idx_thisSite);
            spfWidths_thisSite_ok = nonnans(spfWidths_thisSite);
            if (length(spfWidths_thisSite) >= minNCells_perSite)  &&  ...            
                (length(spfWidths_thisSite_ok) >= minNCells_perSite  || ~require_min_spf_width)

                mean_spfWidth(site_i) = mean(spfWidths_thisSite_ok);
            end
            
        end
        
        n = nnz(all(~isnan([spfPref_scatter(:), mean_spfWidth(:)]), 2));
        [cc_spfStd_meanSpfWidth, p_spfStd_meanSpfWidth] = corr(spfPref_scatter(:), mean_spfWidth(:), 'rows', 'complete');

        
  
%         spf_pref_atEachSite = cellfun(@(idx) [allSpf_pref(idx)], idx_cellsAtEachSite_spf_use, 'un', 0);
%         
% %         circ_std = @(x) sqrt( circVar(ones(size(x)), deg2rad(x)) );
%         std_func = @(x) nanstd( log10(x)); 
% %         circ_std = @(x) std(x);
%         spf_scatter = cellfun(std_func, spf_pref_atEachSite);
%         
%         spf_widths_atEachSite = cellfun(@(idx) [allSpf_width(idx)], idx_cellsAtEachSite_spf_use, 'un', 0);
% 
%         mean_spfWidth = cellfun(@nanmean, spf_widths_atEachSite);        
        
        figure(70); clf;
        plot(mean_spfWidth, spfPref_scatter, 'o');
        ylabel('Spatial Frequency Scatter (std)'); xlabel(sprintf('Spatial Frequency Width (mean)'));
        grat_str = [titleCase(gratingType) ' Gratings'];
        cc_str = sprintf('cc = %.2f. p = %.2g (N = %d)', cc_spfStd_meanSpfWidth, p_spfStd_meanSpfWidth, n);
        title({grat_str , cc_str})
        
        
        
  %}      
        
        
        
        
    end
    
    
    
    if (doSimpleComplexStats) && (subtractSpont == 1)
        
        %% Get Orientation F1/DC stats:        
%         ori_unitF1oDC = [allOriUnits.F1oDC_maxR_avP_sm];
%         ori_pairF1oDC_Wcc = ori_unitF1oDC(Wcc_pairs_ori);       
%         ori_pairF1oDC_Bcc = ori_unitF1oDC(Bcc_pairs_ori);      
%         ori_pairF1oDC_Wrcc = cellfun(@(idxs) ori_unitF1oDC(idxs), Wrcc_pairs_ori, 'un', 0);        
        sc_comments = {};
        
                
        if curPreserveSimpleComplex
            sc_opt.requireSameNumberOfSCpairsInShuffles = 1;
            sc_opt.controlMethod = 'cell_shuffle';
            sc_opt.justRandomizeSClabels = false; 
            sc_opt.calcPairTypeSig = 'indiv';
%             sc_opt.clusterIdxUses = 'brcc';
            sc_opt.clusterIdxUses = 'bcc';
            sc_opt.calcStatRatioCIs = true;
            sc_opt.nBoots = nBoots;
            
        elseif ~curPreserveSimpleComplex
            sc_opt.requireSameNumberOfSCpairsInShuffles = 0;
            sc_opt.controlMethod = 'SC_labelShuffle';
            sc_opt.justRandomizeSClabels = true;
            sc_opt.calcPairTypeSig = 'differences';
            sc_opt.clusterIdxUses = 'brcc';
            sc_opt.calcStatRatioCIs = false;
        end
        
        sc_opt.getBrccPairs = strcmp(sc_opt.clusterIdxUses, 'brcc');
        
        
        %% Get Orientation F1/DC stats:
        ori_unit_simple_tf = double(ori_unitF1oDC >= 1);
        ori_unit_simple_tf(isnan(ori_unitF1oDC)) = nan;
        ori_cell_simple_tf = ori_unit_simple_tf(oriCells_use);
                        
        nSimple_ori_use = nnz(ori_cell_simple_tf == 1 );
        nComplex_ori_use = nnz( ori_cell_simple_tf == 0 );

        ori_nSimp_Wcc = sum( ori_unit_simple_tf(Wcc_pairs_ori), 2); 
        ori_nSimp_Bcc = sum( ori_unit_simple_tf(Bcc_pairs_ori), 2); 
        ori_nSimp_Wrcc = cellfun( @(idxs) sum( ori_unit_simple_tf(idxs), 2), Wrcc_pairs_ori, 'un', 0);
        nSS_ori = nnz(ori_nSimp_Wcc == 2);
        nSC_ori = nnz(ori_nSimp_Wcc == 1);
        nCC_ori = nnz(ori_nSimp_Wcc == 0);
        
        s_ori = sprintf('For %s gratings (orientation) %s: %d simple cells and %d complex cells (total of %d)', gratingType, spont_s, nSimple_ori_use, nComplex_ori_use, nnz( ~isnan(ori_cell_simple_tf) ));
        if showWorking
            fprintf('%s\n', s_ori);            
        end
        sc_comments = [sc_comments, s_ori];        
                                
        %% Get Spatial Frequency F1/DC stats:

        
        spf_unit_simple_tf = double( spf_unitF1oDC >= 1);
        spf_unit_simple_tf (isnan(spf_unitF1oDC)) = nan;
%         spf_complex_tf = spf_unitF1oDC < 1;        
        spf_cell_simple_tf = spf_unit_simple_tf(spfCells_use);
%         spf_cell_complex_tf = spf_complex_tf(spfCells_use);
        
        nSimple_spf_use  = nnz(spf_cell_simple_tf == 1 );
        nComplex_spf_use = nnz(spf_cell_simple_tf == 0 );
        
        spf_nSimp_Wcc = sum( spf_unit_simple_tf(Wcc_pairs_spf), 2); 
        spf_nSimp_Bcc = sum( spf_unit_simple_tf(Bcc_pairs_spf), 2); 
        spf_nSimp_Wrcc = cellfun( @(idxs) sum( spf_unit_simple_tf(idxs), 2), Wrcc_pairs_spf, 'un', 0);        
        nSS_spf = nnz(spf_nSimp_Wcc == 2);
        nSC_spf = nnz(spf_nSimp_Wcc == 1);
        nCC_spf = nnz(spf_nSimp_Wcc == 0);
                
        s_spf = sprintf('For %s gratings (spatial frequency) %s: %d simple cells and %d complex cells (total of %d)\n', gratingType, spont_s, nSimple_spf_use, nComplex_spf_use, nnz(~isnan(spf_cell_simple_tf)));
        if showWorking
            fprintf('%s\n', s_spf);            
        end
        sc_comments = [sc_comments, ' ', s_spf];
        
        %%
        
        miscSCstats.nSimpComp = struct('n_ori_simp', nSimple_ori_use, 'n_ori_comp', nComplex_ori_use, ...
                                       'n_spf_simp', nSimple_spf_use, 'n_spf_comp', nComplex_spf_use, ...
                                       'n_ori_SSpairs', nSS_ori, 'n_ori_SCpairs', nSC_ori, 'n_ori_CCpairs', nCC_ori, ...
                                       'n_spf_SSpairs', nSS_spf, 'n_spf_SCpairs', nSC_spf, 'n_spf_CCpairs', nCC_spf ......
                                   );
                
        
%         spf_nSimp_Wcc = sum(spf_pairF1oDC_Wcc > 1, 2);
%         spf_nSimp_Bcc = sum(spf_pairF1oDC_Bcc > 1, 2);        
%         spf_nSimp_Wrcc = cellfun( @(f1odcs) sum(f1odcs > 1, 2), spf_pairF1oDC_Wrcc, 'un', 0);

        3; 
        
        
        nResamplesTotal = 10000;
        
%         randInfo_ori = Wrcc_params_ori;
%         randInfo_ori.nUnits = nOriUnits;
%         randInfo_ori.idxMtx = idxMtx_ori;
        randInfo_ori.unit_simple_tf = ori_unit_simple_tf;
        randInfo_ori.nResamplesTotal = nResamplesTotal;
        randInfo_ori.getBrccPairs = sc_opt.getBrccPairs;
        
%         randInfo_spf = Wrcc_params_spf;
%         randInfo_spf.nUnits = nSpfUnits;
%         randInfo_spf.idxMtx = idxMtx_spf;
        randInfo_spf.unit_simple_tf = spf_unit_simple_tf;
        randInfo_spf.nResamplesTotal = nResamplesTotal;
        randInfo_spf.getBrccPairs = sc_opt.getBrccPairs;
        %%
        
        getRandomizedPairsNow = (nResamplesTotal < 1000) && sc_opt.getBrccPairs;
        
        if getRandomizedPairsNow 
            %%
            Wrcc_oo_pairIdxs_saved = Wrcc_oo_pairIdxs;
            Wrcc_ss_pairIdxs_saved = Wrcc_ss_pairIdxs;
            
            [Wrcc_oo_pairIdxs, ori_nSimp_Wrcc,    Brcc_oo_pairIdxs, ori_nSimp_Brcc] = getRandomizedPairs(randInfo_ori, true, nResamplesTotal);
            [Wrcc_ss_pairIdxs, spf_nSimp_Wrcc,    Brcc_ss_pairIdxs, spf_nSimp_Brcc] = getRandomizedPairs(randInfo_spf, true, nResamplesTotal);
            
            doCheck = true;
            if doCheck
                assert(isequal(Wrcc_oo_pairIdxs,   Wrcc_oo_pairIdxs_saved(1:nResamplesTotal)) );                
                assert(isequal(Wrcc_ss_pairIdxs,   Wrcc_ss_pairIdxs_saved(1:nResamplesTotal)) );
            end
                        
        else
            [Brcc_oo_pairIdxs, ori_nSimp_Brcc] = deal([]);
            [Brcc_ss_pairIdxs, spf_nSimp_Brcc] = deal([]);            
        end
        3;
        
        
        %% probability of simple/complex pairing
        spf_pairTypeIdxs = {Wcc_ss_pairIdxs, Wrcc_ss_pairIdxs, Bcc_ss_pairIdxs, Brcc_ss_pairIdxs};        
%         allGC_cells = [allSpfCells.Gid]*10000 + [allSpfCells.cellId];
        allGC_units = [allSpfUnits.Gid]*10000 + [allSpfUnits.cellId];
        
        idx_unitsInPairs = binarySearch(allGC_units, unique(allGC_units(Wcc_pairs_spf(:))), [], 0);
        
        ncells_used = length(idx_unitsInPairs);
        nsimp_used = nnz( spf_unitF1oDC(idx_unitsInPairs) >=1 );
        ncomp_used = nnz( spf_unitF1oDC(idx_unitsInPairs) < 1 );
        
        p_simple = nsimp_used/ncells_used;
        p_complex = 1-p_simple;
        
        frac_ss_exp = p_simple^2;
        frac_cc_exp = p_complex^2;
        frac_sc_exp = 2*p_simple*p_complex;
        
        Frac_ss = @(x) nnz(x == 2)/length(x);
        Frac_sc = @(x) nnz(x == 1)/length(x);
        Frac_cc = @(x) nnz(x == 0)/length(x);
        s_idx = find(strcmp(measures_spf, 'D_F1pair_spf'), 1);
        for pt_i = 1:length(pairTypes)
            allVals = S_spf{ s_idx }.val;
            pairIdxs = spf_pairTypeIdxs{pt_i};
            if ~iscell(pairIdxs)
                %                         vals = allVals ( pairIdxs );
                num_ss{pt_i} = Frac_ss(allVals ( pairIdxs ));
                num_sc{pt_i} = Frac_sc(allVals ( pairIdxs ));
                num_cc{pt_i} = Frac_cc(allVals ( pairIdxs ));
            else
                num_ss{pt_i} = cellfun(@(idxs) Frac_ss( allVals ( idxs ) ), pairIdxs ) ;
                num_sc{pt_i} = cellfun(@(idxs) Frac_sc( allVals ( idxs ) ), pairIdxs ) ;
                num_cc{pt_i} = cellfun(@(idxs) Frac_cc( allVals ( idxs ) ), pairIdxs ) ;
            end
            
        end
        
        
        %             gids = pairdata.gids(pairtypeidxs{pt_i}(idx_nonnans),:);
        %             cids = pairdata.cellids(pairtypeidxs{pt_i}(idx_nonnans),:);
        %             allGC = unique([gids(:), cids(:)], 'rows');
        %             ncl{spont_i} = size(allGC,1);
        %%
        nwcc = length(spf_pairTypeIdxs{1});
        n_sc_obs = num_sc{1}*nwcc;
        n_sc_exp = frac_sc_exp*nwcc;
        
        p_sc_pair = 2*p_simple*p_complex;
        p_sc_obs = binocdf(n_sc_obs, nwcc, p_sc_pair);
                
        sc_pair{1} = sprintf('of a total of %d cells at sites containing >1 cell, %d were simple, %d were complex', ncells_used, nsimp_used, ncomp_used);
        sc_pair{2} = sprintf('expected: %.1f%% ss pairs, %.1f%% sc pairs (~%.1f), %.1f%% cc pairs', frac_ss_exp*100, frac_sc_exp*100, n_sc_exp, frac_cc_exp*100);
        sc_pair{3} = sprintf('observed: %.1f%% (%d) ss pairs, %.1f%% (%d) sc pairs, %.1f%% (%d) cc pairs', num_ss{1}*100, num_ss{1}*nwcc, num_sc{1}*100, num_sc{1}*nwcc, num_cc{1}*100, round(num_cc{1}*nwcc));
        sc_pair{4} = sprintf('probability of observing %d sc pairs (or fewer) instead of (expected) %.1f if random (with p_sc = %.3f): p = %.3g', n_sc_obs, n_sc_exp, p_sc_pair, p_sc_obs);        
        if showWorking            
            cellfun(@(str) fprintf('%s\n', str), sc_pair);
        end
        
        sc_comments = [sc_comments, ' ', sc_pair];
        3;
                

%         
        
        
%         controlMethod = 'bootstrap';
%         controlMethod = 'SC_labelShuffle';
      
        
%         compareSimpleComplexPairingStats('Local Orientation width', 'Dw_ori_loc_si', 60, measures_ori, S_ori, ori_allNSimp, ori_pairTypeIdxs);
        3;
        
        
        if doSimpleComplexStats
        
            if strcmp(sc_opt.controlMethod, 'cell_shuffle') && sc_opt.requireSameNumberOfSCpairsInShuffles
                if ~curPreserveSimpleComplex
                    error('Use datafile where # simple/complex at a site are preserved');
                end
                doCheckNSimpleNComplexPairs = 1;
                if doCheckNSimpleNComplexPairs
                    % make sure Wrcc pairs have same # of simple/complex pairs
                    ori_nSimpMtx = bsxfun(@plus, ori_unit_simple_tf(:), ori_unit_simple_tf(:)');
                    [uSC, sc_count] = uniqueCount( ori_nSimpMtx(Wcc_oo_pairIdxs_M) );
                    [uSC2, sc_count2] = uniqueCount( ori_nSimp_Wcc );
                    assert(isequal(sc_count, sc_count2));

                    for i = 1 : min(Npermutes, nResamplesTotal)
                        [uSC_i, sc_count_i] = uniqueCount( ori_nSimpMtx(Wrcc_oo_pairIdxs_M{i}) );
                        [uSC_i, sc_count_i2] = uniqueCount( ori_nSimp_Wrcc{i} );
                        assert(isequal(sc_count_i, sc_count_i2));
                        assert(isequal(sc_count, sc_count_i));
                    end

                    spf_nSimpMtx = bsxfun(@plus, spf_unit_simple_tf(:), spf_unit_simple_tf(:)');
                    [uSC, sc_count] = uniqueCount( spf_nSimpMtx(Wcc_ss_pairIdxs_M) );
                    [uSC2, sc_count2] = uniqueCount( spf_nSimp_Wcc );
                    assert(isequal(sc_count, sc_count2));
                    for i = 1 : min(Npermutes, nResamplesTotal)
                        [uSC_i, sc_count_i] = uniqueCount( spf_nSimpMtx(Wrcc_ss_pairIdxs_M{i}) );
                        [uSC_i, sc_count_i2] = uniqueCount( spf_nSimp_Wrcc{i} );
                        assert(isequal(sc_count_i, sc_count_i2));
                        assert(isequal(sc_count, sc_count_i));
                    end
                    3;
                end
            end
        
            ori_pairTypeIdxs = {Wcc_oo_pairIdxs, Wrcc_oo_pairIdxs,   Bcc_oo_pairIdxs,   Brcc_oo_pairIdxs};
            ori_allNSimp     = {ori_nSimp_Wcc,   ori_nSimp_Wrcc,     ori_nSimp_Bcc,     ori_nSimp_Brcc};

            spf_pairTypeIdxs = {Wcc_ss_pairIdxs, Wrcc_ss_pairIdxs,   Bcc_ss_pairIdxs,   Brcc_ss_pairIdxs};        
            spf_allNSimp     = {spf_nSimp_Wcc,   spf_nSimp_Wrcc,     spf_nSimp_Bcc,     spf_nSimp_Brcc};

            
            oe_mode = curDegreeOEmode;
            oe_offset = switchh(oe_mode, {'aa', 'oe_same', 'oe_diff'}, [0, 1000, 2000]);
            cell_stat_figIds = [51:2:61];
            pair_stat_figIds = [52:2:62] + (curGratingType == 1)*100 + oe_offset;
            showCellStatFigs = 0;
            showPairStatFigs = 1;
            
            showWorking = 1;
            
            printStats = 1 && showWorking;
            if ~showCellStatFigs || ~showWorking
                cell_stat_figIds(:) = nan; 
            end
            if ~showPairStatFigs || ~showWorking
                pair_stat_figIds(:) = nan;
            end
            %%
            simpCompCells_matFile = getFileName('scCellStats');
            simpCompPairs_matFile = getFileName('scPairStats');

            simpCompCells_S.cmpType = 'degree';
            simpCompCells_S.gratingType = gratingType;
            simpCompCells_S.subtractSpont = subtractSpont;       
            simpCompCells_S.bccType = bccType;
            simpCompCells_S.preserveSC = preserveSC;            
            simpCompCells_S.preserveAB = preserveAB;
            simpCompPairs_S = simpCompCells_S;

            %%
            allStatNames    = {'Global Orientation Width',    'Local Orientation Width',   'Preferred Orientation', 'Direction Selectivity Index', 'Spatial Frequency Tuning Width', 'Preferred Spatial Frequency'};
            allCellMeasureNames = {'w_ori_global',                'w_ori_local',            'ori_pref_deg',             'DSI_global'                'w_spf',                            'f_opt'};
            allDiffMeasureNames = {'Dw_ori_glob',                 'Dw_ori_loc',              'D_ori_pref',              'D_dsi_glob',                  'Dw_spf',                         'D_spf_pref'};                

            do_spf_only = 1;
            if do_spf_only
                allStatNames = allStatNames(5:6);
                allCellMeasureNames = allCellMeasureNames(5:6);
                allDiffMeasureNames = allDiffMeasureNames(5:6);
            end
%             idx_order = [3 1 2 4 5 6];
%             allStatNames = allStatNames(idx_order);
%             allCellMeasureNames = allCellMeasureNames(idx_order);
%             allDiffMeasureNames = allDiffMeasureNames(idx_order);
            
            
            if strcmp(gratingType, 'flashed')
                idx_dsi = find(strcmp(allDiffMeasureNames, 'D_dsi_glob'), 1);
                allStatNames(idx_dsi) = [];
                allCellMeasureNames(idx_dsi) = [];
                allDiffMeasureNames(idx_dsi) = [];
            end

            [cell_stat_labels, pairing_labels] = deal(cell(1, length(allStatNames)));
            for si = 1:length(allStatNames)
                stat_name = allStatNames{si};
                cell_measure_name = allCellMeasureNames{si};
                diff_measure_name = allDiffMeasureNames{si};
                [cell_stat_labels{si}, pairing_labels{si}] = deal(stat_name);
                isOri = isempty(strfind(stat_name, 'Spatial'));
                if isOri
                    cell_simple_tf = ori_cell_simple_tf;                
                    all_diff_measures_i = measures_ori;
                    S_i = S_ori;
                    allCellVals = [allOriCellStats.(cell_measure_name)];
                    allNSimp = ori_allNSimp;
                    pairTypeIdxs = ori_pairTypeIdxs;
                    randInfo = randInfo_ori;
                else
                    cell_simple_tf = spf_cell_simple_tf;
                    all_diff_measures_i = measures_spf;
                    S_i = S_spf;
                    allCellVals = [allSpfCellStats.(cell_measure_name)];
                    allNSimp = spf_allNSimp;
                    pairTypeIdxs = spf_pairTypeIdxs;
                    randInfo = randInfo_spf;
                end

                cell_stats(si)    = compareSimpleComplexCellStats(   stat_name, allCellVals,  cell_simple_tf, cell_stat_figIds(si), printStats); 
                if any(strcmp(sc_opt.controlMethod, {'cell_shuffle', 'SC_labelShuffle'}))
                    pairing_stats(si) = compareSimpleComplexPairingStats(stat_name, diff_measure_name, all_diff_measures_i, S_i, allNSimp, pairTypeIdxs, pair_stat_figIds(si), printStats, randInfo, sc_opt); 
                elseif strcmp(sc_opt.controlMethod, 'bootstrap')
                    pairing_stats(si) = compareSimpleComplexPairingStats_boot(stat_name, diff_measure_name, all_diff_measures_i, S_i, allNSimp, pairTypeIdxs, pair_stat_figIds(si), printStats); 
                else
                    error('unknown SC pairing test');
                end
                   

                simpCompCells_S.(cell_measure_name) = cell_stats(si);
                simpCompPairs_S.(diff_measure_name) = pairing_stats(si);

            end

            simpCompCells_S.allMeasureNames = allCellMeasureNames;        
            simpCompCells_S.columns = fieldnames(cell_stats(1));        
            simpCompCells_S.comments = sc_comments;        
            simpCompCells_S.miscStats = miscSCstats;        

            simpCompPairs_S.allMeasureNames = allDiffMeasureNames;                    
            simpCompPairs_S.columns = fieldnames(pairing_stats(1));        
            saveTofile = 1;
            if saveTofile
                save(simpCompCells_matFile, '-struct', 'simpCompCells_S');             
                save(simpCompPairs_matFile, '-struct', 'simpCompPairs_S');             
            else
                warning('Skipping the saving')
            end
            
    %         %% Preferred Orientation 
    %         % Single cell Preferred Orientation;  % Differences in preferred Orientation
    %         s_i = 3; stat_name = 'Preferred Orientation';  [cell_stat_labels{s_i}, pairing_labels{s_i}] = deal(stat_name);        
    %         cell_stats(s_i) = compareSimpleComplexCellStats(stat_name, [allOriCellStats.ori_pref_deg], ori_cell_simple_tf,  cell_stat_figIds(s_i), printStats);
    %         pairing_stats(s_i) = compareSimpleComplexPairingStats(stat_name, 'D_ori_pref',  measures_ori, S_ori, ori_allNSimp, ori_pairTypeIdxs, pair_stat_figIds(s_i), printStats);


            %%
            3;
            if showWorking
                fprintf('\n\n\n');
                    fprintf('            ^Property             ^----------------------------------Cell Stats -------------------------    ---------------------------- Pairing Stats ------------------------------------- \n');
                    fprintf('        ^================          ^S_mean,     ^C_mean,      ^pT,      ^S_med,       ^C_med,       ^pU       ^medR_all,  ^medR_ss,  ^medR_cc,   ^medR_sc,  ^p(SS~=CC), ^p(SS~=SC), ^p(CC~=SC) \n');        
                for i = 1:length(cell_stats)
                    C = cell_stats(i);
                    P = pairing_stats(i);
                    fprintf('^%30s: ', strrep(cell_stat_labels{i}, ' ', '@') )            
                    fprintf(  '^%10.2f ^%10.2f ^%10.2g ^%10.2f ^%10.2f ^%10.2g ', C.mean_simple, C.mean_complex, C.pval_t, C.median_simple, C.median_complex, C.pval_U);
                    fprintf(  '^%10.2f ^%10.2f ^%10.2f ^%10.2f ^%10.4f ^%10.4f ^%10.4f', P.medianRatio_all, P.medianRatio_ss, P.medianRatio_cc, P.medianRatio_sc, P.medianProb_ss_cc, P.medianProb_ss_sc, P.medianProb_cc_sc);
                    fprintf('\n');
                end
            end


        end
            3;
        

        
        
    end
    
    
    
    if doDepthAndDistancePlots
        %%
        fprintf('\n==== Doing Depth & Distance plots === \n');
        makeFiguresForPaper = 1 && 1; %~showWorking;
        
        penDistStats_matFile = getFileName('penDistStats');
        
        
        makeDistFigs = showWorking;
        penDistStats_S.cmpType = 'degree';
        penDistStats_S.gratingType = gratingType;
        penDistStats_S.subtractSpont = subtractSpont;
        penDistStats_S.bccType = curBccType;
        penDistStats_S.preserveSC = preserveSC;
        penDistStats_S.preserveAB = preserveAB;
        
        nResamplesTotal = 2000;
%                     nResamplesTotal = 50;
        randInfo_ori.nResamplesTotal = nResamplesTotal;
        randInfo_ori.getBrccPairs = 1;
        
        saveStatsToFile = (nResamplesTotal == 2000);
        randInfo_spf.nResamplesTotal = nResamplesTotal;
        randInfo_spf.getBrccPairs = 1;
        
        %             [Wrcc_idxs, nsimp_wrcc,   Brcc_idxs, nsimp_brcc] = getRandomizedPairs(randInfo, loop_i==1, nResamplesPerLoop);
        
        dF1oDC_ori_cc = abs(diff(pairData_ori.F1oDCs_pref, [], 2));
        dF1oDC_ori_wcc = dF1oDC_ori_cc(Wcc_oo_pairIdxs);
        dF1oDC_ori_bcc = dF1oDC_ori_cc(Bcc_oo_pairIdxs);
        
        dF1oDC_spf_cc = abs(diff(pairData_spf.F1oDCs_pref, [], 2));
        dF1oDC_spf_wcc = dF1oDC_spf_cc(Wcc_ss_pairIdxs);
        dF1oDC_spf_bcc = dF1oDC_spf_cc(Bcc_ss_pairIdxs);
        
        stimTypes = {'ori', 'spf'};
        if makeFiguresForPaper
%                             stimTypes = {'ori'};
            fig_idx = 200;
            subSpcM = [0.01 0.015, 0.01];
            subSpcN = [0.02 0.00, 0.02];
        else
            fig_idx = curGratingType*100;
        end
        %             stimTypes = {'spf', 'ori'};
        oriMeasureTypes      = {'Global ori width', 'Local ori width', 'Pref. orientation', 'F1/DC (ori)'};
        oriMeasureFieldNames = {'w_ori_global', 'w_ori_local', 'ori_pref', 'F1oDC_ori'};
        ori_isDeg =           [1, 1, 1, 0];
        ori_fmt =             {'%.1f', '%.1f', '%.1f', '%.2f'};
        oriMeasures_cc  = {dOriW_global_cc,  dOriW_local_cc,  dOriPref_cc,  dF1oDC_ori_cc};
        oriMeasures_bcc = {dOriW_global_bcc, dOriW_local_bcc, dOriPref_bcc, dF1oDC_ori_bcc};
        oriMeasures_wcc = {dOriW_global_wcc, dOriW_local_wcc, dOriPref_wcc, dF1oDC_ori_wcc};
        if strcmp(gratingType, 'drifting')
            oriMeasureTypes = [oriMeasureTypes, 'DSI'];
            oriMeasureFieldNames = [oriMeasureFieldNames, 'DSI'];
            oriMeasures_cc = [ oriMeasures_cc , dDSI_cc];
            oriMeasures_bcc = [oriMeasures_bcc, dDSI_bcc];
            oriMeasures_wcc = [oriMeasures_wcc, dDSI_wcc];
            ori_isDeg = [ori_isDeg, 0];
            ori_fmt = [ori_fmt, '%.2f'];
        end
        %             'SF tuning width', 'Preferred SF', 'F1/DC (SF)'};
        spfMeasureTypes = {'SF tuning width', 'Preferred SF', 'F1/DC (SF)'};
        spfMeasureFieldNames = {'spf_w', 'spf_pref', 'F1oDC_spf'};
        spfMeasures_cc  = {dSpf_w_cc,  dSpfPref_cc,  dF1oDC_spf_cc};
        spfMeasures_bcc = {dSpf_w_bcc, dSpfPref_bcc, dF1oDC_spf_bcc};
        spfMeasures_wcc = {dSpf_w_wcc, dSpfPref_wcc, dF1oDC_spf_wcc};
        spf_isDeg = [0, 0, 0];
        spf_fmt = {'%.2f', '%.2f', '%.2f'};
        
        allMeasureFieldNames = [oriMeasureFieldNames, spfMeasureFieldNames];
        
        set(0,'DefaultFigureWindowStyle','docked')
        
        for ti = 1:length(stimTypes)
            if strcmp(stimTypes{ti}, 'ori')
                measureTypes = oriMeasureTypes;
                measureFieldNames = oriMeasureFieldNames;
                measures_cc = oriMeasures_cc;
                measures_bcc = oriMeasures_bcc;
                measures_wcc = oriMeasures_wcc;
                pairData = pairData_ori;
                Bcc_pairIdxs = Bcc_oo_pairIdxs;
                Wcc_pairIdxs = Wcc_oo_pairIdxs;
                isDeg = ori_isDeg;
                fmt = ori_fmt;
                doSSpairs = 0;
                [Wrcc_idxs_ori, ~,   Brcc_idxs_ori, ~] = getRandomizedPairs(randInfo_ori, 1, nResamplesTotal);
                [Wrcc_idxs, Brcc_idxs] = deal(Wrcc_idxs_ori, Brcc_idxs_ori);
                
                
            elseif strcmp(stimTypes{ti}, 'spf')
                measureTypes = spfMeasureTypes;
                measureFieldNames = spfMeasureFieldNames;
                measures_cc  = spfMeasures_cc;
                measures_bcc = spfMeasures_bcc;
                measures_wcc = spfMeasures_wcc;
                pairData = pairData_spf;
                Bcc_pairIdxs = Bcc_ss_pairIdxs;
                Wcc_pairIdxs = Wcc_ss_pairIdxs;
                
                doSSpairs = (curGratingType == 1) && 0;
                isDeg = spf_isDeg;
                fmt = spf_fmt;
                [Wrcc_idxs_spf, ~,   Brcc_idxs_spf, ~] = getRandomizedPairs(randInfo_spf, 1, nResamplesTotal);
                [Wrcc_idxs, Brcc_idxs] = deal(Wrcc_idxs_spf, Brcc_idxs_spf);
                
            end
            
            
            fprintf('\n');
            %%
            if doSSpairs
                bcc_SC_type = pairData.SCtype_pref(Bcc_pairIdxs);
                wcc_SC_type = pairData.SCtype_pref(Wcc_pairIdxs);
                idx_use_bcc = find(bcc_SC_type == 2);
                idx_use_wcc = find(wcc_SC_type == 2);
                
                Bcc_pairIdxs = Bcc_pairIdxs(idx_use_bcc);
                Wcc_pairIdxs = Wcc_pairIdxs(idx_use_wcc);
                %                     sc_str = '(SS only)';
                sc_str = '';
            else
                sc_str = '';
            end
            
            micronsPerStep = 0.12406;
            
            sameAnimal_bcc = pairData.animalCmp(Bcc_pairIdxs);
            sameHemisphere_bcc = pairData.hemiCmp(Bcc_pairIdxs);
            samePen_bcc = pairData.penetrCmp(Bcc_pairIdxs);
            sameLoc_bcc = pairData.locCmp(Bcc_pairIdxs);
            penDist_bcc = pairData.penDist(Bcc_pairIdxs);
            depthDist_bcc = abs(pairData.depthDist(Bcc_pairIdxs)) * micronsPerStep;
            
            penDist_bcc_nz = penDist_bcc;
            penDist_bcc_nz(penDist_bcc_nz == 0) = nan;
            [n,binId] = histc(penDist_bcc_nz, [0, .1:1:7]);
            
            depthDist_bcc_nz = depthDist_bcc;
            depthDist_bcc_nz(depthDist_bcc_nz == 0) = nan;
            [n,binId] = histc(depthDist_bcc_nz, [0, .1:1:7]);
            
            
            if strcmp(gratingType, 'drifting')
                nPairsEach = 1000;
            elseif strcmp(gratingType, 'flashed')
                nPairsEach = 500;
                if doSSpairs
                    nPairsEach = 20;
                end
            end
            %%
            
            
            
            if makeFiguresForPaper && makeDistFigs
                %%
                title_fsize = 11; % 9
                xlabel_fsize = 12; %9;
                ylabel_fsize = 12; %9;
                legend_fsize = 7.5;
                
                randHiLo_color = .3*[1, 1, 1];
                sliding_color = 0*[1, 1, 1];
                sliding_color_err = .6*[1, 1, 1];
                
                showDistributions = 0 && 1;
                
                showWP = 1 && 1;
                showWA = 0 && 1;
                
                
                sub_col_idx = switchh(gratingType, {'drifting', 'flashed'}, [1, 2]);
                %                     column_measures = {'Pref. orientation', 'Global ori width'};
                
                if strcmp(stimTypes{ti}, 'ori')
                    column_measures = {'Global ori width', 'Local ori width', 'Pref. orientation', 'F1/DC (ori)', 'DSI'};
                    cur_fig_idx = fig_idx;
                    let_offset = 0;
                else
                    column_measures = {'SF tuning width', 'Preferred SF', 'F1/DC (SF)'};
                    figure(fig_idx+1);
                    cur_fig_idx = fig_idx+1;
                    %                         let_offset = 7;
                    let_offset = 0;
                    
                end
                %                     xlim_mode = 'auto';
                xlim_mode = [0, 1600];
                figure(cur_fig_idx);
                if strcmp(gratingType, 'drifting')
                    clf;
                end
                set(cur_fig_idx, 'windowStyle', 'normal');
                fig_LB = [150, 120]; fig_height_row = 350; fig_width_row = 730;
                nRows = length(column_measures);
                set(cur_fig_idx, 'position', [fig_LB,  fig_width_row, fig_height_row*nRows], 'color', 'w')
                drawnow;
                
            else
                sliding_color = 'r';
                sliding_color_err = [.7 0 0];
                randHiLo_color = [0 0 .7];
                
                showDistributions = 0 && 1;
                showWP = 1 && 1;
                showWA = 0 && 1;
                xlim_mode = 'auto';
                
            end
            
            
            nSubplots = showDistributions + showWP + showWA;
            
            for mi = 1:length(measureTypes);
                
                if makeFiguresForPaper && makeDistFigs
                    sub_row_idx = find(strcmp(measureTypes{mi}, column_measures), 1);
                    if isempty(sub_row_idx);
                        continue;
                    end
                end
                
                
                linew = 2;
                dX_wcc = measures_wcc{mi};
                dX_bcc = measures_bcc{mi};
                if doSSpairs
                    dX_wcc = dX_wcc(idx_use_wcc);
                    dX_bcc = dX_bcc(idx_use_bcc);
                end
                
                
                %                     dX_wrcc = cellfun(@(idx) measures_cc{mi}(idx), Wrcc_idxs, 'un', 0);
                %                     dX_brcc = cellfun(@(idx) measures_cc{mi}(idx), Brcc_idxs, 'un', 0);
                
                X_name = measureTypes{mi};
                if ~makeFiguresForPaper && makeDistFigs
                    fig_idx = fig_idx+1;
                    figure(fig_idx); clf;
                    
                end
                
                idx_WA = sameAnimal_bcc==1;
                idx_WP = samePen_bcc==1;
                idx_BP_WA = samePen_bcc==0 & sameAnimal_bcc==1;
                idx_BP_BH_WA = samePen_bcc==0 & sameHemisphere_bcc == 0 & sameAnimal_bcc==1;
                idx_BP = samePen_bcc==0;
                idx_BA = sameAnimal_bcc==0;
                idx_BP_WA_close = idx_BP_WA & penDist_bcc_nz < .5;
                idx_BP_WA_far = idx_BP_WA & penDist_bcc_nz > 2;
                
                %%% look at all the distributions:
                Dists_C     = {nonnans(dX_bcc),             nonnans(dX_bcc(idx_BP_WA)),  nonnans(dX_bcc(idx_WP)),   nonnans(dX_wcc)};
                dists_names = {'BS             ', 'BS-BP-WA ',        'BS-WP      ',   'WS            '};
                i_bs = 1; i_bs_bp_wa = 2; i_bs_wp = 3; i_ws = 4;
                medians = cellfun(@nanmedian, Dists_C);
                bs_median = nanmedian(dX_bcc);
                bs_wp_median = nanmedian(dX_bcc(idx_WP));
                bs_bp_wa_median = nanmedian(dX_bcc(idx_BP_WA));
                ws_median = nanmedian(dX_wcc);
                
                N = cellfun(@length, Dists_C);
                nBS = N(1); nWABP= N(2);  nWP = N(3); nWA = nWABP+nWP;
%                 fprintf('N (BS) = %d. N(WA) = %d. N(WA-BP) = %d. R(BS/WA) = %.1f. R(BS/WA-BP) = %.1f)\n', nBS, nWA, nWABP, nBS/nWA, nBS/nWABP);
                deg_str = iff(isDeg(mi), '�', '');
                median_strs = arrayfun(@(m) sprintf([fmt{mi} deg_str ], m), medians, 'un', 0);
                leg_strs = arrayfun(@(nm, n, s) sprintf(['%s (N = %d, median = %s) .'], nm{1}, n, s{1}), dists_names, N, median_strs, 'un', 0);
                
                
                if showDistributions && makeDistFigs
                    subplotGap(1,nSubplots,1);
                    %                     Dists_C = {dX_bcc, dX_bcc(idx_WA), dX_bcc(idx_WP), dX_bcc(idx_BP_WA), dX_wcc};
                    %                     dists_names = {'BS              ',  'BS-WA       ', 'BS-WP       ', 'BS-BP-WA ', 'WS              '};
                    %                     Dists_C = {dX_bcc, dX_bcc(idx_BP), dX_bcc(idx_WA), dX_bcc(idx_WP), dX_bcc(idx_BP_WA), dX_wcc};
                    %                     dists_names = {'BS              ',  'BS-BP       ', 'BS-WA       ', 'BS-WP       ', 'BS-BP-WA ', 'WS              '};
                    
                    
                    h_hists = hist2(Dists_C, 20, 'line', 'norm');
                    mks = 'sovs.';
                    for i = 1:length(h_hists)
                        set(h_hists(i), 'marker', mks(i), 'linewidth', linew)
                    end
                    legend(leg_strs, 'location', 'NE', 'fontsize', 8);
                    
                    c_order = get(gca, 'colorOrder');
                    for i = 1:length(N)
                        drawVerticalLine(medians(i), 'color', c_order(i,:), 'linestyle', ':', 'linewidth', 2);
                    end
                    xlabel(sprintf('Difference in %s', X_name), 'fontsize', 11); ylabel('Normalized count', 'fontsize', 11);
                    title(sprintf('%s (%s  gratings) %s', X_name, gratingType, sc_str), 'fontsize', 11);
                end
                3;
                
                %                     allLocData = [allOriUnits.locData];
                %                     allLocIds = [allLocData.LocId];
                
                %                     allLocIds_nan = [allLocIds(Bcc_pairs_ori(idx_nanDist,1))];
                
                %%
                %
                %                     dMax = prctile(depthDist_bcc, 95);
                %                     nSteps = 100;
                % %                     depthDist_bcc_sorted = sort(nonnans(depthDist_bcc));
                %                     nSize = 300;
                %                     dLims = linspace(0, dMax, nSteps);
                %                     depth_medians = zeros(1, nSteps);
                %                     for i = 1:nSteps
                %                         depth_medians(i) = nanmedian(  dX_bcc(samePen_bcc==1 & depthDist_bcc > dLims(i) ) );
                %
                %                     end
                
                
                
                
                %% look at within-penetration distributions
                
                
                %                     figure(402);
                %                     plot(dLims, depth_medians, 'o-');
                %                     opt.nPairsEach = nPairsEach;
                %                     nPairsEach = 250;
                opt.nPairsEach = nPairsEach;
                opt.overlap_f = 5;
                opt.bs_median = bs_median;
                opt.sigTestTail = 'both';
                %                     opt.sigTestTail = 'left';
                
                
                %                     dX_wrcc = cellfun(@(idx) measures_cc{mi}(idx), Wrcc_idxs, 'un', 0);
                
                %                     opt.allDistances_sorted = sort(nonnans(depthDist_bcc));
                idx_WP = find(samePen_bcc==1);
                depthDist_bcc_WP = depthDist_bcc(idx_WP);
                dX_bcc_WP = dX_bcc(idx_WP);
                dX_brcc_WP = cellfun(@(idx) measures_cc{mi}(idx(idx_WP)), Brcc_idxs, 'un', 0);
                doCheck = 0;
                if doCheck
                    %%
                    dX_brcc = cellfun(@(idx) measures_cc{mi}(idx), Brcc_idxs, 'un', 0);
                    dX_brcc_WP2  = cellfun(@(dX) dX(idx_WP), dX_brcc, 'un', 0);
                    assert(isequaln(dX_brcc_WP, dX_brcc_WP2));
                end
                
                
                
                idx_BP_WA = find(samePen_bcc==0 & sameAnimal_bcc == 1);
                penDist_bcc_BP_WA = penDist_bcc_nz(idx_BP_WA);
                dX_bcc_BP_WA = dX_bcc(idx_BP_WA);
                dX_brcc_BP_WA = cellfun(@(idx) measures_cc{mi}(idx(idx_BP_WA)), Brcc_idxs, 'un', 0);
                
                if doCheck
                    dX_brcc = cellfun(@(idx) measures_cc{mi}(idx), Brcc_idxs, 'un', 0);
                    dX_brcc_BP_WA2 = cellfun(@(dX) dX(idx_BP_WA), dX_brcc, 'un', 0);
                    assert(isequaln(dX_brcc_BP_WA, dX_brcc_BP_WA2));
                end
                %%
                [depth_medians, depth_rand_lo, depth_rand_hi] = deal([]);
                %                     if showWP
                opt.overlap_f = 3;
                [xs_depth, x_depth_lo, x_depth_hi, depth_medians, depth_median_lo, depth_median_hi, depth_pvalues, depth_rand_med, depth_rand_lo, depth_rand_hi] = ...
                    getSlidingWindowMediansVsDist(depthDist_bcc_WP, dX_bcc_WP, dX_brcc_WP, opt);
                %                     end
                %%
                [horiz_medians, horiz_rand_lo, horiz_rand_hi] = deal([]);
                if showWA
                    opt.overlap_f = 3.5;
                    [xs_horiz, x_horiz_lo, x_horiz_hi, horiz_medians, horiz_medians_lo, horiz_medians_hi, horiz_pvalues, horiz_rand_med, horiz_rand_lo, horiz_rand_hi] = ...
                        getSlidingWindowMediansVsDist(penDist_bcc_BP_WA, dX_bcc_BP_WA, dX_brcc_BP_WA, opt);
                end
                
                %%
                %                     clf;
                pval_th = 0.01;
                if strcmp(X_name, 'Global ori width') && strcmp(gratingType, 'drifting')
                    pval_th = 0.015;
                elseif strcmp(X_name, 'F1/DC (ori)') && strcmp(gratingType, 'flashed')
                    pval_th = 0.019;
                end
                
                idx_notSig = find(depth_pvalues > pval_th, 1);
                if idx_notSig > 1
                    %                         drawHorizontalLine(depth_medians_sm(idx_notSig), 'linestyle', ':', 'color', grey_col, 'linewidth', linew);
                    x_notSig = (xs_depth(idx_notSig) + xs_depth(idx_notSig-1))/2;
                    x_notSig = roundToNearest(x_notSig, 5);
                    x_notSig_str = num2str(x_notSig);
                    %                             x_notSig_str = sprintf('%.0f', x_notSig);
                else
                    x_notSig = nan;
                    x_notSig_str = ['<' num2str(roundToNearest( xs_depth(1), 5))];
                end
                
                fprintf('(%s gratings) : %s : "%s" [%.1f]\n', gratingType, X_name, x_notSig_str, x_notSig);
                penDistStats_S.(measureFieldNames{mi}).dist = x_notSig_str;
                
                if showWP && makeDistFigs
                    
                    %%% Plot
                    %                     clf;
                    %                     figure(fig_idx + 100); clf;
                    if makeFiguresForPaper
                        %                             subplotGap(2,2, sub_row_idx, sub_col_idx, subSpcM, subSpcN);
                        [subSpcM_use, subSpcN_use, nCols_use] = deal(subSpcM, subSpcN, 2);
                        if strcmp(X_name, 'DSI')
                            subSpcN_use([1, 3]) = .25;
                            nCols_use = 1;
                        end
                        subplotGap(nRows,nCols_use, sub_row_idx, sub_col_idx, subSpcM_use, subSpcN_use);
                        sub_idx = sub2ind([2, nRows], sub_col_idx, sub_row_idx);
                        addSubplotLetter(nRows, nCols_use, sub_row_idx, sub_col_idx, subSpcM_use, subSpcN_use, char('A'+sub_idx-1+let_offset), [], 'fontname', subplotLetFont);
                    else
                        subplotGap(1,nSubplots, showDistributions+1);
                    end
                    %                     h = ploterr(xs_depth, depth_medians, {x_depth_lo, x_depth_hi}, {depth_median_lo, depth_median_hi}, 'hhxy', .7); hold on;
                    
                    if ~makeFiguresForPaper
                        %                             h_rand = plot(xs_depth, depth_rand_med, 'linewidth', 1, 'color', 'b', 'marker', 's', 'markersize', 6); hold on;
                    end
                    line_w_hiLo = iff(makeFiguresForPaper, 1, 2);
                    h_rand_hiLo = plot(xs_depth, [depth_rand_lo', depth_rand_hi'], 'linewidth', line_w_hiLo, 'color', randHiLo_color, 'linestyle', '-'); hold on;
                    
                    %                     h_rand = ploterr(xs_depth, depth_rand_med, [], {depth_rand_lo, depth_rand_hi}, 'hhxy', .7); hold on;
                    %                     set(h_rand, 'linewidth', 2, 'color', 'b') ;
                    %                     set(h_rand(1), 'marker', 's', 'markersize', 6)
                    %                     set(h_rand(2), 'color', [0 0 .7]) ;
                    
                    h_sliding = ploterr(xs_depth, depth_medians, {x_depth_lo, x_depth_hi}, [], 'hhxy', 1); hold on;
                    set(h_sliding, 'linewidth', 1.5, 'color', sliding_color) ;
                    set(h_sliding(1), 'marker', 'o', 'markersize', 6)
                    set(h_sliding(2), 'color', sliding_color_err) ;
                    
                    %                     h_hiLo = plot(xs_depth, [depth_median_lo', depth_median_hi'], 'color', [.7 0 0], 'linestyle', ':', 'linewidth', 2);
                    
                    
                    
                    %                     , 'o-', 'linewidth', linew, 'color', 'r'); hold on;
                    
                    ylims1 = lims([depth_medians, medians, depth_rand_lo, depth_rand_hi], [.5, .3]);
                    ylims2 = lims([horiz_medians, medians, horiz_rand_lo, horiz_rand_hi], .05);
                    ylims = lims([ylims1, ylims2]);
                    ylim(ylims);
                    %                     x_max = maxElements([x_depth_lo, x_depth_hi], 5);
                    if ischar(xlim_mode) && strcmp(xlim_mode, 'auto')
                        xlims = lims(xs_depth, .1);
                        xlim([0, xlims(2)]);
                    else
                        xlim(xlim_mode);
                    end
                    %                     error
                    %                     h_ax(1) = gca;
                    yval_stars = ylims(1) + diff(ylims)*.9;
                    showSmoothed = 0;
                    grey_col = .5*[1, 1, 1];
                    if showSmoothed
                        
                        depth_medians_sm = gaussSmooth(depth_medians, 3);
                        height_max = min(max(depth_medians_sm), bs_median);
                        h = (height_max - ws_median)/2 + ws_median;
                        plot(xs_depth, depth_medians_sm, '.:', 'linewidth', linew, 'color', grey_col)
                        idx_half = find(depth_medians_sm > h, 1);
                        drawHorizontalLine(depth_medians_sm(idx_half), 'linestyle', ':', 'color', grey_col, 'linewidth', linew);
                        x_half = xs_depth(idx_half);
                        h_half = drawVerticalLine(x_half, 'linestyle', ':', 'color', grey_col, 'linewidth', linew);
                    end
                    
                    if idx_notSig > 1
                        %                         drawHorizontalLine(depth_medians_sm(idx_notSig), 'linestyle', ':', 'color', grey_col, 'linewidth', linew);
                        h_notSig = drawVerticalLine(x_notSig, 'linestyle', ':', 'color', grey_col, 'linewidth', linew);
                    end
                    
                    
                    if  ~makeFiguresForPaper
                        h_B = drawHorizontalLine(bs_median, 'linestyle', '--', 'color', 'b', 'linewidth', linew);
                        h_BW = drawHorizontalLine(bs_wp_median, 'linestyle', '--', 'color', 'r', 'linewidth', linew);
                        h_W = drawHorizontalLine(ws_median, 'linestyle', '--', 'color', [0, .75, .75], 'linewidth', linew);
                        h_lines_legend = [h_B, h_BW, h_W, h_sliding(1), h_rand_hiLo(1)];
                        legend_str_lines = {sprintf('BS median (%s)', median_strs{i_bs}), sprintf('BS-WP median (%s)', median_strs{i_bs_wp}), sprintf('WS median (%s)', median_strs{i_ws}), 'BS-WP(window)', 'BS-WP(window)-rand'};
                    else
                        showBS = 1;
                        h_W = drawHorizontalLine(ws_median, 'linestyle', '--', 'color', 'k', 'linewidth', linew);
                        h_lines_legend = [h_W, h_sliding(1), h_rand_hiLo(1)];
                        legend_str_lines = {'WS median', 'BS-WP (observed)', 'BS-WP (control)'};
                        if showBS
                            h_B = drawHorizontalLine(bs_median, 'linestyle', ':', 'color', .3*[1, 1, 1], 'linewidth', linew);
                            h_lines_legend = [h_B, h_lines_legend];
                            legend_str_lines = ['BS median', legend_str_lines];
                        end
                        
                    end
                    
                    show05 = 0;
                    show01 = 1;
                    alsoShow001 = 0;
                    if alsoShow001
                        idx_lt_001 = find( depth_pvalues < .001 );
                        idx_lt_01 = find( depth_pvalues >= .001 & depth_pvalues < .01 );
                    else
                        idx_lt_01 = find( depth_pvalues < pval_th );
                    end
                    idx_lt_05 = find( depth_pvalues >=  .01 & depth_pvalues < .05 );
                    
                    %                         idx_lt_001 = find( depth_pvalues < .001 );
                    %                         idx_lt_01 = find( depth_pvalues >= .001 & depth_pvalues < .01 );
                    %                         if makeFiguresForPaper
                    yvals_stars = @(idx) yval_stars * ones(size(idx));
                    %                         else
                    %                             yvals_stars = @(idx) depth_medians(idx);
                    %                         end
                    
                    
                    if alsoShow001 && ~isempty(idx_lt_001), h_sig_depth(1) = plot(xs_depth(idx_lt_001), yvals_stars(idx_lt_001), 'ks', 'markersize', 6, 'linewidth', 1); end
                    if show01      && ~isempty(idx_lt_01),  h_sig_depth(2) = plot(xs_depth(idx_lt_01), yvals_stars(idx_lt_01), 'k*', 'markersize', 6, 'linewidth', 1); end
                    if show05      && ~isempty(idx_lt_05),  h_sig_depth(3) = plot(xs_depth(idx_lt_05), yvals_stars(idx_lt_05), 'ko', 'markersize', 6, 'linewidth', 1); end
                    3;
                    if ~makeFiguresForPaper || ((cur_fig_idx == fig_idx) && sub_row_idx == 1 && sub_col_idx == 2) || ((cur_fig_idx == fig_idx+1) && sub_row_idx == 1 && sub_col_idx == 1)
                        legend(h_lines_legend, legend_str_lines, 'location', 'SE', 'fontsize', 8);
                    end
                    
                    %                         nPairs_str = iff(~makeFiguresForPaper, sprintf(' of %d pairs in sliding window', nPairsEach), '');
                    nPairs_str = '';
                    xlabel(sprintf('Median distance between sites (\\mum)%s', nPairs_str), 'fontsize', xlabel_fsize);
                    X_name_lower = X_name;
                    if ~any(strncmp(X_name, {'DSI', 'SF', 'F1'}, 2))
                        X_name_lower(1) = lower(X_name(1));
                    end
                    ylabel(sprintf('Median diff in %s', X_name_lower), 'fontsize', ylabel_fsize);
                    if ~makeFiguresForPaper
                        %                             title(sprintf('%s (%s gratings) : BS-WP distribution. [h = %.0f   \\mum] %s', X_name, gratingType, x_notSig, sc_str), 'fontsize', 9);
                    else
                        %                             title({X_name, sprintf('(%s gratings)', gratingType)}, 'fontsize', title_fsize);
                        %                             title(sprintf('%s (%s gratings)', X_name, gratingType), 'fontsize', title_fsize);
                    end
                    title(sprintf('%s (%s gratings)', X_name, gratingType), 'fontsize', title_fsize);
                    
                    
                end
                3;
                %%
                %                     figure(fig_idx + 200); clf;
                
                if showWA && makeDistFigs
                    subplot(1,nSubplots,showDistributions+showWP + 1);
                    
                    h_rand = plot(xs_horiz, horiz_rand_med, 'linewidth', 2, 'color', 'b', 'marker', 's', 'markersize', 6); hold on;
                    h_rand_hiLo = plot(xs_horiz, [horiz_rand_lo', horiz_rand_hi'], 'linewidth', 2, 'color', [0 0 .7], 'linestyle', ':');
                    
                    
                    %                      h_rand = ploterr(xs_horiz, horiz_rand_med, [], {horiz_rand_lo, horiz_rand_hi}, 'hhxy', .7); hold on;
                    %                      h_rand = plot(xs_horiz, horiz_rand_med, [], {horiz_rand_lo, horiz_rand_hi}, 'hhxy', .7); hold on;
                    %
                    %                     set(h_rand, 'linewidth', 2, 'color', 'b') ;
                    %                     set(h_rand(1), 'marker', 's', 'markersize', 6)
                    %                     set(h_rand(2), 'color', [0 0 .7]) ;
                    
                    %                     plot(xs_horiz, horiz_medians, 'o-', 'linewidth', linew, 'color', [0 .75 0]); hold on;
                    h_sliding = ploterr(xs_horiz, horiz_medians, {x_horiz_lo, x_horiz_hi}, [], 'hhxy', 1.2); hold on;
                    %                     h = ploterr(xs_horiz, horiz_medians, {x_horiz_lo, x_horiz_hi}, {horiz_medians_lo, horiz_medians_hi}, 'hhxy', 1.2); hold on;
                    set(h_sliding, 'linewidth', 2, 'color', [0 .75 0]) ;
                    set(h_sliding(1), 'marker', 'o', 'markersize', 6)
                    set(h_sliding(2), 'color', [0 .5 0]) ;
                    
                    %                     h_hiLo = plot(xs_horiz, [horiz_medians_lo', horiz_medians_hi'], 'color', [0 .5 0], 'linestyle', ':', 'linewidth', 2);
                    
                    
                    
                    
                    %                     h_ax(2) = gca;
                    ylim(ylims);
                    %                     set(h_ax, 'ylim', lims([ylims1, ylims2]));
                    
                    xlims = lims(xs_horiz, .1);
                    xlim(xlims);
                    
                    %                     xlims = lims(xs_horiz, .1);
                    %                     xlim(xlims);
                    
                    doHalfMax = 0;
                    if doHalfMax
                        horiz_medians_sm = gaussSmooth(horiz_medians, 10);
                        height_max = min(max(horiz_medians_sm), bs_median);
                        h = (height_max - ws_median)/2 + ws_median;
                        plot(xs_horiz, horiz_medians_sm, '.:', 'color', grey_col)
                        idx_half = find(horiz_medians_sm > h, 1);
                        
                        drawHorizontalLine(horiz_medians_sm(idx_half), 'linestyle', ':', 'color', grey_col);
                        x_half = xs_horiz(idx_half);
                        drawVerticalLine(x_half, 'linestyle', ':', 'color', grey_col);
                        x_half_str = sprintf('[h = %.0f      \\mum]')
                    else
                        x_half_str = '';
                    end
                    
                    
                    h_B = drawHorizontalLine(bs_median, 'linestyle', '--', 'color', 'b', 'linewidth', linew);
                    h_BP_WA = drawHorizontalLine(bs_bp_wa_median, 'linestyle', '--', 'color', [0 .75 0], 'linewidth', linew);
                    h_W = drawHorizontalLine(ws_median, 'linestyle', '--', 'color', [0, .75, .75], 'linewidth', linew);
                    
                    idx_lt_001 = find( horiz_pvalues < .001 );
                    idx_lt_01 = find( horiz_pvalues >= .001 & horiz_pvalues < .01 );
                    idx_lt_05 = find( horiz_pvalues >=  .01 & horiz_pvalues < .05 );
                    
                    if ~isempty(idx_lt_001), h_sig_horiz(1) = plot(xs_horiz(idx_lt_001), horiz_medians(idx_lt_001), 'k*', 'markersize', 8, 'linewidth', 2); end
                    if ~isempty(idx_lt_01), h_sig_horiz(2) = plot(xs_horiz(idx_lt_01), horiz_medians(idx_lt_01), 'ks', 'markersize', 8, 'linewidth', 2); end
                    if ~isempty(idx_lt_05), h_sig_horiz(3) = plot(xs_horiz(idx_lt_05), horiz_medians(idx_lt_05), 'ko', 'markersize', 8, 'linewidth', 2); end
                    
                    legend([h_B, h_BP_WA, h_W, h_sliding(1), h_rand], {sprintf('BS median (%s)', median_strs{i_bs} ), ...
                        sprintf('BS-BP-WA median (%s)', median_strs{i_bs_bp_wa} ), ...
                        sprintf('WS median (%s)', median_strs{i_ws} ), ...
                        'BS-BP-WA(window)', 'BS-BP-WA(window)-rand'}, 'location', 'SE', 'fontsize', 8);
                    
                    xlabel(sprintf('Median AP-ML distance (mm) of %d pairs in sliding window', nPairsEach), 'fontsize', 9);
                    ylabel(sprintf('Median diff in %s', X_name));
                    title(sprintf('%s (%s gratings) : BS-BP-WA distribution. %s %s', X_name, gratingType, x_half_str, sc_str), 'fontsize', 9);
                    
                end
                3;
                
            end
            
        end
        
        
        penDistStats_S.columns = {'dist'};
        penDistStats_S.allMeasureNames = allMeasureFieldNames;
        penDistStats_S.nResamples = nResamplesTotal;
        if saveStatsToFile
            save(penDistStats_matFile, '-struct', 'penDistStats_S');
        end
        if strcmp(gratingType, 'flashed') && makeDistFigs && saveStatsToFile
            keyboard;
        end
        
        if 0
            %%
            fig_file_name_a = sprintf('%sFigure13_clusteringExtent_ori.pdf', figureFolder);
            export_fig(fig_idx, fig_format, fig_file_name_a);
            
%             fig_file_name_b = sprintf('%sFigure14_clusteringExtent_sf.pdf', figureFolder);
%             export_fig(fig_idx+1, fig_format, fig_file_name_b);
            
            
        end
        
        
        
        %%
        
        %%
        
        
        3;
        
    end

    
    if doPairStats
        pairDiffs_S.miscStats = miscPairStats;
        save(pairDiffs_matFile, '-struct', 'pairDiffs_S');
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



function [allOriUnits, pairData_ori, S_ori, ...
    pairTypes, measures_ori, pairIdxs_ori, pairIdxList_ori, idxMtx_ori, ...
    Wcc_oo_pairIdxs, Wcc_oo_pairIdxs_M, Wcc_pairs_ori, ...
    Wrcc_oo_pairIdxs, Wrcc_oo_pairIdxs_M, Wrcc_pairs_ori, ...
    Bcc_oo_pairIdxs, Bcc_oo_pairIdxs_M, Bcc_pairs_ori, ...    
    oriCells_use, allOriUnitStats, Wrcc_details, ...
    idx_cellsAtEachSite_full, idx_pairsAtEachSite_full, miscStats] ...
        = loadDegreeOfTuningData(osType, opt, miscStats) %limitToNPermutes, showWorking)

    file_ext = ['_' osType];

    % load ORI data
    osps_file = getFileName('osps', file_ext);
    pairs_file = getFileName('pairs', file_ext);
    cmps_file = getFileName('comparisons', file_ext);
    
    if 0
        file_ext = '_ori';
        osps_file = getFileName('osps', file_ext)
        pairs_file = getFileName('pairs', file_ext)
        cmps_file = getFileName('comparisons', file_ext)
    end
    
    if opt.showWorking,   tic; fprintf('Loading OSP file %s ...', extractFileName(osps_file));  end
    S1_ori = load(osps_file); 
    if opt.showWorking,   toc; end
        
    allOriUnits = S1_ori.allCells;        
    nUnits_ori = length(allOriUnits);    

    if opt.showWorking,   tic; fprintf('Loading pairs file %s ...', extractFileName(pairs_file)); end
    S2_ori = load(pairs_file); 
    if opt.showWorking,   toc; end
    
    if opt.showWorking,   tic; fprintf('Loading comparisons file %s ...', extractFileName(cmps_file)); end
    S3_ori = load(cmps_file); 
    if opt.showWorking,   toc; end
    
    [pairData_ori, S_ori, pairTypes, measures_ori] = deal(S3_ori.pairData, S3_ori.allStatsC, S3_ori.pairTypes, S3_ori.measures);    
    [pairIdxs_ori, pairIdxList_ori, idxMtx_ori] = useSubsetOfPairIdxs(S2_ori, pairTypes, nUnits_ori);
    Wrcc_details = S2_ori.Wrcc_details;

    
    removeOddEvenStats = 1;
    if removeOddEvenStats
        for i = 1:length(allOriUnits)
            ts = allOriUnits(i).tuningStats;
            fn = fieldnames(ts);
            fn_remove_idx = find(cellfun(@(s) ~isempty(strfind(s, 'odd')) || ~isempty(strfind(s, 'even')), fn));
            for j = fn_remove_idx(:)'
                ts = rmfield(ts, fn{j});
            end
            allOriUnits(i).tuningStats = ts;
        end    
        
    end

    suffix_str = iff(curSubtractSpont, 'ss', 'si');    
    switch osType
        case 'ori', statsField = ['oriStats_' suffix_str];  stimType_fld = 'Grating:Orientation';
        case 'spf', statsField = ['spfStats_' suffix_str];  stimType_fld = 'Grating:Spatial Freq';            
    end
    3;
    
    %%
    % Choose cells (cellID > 0) that are OK
    oriIsCell = [allOriUnits.cellId] > 0;
    
    allID = nestedFields(allOriUnits, 'spkFeatures', 'IsolationDistance');
    unitWellIsolated  = allID > opt.criteria.minID;
%     unitWellIsolated(oriIsCell)
    
%     unitPassedCriterion = nestedFields(allOriUnits, 'tuningStats', statsField, 'cellOK');
    unitPassedCriterion = arrayfun(@(s) s.tuningStats.(statsField).cellOK, allOriUnits(:)');

        %%
    
%     allID_cell = allID(oriIsCell);
    oriCells_use = oriIsCell & unitWellIsolated & unitPassedCriterion;
    oriCells_use_idx = find(oriCells_use);
    allCellsUsed = all(oriCells_use(oriIsCell));
   
%     cellsUsable_frac = [nnz(idx_cellOK) / nnz(oriIsCell)];

    %  reindexing
    Wcc_oo_pairIdxs_M = pairIdxs_ori{ find(strcmp(pairTypes, 'Wcc'), 1) };
    Wcc_oo_pairIdxs = idxMtx_ori( Wcc_oo_pairIdxs_M );
    Wcc_pairs_ori = ind2subV([nUnits_ori, nUnits_ori], Wcc_oo_pairIdxs_M);
    if ~allCellsUsed
        wcc_pairs_use = all(binarySearch(oriCells_use_idx, Wcc_pairs_ori, 1, 0), 2); % both cells in a pair are 'good' cells.
        Wcc_pairs_ori = Wcc_pairs_ori(wcc_pairs_use,:);
        Wcc_oo_pairIdxs = Wcc_oo_pairIdxs(wcc_pairs_use);
        Wcc_oo_pairIdxs_M = Wcc_oo_pairIdxs_M(wcc_pairs_use);
    end
    
    if any(strcmp(pairTypes, 'Wrcc'))
        if ~isempty(opt.limitToNPermutes)
            i_wrcc = find(strcmp(pairTypes, 'Wrcc'), 1);
            pairIdxs_ori{ i_wrcc } = pairIdxs_ori { i_wrcc }(1:opt.limitToNPermutes);
        end

        Wrcc_oo_pairIdxs_M = pairIdxs_ori{ find(strcmp(pairTypes, 'Wrcc'), 1) };        
        Wrcc_oo_pairIdxs =  cellfun(@(idxs_M) idxMtx_ori( idxs_M ), Wrcc_oo_pairIdxs_M, 'un', 0);                
        Wrcc_pairs_ori = cellfun(@(idxs) ind2subV([nUnits_ori, nUnits_ori], idxs), Wrcc_oo_pairIdxs_M, 'un', 0);     
        %%
        if ~allCellsUsed
            for j = 1:length(Wrcc_pairs_ori)
                wrcc_pairs_use_j = all(binarySearch(oriCells_use_idx, Wrcc_pairs_ori{j}, 1, 0), 2); % both cells in a pair are 'good' cells.
                Wrcc_pairs_ori{j} = Wrcc_pairs_ori{j}(wrcc_pairs_use_j,:);
                Wrcc_oo_pairIdxs{j} = Wrcc_oo_pairIdxs{j}(wrcc_pairs_use_j);
                Wrcc_oo_pairIdxs_M{j} = Wrcc_oo_pairIdxs_M{j}(wrcc_pairs_use_j);
            end
        end
        
    end
    Bcc_oo_pairIdxs_M  = pairIdxs_ori{ find(strcmp(pairTypes, 'Bcc'), 1) };
    Bcc_oo_pairIdxs = idxMtx_ori( Bcc_oo_pairIdxs_M );
    Bcc_pairs_ori = ind2subV([nUnits_ori, nUnits_ori], Bcc_oo_pairIdxs_M);
    if ~allCellsUsed
        bcc_pairs_use = all(binarySearch(oriCells_use_idx, Bcc_pairs_ori, 1, 0), 2); % both cells in a pair are 'good' cells.
        Bcc_pairs_ori = Bcc_pairs_ori(bcc_pairs_use,:);
        Bcc_oo_pairIdxs = Bcc_oo_pairIdxs(bcc_pairs_use);
        Bcc_oo_pairIdxs_M = Bcc_oo_pairIdxs_M(bcc_pairs_use);
    end

%%
    % Choose cells (cellID > 0) that are OK

%     switch curGratingType(''), 
%         case 'flashed',
%         case 'drifting',
%             oriCells_use = oriIsCell & cellIsOK; % strncmp({allOriUnits.stimType}, stimType_fld, 19);            
%     end
    
%     allOriUnitTuningStats = [allOriUnits.tuningStats];    
    
    allOriUnitStats = nestedFields(allOriUnits, {'tuningStats', statsField}, 1);
    for i = 1:length(allOriUnitStats)
        if ~isfield(allOriUnitStats(i), 'error_jack') || isequaln(allOriUnitStats(i).error_jack, nan)
            allOriUnitStats(i).error_jack = struct;
        end
    end
       
    assert(all([allOriUnitStats(oriCells_use).cellOK]));        

        
    printNCellsPerSiteStats = 1;
    %%
    if printNCellsPerSiteStats       
        
        
        %%
%         GC = [allOriUnits.Gid]*1000 + [allOriUnits.cellId];
%         idx_cells = [allOriUnits.cellId] > 0;
        allCellGids = [allOriUnits(oriCells_use).Gid];
%         allCellGids2 = [allOriUnits(idx_cells).Gid];
        
%         GC_use = GC(oriCells_use);
%         GC_cell = GC(idx_cells);
%         setdiff(GC_cell, GC_use);
        %%
        
        [uGid, gidCount] = uniqueCount(allCellGids);
        nSites = length(uGid);
        nCells = length(allCellGids);
        nCells_per_site = mean(gidCount); assert(nCells_per_site == nCells/nSites);
        nSites_1cell = nnz(gidCount == 1);
        nPairs = length(S2_ori.Wcc_idxs);
        %%
        cell_s{1} = sprintf('**** %s (%s gratings) ***', osType, curGratingType(''));
        cell_s{2} = sprintf('Recorded from %d cells from %d sites (mean %.2f cells per site).', nCells, nSites, nCells_per_site);
        cell_s{3} = sprintf('%d sites had 1 cell. Thus %d cells yielded a total of %d pairs',nSites_1cell, nCells-nSites_1cell, nPairs  );
%         cell_s{4} = sprintf('%d / %d cells (%.1f%%) had an isolation distance of at least %d.', nnz(allID_cell > 10), length(allID_cell), nnz(allID_cell > 10)/length(allID_cell)*100, opt.criteria.minID );

        cellStats = struct('numSites', nSites, 'nCellsTot', nCells, 'nSites_1Cell', nSites_1cell, 'nPairs', nPairs);
        miscStats.([osType 'CellsStats']) = cellStats;
        
        if opt.showWorking
            cellfun(@(str) fprintf('%s\n', str), cell_s);
        end
%         comments = cell_s;
            
    end    
        
%             subtractSpont = any(cellfun(@(s) ~isempty(strfind(s, '_ss')), measures)) && ...
%                 isfield(allOriCells(1).stats.tuningStats, 'oriStats_ss');
% 
%             if subtractSpont
%                 allOriStats_ss = nestedFields(allOriCells, {'stats', 'tuningStats', 'oriStats_ss'});
%             end            
        
        

    % Double check that the renumbering is all ok.
    checkRenumbering = 1;
    if checkRenumbering
        allOriUnitGids = [allOriUnits.Gid];
        assert( all(allOriUnitGids(Wcc_pairs_ori(:,1)) == allOriUnitGids(Wcc_pairs_ori(:,2))));
        assert( all(allOriUnitGids(Bcc_pairs_ori(:,1)) ~= allOriUnitGids(Bcc_pairs_ori(:,2))));
    end    

    
    if any(strcmp(opt.bccType, {'Pen', 'Animal'}))
        %%                         
        locData = [allOriUnits.locData]; 
        if strcmp(opt.bccType, 'Pen')
            nSitesAtPen  = [locData.nSitesAtPen]';
                %         BccFracAtPen = [locData.BccFracAtPen];
            cell_nSites_ok = nSitesAtPen >= opt.minNSitesPerPen;
            
            %% make sure match nSitesAtPen:
            Wrcc_pair_nsites_match_Npen_ok_C = cellfun(@(idx) (nSitesAtPen(idx(:,1)) == nSitesAtPen(idx(:,2)) ), Wrcc_pairs_ori, 'un', 0);
            assert(all(cellfun(@all, Wrcc_pair_nsites_match_Npen_ok_C)))

            %% make sure each permutation has same number
            Wrcc_pair_nsites_Npen_C = cellfun(@(idx) nSitesAtPen(idx(:,1))', Wrcc_pairs_ori, 'un', 0);
            sort_vals1 = sort(Wrcc_pair_nsites_Npen_C{1});
            assert( all( cellfun(@(x) isequal(sort(x), sort_vals1), Wrcc_pair_nsites_Npen_C) )  );

        elseif strcmp(opt.bccType, 'Animal')            
            nSitesInAnimal = [locData.nSitesInAnimal]';
    %         BccFracAtPen = [locData.BccFracAtPen];
            cell_nSites_ok = nSitesInAnimal >= opt.minNSitesPerAnimal;
            
            
             %% make sure match nSitesInAnimal:
            Wrcc_pair_nsites_match_Npen_ok_C = cellfun(@(idx) (nSitesInAnimal(idx(:,1)) == nSitesInAnimal(idx(:,2)) ), Wrcc_pairs_ori, 'un', 0);
            assert(all(cellfun(@all, Wrcc_pair_nsites_match_Npen_ok_C)))

            %% make sure each permutation has same number
            Wrcc_pair_nsites_Npen_C = cellfun(@(idx) nSitesInAnimal(idx(:,1))', Wrcc_pairs_ori, 'un', 0);
            sort_vals1 = sort(Wrcc_pair_nsites_Npen_C{1});
            assert( all( cellfun(@(x) isequal(sort(x), sort_vals1), Wrcc_pair_nsites_Npen_C) )  );
        end
        

        %%
        
        Wrcc_pair_nSites_ok_C = cellfun(@(idx) (cell_nSites_ok(idx(:,1)) | cell_nSites_ok(idx(:,2)) ) , Wrcc_pairs_ori, 'un', 0);                
        
%         pair_nSitesPerPen_ok = pair_nSitesPerPen_ok_C{1}(:);        
        
        
        %%
        if opt.excludeSameGroupPairs
            sameGroup = pairData_ori.Gids(:,1) == pairData_ori.Gids(:,2);
            
            Wrcc_idx_diffSites = cellfun(@(idx) ~sameGroup(idx), Wrcc_oo_pairIdxs, 'un', 0);            
            Wrcc_idx_diffSites_and_nSitesPerPenOK = cellfun(@and, Wrcc_idx_diffSites, Wrcc_pair_nSites_ok_C, 'un', 0);
            Wrcc_idx_ok = Wrcc_idx_diffSites_and_nSitesPerPenOK;       
        else
            Wrcc_idx_ok = Wrcc_pair_nSites_ok_C; %repmat(cellfun(@(tf) tf & pair_nSitesPerPen_ok, Wrcc_idx_diffSites, 'un', 0);            
        end
        
        Wrcc_oo_pairIdxs = cellfun(@(idx, tf) idx(tf), Wrcc_oo_pairIdxs, Wrcc_idx_ok, 'un', 0);                                                
                
        3;
    end    
    
    
    %%
    allUnitGids = [allOriUnits.Gid];
    uUnitGids_all = unique(allUnitGids);
    allUnitGids_cellsOnly = allUnitGids;
    allUnitGids_cellsOnly(~oriCells_use) = nan;
       
    [uUnitGids, idx_cellsAtEachSite] = uniqueList(allUnitGids_cellsOnly);
    
    idx_rm = isnan(uUnitGids);
    uUnitGids(idx_rm) = [];
    idx_cellsAtEachSite(idx_rm) = [];
    
    idx_unit_gids = binarySearch(uUnitGids_all, uUnitGids);
    idx_cellsAtEachSite_full = cell(size(uUnitGids_all));
    idx_cellsAtEachSite_full(idx_unit_gids) = idx_cellsAtEachSite;
    
    nCellsAtEachSite = cellfun(@length, idx_cellsAtEachSite_full);
    
    nPairsAtEachSite_expected = zeros(size(nCellsAtEachSite));
    for i = 1:length(nCellsAtEachSite)
        if nCellsAtEachSite(i) > 1
            nPairsAtEachSite_expected(i) = nchoosek(nCellsAtEachSite(i), 2);
        end
    end
   
        %%
    
    [uPairUnitGids, idx_pairsAtEachSite] = uniqueList(allUnitGids(Wcc_pairs_ori(:,1)));
        
    idx_unitpair_gids = binarySearch(uUnitGids_all, uPairUnitGids);
    idx_pairsAtEachSite_full = cell(size(uUnitGids_all));
    idx_pairsAtEachSite_full(idx_unitpair_gids) = idx_pairsAtEachSite;
    
    nPairsAtEachSite = cellfun(@length, idx_pairsAtEachSite_full);
    
    assert(isequal(nPairsAtEachSite, nPairsAtEachSite_expected));

    3;
    
end
    

function fn = extractFileName(full_filename)
    [pth, file_name, ext] = fileparts(full_filename);
    fn = [file_name ext];
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
%      x_median = nanmedian(x);
%      x_p25 = prctile(x,25);
%      x_p75 = prctile(x,75);
          
     x_prctiles = prctile(x, [25, 50, 75]);
        x_p25 = x_prctiles(1);
        x_median = x_prctiles(2);
        x_p75 = x_prctiles(3);
          
     n = nnz(~isnan(x));    
end

function [x_mean, x_median, ks_stat, n] = getMeanMedianKS(x, ctrl_dist, pairIdxs)
    doKS = (nargin >= 2) && ~isempty(ctrl_dist);

    
    if (nargin < 3) && isvector(x)
        x_mean = nanmean(x);
        x_median = nanmedian(x);        
        n = nnz(~isnan(x));
        if doKS
            ks_stat = getKSstat(x, ctrl_dist);
        else
            ks_stat = 0;
        end

    else
            
                
        x = single(x);
        if nargin == 3 % input pairIdxs
            nEachSet = cellfun(@length, pairIdxs);
            
            if any(nEachSet ~= nEachSet(1)) 
               [x_mean, x_median, ks_stat, n] = cellfun(@(idx) getMeanMedianKS(tovector( x(idx,:,:)), ctrl_dist), pairIdxs);                
               return;
            end
            
            pairIdxs_cat = [pairIdxs{:}];
            X_all = x(pairIdxs_cat);
            if size(x,3)>1
                %%
                x_2 = x(:,:,2);
                X_all = cat(1, X_all, x_2(pairIdxs_cat));
            end
        elseif ismatrix(x)
            X_all = x;
        end
        nPermutes = size(X_all, 2);
        anyNans = any(isnan(X_all(:)));

        x_mean = nanmean(X_all, 1);
        x_median = nanmedian(X_all, 1);        
        n = sum(~isnan(X_all), 1);
        if doKS
            useFastKS = ~anyNans;            
            
            if useFastKS
                divideIntoGroups = 20;
                if divideIntoGroups == 0
                    ks_stat = getKSstat_Mult(X_all, ctrl_dist);        
                else
                    ks_stat = zeros(1, nPermutes);
                    idxs = arrayfun(@(i1) i1 + [1:(nPermutes/divideIntoGroups)], ...
                        (nPermutes/divideIntoGroups)*[0:divideIntoGroups-1], 'un', 0);
%                     assert(isequal([idxs{:}], 1:nPermutes));
                    for i = 1:divideIntoGroups
                        ks_stat(idxs{i}) = getKSstat_Mult(X_all(:, idxs{i}), ctrl_dist);
                    end
                end
                
                
            else     
                ks_stat = zeros(1, nPermutes);
                for j = 1:nPermutes         
                    ks_stat(j) = getKSstat(nonnans(X_all(:,j)), ctrl_dist);
%                     ks_stat(j) = getKSstat(X_all(:,j), ctrl_dist);
%                     if anyNans
%                         ks_stat(j) = getKSstat(nonnans(X_all(:,j)), ctrl_dist);
%                     else
%                         ks_stat(j) = getKSstat_Mult(X_all, ctrl_dist);
%                     end
                end 
            end
            
            
        else
            ks_stat = nan(1, nPermutes);
        end
        
        
%          pairIdxs_cat = [pairIdxs{:}];
%          X_all = x(pairIdxs_cat);
%          x_mean = nanmean(X_all, 1);
%          x_median = nanmedian(X_all, 1);
%          n = sum(~isnan(X_all), 1);    
%          if doKS
%             ks_stat = getKSstat_Mult(X_all, ctrl_dist);
%          end

         3;
         
     
     end
     
end







function pos_inset = insetAxesPosition(pos, relDist)
    [L, B, W, H] = dealV(relDist);
%     [.5, .5, .05, .05]

    pos_inset = [pos(1) + pos(3)*L;
                 pos(2) + pos(4)*B;
                 pos(3)*W;
                 pos(4)*H];
end

function cell_stats = compareSimpleComplexCellStats(pname, p, idx_simple, fig_id, printStats)
%%
    p_simple = p(idx_simple == 1);
    p_complex = p(idx_simple == 0);
    mean_simple = nanmean(p_simple);   mean_complex = nanmean(p_complex);
    median_simple = nanmedian(p_simple);   median_complex = nanmedian(p_complex);
    [~, pval_ks] = kstest2(p_simple, p_complex); % ks test:
    [~, pval_t] = ttest2(p_simple, p_complex);
    pval_U = ranksum(p_simple, p_complex);
    
    cell_stats = struct('median_simple', median_simple, 'median_complex', median_complex, 'pval_U', pval_U, ...
        'mean_simple', mean_simple, 'mean_complex', mean_complex, 'pval_t', pval_t );
        3;

    if printStats
        fprintf(' ** %s  (N simple = %d. N complex = %d. N tot = %d): \n', pname, nnz(idx_simple == 1), nnz(idx_simple == 0), nnz( ~isnan(idx_simple) ));
        fprintf('    Simple Cells : Mean : %.3g.  Std: %.3g.  Median %.3g\n', mean_simple, nanstd(p_simple), median_simple )
        fprintf('    Complex Cells: Mean : %.3g.  Std: %.3g.  Median %.3g\n', mean_complex, nanstd(p_complex), median_complex )
        fprintf('    T-test : p = %.3g, U-test: p = %.3g.  KS test p = %.3g \n', pval_t, pval_U, pval_ks);    
    end
        
    if (nargin >= 4) && ~isnan(fig_id)
        %%
        figure(fig_id); clf;
        h = hist2({p_simple, p_complex}, 20, 'norm', 'line');
        set(h(2), 'color', 'r');
        set(h, 'linewidth', 2)
        gratingStr = curGratingType('');
        title({sprintf('\\bf %s (%s gratings)\\rm', pname, gratingStr), ...
               sprintf('Mean: Simple = %.3g, Complex = %.3g', mean_simple, mean_complex ), ...
               sprintf('Median: Simple = %.3g, Complex = %.3g', median_simple, median_complex ), ...
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


function pair_stats = compareSimpleComplexPairingStats(fig_name, measure, allMeasures, S, allNSimp, allPairIdxs, fig_id, printStats, randInfo, sc_opt)
    ms_idx = find(strcmp(measure, allMeasures), 1);
    allVals = S{ms_idx}.val;   
    
    haveOddEven = size(allVals, 3) == 2;
    if haveOddEven
        allVals = squeeze(allVals);
    end
    
    Wcc_idxs = allPairIdxs{1};
    Bcc_idxs = allPairIdxs{3};

    nsimp_wcc = allNSimp{1};
    nsimp_bcc = allNSimp{3};
    
    dval_wcc = allVals(Wcc_idxs,:);
    dval_bcc = allVals(Bcc_idxs,:);

    idx_ss_pair_wcc = nsimp_wcc == 2;
    idx_sc_pair_wcc = nsimp_wcc == 1;
    idx_cc_pair_wcc = nsimp_wcc == 0;

    idx_ss_pair_bcc = nsimp_bcc == 2;
    idx_sc_pair_bcc = nsimp_bcc == 1;
    idx_cc_pair_bcc = nsimp_bcc == 0;
    
    %% get OBSERVED SS/SC/CC distributions.
    dval_wcc_ss = nonnans( dval_wcc(idx_ss_pair_wcc, :) );
    dval_wcc_sc = nonnans( dval_wcc(idx_sc_pair_wcc, :) );
    dval_wcc_cc = nonnans( dval_wcc(idx_cc_pair_wcc, :) );
    dval_wcc_all = nonnans( dval_wcc );

    dval_bcc_ss = nonnans( dval_bcc(idx_ss_pair_bcc, :) );
    dval_bcc_sc = nonnans( dval_bcc(idx_sc_pair_bcc, :) );
    dval_bcc_cc = nonnans( dval_bcc(idx_cc_pair_bcc, :) );
    dval_bcc_all = nonnans( dval_bcc );

    nw_ss = length(dval_wcc_ss);
    nw_sc = length(dval_wcc_sc);
    nw_cc = length(dval_wcc_cc);
    nw_all = length(dval_wcc_all);

    doStandard = 0;
    doClusterIndex = 1; 
    % standard analysis
    
    utest_str = {};
    
   
    %         clusterIdxUses = 'wrcc_mean';
    %         clusterIdxUses = 'wrcc_median';
    
%     ops = {@median};           stat_names = {'median'};
    ops = {@median, @mean};    stat_names = {'median', 'mean'};
    
    %     printStats = 1;
    if printStats
        fprintf(' \n *** %s\n', measure)
    end
    
    gauss1stddev_pct = erf(1/sqrt(2)) * 100; % integral(@(x) gaussian(x, 0, 1), -1, 1) * 100;

    pair_stats = struct;
    
    for stat_i = 1:length(stat_names)
        op = ops{stat_i};
        stat_name = stat_names{stat_i};
        
        % calculate OBSERVED SS/SC/CC median/mean ratios..
        
        wcc_stat_ss = op(dval_wcc_ss); 
        wcc_stat_sc = op(dval_wcc_sc); 
        wcc_stat_cc = op(dval_wcc_cc); 
        wcc_stat_all = op(dval_wcc_all);
        
        bcc_stat_ss = op(dval_bcc_ss); 
        bcc_stat_sc = op(dval_bcc_sc); 
        bcc_stat_cc = op(dval_bcc_cc); 
        bcc_stat_all = op(dval_bcc_all); 
        
        stat_ratio_ss = bcc_stat_ss / wcc_stat_ss;
        stat_ratio_sc = bcc_stat_sc / wcc_stat_sc;
        stat_ratio_cc = bcc_stat_cc / wcc_stat_cc;
        stat_ratio_all = bcc_stat_all / wcc_stat_all;
        
        % calculate OBSERVED SS/SC/CC ratios DIFFERENCES.
        dStatRatio_ss_cc_obs = stat_ratio_ss - stat_ratio_cc;
        dStatRatio_ss_sc_obs = stat_ratio_ss - stat_ratio_sc;
        dStatRatio_cc_sc_obs = stat_ratio_cc - stat_ratio_sc;

      
        
        if doStandard
            pu_ss_sc = ranksum(dval_wcc_ss, dval_wcc_sc);
            pu_ss_cc = ranksum(dval_wcc_ss, dval_wcc_cc);
            pu_sc_cc = ranksum(dval_wcc_cc, dval_wcc_sc);
            
            [~, pt_ss_sc] = ttest2(dval_wcc_ss, dval_wcc_sc);
            [~, pt_ss_cc] = ttest2(dval_wcc_ss, dval_wcc_cc);
            [~, pt_sc_cc] = ttest2(dval_wcc_cc, dval_wcc_sc);
            
            [~, pks_ss_sc] = kstest2(dval_wcc_ss, dval_wcc_sc);
            [~, pks_ss_cc] = kstest2(dval_wcc_ss, dval_wcc_cc);
            [~, pks_sc_cc] = kstest2(dval_wcc_cc, dval_wcc_sc);
            
            if printStats
                fprintf('Differences in %s of %s \n', stat_name, fig_name);
                fprintf('  SS : %.3g. SC : %.3g. CC: %.3g.\n', wcc_stat_ss, wcc_stat_sc, wcc_stat_cc);
                fprintf('  SS/SC : U-test: p = %6.2g.  T-test p = %6.2g. KS-test: p = %6.g\n', pu_ss_sc, pt_ss_sc, pks_ss_sc);
                fprintf('  SC/CC : U-test: p = %6.2g.  T-test p = %6.2g. KS-test: p = %6.g\n', pu_sc_cc, pt_sc_cc, pks_sc_cc);
                fprintf('  SS/CC : U-test: p = %6.2g.  T-test p = %6.2g. KS-test: p = %6.g\n\n', pu_ss_cc, pt_ss_cc, pks_ss_cc);
                
                utest_str = {sprintf('Md: SS= %.3g, SC = %.3g, CC = %.3g', wcc_stat_ss, wcc_stat_sc, wcc_stat_cc), ...
                    sprintf('p_U: SS/SC=%.3g,   SC/CC = %.3g,   SS/CC = %.3g', pu_ss_sc, pu_sc_cc, pu_ss_cc)};
            end
            
%              ...'median_SS', median_SS, 'median_SC', median_SC, 'median_CC', median_CC, ...
%                             ...'mean_SS', mean_SS, 'mean_SC', mean_SC, 'mean_CC', mean_CC, ...
%                             ...'pU_ss_sc', pu_ss_sc, 'pU_ss_cc', pu_ss_cc, 'pU_sc_cc', pu_sc_cc, ...
%                             ...'pT_ss_sc', pt_ss_sc, 'pT_ss_cc', pt_ss_cc, 'pT_sc_cc', pt_sc_cc, ...
%                             ...'pKS_ss_sc', pks_ss_sc, 'pKS_ss_cc', pks_ss_cc, 'pKS_sc_cc', pks_sc_cc); 
            
        end
        
        
        if doClusterIndex   % clustering index approach.

            nResamplesTotal = randInfo.nResamplesTotal;
            haveAllRandSamplesAlready = ~isempty(allPairIdxs{4}) || (~sc_opt.getBrccPairs && length(allPairIdxs{2}) == nResamplesTotal);

            randInfo.doProgressBar = 0;
            
            nResamplesPerLoop = min(200, nResamplesTotal);
            if haveAllRandSamplesAlready
                nResamplesPerLoop = nResamplesTotal;
            end
            nLoops = nResamplesTotal / nResamplesPerLoop;
            
            all_loop_idxs = arrayfun(@(loop_i)  [1:nResamplesPerLoop] + (loop_i-1)*nResamplesPerLoop,  1:nLoops, 'un', 0);
            
            [wrcc_stat_ss, wrcc_stat_sc, wrcc_stat_cc, wrcc_stat_all,   ...
             brcc_stat_ss, brcc_stat_sc, brcc_stat_cc, brcc_stat_all] = deal(zeros(1, nResamplesTotal));
            %%
            if nLoops > 1
                fprintf('\n[%s]', measure);
                progressBar('init', nLoops, 30);
            end
            for loop_i = 1:nLoops
                loop_idxs = all_loop_idxs{loop_i};
                
                if haveAllRandSamplesAlready
                    Wrcc_idxs = allPairIdxs{2}(loop_idxs);  nsimp_wrcc = allNSimp{2}(loop_idxs);

                    if sc_opt.getBrccPairs
                        Brcc_idxs = allPairIdxs{4}(loop_idxs);
                        nsimp_brcc = allNSimp{4}(loop_idxs);
                    end
                else
                    [Wrcc_idxs, nsimp_wrcc,   Brcc_idxs, nsimp_brcc] = getRandomizedPairs(randInfo, loop_i==1, nResamplesPerLoop);                
                end
                
                
                idx_ss_pair_wrcc = cellfun(@(n) find(n == 2), nsimp_wrcc, 'un', 0);
                idx_sc_pair_wrcc = cellfun(@(n) find(n == 1), nsimp_wrcc, 'un', 0);
                idx_cc_pair_wrcc = cellfun(@(n) find(n == 0), nsimp_wrcc, 'un', 0);
                %%
                if sc_opt.getBrccPairs
                    idx_ss_pair_brcc = cellfun(@(n) find(n == 2), nsimp_brcc, 'un', 0);
                    idx_sc_pair_brcc = cellfun(@(n) find(n == 1), nsimp_brcc, 'un', 0);
                    idx_cc_pair_brcc = cellfun(@(n) find(n == 0), nsimp_brcc, 'un', 0);
                end          
            
%                 justRandomizeSClabels = 1  && 1;
                if sc_opt.justRandomizeSClabels
                    % use the original within-site distribution, just index according to randomized labels
                    dval_wrcc_ss = cellfun(@(tf) nonnans( dval_wcc(tf) ), idx_ss_pair_wrcc, 'un', 0);
                    dval_wrcc_sc = cellfun(@(tf) nonnans( dval_wcc(tf) ), idx_cc_pair_wrcc, 'un', 0);
                    dval_wrcc_cc = cellfun(@(tf) nonnans( dval_wcc(tf) ), idx_sc_pair_wrcc, 'un', 0);
                    dval_wrcc_all(1:nResamplesPerLoop) = { nonnans( dval_wcc ) };

                    if sc_opt.getBrccPairs
                        dval_brcc_ss = cellfun(@(tf) nonnans( dval_bcc(tf) ), idx_ss_pair_brcc, 'un', 0);
                        dval_brcc_sc = cellfun(@(tf) nonnans( dval_bcc(tf) ), idx_cc_pair_brcc, 'un', 0);
                        dval_brcc_cc = cellfun(@(tf) nonnans( dval_bcc(tf) ), idx_sc_pair_brcc, 'un', 0);
                        dval_brcc_all(1:nResamplesPerLoop) = { nonnans( dval_bcc ) };
                    end
                else
                    %%
                    dval_wrcc = cellfun(@(idx) allVals(idx,:), Wrcc_idxs, 'un', 0);
                    
                    % randomize entire cells (values with S/C labels)
                    dval_wrcc_ss = cellfun(@(x,idx) nonnans( x(idx,:) ), dval_wrcc, idx_ss_pair_wrcc, 'un', 0);
                    dval_wrcc_sc = cellfun(@(x,idx) nonnans( x(idx,:) ), dval_wrcc, idx_cc_pair_wrcc, 'un', 0);
                    dval_wrcc_cc = cellfun(@(x,idx) nonnans( x(idx,:) ), dval_wrcc, idx_sc_pair_wrcc, 'un', 0);
                    dval_wrcc_all = cellfun(@(idxs) nonnans( allVals(idxs,:) ), Wrcc_idxs, 'un', 0);
                    
                    if sc_opt.getBrccPairs
                        dval_brcc = cellfun(@(idx) allVals(idx), Brcc_idxs, 'un', 0);
                        
                        dval_brcc_ss = cellfun(@(x,idx) nonnans( x(idx,:) ), dval_brcc, idx_ss_pair_brcc, 'un', 0);
                        dval_brcc_sc = cellfun(@(x,idx) nonnans( x(idx,:) ), dval_brcc, idx_cc_pair_brcc, 'un', 0);
                        dval_brcc_cc = cellfun(@(x,idx) nonnans( x(idx,:) ), dval_brcc, idx_sc_pair_brcc, 'un', 0);
                        dval_brcc_all = cellfun(@(idxs) nonnans( allVals(idxs,:) ), Brcc_idxs, 'un', 0); 
                    end

                end
            

                % calculate RANDOMIZED SS/SC/CC distributions..
                wrcc_stat_ss(loop_idxs) = cellfun(op, dval_wrcc_ss);
                wrcc_stat_sc(loop_idxs) = cellfun(op, dval_wrcc_sc);
                wrcc_stat_cc(loop_idxs) = cellfun(op, dval_wrcc_cc);
                wrcc_stat_all(loop_idxs)= cellfun(op, dval_wrcc_all);

                switch sc_opt.clusterIdxUses
                    case 'bcc',       
                        brcc_stat_ss(loop_idxs)  = bcc_stat_ss;      
                        brcc_stat_sc(loop_idxs)  = bcc_stat_sc;      
                        brcc_stat_cc(loop_idxs)  = bcc_stat_cc;      
                        brcc_stat_all(loop_idxs) = bcc_stat_all;         
                    case 'brcc',       
                        brcc_stat_ss(loop_idxs) = cellfun(op, dval_brcc_ss);
                        brcc_stat_sc(loop_idxs) = cellfun(op, dval_brcc_sc);
                        brcc_stat_cc(loop_idxs) = cellfun(op, dval_brcc_cc);
                        brcc_stat_all(loop_idxs)= cellfun(op, dval_brcc_all);
                end
                
                if nLoops > 1
                    progressBar(loop_i);
                end
                
            end
            if nLoops > 1
                progressBar('done');
            end

            % calculate RANDOMIZED median/mean ratios:
            stat_ratio_all_rand = brcc_stat_all ./ wrcc_stat_all;
            stat_ratio_ss_rand  = brcc_stat_ss ./ wrcc_stat_ss;
            stat_ratio_sc_rand  = brcc_stat_sc ./ wrcc_stat_sc;
            stat_ratio_cc_rand  = brcc_stat_cc ./ wrcc_stat_cc;
            
            
            pair_stats.(sprintf('%sRatio_all', stat_name))  = stat_ratio_all;
            pair_stats.(sprintf('%sRatio_ss', stat_name))   = stat_ratio_ss;
            pair_stats.(sprintf('%sRatio_sc', stat_name))   = stat_ratio_sc;
            pair_stats.(sprintf('%sRatio_cc', stat_name))   = stat_ratio_cc;

            
            pair_stats.(sprintf('%sRatio_all__rand', stat_name))  = stat_ratio_all_rand;
            pair_stats.(sprintf('%sRatio_ss__rand', stat_name))   = stat_ratio_ss_rand;
            pair_stats.(sprintf('%sRatio_sc__rand', stat_name))   = stat_ratio_sc_rand;
            pair_stats.(sprintf('%sRatio_cc__rand', stat_name))   = stat_ratio_cc_rand;
            
            pair_stats.(sprintf('tmp_bcc%s_all', stat_name))  = bcc_stat_all;
            pair_stats.(sprintf('tmp_bcc%s_ss', stat_name))   = bcc_stat_ss;
            pair_stats.(sprintf('tmp_bcc%s_sc', stat_name))   = bcc_stat_sc;
            pair_stats.(sprintf('tmp_bcc%s_cc', stat_name))   = bcc_stat_cc;
            pair_stats.(sprintf('tmp_wcc%s_all', stat_name))  = wcc_stat_all;
            pair_stats.(sprintf('tmp_wcc%s_ss', stat_name))   = wcc_stat_ss;
            pair_stats.(sprintf('tmp_wcc%s_sc', stat_name))   = wcc_stat_sc;
            pair_stats.(sprintf('tmp_wcc%s_cc', stat_name))   = wcc_stat_cc;

            
            statRatio_str = sprintf('  %s ratios: All pairs: %.3f (%d).  SS: %.3f (%d)   SC: %.3f (%d)   CC: %.3f (%d)', ...
                stat_name, stat_ratio_all, nw_all, stat_ratio_ss, nw_ss, stat_ratio_sc, nw_sc, stat_ratio_cc, nw_cc);
           
            
            switch sc_opt.calcPairTypeSig
                case 'differences'
            
                % calculate RANDOMIZED SS/SC/CC ratios DIFFERENCES.
                dStatRatio_ss_cc_rand = stat_ratio_ss_rand - stat_ratio_cc_rand;
                dStatRatio_ss_sc_rand = stat_ratio_ss_rand - stat_ratio_sc_rand;
                dStatRatio_cc_sc_rand = stat_ratio_cc_rand - stat_ratio_sc_rand;
                %%

                % compare OBSERVED and RANDOMIZED differences, and get p-values
                [p_stat_ss_cc, stat_ss_cc_sgn] = getRandomizedProb(dStatRatio_ss_cc_obs, dStatRatio_ss_cc_rand);
                [p_stat_ss_sc, stat_ss_sc_sgn] = getRandomizedProb(dStatRatio_ss_sc_obs, dStatRatio_ss_sc_rand);
                [p_stat_cc_sc, stat_cc_sc_sgn] = getRandomizedProb(dStatRatio_cc_sc_obs, dStatRatio_cc_sc_rand);
                stat_ss_cc_str = getSCstr('SS', 'CC', stat_ss_cc_sgn, p_stat_ss_cc);
                stat_ss_sc_str = getSCstr('SS', 'SC', stat_ss_sc_sgn, p_stat_ss_sc);
                stat_cc_sc_str = getSCstr('CC', 'SC', stat_cc_sc_sgn, p_stat_cc_sc);
                
                statProbs_str = sprintf('  %s probs: %s  %s  %s', stat_name, stat_ss_cc_str, stat_ss_sc_str, stat_cc_sc_str);

                pair_stats.(sprintf('%sProb_ss_cc', stat_name)) = p_stat_ss_cc;
                pair_stats.(sprintf('%sProb_ss_sc', stat_name)) = p_stat_ss_sc;
                pair_stats.(sprintf('%sProb_cc_sc', stat_name)) = p_stat_cc_sc;
             case 'indiv'
                statProb_ss = getRandomizedProb(stat_ratio_ss, stat_ratio_ss_rand, 'right');
                statProb_sc = getRandomizedProb(stat_ratio_sc, stat_ratio_sc_rand, 'right');
                statProb_cc = getRandomizedProb(stat_ratio_cc, stat_ratio_cc_rand, 'right');

                pair_stats.(sprintf('%sProb_ss', stat_name)) = statProb_ss;
                pair_stats.(sprintf('%sProb_sc', stat_name)) = statProb_sc;
                pair_stats.(sprintf('%sProb_cc', stat_name)) = statProb_cc;

                statProbs_str = sprintf('  %s probs: SS: p = %.1g.  SC: p = %.1g   CC: p = %.1g.', stat_name, statProb_ss, statProb_sc, statProb_cc);

            end

            
            
            if sc_opt.calcStatRatioCIs
                nBoots = sc_opt.nBoots;
                stats_wcc_ss_boot = bootstrp(nBoots, op, dval_wcc_ss);  stats_bcc_ss_boot = bootstrp(nBoots, op, dval_bcc_ss);
                stats_wcc_sc_boot = bootstrp(nBoots, op, dval_wcc_sc);  stats_bcc_sc_boot = bootstrp(nBoots, op, dval_bcc_sc);
                stats_wcc_cc_boot = bootstrp(nBoots, op, dval_wcc_cc);  stats_bcc_cc_boot = bootstrp(nBoots, op, dval_bcc_cc);
                
                statRatios_ss_boot = stats_bcc_ss_boot ./ stats_wcc_ss_boot;
                statRatios_sc_boot = stats_bcc_sc_boot ./ stats_wcc_sc_boot;
                statRatios_cc_boot = stats_bcc_cc_boot ./ stats_wcc_cc_boot;
                
                statRatio_ss_ci_68 = getBootCI(statRatios_ss_boot, gauss1stddev_pct);
                statRatio_sc_ci_68 = getBootCI(statRatios_sc_boot, gauss1stddev_pct);
                statRatio_cc_ci_68 = getBootCI(statRatios_cc_boot, gauss1stddev_pct);
                
                pair_stats.(sprintf('%sRatio_ss_lo', stat_name)) = statRatio_ss_ci_68(1); pair_stats.(sprintf('%sRatio_ss_hi', stat_name)) = statRatio_ss_ci_68(2);
                pair_stats.(sprintf('%sRatio_sc_lo', stat_name)) = statRatio_sc_ci_68(1); pair_stats.(sprintf('%sRatio_sc_hi', stat_name)) = statRatio_sc_ci_68(2);
                pair_stats.(sprintf('%sRatio_cc_lo', stat_name)) = statRatio_cc_ci_68(1); pair_stats.(sprintf('%sRatio_cc_hi', stat_name)) = statRatio_cc_ci_68(2);
                
            end

            
            
            
            
        end
        
        
        
        if printStats
            fprintf(' * Clustering indices\n');
            fprintf('%s\n', statRatio_str);
            fprintf('%s\n', statProbs_str);
        end
        clusterIdx_str = {statRatio_str, statProbs_str};
        
        
        addNPairs = 0;
        if addNPairs
            pair_stats.N_ss = length(dval_wcc_ss);
            pair_stats.N_sc = length(dval_wcc_sc);
            pair_stats.N_cc = length(dval_wcc_cc);                        
            
        end
        
%         pair_stats = struct('medianRatio_all', stat_bw_ratio_all, 'medianRatio_ss', stat_bw_ratio_ss, 'medianRatio_cc', stat_bw_ratio_cc, 'medianRatio_sc', stat_bw_ratio_sc, ...
%                             'medianProb_ss_cc', p_stat_ss_cc, 'medianProb_ss_sc', p_stat_ss_sc, 'medianProb_cc_sc', p_stat_cc_sc, ...
%                             'meanRatio_all', mean_bw_ratio_all, 'meanRatio_ss', mean_bw_ratio_ss,  'meanRatio_cc', mean_bw_ratio_cc, 'meanRatio_sc', mean_bw_ratio_sc, ...
%                             'meanProb_ss_cc', p_mean_ss_cc, 'meanProb_ss_sc', p_mean_ss_sc, 'meanProb_cc_sc', p_mean_cc_sc ...
%                             ...'median_SS', median_SS, 'median_SC', median_SC, 'median_CC', median_CC, ...
%                             ...'mean_SS', mean_SS, 'mean_SC', mean_SC, 'mean_CC', mean_CC, ...
%                             ...'pU_ss_sc', pu_ss_sc, 'pU_ss_cc', pu_ss_cc, 'pU_sc_cc', pu_sc_cc, ...
%                             ...'pT_ss_sc', pt_ss_sc, 'pT_ss_cc', pt_ss_cc, 'pT_sc_cc', pt_sc_cc, ...
%                             ...'pKS_ss_sc', pks_ss_sc, 'pKS_ss_cc', pks_ss_cc, 'pKS_sc_cc', pks_sc_cc); 
%                             );

   
    3;

        if exist('fig_id', 'var') && ~isnan(fig_id);
            %%
             oe_str = strrep(curDegreeOEmode(''), '_', '-');
            figure(fig_id + (stat_i-1)*20 ); clf;
            if doStandard && ~doClusterIndex
                h = hist2({dval_wcc_ss, dval_wcc_sc, dval_wcc_cc}, 15, 'line', 'norm');
                set(h(1), 'color', 'b');
                set(h(2), 'color', 'm');
                set(h(3), 'color', 'r');
                title(fig_name); 
                legend('SS', 'SC', 'CC');                 

                set(h, 'linewidth', 2)
                gratingStr = curGratingType('');
               
                title([{sprintf('%s (%s gratings) [%s]', fig_name, gratingStr, oe_str)}, utest_str, clusterIdx_str], 'fontsize', 10) 
                ylims = ylim;
        %         if pval_U < .05
                line(median(dval_wcc_ss)*[1, 1],  [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'b', 'linestyle', '--');
                line(median(dval_wcc_sc)*[1, 1],  [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'm', 'linestyle', '--');
                line(median(dval_wcc_cc)*[1, 1],  [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'r', 'linestyle', '--');

                line(mean(dval_wcc_ss)*[1, 1],  [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'b', 'linestyle', ':');
                line(mean(dval_wcc_sc)*[1, 1],  [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'm', 'linestyle', ':');
                line(mean(dval_wcc_cc)*[1, 1],  [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'r', 'linestyle', ':');
            elseif doClusterIndex
    %%
                nBin = 30;
                showOrigDists = 1;
                showDifferences = 1;
                nSubM = showOrigDists + showDifferences;
                prct_lims = 0.02;
                clf;
                colors_C = {'g', 'b', 'r'};
                gratingStr = curGratingType('');        

                switch sc_opt.calcPairTypeSig
                    case 'indiv'
                        plotDData_rnd = {stat_ratio_ss_rand, stat_ratio_sc_rand, stat_ratio_cc_rand};
                        plotDData_obs = [stat_ratio_ss,      stat_ratio_sc,      stat_ratio_cc];
                        diff_str = '';
                        leg_strs = {'SS', 'SC', 'CC'};

                    case 'differences'
                        plotDData_rnd = {dStatRatio_ss_cc_rand, dStatRatio_ss_sc_rand, dStatRatio_cc_sc_rand};
                        plotDData_obs = [dStatRatio_ss_cc_obs,  dStatRatio_ss_sc_obs,  dStatRatio_cc_sc_obs];
                        
                        diff_str = 'Differences in ';
                        leg_strs = {'SS-CC', 'SS-SC', 'CC-SC'};
                end
                
                title_str = {statRatio_str, statProbs_str};

                dxlims = lims( [prctile([plotDData_rnd{:}], [prct_lims, 100-prct_lims]), plotDData_obs], .05 );

                cur_row_idx = 1;
%                 if stat_i == 1
                    title_str = [title_str, ['\bf ' diff_str fig_name ' (Cluster Indices) (' gratingStr ' [' oe_str ']) \rm' ]];
%                 end
                
                if showDifferences
%                 subplotGap(length(stat_names),2,stat_i, 1);
                    subplotGap(nSubM, 1, cur_row_idx, 1);

                    h = hist2(plotDData_rnd, nBin, 'line', 'norm');
    %                 xlims = xlim;            
                    set(h, {'color'}, colors_C');
                    set(h, 'linewidth', 2);
    %                 if stat_i == 1
                        legend(leg_strs, 'location', 'NE');        
    %                 end                

                    uistack(h(1), 'top');
                    title(title_str);

                    ylims = ylim;
                    for j = 1:3
                        line(plotDData_obs(j)*[1, 1],  ylims, 'color', colors_C{j}, 'linestyle', '-', 'linewidth', 2);
                    end
    %                 xlims = lims([xlims(:)', plotDData_obs], .05);
                    xlim(dxlims);
                    xlabel(sprintf('Difference in %s ratio', stat_name))

                    cur_row_idx = cur_row_idx + 1;
   
                end
                
                if showOrigDists
%                 h_ax2 = subplotGap(length(stat_names),2,stat_i, 2);
                    h_ax2 = subplotGap(nSubM, 1, cur_row_idx, 1);
    %                     h_ax2 = axes('position', [0.6, 0.5, .3, .3]);


                    plotData_rnd = {stat_ratio_ss_rand, stat_ratio_sc_rand, stat_ratio_cc_rand};
                    plotData_obs = [stat_ratio_ss, stat_ratio_sc, stat_ratio_cc];

                    xlims = lims( [prctile([plotData_rnd{:}], [prct_lims, 100-prct_lims]), plotData_obs], .05 );

                    h_hist2 = hist2(plotData_rnd, nBin, 'line', 'norm');
                    colors_C2 = {'m', 'k', 'c'};
                    set(h_hist2, {'color'}, colors_C2');
                    set(h_hist2, 'linewidth', 2);
    %                 set(h_hist2(4), 'linestyle', ':')
                    set(h_ax2, 'xlim', xlims);

                    ylims = ylim;
                    for j = 1:3
                        line(plotData_obs(j)*[1, 1],  ylims, 'color', colors_C2{j}, 'linestyle', '-', 'linewidth', 2);
                    end
                    h_line_all = drawVerticalLine(stat_ratio_all, 'color', 'r', 'linewidth', 3, 'linestyle', '--');
                    xlabel(sprintf('%s ratio', stat_name))
                
                    if cur_row_idx == 1
                        title(title_str);
                    end
%                 if stat_i == 1
                    legend([h_hist2; h_line_all], {'SS', 'SC', 'CC', '(all)'}, 'location', 'NE');        
%                 end    
                end                 

3;

            end

        end    
    end
                                                       
                                                   
                                                   

end


function pair_stats = compareSimpleComplexPairingStats_boot(fig_name, measure, allMeasures, S, allNSimp, allPairIdxs, fig_id, printStats)
         
    ms_idx = find(strcmp(measure, allMeasures), 1);
%     medianRatio_all: 1.1551
%       medianRatio_ss: 1.0249
%       medianRatio_cc: 1.0382
%       medianRatio_sc: 1.3725
%     medianProb_ss_cc: 0.8436
%     medianProb_ss_sc: 0.1157
%     medianProb_cc_sc: 0.0346
%        meanRatio_all: 1.1556
%         meanRatio_ss: 1.1109
%         meanRatio_cc: 1.0319
%         meanRatio_sc: 1.3806
%       meanProb_ss_cc: 0.6344
%       meanProb_ss_sc: 0.0695
%       meanProb_cc_sc: 0.0017
    

%%
    wcc_idx = 1;
    bcc_idx = 3;

    allVals = S{ms_idx}.val;
    idx_wcc_all = allPairIdxs{wcc_idx};
    idx_wcc_ss = idx_wcc_all(  allNSimp{wcc_idx} == 2  );
    idx_wcc_sc = idx_wcc_all(  allNSimp{wcc_idx} == 1  );
    idx_wcc_cc = idx_wcc_all(  allNSimp{wcc_idx} == 0  );

    idx_bcc_all = allPairIdxs{bcc_idx};
    idx_bcc_ss =  idx_bcc_all( allNSimp{bcc_idx} == 2 );
    idx_bcc_sc =  idx_bcc_all( allNSimp{bcc_idx} == 1 );
    idx_bcc_cc =  idx_bcc_all( allNSimp{bcc_idx} == 0 );

    vals_w_all = nonnans( allVals(idx_wcc_all) );
    vals_w_ss =  nonnans( allVals(idx_wcc_ss) );
    vals_w_sc = nonnans( allVals(idx_wcc_sc) );
    vals_w_cc = nonnans( allVals(idx_wcc_cc) );
    
    vals_b_all = nonnans( allVals(idx_bcc_all) );
    vals_b_ss = nonnans( allVals(idx_bcc_ss) );
    vals_b_sc = nonnans( allVals(idx_bcc_sc) );
    vals_b_cc = nonnans( allVals(idx_bcc_cc) );
%%
%     ops = {@median, @mean};
%     stat_names = {'median', 'mean'};

    ops = {@median};
    stat_names = {'median'};

%     printStats = 1;
    if printStats
        fprintf(' \n *** %s\n', measure)
    end
    
    gauss1stddev_pct = integral(@(x) gaussian(x, 0, 1), -1, 1) * 100;
    
    for i = 1:length(stat_names)
        op = ops{i};
        stat_name = stat_names{i};

        stats_w_ss_boot = bootstrp(nBoots, op, vals_w_ss); fprintf('.')
        stats_w_sc_boot = bootstrp(nBoots, op, vals_w_sc); fprintf('.')
        stats_w_cc_boot = bootstrp(nBoots, op, vals_w_cc); fprintf('.')

        stats_b_ss_boot = bootstrp(nBoots, op, vals_b_ss); fprintf('.')
        stats_b_sc_boot = bootstrp(nBoots, op, vals_b_sc); fprintf('.')
        stats_b_cc_boot = bootstrp(nBoots, op, vals_b_cc); fprintf('.')

        statRatio_all = op(vals_b_all) / op(vals_w_all);
        statRatio_ss = op(vals_b_ss) / op(vals_w_ss);
        statRatio_sc = op(vals_b_sc) / op(vals_w_sc);
        statRatio_cc = op(vals_b_cc) / op(vals_w_cc);
        
        
        statRatios_ss_boot = stats_b_ss_boot ./ stats_w_ss_boot;
        statRatios_sc_boot = stats_b_sc_boot ./ stats_w_sc_boot;
        statRatios_cc_boot = stats_b_cc_boot ./ stats_w_cc_boot;

        statRatios_ss_boot_ci = getBootCI(statRatios_ss_boot, gauss1stddev_pct);
        statRatios_sc_boot_ci = getBootCI(statRatios_sc_boot, gauss1stddev_pct);
        statRatios_cc_boot_ci = getBootCI(statRatios_cc_boot, gauss1stddev_pct);
                
        
        dStatRatio_ss_cc_obs = statRatio_ss - statRatio_cc; dStatRatio_ss_cc_rand = statRatios_ss_boot - statRatios_cc_boot;
        dStatRatio_ss_sc_obs = statRatio_ss - statRatio_sc; dStatRatio_ss_sc_rand = statRatios_ss_boot - statRatios_sc_boot;
        dStatRatio_cc_sc_obs = statRatio_cc - statRatio_sc; dStatRatio_cc_sc_rand = statRatios_cc_boot - statRatios_sc_boot;
%%
        [p_stat_ss_cc3, med_ss_cc_sgn] = getRandomizedProb(abs(dStatRatio_ss_cc_obs), abs(dStatRatio_ss_cc_rand), 'right');
        [p_stat_ss_sc3, med_ss_sc_sgn] = getRandomizedProb(abs(dStatRatio_ss_sc_obs), abs(dStatRatio_ss_sc_rand), 'right');
        [p_stat_cc_sc3, med_cc_sc_sgn] = getRandomizedProb(abs(dStatRatio_cc_sc_obs), abs(dStatRatio_cc_sc_rand), 'right');

%         [p_stat_ss_cc, med_ss_cc_sgn] = getRandomizedProb(dStatRatio_ss_cc_obs), abs(dStatRatio_ss_cc_rand), 'right');
%         [p_stat_ss_sc, med_ss_sc_sgn] = getRandomizedProb(dStatRatio_ss_sc_obs), abs(dStatRatio_ss_sc_rand), 'right');
%         [p_stat_cc_sc, med_cc_sc_sgn] = getRandomizedProb(dStatRatio_cc_sc_obs), abs(dStatRatio_cc_sc_rand), 'right');

%         [p_stat_ss_cc, med_ss_cc_sgn] = getRandomizedProb(dStatRatio_ss_cc_obs, dStatRatio_ss_cc_rand);
%         [p_stat_ss_sc, med_ss_sc_sgn] = getRandomizedProb(dStatRatio_ss_sc_obs, dStatRatio_ss_sc_rand);
%         [p_stat_cc_sc, med_cc_sc_sgn] = getRandomizedProb(dStatRatio_cc_sc_obs, dStatRatio_cc_sc_rand);
%%
        p_stat_ss_cc = getFractionOfDifferencesWithSameSignAs(statRatios_ss_boot, statRatios_cc_boot, dStatRatio_ss_cc_obs, 1);
        p_stat_ss_sc = getFractionOfDifferencesWithSameSignAs(statRatios_ss_boot, statRatios_sc_boot, dStatRatio_ss_cc_obs, 1);
        p_stat_cc_sc = getFractionOfDifferencesWithSameSignAs(statRatios_cc_boot, statRatios_sc_boot, dStatRatio_ss_cc_obs, 1);
        
%         p_stat_ss_cc2 = getFractionOfDifferencesWithSameSignAs(statRatios_ss_boot, statRatios_cc_boot, dStatRatio_ss_cc_obs, 1);
%         p_stat_ss_sc2 = getFractionOfDifferencesWithSameSignAs(statRatios_ss_boot, statRatios_sc_boot, dStatRatio_ss_cc_obs, 1);
%         p_stat_cc_sc2 = getFractionOfDifferencesWithSameSignAs(statRatios_cc_boot, statRatios_sc_boot, dStatRatio_ss_cc_obs, 1);        
%         
%         p_stat_ss_cc3 = nnz( (dStatRatio_ss_cc_rand) * -sign(dStatRatio_ss_cc_obs) > 0)/nBoots;
%         p_stat_ss_sc3 = nnz( (dStatRatio_ss_sc_rand) * -sign(dStatRatio_ss_sc_obs) > 0)/nBoots;
%         p_stat_cc_sc3 = nnz( (dStatRatio_cc_sc_rand) * -sign(dStatRatio_ss_sc_obs) > 0)/nBoots;
        
        pU_stat_ss_cc = ranksum(statRatios_ss_boot, statRatios_cc_boot);
        pU_stat_ss_sc = ranksum(statRatios_ss_boot, statRatios_sc_boot);
        pU_stat_cc_sc = ranksum(statRatios_cc_boot, statRatios_sc_boot);


        if printStats
            %%
            med_ss_cc_str = getSCstr('SS', 'CC', med_ss_cc_sgn, p_stat_ss_cc);
            med_ss_sc_str = getSCstr('SS', 'SC', med_ss_sc_sgn, p_stat_ss_sc);
            med_cc_sc_str = getSCstr('CC', 'SC', med_cc_sc_sgn, p_stat_cc_sc);

            statRatio_str = sprintf('  %s ratios (observed): All pairs: %.3f.  SS: %.3f   SC: %.3f   CC: %.3f', stat_name,  statRatio_all, statRatio_ss, statRatio_sc, statRatio_cc );
            statProbs_str = sprintf('  %s probs: %s  %s  %s', stat_name, med_ss_cc_str, med_ss_sc_str, med_cc_sc_str);

            statProbsU_str = sprintf('  %s probs (U-test): SS-CC: %.1g.  SS-SC: %.1g.  CC-SC: %.1g. ', stat_name, pU_stat_ss_cc, pU_stat_ss_sc, pU_stat_cc_sc);
            
            if any([p_stat_ss_cc, p_stat_ss_sc, p_stat_cc_sc] < .1)
                beep;
            end
            fprintf(' * Clustering indices\n');
            fprintf('%s\n', statRatio_str);
            fprintf('%s\n', statProbs_str);
            fprintf('%s\n', statProbsU_str);
        end        
        
        pair_stats.(sprintf('%sRatio_all', stat_name))  = statRatio_all;
        
        pair_stats.(sprintf('%sRatio_ss', stat_name))   = statRatio_ss;
        pair_stats.(sprintf('%sRatio_sc', stat_name))   = statRatio_sc;
        pair_stats.(sprintf('%sRatio_cc', stat_name))   = statRatio_cc;
        
        pair_stats.(sprintf('%sRatio_ssLO', stat_name))   = statRatios_ss_boot_ci(1);
        pair_stats.(sprintf('%sRatio_ssHI', stat_name))   = statRatios_ss_boot_ci(2);

        pair_stats.(sprintf('%sRatio_scLO', stat_name))   = statRatios_sc_boot_ci(1);
        pair_stats.(sprintf('%sRatio_scHI', stat_name))   = statRatios_sc_boot_ci(2);

        pair_stats.(sprintf('%sRatio_ccLO', stat_name))   = statRatios_sc_boot_ci(1);
        pair_stats.(sprintf('%sRatio_ccHI', stat_name))   = statRatios_cc_boot_ci(2);
        
%         pair_stats.(sprintf('%sProb_ss_cc', stat_name)) = p_stat_ss_cc;
%         pair_stats.(sprintf('%sProb_ss_sc', stat_name)) = p_stat_ss_sc;
%         pair_stats.(sprintf('%sProb_cc_sc', stat_name)) = p_stat_cc_sc;
        
        if strcmp(stat_name, 'median') && ~isnan(fig_id)
            %%
            figure(fig_id)
            h = hist2({statRatios_ss_boot, statRatios_sc_boot, statRatios_cc_boot}, 40, 'stairs');
            set(h, 'linewidth', 2)
            tit_str = [fig_name ' (' titleCase( curGratingType(''))  ' Gratings)'  ];
            allTitle_str = {tit_str, statRatio_str, statProbs_str, statProbsU_str};
            legend({'SS', 'SC', 'CC'}); xlabel('Median Ratios'); ylabel('Count'); 
            title(allTitle_str, 'interpreter', 'none')
%             subplot(1,3,1); hist(dStatRatio_ss_cc_rand, 30);
%             subplot(1,3,2); hist(dStatRatio_ss_sc_rand, 30);
%             subplot(1,3,3); hist(dStatRatio_cc_sc_rand, 30);
            
        end
        3;
    end
    3;

%%    
%         pair_stats = struct('medianRatio_all', median_bw_ratio_all, 'medianRatio_ss', median_bw_ratio_ss, 'medianRatio_cc', median_bw_ratio_cc, 'medianRatio_sc', median_bw_ratio_sc, ...
%                             'medianProb_ss_cc', p_median_ss_cc, 'medianProb_ss_sc', p_median_ss_sc, 'medianProb_cc_sc', p_median_cc_sc, ...
%                             'meanRatio_all', mean_bw_ratio_all, 'meanRatio_ss', mean_bw_ratio_ss,  'meanRatio_cc', mean_bw_ratio_cc, 'meanRatio_sc', mean_bw_ratio_sc, ...
%                             'meanProb_ss_cc', p_mean_ss_cc, 'meanProb_ss_sc', p_mean_ss_sc, 'meanProb_cc_sc', p_mean_cc_sc ...
end




% function pair_stats = compareSimpleComplexPairingStats_main(dval_wcc, dval_bcc, dval_wrcc, nsimp_wcc, nsimp_bcc, nsimp_wrcc, dname, fig_id, printStats)
%     3;
%     
% end






function str = getSCstr(s1, s2, sgn, pval)
%     if sgn == -1
%         [s1, s2] = deal(s2, s1);
%     end    
%     str = sprintf('%s>%s, p = %.4f', s1, s2, pval);
    symb = iff(sgn == 1, '>', '<');
    str = sprintf('%s%s%s, p = %.4f', s1, symb, s2, pval);
end


function [pairIdxs, pairIdxList, idxMtx] = useSubsetOfPairIdxs(allPT, pairTypes, nUnits)
    %%
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
    tf_use = ~isnan(x1) & ~isnan(x2) & ~any(isnan(pairF1oDCs), 2);
    x1 = x1(tf_use);
    x2 = x2(tf_use);
    pairF1oDCs = pairF1oDCs(tf_use, :);

    nsimp = sum(pairF1oDCs > 1, 2);
    idx_ss_pair = nsimp == 2;
    idx_cc_pair = nsimp == 0;
    idx_sc_pair = nsimp == 1;
    cc_str = switchh(corrType, {'pearson', 'spearman'}, {'cc', '\\rho'});
    
    x1 = x1(:);
    x2 = x2(:);
    [v_all.cc, v_all.p] = corr(x1, x2, 'type', corrType);
    [v_ss.cc, v_ss.p] = corr(x1(idx_ss_pair), x2(idx_ss_pair), 'type', corrType);
    [v_sc.cc, v_sc.p] = corr(x1(idx_sc_pair), x2(idx_sc_pair), 'type', corrType);
    [v_cc.cc, v_cc.p] = corr(x1(idx_cc_pair), x2(idx_cc_pair), 'type', corrType);
    str1 = sprintf('All (%d) : %s = %.2f. p = %.2g', nnz(x1), cc_str, v_all.cc, v_all.p);
    str2 = sprintf('S/S (%d) : %s = %.2f. p = %.2g', nnz(idx_ss_pair), cc_str, v_ss.cc, v_ss.p);
    str3 = sprintf('C/C (%d) : %s = %.2f. p = %.2g', nnz(idx_cc_pair), cc_str, v_cc.cc, v_cc.p);
    str4 = sprintf('S/C (%d) : %s = %.2f. p = %.2g', nnz(idx_sc_pair), cc_str, v_sc.cc, v_sc.p);
    str = {str1, str2, str3, str4};

end

function varargout = getPairDiffs(measure, allMeasures, S, varargin)
    ms_idx = find(strcmp(measure, allMeasures), 1);
    for i = 1:length(varargin)
        if ~iscell(varargin{i})
            varargout{i} = S{ms_idx}.val(varargin{i});
        else
            varargout{i} = cellfun(@(idxs) S{ms_idx}.val(idxs), varargin{i}, 'un', 0); 
        end
    end                    
end

function v = tovector(x)
    v = x(:);
end


function p = getFractionOfDifferencesWithSameSignAs(a, b, ref, pair_flag)
    doAllPairs = exist('pair_flag', 'var') && isequal(pair_flag, 1);
    if doAllPairs
        diffs = bsxfun(@minus, a(:), b(:)');
    else
        diffs = a(:)-b(:);
    end
    n = numel(diffs);
    if ref > 0
        p = nnz( diffs(:) < 0)/n;
    elseif ref < 0
        p = nnz( diffs(:) > 0)/n;
    end

end


function ci = getBootCI(xs_boot, ci_margin_pct)
%     xs_boot = sort(xs_boot);
    assert(ibetween(ci_margin_pct, [1, 100]));
    
    halfMargin = (100 - ci_margin_pct)/2;
    ci = prctile(xs_boot, [halfMargin, 100-halfMargin]);
    
%     idx_central = indmin(abs(xs_boot - mean(xs_boot) ));
%     nMargin = round(length(xs_boot) * ci_margin/2);
%     
%     ci = xs_boot([idx_central - nMargin, idx_central + nMargin]);
    
%     a = 
    
    
end




%{


function pair_stats = compareSimpleComplexPairingStats(fig_name, measure, allMeasures, S, allNSimp, allPairIdxs, fig_id, printStats)
    ms_idx = find(strcmp(measure, allMeasures), 1);

    Wcc_idxs = allPairIdxs{1};
    Wrcc_idxs = allPairIdxs{2};
    Bcc_idxs = allPairIdxs{3};

    nsimp_wcc = allNSimp{1};
    nsimp_wrcc = allNSimp{2};
    nsimp_bcc = allNSimp{3};

    dval_wcc = S{ms_idx}.val(Wcc_idxs);
    dval_bcc = S{ms_idx}.val(Bcc_idxs);
    dval_wrcc = cellfun(@(idxs) S{ms_idx}.val(idxs), Wrcc_idxs, 'un', 0);

%     pair_stats = compareSimpleComplexPairingStats_main(Wcc_pairDiffs, Bcc_pairDiffs, Wrcc_pairDiffs, ...
%                                                        Wcc_nsimp,     Bcc_nsimp,     Wrcc_nsimp,     fig_name, fig_id, printStats);

%        [dval_wcc,      dval_bcc,      dval_wrcc,      nsimp_wcc,     nsimp_bcc,    nsimp_wrcc,      dname,    fig_id, printStats] = deal(...
%          Wcc_pairDiffs, Bcc_pairDiffs, Wrcc_pairDiffs, Wcc_nsimp,     Bcc_nsimp,     Wrcc_nsimp,     fig_name, fig_id, printStats);
%                                                    
    %%
    idx_ss_pair_wcc = nsimp_wcc == 2;
    idx_cc_pair_wcc = nsimp_wcc == 0;
    idx_sc_pair_wcc = nsimp_wcc == 1;

    idx_ss_pair_bcc = nsimp_bcc == 2;
    idx_cc_pair_bcc = nsimp_bcc == 0;
    idx_sc_pair_bcc = nsimp_bcc == 1;
%%
    idx_ss_pair_wrcc = cellfun(@(n) (n == 2), nsimp_wrcc, 'un', 0);
    idx_cc_pair_wrcc = cellfun(@(n) (n == 0), nsimp_wrcc, 'un', 0);
    idx_sc_pair_wrcc = cellfun(@(n) (n == 1), nsimp_wrcc, 'un', 0);
  %%  
    
    dval_ss_wcc = nonnans( dval_wcc(idx_ss_pair_wcc) );
    dval_sc_wcc = nonnans( dval_wcc(idx_sc_pair_wcc) );
    dval_cc_wcc = nonnans( dval_wcc(idx_cc_pair_wcc) );

    dval_ss_bcc = nonnans( dval_bcc(idx_ss_pair_bcc) );
    dval_sc_bcc = nonnans( dval_bcc(idx_sc_pair_bcc) );
    dval_cc_bcc = nonnans( dval_bcc(idx_cc_pair_bcc) );

    justRandomizeSClabels = 1;
    if justRandomizeSClabels
        % use the original within-site distribution, just index according to randomized labels
        dval_ss_wrcc = cellfun(@(tf) nonnans( dval_wcc(tf) ), idx_ss_pair_wrcc, 'un', 0);  
        dval_sc_wrcc = cellfun(@(tf) nonnans( dval_wcc(tf) ), idx_cc_pair_wrcc, 'un', 0);
        dval_cc_wrcc = cellfun(@(tf) nonnans( dval_wcc(tf) ), idx_sc_pair_wrcc, 'un', 0);
        
    else
        % randomize entire cells (values with S/C labels)
        dval_ss_wrcc = cellfun(@(a,tf) nonnans( a(tf) ), dval_wrcc, idx_ss_pair_wrcc, 'un', 0);
        dval_sc_wrcc = cellfun(@(a,tf) nonnans( a(tf) ), dval_wrcc, idx_cc_pair_wrcc, 'un', 0);
        dval_cc_wrcc = cellfun(@(a,tf) nonnans( a(tf) ), dval_wrcc, idx_sc_pair_wrcc, 'un', 0);
        
    end
    

    doStandard = 1;
    doClusterIndex = 1; 
    % standard analysis
    
    utest_str = {};
    clusterIdx_str = {};
    if doStandard
        median_SS = median(dval_ss_wcc); mean_SS = mean(dval_ss_wcc);
        median_SC = median(dval_sc_wcc); mean_SC = mean(dval_sc_wcc);
        median_CC = median(dval_cc_wcc); mean_CC = mean(dval_cc_wcc);        
        
        pu_ss_sc = ranksum(dval_ss_wcc, dval_sc_wcc);
        pu_ss_cc = ranksum(dval_ss_wcc, dval_cc_wcc);
        pu_sc_cc = ranksum(dval_cc_wcc, dval_sc_wcc);

        [~, pt_ss_sc] = ttest2(dval_ss_wcc, dval_sc_wcc);
        [~, pt_ss_cc] = ttest2(dval_ss_wcc, dval_cc_wcc);
        [~, pt_sc_cc] = ttest2(dval_cc_wcc, dval_sc_wcc);

        [~, pks_ss_sc] = kstest2(dval_ss_wcc, dval_sc_wcc);
        [~, pks_ss_cc] = kstest2(dval_ss_wcc, dval_cc_wcc);
        [~, pks_sc_cc] = kstest2(dval_cc_wcc, dval_sc_wcc);

        if printStats
            fprintf('Differences in %s \n', fig_name);
            fprintf('  SS : median %.3g, mean %.3g, rms = %.3g\n', median_SS, mean_SS, rms(dval_ss_wcc));
            fprintf('  SC : median %.3g, mean %.3g, rms = %.3g\n', median_SC, mean_SC, rms(dval_sc_wcc));
            fprintf('  CC : median %.3g, mean %.3g, rms = %.3g\n', median_CC, mean_CC, rms(dval_cc_wcc));
            fprintf('  SS/SC : U-test: p = %6.2g.  T-test p = %6.2g. KS-test: p = %6.g\n', pu_ss_sc, pt_ss_sc, pks_ss_sc);
            fprintf('  SC/CC : U-test: p = %6.2g.  T-test p = %6.2g. KS-test: p = %6.g\n', pu_sc_cc, pt_sc_cc, pks_sc_cc);
            fprintf('  SS/CC : U-test: p = %6.2g.  T-test p = %6.2g. KS-test: p = %6.g\n\n', pu_ss_cc, pt_ss_cc, pks_ss_cc);

            utest_str = {sprintf('Md: SS= %.3g, SC = %.3g, CC = %.3g', median(dval_ss_wcc), median(dval_sc_wcc), median(dval_cc_wcc) ), ...
                         sprintf('p_U: SS/SC=%.3g,   SC/CC = %.3g,   SS/CC = %.3g', pu_ss_sc, pu_sc_cc, pu_ss_cc)};
        end
    end
    
    if doClusterIndex   % clustering index approach.
        clusterIdxUses = 'bcc';
%         clusterIdxUses = 'wrcc_mean';
%         clusterIdxUses = 'wrcc_median';

        wrcc_medians_ss = cellfun(@median, dval_ss_wrcc);
        wrcc_medians_cc = cellfun(@median, dval_cc_wrcc);
        wrcc_medians_sc = cellfun(@median, dval_sc_wrcc);
        wrcc_medians_all= cellfun(@median, dval_wrcc);

        wrcc_means_ss = cellfun(@mean, dval_ss_wrcc);
        wrcc_means_cc = cellfun(@mean, dval_cc_wrcc);
        wrcc_means_sc = cellfun(@mean, dval_sc_wrcc);
        wrcc_means_all= cellfun(@mean, dval_wrcc);

        switch clusterIdxUses
            case 'bcc',       
                ctrl_median_ss  = median(dval_ss_bcc);      ctrl_mean_ss  = mean(dval_ss_bcc);   
                ctrl_median_sc  = median(dval_sc_bcc);      ctrl_mean_sc  = mean(dval_sc_bcc);
                ctrl_median_cc  = median(dval_cc_bcc);      ctrl_mean_cc  = mean(dval_cc_bcc);
                ctrl_median_all = median(dval_bcc);         ctrl_mean_all = mean(dval_bcc);
            case 'wrcc_mean', 
                ctrl_median_ss  = mean(wrcc_medians_ss);    ctrl_mean_ss  = mean(wrcc_means_ss);
                ctrl_median_sc  = mean(wrcc_medians_sc);    ctrl_mean_sc  = mean(wrcc_means_sc);
                ctrl_median_cc  = mean(wrcc_medians_cc);    ctrl_mean_cc  = mean(wrcc_means_cc);
                ctrl_median_all = mean(wrcc_medians_all);   ctrl_mean_all = mean(wrcc_means_all);
            case 'wrcc_median' 
                ctrl_median_ss  = median(wrcc_medians_ss);  ctrl_mean_ss  = median(wrcc_means_ss);
                ctrl_median_sc  = median(wrcc_medians_sc);  ctrl_mean_sc  = median(wrcc_means_sc);
                ctrl_median_cc  = median(wrcc_medians_cc);  ctrl_mean_cc  = median(wrcc_means_cc);        
                ctrl_median_all = median(wrcc_medians_all); ctrl_mean_all = median(wrcc_medians_all);
        end
        
        
        median_bw_ratio_ss = ctrl_median_ss/median(dval_ss_wcc);
        median_bw_ratio_sc = ctrl_median_sc/median(dval_sc_wcc);
        median_bw_ratio_cc = ctrl_median_cc/median(dval_cc_wcc);
        median_bw_ratio_all = ctrl_median_all/median(dval_wcc);

        mean_bw_ratio_ss = ctrl_mean_ss/mean(dval_ss_wcc);
        mean_bw_ratio_sc = ctrl_mean_sc/mean(dval_sc_wcc);
        mean_bw_ratio_cc = ctrl_mean_cc/mean(dval_cc_wcc);
        mean_bw_ratio_all = ctrl_mean_all/mean(dval_wcc);
        
        
        mean_br_ratio_ss = ctrl_mean_ss./wrcc_means_ss;
        mean_br_ratio_sc = ctrl_mean_sc./wrcc_means_sc;
        mean_br_ratio_cc = ctrl_mean_cc./wrcc_means_cc;

        median_br_ratio_ss = ctrl_median_ss./wrcc_medians_ss;
        median_br_ratio_sc = ctrl_median_sc./wrcc_medians_sc;
        median_br_ratio_cc = ctrl_median_cc./wrcc_medians_cc;
        
        
%         p_median_ss_cc = getRandomizedProb(wrcc_medians_ss, wrcc_medians_cc);
%         p_median_ss_sc = getRandomizedProb(wrcc_medians_ss, wrcc_medians_sc);
%         p_median_cc_sc = getRandomizedProb(wrcc_medians_cc, wrcc_medians_sc);
% 
%         p_mean_ss_cc = getRandomizedProb(wrcc_means_ss, wrcc_means_cc);
%         p_mean_ss_sc = getRandomizedProb(wrcc_means_ss, wrcc_means_sc);
%         p_mean_cc_sc = getRandomizedProb(wrcc_means_cc, wrcc_means_sc);
        
        
        dMedRatio_ss_cc_obs = median_bw_ratio_ss-median_bw_ratio_cc; dMedRatio_ss_cc_rand = median_br_ratio_ss-median_br_ratio_cc;
        dMedRatio_ss_sc_obs = median_bw_ratio_ss-median_bw_ratio_sc; dMedRatio_ss_sc_rand = median_br_ratio_ss-median_br_ratio_sc;
        dMedRatio_cc_sc_obs = median_bw_ratio_cc-median_bw_ratio_sc; dMedRatio_cc_sc_rand = median_br_ratio_cc-median_br_ratio_sc;

        dMeanRatio_ss_cc_obs = mean_bw_ratio_ss-mean_bw_ratio_cc; dMeanRatio_ss_cc_rand = mean_br_ratio_ss-mean_br_ratio_cc;
        dMeanRatio_ss_sc_obs = mean_bw_ratio_ss-mean_bw_ratio_sc; dMeanRatio_ss_sc_rand = mean_br_ratio_ss-mean_br_ratio_sc;
        dMeanRatio_cc_sc_obs = mean_bw_ratio_cc-mean_bw_ratio_sc; dMeanRatio_cc_sc_rand = mean_br_ratio_cc-mean_br_ratio_sc;
        %%
        [p_median_ss_cc, med_ss_cc_sgn] = getRandomizedProb(dMedRatio_ss_cc_obs, dMedRatio_ss_cc_rand);
        [p_median_ss_sc, med_ss_sc_sgn] = getRandomizedProb(dMedRatio_ss_sc_obs, dMedRatio_ss_sc_rand);
        [p_median_cc_sc, med_cc_sc_sgn] = getRandomizedProb(dMedRatio_cc_sc_obs, dMedRatio_cc_sc_rand);
        med_ss_cc_str = getSCstr('SS', 'CC', med_ss_cc_sgn, p_median_ss_cc);
        med_ss_sc_str = getSCstr('SS', 'SC', med_ss_sc_sgn, p_median_ss_sc);
        med_cc_sc_str = getSCstr('CC', 'SC', med_cc_sc_sgn, p_median_cc_sc);
        
        [p_mean_ss_cc, mean_ss_cc_sgn] = getRandomizedProb(dMeanRatio_ss_cc_obs, dMeanRatio_ss_cc_rand);
        [p_mean_ss_sc, mean_ss_sc_sgn] = getRandomizedProb(dMeanRatio_ss_sc_obs, dMeanRatio_ss_sc_rand);
        [p_mean_cc_sc, mean_cc_sc_sgn] = getRandomizedProb(dMeanRatio_cc_sc_obs, dMeanRatio_cc_sc_rand);
        mean_ss_cc_str = getSCstr('SS', 'CC', mean_ss_cc_sgn, p_mean_ss_cc);
        mean_ss_sc_str = getSCstr('SS', 'SC', mean_ss_sc_sgn, p_mean_ss_sc);
        mean_cc_sc_str = getSCstr('CC', 'SC', mean_cc_sc_sgn, p_mean_cc_sc);
        
        %%
        medianRatio_str = sprintf('  Median ratios: All pairs: %.3f.  SS: %.3f   SC: %.3f   CC: %.3f', median_bw_ratio_all, median_bw_ratio_ss, median_bw_ratio_sc, median_bw_ratio_cc);
        medianProbs_str = sprintf('  Median probs: %s  %s  %s', med_ss_cc_str, med_ss_sc_str, med_cc_sc_str);
        meanRatio_str   = sprintf('  Mean ratios:  All pairs: %.3f.  SS: %.3f   SC: %.3f   CC: %.3f', mean_bw_ratio_all, mean_bw_ratio_ss, mean_bw_ratio_sc, mean_bw_ratio_cc);
        meanProbs_str   = sprintf('  Mean probs: %s  %s  %s', mean_ss_cc_str, mean_ss_sc_str, mean_cc_sc_str);
        clusterIdx_str = {medianRatio_str, medianProbs_str, meanRatio_str, meanProbs_str};
            
        if printStats
            fprintf(' * Clustering indices\n');
            fprintf('%s\n', medianRatio_str);
            fprintf('%s\n', medianProbs_str);
            fprintf('%s\n', meanRatio_str);
            fprintf('%s\n', meanProbs_str);            
        end
        pair_stats = struct('medianRatio_all', median_bw_ratio_all, 'medianRatio_ss', median_bw_ratio_ss, 'medianRatio_cc', median_bw_ratio_cc, 'medianRatio_sc', median_bw_ratio_sc, ...
                            'medianProb_ss_cc', p_median_ss_cc, 'medianProb_ss_sc', p_median_ss_sc, 'medianProb_cc_sc', p_median_cc_sc, ...
                            'meanRatio_all', mean_bw_ratio_all, 'meanRatio_ss', mean_bw_ratio_ss,  'meanRatio_cc', mean_bw_ratio_cc, 'meanRatio_sc', mean_bw_ratio_sc, ...
                            'meanProb_ss_cc', p_mean_ss_cc, 'meanProb_ss_sc', p_mean_ss_sc, 'meanProb_cc_sc', p_mean_cc_sc ...
                            ...'median_SS', median_SS, 'median_SC', median_SC, 'median_CC', median_CC, ...
                            ...'mean_SS', mean_SS, 'mean_SC', mean_SC, 'mean_CC', mean_CC, ...
                            ...'pU_ss_sc', pu_ss_sc, 'pU_ss_cc', pu_ss_cc, 'pU_sc_cc', pu_sc_cc, ...
                            ...'pT_ss_sc', pt_ss_sc, 'pT_ss_cc', pt_ss_cc, 'pT_sc_cc', pt_sc_cc, ...
                            ...'pKS_ss_sc', pks_ss_sc, 'pKS_ss_cc', pks_ss_cc, 'pKS_sc_cc', pks_sc_cc); 
                            );
                            
        3;
%          wrcc_medians_ss = cellfun(@median, dval_ss_wrcc);
%         wrcc_medians_cc = cellfun(@median, dval_cc_wrcc);
%         wrcc_medians_sc = cellfun(@median, dval_sc_wrcc);
%         wrcc_medians_all= cellfun(@median, dval_wrcc);
% 
%         wrcc_means_ss = cellfun(@mean, dval_ss_wrcc);
%         wrcc_means_cc = cellfun(@mean, dval_cc_wrcc);
%         wrcc_means_sc = cellfun(@mean, dval_sc_wrcc);
%         wrcc_means_all= cellfun(@mean, dval_wrcc);
    end    
    3;
    
    if exist('fig_id', 'var') && ~isnan(fig_id);
        %%
        figure(fig_id); clf;
        if doStandard && ~doClusterIndex
            h = hist2({dval_ss_wcc, dval_sc_wcc, dval_cc_wcc}, 15, 'line', 'norm');
            set(h(1), 'color', 'b');
            set(h(2), 'color', 'm');
            set(h(3), 'color', 'r');
            title(fig_name); 
            legend('SS', 'SC', 'CC');                 

            set(h, 'linewidth', 2)
            gratingStr = curGratingType('');
            title([{sprintf('%s (%s gratings)', fig_name, gratingStr)}, utest_str, clusterIdx_str], 'fontsize', 10) 
            ylims = ylim;
    %         if pval_U < .05
            line(median(dval_ss_wcc)*[1, 1],  [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'b', 'linestyle', '--');
            line(median(dval_sc_wcc)*[1, 1],  [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'm', 'linestyle', '--');
            line(median(dval_cc_wcc)*[1, 1],  [ylims(1), ylims(1)+diff(ylims)/2], 'color', 'r', 'linestyle', '--');

            line(mean(dval_ss_wcc)*[1, 1],  [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'b', 'linestyle', ':');
            line(mean(dval_sc_wcc)*[1, 1],  [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'm', 'linestyle', ':');
            line(mean(dval_cc_wcc)*[1, 1],  [ylims(1)+diff(ylims)/2, ylims(2)], 'color', 'r', 'linestyle', ':');
        else
%%
            clf;
            colors_C = {'g', 'b', 'r'};
            gratingStr = curGratingType('');        
            for plot_i = 1:2
                switch plot_i
                    case 1, plotData_rnd = {dMedRatio_ss_cc_rand, dMedRatio_ss_sc_rand, dMedRatio_cc_sc_rand};
                            plotData_obs = [dMedRatio_ss_cc_obs,  dMedRatio_ss_sc_obs,  dMedRatio_cc_sc_obs];
                            title_str = {['\bf Differences in ' fig_name ' (Cluster Indices) (' gratingStr ') \rm' ], medianRatio_str,medianProbs_str};
                            
                    case 2, plotData_rnd = {dMeanRatio_ss_cc_rand, dMeanRatio_ss_sc_rand, dMeanRatio_cc_sc_rand};
                            plotData_obs = [dMeanRatio_ss_cc_obs,  dMeanRatio_ss_sc_obs,  dMeanRatio_cc_sc_obs];
                            title_str = {meanRatio_str, meanProbs_str};
                end
                subplotGap(2,1,plot_i);

                h = hist2(plotData_rnd, 30, 'line', 'norm');
                xlims = xlim;            
                set(h, {'color'}, colors_C');
                set(h, 'linewidth', 2);
                if plot_i == 1
                    legend('SS-CC', 'SS-SC', 'CC-SC');        
                end                
                
                uistack(h(1), 'top');
                title(title_str);

                ylims = ylim;
                for j = 1:3
                    line(plotData_obs(j)*[1, 1],  ylims, 'color', colors_C{j}, 'linestyle', '-');
                end
                xlims = lims([xlims(:)', plotData_obs], .05);
                xlim(xlims);
                
                3;
            end
        end
        
    end    

                                                   
                                                   
                                                   

end

%}

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
        spf_cellF1oDC = [allSpfUnitStats_si.F1oDC];
        spf_pairF1oDC_Wcc = spf_cellF1oDC(Wcc_pairs_spf);
        spf_pairF1oDC_Bcc = spf_cellF1oDC(Bcc_pairs_spf);
        
        spf_nSimp_Wcc = sum(spf_pairF1oDC_Wcc > 1, 2);
        spf_nSimp_Bcc = sum(spf_pairF1oDC_Bcc > 1, 2);
    else
        [Wcc_ss_pairIdxs, Wrcc_ss_pairIdxs, Bcc_ss_pairIdxs, spf_cellF1oDC, spf_pairF1oDC_Wcc] = deal([]);
    end
    
%}

            %{
            subplotGap(1,3,1); h1 = plot(Dori_pref_MU(idx_mark0), spkAmps(idx_mark0), sym1, ...
                                         Dori_pref_MU(idx_mark1), spkAmps(idx_mark1), sym2, ...
                                         Dori_pref_MU(idx_mark2), spkAmps(idx_mark2), sym3);              
            xlabel(' '); %xlabel('Difference from multiunits, degrees'); 
            ylabel('spike amplitude, mV');        
            xlim([0 90]); set(gca, 'xtick', [0:15:90]);
            pvalU_spkAmp = ranksum(spkAmps(~idx_mark0), spkAmps( idx_mark0) );
            [h, pvalT_spkAmp] = ttest2(spkAmps(~idx_mark0), spkAmps( idx_mark0) );
            title({' ', sprintf('p_U = %.2g, p_T = %.2g', pvalU_spkAmp, pvalT_spkAmp)});
            med_amps_marked   = nanmedian( spkAmps( ~idx_mark0) );
            med_amps_unmarked = nanmedian( spkAmps( idx_mark0) );
            fprintf('Median amp of of outliers: %.2f. Median amp of non-outliers: %.2f\n', ...
                med_amps_marked, med_amps_unmarked);
            
            subplotGap(1,3,2); h2 = plot(Dori_pref_MU(idx_mark0),w_ori_global(idx_mark0), sym1, ...
                                         Dori_pref_MU(idx_mark1),w_ori_global(idx_mark1), sym2, ...
                                         Dori_pref_MU(idx_mark2),w_ori_global(idx_mark2), sym3);             
            xlabel('Difference from multiunits, degrees'); 
            ylabel('w_{ORI}^{Global}');
            xlim([0 90]); set(gca, 'xtick', [0:15:90]);                                  
            pvalU_oriGlobal = ranksum(w_ori_global(~idx_mark0), w_ori_global( idx_mark0) );
            [h, pvalT_oriGlobal] = ttest2(w_ori_global(~idx_mark0), w_ori_global( idx_mark0) );
            title({gratingStr, sprintf('p_U = %.2g, p_T = %.2g', pvalU_oriGlobal, pvalT_oriGlobal)});
            med_w_ori_glob_marked   = nanmedian( w_ori_global( ~idx_mark0) );
            med_w_ori_glob_unmarked = nanmedian( w_ori_global( idx_mark0) );
            fprintf('Median w_ori_global of of outliers: %.2f. Median w_ori_global of non-outliers: %.2f\n',...
                med_w_ori_glob_marked, med_w_ori_glob_unmarked);
            
%             w_ori_local = F1oDCs;
            subplotGap(1,3,3); h3 = plot(Dori_pref_MU(idx_mark0),w_ori_local(idx_mark0), sym1, ...
                                         Dori_pref_MU(idx_mark1),w_ori_local(idx_mark1), sym2, ...
                                         Dori_pref_MU(idx_mark2),w_ori_local(idx_mark2), sym3);         
            xlabel(' '); %xlabel('Difference from multiunits, degrees'); 
            ylabel('w_{ORI}^{Local}');
            xlim([0 90]); set(gca, 'xtick', [0:15:90]);                                  
            pvalU_oriLocal = ranksum( nonnans( w_ori_local(~idx_mark0)), nonnans( w_ori_local( idx_mark0) ));
            [h, pvalT_oriLocal] = ttest2(w_ori_local(~idx_mark0), w_ori_local( idx_mark0) );
            title({' ', sprintf('p_U = %.2g, p_T = %.2g', pvalU_oriLocal, pvalT_oriLocal)});            
            med_w_ori_loc_marked   = nanmedian( w_ori_local( ~idx_mark0) );
            med_w_ori_loc_unmarked = nanmedian( w_ori_local( idx_mark0) );
            fprintf('Median amp of of outliers: %.2f. Median amp of non-outliers: %.2f\n', ...
                med_w_ori_loc_marked, med_w_ori_loc_unmarked);
            3;
            
%             subplotGap(1,4,); h3 = plot(Dori_pref_MU(~idx_mark),w_ori_local(~idx_mark), 'b+', ...
%                                       Dori_pref_MU( idx_mark),w_ori_local( idx_mark), 'r+');         
%             xlabel(' '); %xlabel('Difference from multiunits, degrees'); 
%             ylabel('w_{ORI}^{Local}');
%             xlim([0 90]); set(gca, 'xtick', [0:15:90]);                                  
%             pvalU_oriLocal = ranksum( nonnans( w_ori_local(~idx_mark)), nonnans( w_ori_local( idx_mark) ));
%             [h, pvalT_oriLocal] = ttest2(w_ori_local(~idx_mark), w_ori_local( idx_mark) );
%             title({' ', sprintf('p_U = %.2g, p_T = %.2g', pvalU_oriLocal, pvalT_oriLocal)});            
%             med_w_ori_loc_marked   = nanmedian( w_ori_local( idx_mark) );
%             med_w_ori_loc_unmarked = nanmedian( w_ori_local( ~idx_mark) );
%             fprintf('Median amp of of outliers: %.2f. Median amp of non-outliers: %.2f\n', ...
%                 med_w_ori_loc_marked, med_w_ori_loc_unmarked);
            
            set([h1, h2, h3], 'markersize',5);
            %}


        %{
        %%%% SPATIAL FREQUENCY MEAUSURES
        %% Preferred spatial frequency
        % A. Does S/C affect preferred spatial frequency?
        s_i = 1; stat_name = 'Preferred Spatial Frequency';  
        cell_stats(s_i) = compareSimpleComplexCellStats(stat_name, [allSpfCellStats.f_opt], spf_cell_simple_tf, cell_stat_figIds(s_i), printStats);
        pairing_stats(s_i) = compareSimpleComplexPairingStats(stat_name, 'D_spf_pref', measures_spf, S_spf, spf_allNSimp, spf_pairTypeIdxs, pair_stat_figIds(s_i), printStats);
        
        
        %% Spatial frequency tuning width
        % A. Does S/C affect single cell tuning?   % B. Does SS/SC/CC affect differences in spatial frequency tuning width?
        s_i = 2; stat_name = 'Spatial Frequency Tuning Width';  [cell_stat_labels{s_i}, pairing_labels{s_i}] = deal(stat_name);        
        cell_stats(s_i) = compareSimpleComplexCellStats(stat_name, [allSpfCellStats.w_spf], spf_cell_simple_tf, cell_stat_figIds(s_i), printStats);        
        pairing_stats(s_i) = compareSimpleComplexPairingStats(stat_name, 'Dw_spf', measures_spf, S_spf, spf_allNSimp, spf_pairTypeIdxs, pair_stat_figIds(s_i), printStats);
        
        %%%% ORIENTATION MEAUSURES    
        %% Preferred Orientation 
        % Single cell Preferred Orientation;  % Differences in preferred Orientation
        s_i = 3; stat_name = 'Preferred Orientation';  [cell_stat_labels{s_i}, pairing_labels{s_i}] = deal(stat_name);        
        cell_stats(s_i) = compareSimpleComplexCellStats(stat_name, [allOriCellStats.ori_pref_deg], ori_cell_simple_tf,  cell_stat_figIds(s_i), printStats);
        pairing_stats(s_i) = compareSimpleComplexPairingStats(stat_name, 'D_ori_pref',  measures_ori, S_ori, ori_allNSimp, ori_pairTypeIdxs, pair_stat_figIds(s_i), printStats);


        %% Orientation Tuning Width
        % Global orientation tuning width?
        s_i = 4; stat_name = 'Global Orientation Width';  [cell_stat_labels{s_i}, pairing_labels{s_i}] = deal(stat_name);        
        cell_stats(s_i) = compareSimpleComplexCellStats( stat_name, [allOriCellStats.w_ori_global], ori_cell_simple_tf, cell_stat_figIds(s_i), printStats);                        
        pairing_stats(s_i) = compareSimpleComplexPairingStats(stat_name, 'Dw_ori_glob', measures_ori, S_ori, ori_allNSimp, ori_pairTypeIdxs, pair_stat_figIds(s_i), printStats);
        
        % Local orientation tuning width?
        s_i = 5; stat_name = 'Local Orientation Width';  [cell_stat_labels{s_i}, pairing_labels{s_i}] = deal(stat_name);        
        cell_stats(s_i) =  compareSimpleComplexCellStats(stat_name , [allOriCellStats.w_ori_local],  ori_cell_simple_tf, cell_stat_figIds(s_i), printStats);
        pairing_stats(s_i) = compareSimpleComplexPairingStats(stat_name, 'Dw_ori_loc', measures_ori, S_ori, ori_allNSimp, ori_pairTypeIdxs, pair_stat_figIds(s_i), printStats);


        if strcmp(gratingType, 'drifting')
            %% DSI
            % Single cell DSI
            s_i = 6; stat_name = 'Direction Selectivity Index';  [cell_stat_labels{s_i}, pairing_labels{s_i}] = deal(stat_name);        
            cell_stats(s_i) =  compareSimpleComplexCellStats(stat_name, [allOriCellStats.DSI_global], ori_cell_simple_tf, cell_stat_figIds(s_i), printStats);            
            pairing_stats(s_i) = compareSimpleComplexPairingStats(stat_name, 'D_dsi_glob_si', measures_ori, S_ori, ori_allNSimp, ori_pairTypeIdxs, pair_stat_figIds(s_i), printStats);                
        end
        %}


%{
from fig 4a -- for comparison with Albus 1975
      binVals_cc_all_frac = binVals_cc_all/sum(binVals_cc_all)*100;
            h_bar4a = bar(binC, binVals_cc_all_frac, 1, 'stacked');
%             h_bar4a = bar(binC, [binVals_cc_norm, (binVals_cc_1outlier + binVals_cc_2outliers)], 1, 'stacked');
            set(h_bar4a(1), 'facecolor', 'w');
%             set(h_bar4a(2), 'facecolor', bar_out);
            
            xlim(binE([1, end]));
%             set(gca, 'xtick', [0:15:90]);
            set(gca, 'xtick', [0:30:90]);
%             title(add_fg('Pref. Orientation : Pairwise differences'), 'fontsize', title_fsize);
            title({'dOri', '(drifting)'});
            xlabel('d(Pref Ori)');
            ylabel('% of cell-cell pairs');
            ylim([0 70]);
%             legend({'Typical cells', 'Outlier cells'})
%}


%{
from supp figure 5 -- outlier plots using two kinds of outliers
for var_i = 1:nVariables
%                 subplotGap(2,nPlots/2,sub_i);
    Y_i = Y_vals{var_i};
    sub_i = find(idx_plot == var_i);

    pvalU_i12 = ranksum_nonnans(Y_i(idx_mark0), Y_i( ~idx_mark0) );
    pvalU_i1 = ranksum_nonnans(Y_i(idx_mark0), Y_i( idx_mark1) );
    pvalU_i2 = ranksum_nonnans(Y_i(idx_mark0), Y_i( idx_mark2) );

    [h, pvalT_i12] = ttest2(Y_i(idx_mark0), Y_i(~idx_mark0) );
    [h, pvalT_i1] = ttest2(Y_i(idx_mark0), Y_i(idx_mark1) );
    [h, pvalT_i2] = ttest2(Y_i(idx_mark0), Y_i(idx_mark2) );

    title_str1 = sprintf('(1) p_U = %.2g, p_T = %.2g', pvalU_i1, pvalT_i1);
    title_str2 = sprintf('(2) p_U = %.2g, p_T = %.2g', pvalU_i2, pvalT_i2);
    title_str12 = sprintf('(1+2) p_U = %.2g, p_T = %.2g', pvalU_i12, pvalT_i12);
    title_str_now = sprintf('p_U = %.2g', pvalU_i12);

    gratingType_str = sprintf('(%s gratings)', titleCase(gratingType) );

%                 title({title_str0, title_str1, title_str2, title_str12});

    med_amps_0 = nanmedian( Y_i( idx_mark0) );
    med_amps_1 = nanmedian( Y_i( idx_mark1) );
    med_amps_2 = nanmedian( Y_i( idx_mark2) );
    med_amps_12 = nanmedian( Y_i( ~idx_mark0) );
    if separateWideOri
        fprintf('%s (Median) : Outliers #1 (<35): %.2f. Outliers #2 (>35): %.2f. Outliers(1+2): %.2f.  Non-outliers: %.2f...  pU(1): %.2g. pU(2): %.2g. pU(1+2): %.2g \n', ...
            Y_name_short{var_i}, med_amps_1, med_amps_2, med_amps_12, med_amps_0, pvalU_i1, pvalU_i2, pvalU_i12);
    else
        fprintf('%s (Median) : Outliers: %.2f.  Non-outliers: %.2f...  pU(1): pU = %.2g \n', ...
            Y_name_short{var_i}, med_amps_12,  med_amps_0, pvalU_i12);
    end


    if ~isempty(sub_i)
            subplotGap(1,nPlots,sub_i);


%                 h5{sub_i} = plot(Dori_pref_MU(idx_mark0), Y_i(idx_mark0), sym1, ...
%                                  Dori_pref_MU(idx_mark1), Y_i(idx_mark1), sym2, ...
%                                  Dori_pref_MU(idx_mark2), Y_i(idx_mark2), sym3);  %#ok<AGROW>
            h5{sub_i} = plot(Dori_pref_MU(idx_mark0), Y_i(idx_mark0), sym1, ...
                             Dori_pref_MU(idx_mark_orig), Y_i(idx_mark_orig), sym2 );

        if sub_i == 1 % floor((nPlots+1)/2)
            xlabel('Difference from multiunits, degrees');
            title_str0 = gratingStr;
        else
            xlabel(' ');
            title_str0 = ' ';
        end
        ylabel(Y_labels{var_i});

        xlim([0 90]); set(gca, 'xtick', [0:15:90]);
        title({Y_name_short{var_i}, gratingType_str, title_str_now});

        if strcmp(Y_name_short{var_i}, 'PtP width')
            ylim([0.1, 0.7]);
        end

        if sub_i == nPlots
            %  legend({'Other cells', 'Outliers #1', 'Outliers #2'}, 'location', 'SE')
            legend({'Typical cells', 'Outliers'}, 'location', 'NE')
        end
        if strcmp(Y_labels{sub_i}, 'w_{ORI}^{Global}')
            % line([45 45], [0 50], 'linestyle', ':', 'color', 'k');
            % line([45 90], [35 35], 'linestyle', ':', 'color', 'k');
        end
        if addSubplotLetters
            addSubplotLetter(1, nPlots, 1, sub_i, char('A'+(grating_offset/2)*nPlots+sub_i-1));
        end
    end



end
%}


%{
if isOri(xi) && ori_separateNormAndOutliers
                        if separateLowerPctile
                            plot(X(idx_noOutliers   & ~idx_bothLarge), Y(idx_noOutliers   & ~idx_bothLarge), 'bo', 'markersize', 3);
                            plot(X(idx_withOutliers & ~idx_bothLarge), Y(idx_withOutliers & ~idx_bothLarge), 'ro', 'markersize', 3);
                            plot(X(idx_noOutliers   & idx_bothLarge),  Y(idx_noOutliers   & idx_bothLarge), 'b+', 'markersize', 3);
                            plot(X(idx_withOutliers & idx_bothLarge),  Y(idx_withOutliers & idx_bothLarge), 'r+', 'markersize', 3);
                        else
                            plot(X(idx_noOutliers), Y(idx_noOutliers), 'bo', 'markersize', 3);
                            plot(X(idx_withOutliers), Y(idx_withOutliers), 'ro', 'markersize', 3);
                        end
                    else
                        if separateLowerPctile
                            plot(X(~idx_bothLarge), Y(~idx_bothLarge), 'bo', 'markersize', 3);
                            plot(X(idx_bothLarge), Y(idx_bothLarge), 'b+', 'markersize', 3);
                        else
                            plot(X, Y, 'bo', 'markersize', 3);
                        end
                        
                    end
%}


%{

  switch clusterIdxUses
                case 'bcc',       
                    ctrl_stat_ss  = bcc_stat_ss;      
                    ctrl_stat_sc  = bcc_stat_sc;      
                    ctrl_stat_cc  = bcc_stat_cc;      
                    ctrl_stat_all = bcc_stat_all;         
                case 'brcc',       
                    ctrl_stat_ss = cellfun(op, dval_ss_brcc);
                    ctrl_stat_cc = cellfun(op, dval_cc_brcc);
                    ctrl_stat_sc = cellfun(op, dval_sc_brcc);
                    ctrl_stat_all= cellfun(op, dval_brcc);

                case 'wrcc_mean', 
                    ctrl_stat_ss  = mean(wrcc_stat_ss);    
                    ctrl_stat_sc  = mean(wrcc_stat_sc);    
                    ctrl_stat_cc  = mean(wrcc_stat_cc);    
                    ctrl_stat_all = mean(wrcc_stat_all);   
                case 'wrcc_median' 
                    ctrl_stat_ss  = median(wrcc_stat_ss);  
                    ctrl_stat_sc  = median(wrcc_stat_sc);  
                    ctrl_stat_cc  = median(wrcc_stat_cc);  
                    ctrl_stat_all = median(wrcc_stat_all); 
            end


%}
