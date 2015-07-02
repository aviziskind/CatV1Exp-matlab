function printPhaseTuningComparisons(suppressOutput_flag, criteria_input)
%     doFigs = [1:12];
%     doFigs = [12];
    curCmpType('phase');

    doPrintStats = 1        && 1;
    doPrintSigTests = 1     && 1;
    doPlots = 1             && 1;
    
        
    
    useShuffledPhaseControl = 1;

    plotInColor = 0;
    useCorrectProbs = 1;
    
    doCCvsTetrodeDepthPlots = 0;
    
    saveControlDistribs = 0;
        
    gratingType = curGratingType('');
    if useShuffledPhaseControl
        curPairTypes('Wcc', 'Wscc');            
%         curPairTypes('Wcc', 'Wcm', 'Bcc', 'Bcm', 'Bmm', 'Wscc')
    else
        curPairTypes('Wcc', 'Wrcc', 'Bcc');                   
%         curPairTypes('Wcc', 'Wrcc', 'Bcc');                   
    end
    [sorting_str, cmp_str] = curStatus;
    
    if nargin == 0 || isempty(suppressOutput_flag)
        suppressOutput_flag = 0;
    end
    
    showWorking = all(suppressOutput_flag == 0);
    if ~showWorking
        [doPrintStats, doPrintSigTests, doPlots] = deal(0);
        doPrintStats = any(suppressOutput_flag == 1);
        doPrintSigTests = any(suppressOutput_flag == 2);
    end
    
    
%     curPairTypes('Wcc', 'Wrcc', 'Bcc');            
    
    ospDatafile = getFileName('osps');    
    pairDatafile= getFileName('pairs');
    cmpDatafile = getFileName('comparisons');    
%     statsDatafile = getFileName('controls');
        
    fprintf('Loading ... '); tic;    
    fprintf(' cells : %s \n', removePathFromFilename(ospDatafile) );
    S1 = load(ospDatafile);        
    allCells = S1.allCells;
    clear S1;
    nUnits = length(allCells);    

    fprintf(' pairs : %s \n', removePathFromFilename(pairDatafile) );
    S2 = load(pairDatafile);            
    fprintf(' cmps  : %s ', removePathFromFilename(cmpDatafile) );
    S3 = load(cmpDatafile); 
    [pairData, S, pairTypes, measures, locations, cmp_opts] = deal(S3.pairData, S3.allStatsC, S3.pairTypes, S3.measures, S3.locations, S3.opt);
    [pairIdxs, pairIdxList, idxMtx] = useSubsetOfPairIdxs(S2, pairTypes, nUnits);    
    randInfo = S2.Wrcc_details;
    
    clear S2;
    clear S3;
        
%     measures = measures(~strncmp(measures, 'MID', 3));
    
    nMeasures = length(measures);
    nLocations = length(locations);
    
%     allPairIdxs = {Wcc_pairIdxs, Wrcc_pairIdxs, Bcc_pairIdxs, Wcm_pairIdxs, Wrcm_pairIdxs, Bcm_pairIdxs};
%     allPairIdxs_list = {Wcc_pairIdxs, unique([Wrcc_pairIdxs{:}]), Bcc_pairIdxs, Wcm_pairIdxs, unique([Wrcm_pairIdxs{:}]), Bcm_pairIdxs};
%     allPairTypes = {'Wcc', 'Wrcc', 'Bcc',  'Wcm', 'Wrcm', 'Bcm'};
%     pairTypes = pairTypes(ord(cellfun(@(s) find(strcmp(s, allPairTypes)), pairTypes)));
    
    % 
    limitToNPermutes = [10000];

    Wcc_pairIdxs_M  = pairIdxs{ find(strcmp(pairTypes, 'Wcc'), 1) };
    Wcc_pairs = ind2subV([nUnits, nUnits], Wcc_pairIdxs_M);
    Wcc_pairIdxs = idxMtx(Wcc_pairIdxs_M);    

    frac_pct = @(i,n) sprintf('%3d/%3d (%.1f%%)', i, n, i/n*100);
    toc;
    
    if ~useShuffledPhaseControl

        if any(strcmp(pairTypes, 'Wrcc'))
            Wrcc_pairIdxs_M = pairIdxs{ find(strcmp(pairTypes, 'Wrcc'), 1) };
            if ~isempty(limitToNPermutes)
                Wrcc_pairIdxs_M = Wrcc_pairIdxs_M(1:limitToNPermutes); %%#ok<MSNU> #ok<BDSCI>
            end        
            Wrcc_pairs = cellfun(@(idxs) ind2subV([nUnits, nUnits], idxs), Wrcc_pairIdxs_M, 'un', 0);
            Wrcc_pairIdxs = cellfun(@(idxs) idxMtx(idxs),  Wrcc_pairIdxs_M, 'un', 0);
        end
        
    else
        
        if ~isempty(limitToNPermutes)
            for i = 1:numel(S)
                if isfield(S{i}, 'shuffVal') &&  (size(S{i}.shuffVal, 2) > limitToNPermutes)
                    S{i}.shuffVal = S{i}.shuffVal(:,1:limitToNPermutes); %%#ok<MSNU> #ok<BDSCI>
                end
            end            
        end
        
    end
    
    if any(strcmp(pairTypes, 'Bcc'))        
        Bcc_pairIdxs_M  = pairIdxs{ find(strcmp(pairTypes, 'Bcc'), 1) };
        Bcc_pairs = ind2subV([nUnits, nUnits], Bcc_pairIdxs_M);
        Bcc_pairIdxs = idxMtx(Bcc_pairIdxs_M);
    end
    W_pairIdxs = find( pairData.Gids(:,1) == pairData.Gids(:,2) );
    B_pairIdxs = find( pairData.Gids(:,1) ~= pairData.Gids(:,2) );
    
    useGlobalCriteria = 0;
            
    useSpecificCriteria = 1;
    haveCriteriaInput = nargin >= 2 && ~isempty(criteria_input);
%     criteria.MID_cc.minF1oDC = 0;
    if haveCriteriaInput
        criteria = criteria_input;
        
    elseif useSpecificCriteria
        %%
                
            %%
                    criteria = struct;
%                     nPhase_S = 4;                                         

                    % 1. CELL SELECTION CRITERIA
%                     criteria.n_phases = struct('op', @eq, 'value', 8);
%                     criteria.n_phases = struct('op', @eq, 'value', 60);
                    GLF_overlap = 0; criteria.GLF_overlap = struct('op', @gt, 'value', GLF_overlap);
%                     criteria.negAmps_overlap = struct('op', @gt, 'value', amps_overlap);
                    minID = 10;   criteria.minID = struct('op', @gt, 'value', minID);

                    % PAIR SELECTION CRITERIA
%                     minFrac = 0.25;     criteria.loc_minFracOfMaxes     = minFrac;                    
%                     minF1oDC_cmp = 0;   criteria.loc_minF1oDC_cmp       = minF1oDC_cmp;

                    % simple/complex pairing
                    SCpref_type = 'maxR_avP_sm';
%                     SCpref_type = 'maxP_sm';
                    SC_f1odc_th = 0.25;
                    criteria.(['SCtype_pref_' SCpref_type]) = struct('op', @eq, 'value', 2);
                    criteria.(['F1oDC_' SCpref_type '_maxJackStd']) = struct('op', @lt, 'value', SC_f1odc_th);
                    
                    
%                     criteria.SCtype_pref = struct('op', @eq, 'value', 2);
                    
                    % spatial frequency tuning similarity 
%                     criteria.D_spf_pref = struct('op', @gt, 'value', 1);
%                     criteria.Dw_spf = struct('op', @gt, 'value', 1);

                    % VALUE SELECTION CRITERIA - odd/even reproducibility
                    if strcmp(gratingType, 'flashed')
                        criteria.MID_cc.min_rsqr_oe = 0.25;            
                        criteria.MID_fit_cc.min_rsqr_oe = 0.25;    
                        criteria.MID_fit_cc.min_rsqr_fit = 0.25;        
                    end


                    % VALUE SELECTION CRITERIA - jackknive std err limits
                    doJackStdLimits = 0;
                    if doJackStdLimits
                        jackStd_cc_max = .2;
                        jackStd_dphi_max = 15;
                        jackStd_RF_cc = .1;

                        jackStdCrit.cc     = jackStd_cc_max;
                        jackStdCrit.dphi   = jackStd_dphi_max; 
                        jackStdCrit.MID_cc       = jackStd_RF_cc; 
                        jackStdCrit.MID_ovlp_cc  = jackStd_RF_cc; 
                        jackStdCrit.MID_fit_cc   = jackStd_RF_cc; 
                        jackStdCrit.STA_cc       = jackStd_RF_cc; 
                    else
                        jackStdCrit = struct;
                    end
                    
    else
        %%
        criteria = struct;
    end
    
    if useGlobalCriteria
        
%         if useCriteria_F1oDC
%             pair_ok = pair_ok & (pairData.minF1oDC_pref > minF1oDC);
%             criteria_str = [criteria_str sprintf('F1/DC > %.1f', minF1oDC)];
%         end
%         if useCriteria_minRsqr
%             pair_ok = pair_ok & (pairData.minF1oDC_pref > minF1oDC);
%             criteria_str = [criteria_str sprintf('r^2 > %.1f', minRsqr)];
%         end            
    %         pair_ok = pairData.nphases == 40;
%         pair_ok = pairData.minF1oDC_pref > 1;   
%         pair_ok = pairData.maxF1oDC_pref < 1;
        

        
        wcc_ok = pair_ok(Wcc_pairIdxs);
        if ~useShuffledPhaseControl
            bcc_ok = pair_ok(Bcc_pairIdxs);
        end        
        
    else
        
        criteria_str = '';
        wcc_ok = true(size(Wcc_pairIdxs));        
        
        if ~useShuffledPhaseControl
            bcc_ok = true(size(Bcc_pairIdxs));        
        end
    end
    
    
    
    
    % Get Ori_si, Ori_ss, and Spf_si stats
    cells_tf = [allCells.cellId] > 0;
        
    cells_idx = find(cells_tf);    
    allCells = allCells( cells_tf);

    
    

    printNCellsPerSiteStats = 1 && ~any(suppressOutput_flag);
        
        %%        
    allCellGids = [allCells.Gid];
    [uGid, gidCount] = uniqueCount(allCellGids);
    nSites = length(uGid);
    nCells = length(allCellGids);
    nCells_per_site = mean(gidCount);
    nSites_1cell = nnz(gidCount == 1);
    s1 = sprintf('**** Phase tuning (%s gratings) ***', curGratingType(''));
    s2 = sprintf('Recorded from %d cells from %d sites (mean %.2f cells per site).\n', nCells, nSites, nCells_per_site);
    s3 = sprintf('%d sites had 1 cell. Thus %d cells yielded a total of %d pairs\n',nSites_1cell, nCells-nSites_1cell, length(Wcc_pairIdxs)  );        

    if showWorking    
%         fprintf('%s\n%s\n%s\n', s1, s2, s3)
    end

    cellsPerSite_comments = {s1, s2, s3};        
        
    
    
    
    
%     allOriStats_si = nestedFields(allOriCells, {'stats', 'tuningStats', 'oriStats_si'});
%     allSpfStats_si = nestedFields(allSpfCells, {'stats', 'tuningStats', 'spfStats_si'});
%     allOriSpkFeatures = [allOriCells.spkFeatures];
    
    allGids = [allCells.Gid];

    % Double check that the renumbering is all ok.
    
    
    3;
    PhaseTCMeasures = {'cc', 'rho', 'dphi', 'dF1'};
    RF_Measures = {'STA_cc', 'MID_cc', 'MID_ovlp_cc', 'MID_fit_cc', 'dPh_rel'};


    if ~useShuffledPhaseControl
        Npermutes = length(Wrcc_pairIdxs);
    else        
        Npermutes = cmp_opts.nPhaseShuffles;
    end
    
    
    
    %%
    if ~useShuffledPhaseControl    
        pairTypeIdxs = {Wcc_pairIdxs(wcc_ok), Bcc_pairIdxs(bcc_ok)};
        pairs_C      = {Wcc_pairs(wcc_ok, :), Bcc_pairs(bcc_ok,:)};
        pairTypes_str = {'Wcc', 'Bcc'};                            

    else
        pairTypeIdxs = {Wcc_pairIdxs(wcc_ok,:),  Wcc_pairIdxs(wcc_ok,:)};
        pairs_C      = {Wcc_pairs(wcc_ok,:), Wcc_pairs(wcc_ok,:)};
        pairTypes_str = {'Wcc', 'Wscc'};                            

    end
    
    allPossPhaseMeasures = {'cc', 'dphi', 'MID_cc', 'MID_fit_cc', 'MID_ovlp_cc', 'STA_cc'};
    %%
    doExtraFigures = 0;
    if doExtraFigures
        %%
       figure(301); clf;
       minFracR = pairData.loc_minFracOfMaxes(:,1);
       hist( minFracR, 25 );              
       xlabel('min FracR = Min(cells firing rate, relative to max)'); ylabel('Pair Count');
       
       figure(302); clf;
       minPtcOEcorr = pairData.loc_ptcOEcorr(:,1);
       hist( minPtcOEcorr, 25);
       xlabel('min O/E corr = Min(Same-cell CC between odd/even tuning curves)'); ylabel('Pair Count');        
        
    end
    
    if doPrintStats
        %%
        pairDiffs_matFile = getFileName('pairDiffs');
        pairDiffs_S.cmpType = 'phase';
        pairDiffs_S.gratingType = gratingType;
        pairDiffs_S.timeWindow = curTimeWindow;
        
        
        if showWorking
            fprintf('********************* TABLE 2: STATISTICS FOR PAIRWISE DIFFERENCES **********************\n');
            fprintf(' (%s) \n', cmp_str);
            fprintf('    Parameter       |    Mean     |      Std    |    Median   |     P25     |     P75      | N\n');        
        end
%         Wcc_vals = S{ms_idx}.val(Wcc_oo_pairIdxs);
%         Bcc_vals = S{ms_idx}.val(Bcc_oo_pairIdxs);       

        
        
        allMeasureNames = cell(1, length(measures)*length(locations)); all_ms_idx = 1;
    
        for loc_i = 1:length(locations)
            if showWorking
                fprintf('** %s ** ', locations{loc_i});
            end
        
            for ms_i = 1:length(measures)            

                measure_i = measures{ms_i};    
                                
                isRFmeasure = any(strcmp(measure_i, RF_Measures));
                if isRFmeasure && (loc_i > 1)
                    continue;  % RF measures don't have a location ==> are always put in loc = 1.
                end
                Si = S{loc_i, ms_i};
%                 mu_measure = ~isempty(strfind(measure, 'MU'));
%                 measure = strrep(measure, '_MU', '');


                for pt_i = 1:length(pairTypes_str)
%                     [vals_mean, vals_std, vals_median, vals_P25, vals_P75, Npr, Ncl] ;%= deal(zeros(1, nSpont));
                    %%
                    pt_idxs = pairTypeIdxs{pt_i};
                    
                    tf_pair_passedCrit = getPairIdxsThatPassCriteria(criteria, measure_i, loc_i, pt_idxs, pairData, allPossPhaseMeasures);
                    tf_vals_nonnan = ~any( isnan(Si.val(  pt_idxs, 1, :)),3);
                    
                    idx_use = find( tf_pair_passedCrit & tf_vals_nonnan );
                    
                    
                    if (ms_i == 1) && (pt_i == 1) && showWorking
                        fprintf(' (pairs that satisfied criteria : %s ) \n', frac_pct( nnz(tf_pair_passedCrit), length(tf_pair_passedCrit) ) );
                       3; 
                        
                    end
                    
                    if ~strcmp(pairTypes_str{pt_i}, 'Wscc')
                        vals_use = Si.val(  pt_idxs(idx_use), 1, :  );
                    else
                        vals_use = vertcat( Si.shuffVal {  pt_idxs(idx_use) } );
                    end
                    
%                     tf_nan_vals = isnan(vals_passedCrit);
%                     idx_nonnan_vals = ~any(any(tf_nan_vals, 2),3); 
%                         assert(isequal(any(tf_nan_vals, 2), all(tf_nan_vals, 2)) );                                    
%                         
%                     idx_pair_passedCrit = find(tf_pair_passedCrit);
%                     idx_pair_use = idx_pair_passedCrit(idx_nonnan_vals);                                                                
%                     vals_use = vals_passedCrit(idx_nonnan_vals, :, :);
                    
                    if size(vals_use, 3) > 1
                        vals_use = [vals_use(:,:,1); vals_use(:,:,2)];
                    end
                    
                    Npr = nnz(idx_use);
                    
%                     if pt_i == 1
%                         tf_Wcc_vals_nonnan = ~isnan(Si.val(  pt_idxs ));   
%                         tf_use = tf_pair_passedCrit & tf_Wcc_vals_nonnan;
%                         vals_use2 = Si.val(  pt_idxs(tf_use)  );
%                         assert(isequal(vals_use, vals_use2));                        
%                     end
                    
                    [vals_mean, vals_std, vals_median, ...
                        vals_P25, vals_P75, ~] = getMeanStdStats(vals_use(:));                      
                    
                    Ncl = length(unique(pairs_C{pt_i}(idx_use, :)) );      
                    if strcmp(measure_i, 'MID_cc')
                        3;
                    end                    
                    
                    if useShuffledPhaseControl
                        if strcmp(pairTypes{pt_i}, 'Wcc')
                            Npr_Wcc = Npr;
                        elseif strcmp(pairTypes{pt_i}, 'Wscc')
                            Npr = Npr_Wcc;
                        end
                    end                                        
                    
                    w = 8; 
                    fmt = ['%' num2str(w) '.3f'];

                    pairTmp_str = ['| ' fmt ' '];                    
                    nPairs_str = sprintf('%d Pr; ', Npr);
                    
                    nCells_str = sprintf('%d Cl', Ncl);  % RF_Measures = {'STA_cc', 'MID_cc', 'MID_ovlp_cc', 'MID_fit_cc', 'dPh_rel'};

                    nMU_str = '';
                    nCells_str = [nCells_str nMU_str];  %#ok<AGROW>
                    
                    str_template = ['%14s (' pairTypes_str{pt_i} ')' repmat( pairTmp_str, 1, 5) ' | ' nPairs_str '%s \n'];        

                    if showWorking
                        fprintf(str_template, ...
                            measure_i, vals_mean, vals_std, vals_median, vals_P25, vals_P75, nCells_str);
                    end
                    
                    measure_fld = [measure_i '_' pairTypes_str{pt_i}];
                    if ~isRFmeasure
                        measure_fld = [locations{loc_i} '_' measure_fld]; %#ok<AGROW>
                    end
                    
                    pairDiffs_S.(measure_fld) = struct('mean', vals_mean, 'std', vals_std, 'vals_median', vals_median, 'vals_P25', vals_P25, 'vals_P75', vals_P75, 'N', [nPairs_str  nCells_str], ...
                        'N_pairs', Npr, 'N_cells', Ncl); 
                    allMeasureNames{all_ms_idx} = measure_fld; all_ms_idx = all_ms_idx+1;
                                        
                end        
                if showWorking
                    fprintf('\n');        
                end
            end

        end
        
        
        pairDiffs_S.allMeasureNames_orig = measures;
        pairDiffs_S.allMeasureNames = allMeasureNames(1:all_ms_idx-1);        
        pairDiffs_S.columns = fieldnames(pairDiffs_S.(measure_fld));
        pairDiffs_S.minRoe = curMinR_oe;
        pairDiffs_S.comments = cellsPerSite_comments;
        save(pairDiffs_matFile, '-struct', 'pairDiffs_S');             

        if showWorking        
            fprintf('*******************************************************************************************\n\n\n');
        end
        3;
        
        
    end
        
    
    
    %%
    if ~useShuffledPhaseControl
        pairTypeIdxs = {Wcc_pairIdxs(wcc_ok), Wrcc_pairIdxs, Bcc_pairIdxs(bcc_ok,:)};
        Wcc_idx = find(strcmp(pairTypes, 'Wcc'), 1);
        Bcc_idx = find(strcmp(pairTypes, 'Bcc'), 1);
        rand_idx = find(strcmp(pairTypes, 'Wrcc'), 1);
        
        % pairTypeIdxs = {Wcc_idxs(wcc_ok,:), Bcc_idxs(bcc_ok,:)};
        
    else
        pairTypeIdxs = {Wcc_pairIdxs(wcc_ok,:), Wcc_pairIdxs(wcc_ok,:)};
        
        Wcc_idx = find(strcmp(pairTypes, 'Wcc'), 1);
        rand_idx = find(strcmp(pairTypes, 'Wscc'), 1);
    end
    
    %%
    if doPrintSigTests
        %%
        pairStats_matFile = getFileName('pairStats');
        pairStats_S.cmpType = 'phase';
        pairStats_S.gratingType = gratingType;
        pairStats_S.timeWindow = curTimeWindow;

        doChiSqrForDphi = 0;
        
        doRandKStest = (useShuffledPhaseControl && (Npermutes <= 10000)) || (~useShuffledPhaseControl);
        if showWorking
    %       fprintf('\n\n ************************** SIGNIFICANCE TESTS (npermute = %d) **************************** \n', Npermutes);
    %       fprintf('  Parameter    | Median Prob(1s)| Median Prob(2s)| Mean Ratio (1s)| Mean Prob (2s) |   KS Stat   |   KS prob   |\n') 

            fprintf('\n\n ************************** SIGNIFICANCE TESTS (Npermute = %d) ********************************* \n', Npermutes);
            fprintf(' ** Using %s as a control\n', iff(useShuffledPhaseControl, 'Phase-shuffling', 'Between-site distribution'))
            fprintf(' *** %s *** \n', cmp_str);
            %                             1              2             3              4           5             6            7               8              9              10                              
            %                          
            %                             
            fprintf('  Parameter    | Median      : Median Prob |  Mean       :  Mean Prob  |   KS_stat   :   KS prob   | Median(W)    :Mean Prob(T) |KS prob(std)  |KS stat(std)  | other ...\n') 
        end
        
        if ~isempty(criteria_str)
            fprintf(' Criteria %s \n', criteria_str);        
        end
        if ~useShuffledPhaseControl    
            pairTypeIdxs = {Wcc_pairIdxs(wcc_ok), Wrcc_pairIdxs, Bcc_pairIdxs(bcc_ok,:)};                
            Wcc_idx = find(strcmp(pairTypes, 'Wcc'), 1);
            Bcc_idx = find(strcmp(pairTypes, 'Bcc'), 1);
            rand_idx = find(strcmp(pairTypes, 'Wrcc'), 1);
                                    
%             pairTypeIdxs = {Wcc_idxs(wcc_ok,:), Bcc_idxs(bcc_ok,:)};
            
        else
            pairTypeIdxs = {Wcc_pairIdxs(wcc_ok,:), Wcc_pairIdxs(wcc_ok,:)};            

            Wcc_idx = find(strcmp(pairTypes, 'Wcc'), 1);
            rand_idx = find(strcmp(pairTypes, 'Wscc'), 1);            
        end

        
        allMeasureNames = cell(1, length(measures)*length(locations)); all_ms_idx = 1;
        
        for loc_i = 1:length(locations)
            if showWorking
                fprintf('** %s **\n', locations{loc_i});
            end
        
%             figure(loc_i); clf;
            for ms_i = 1:length(measures)
%                 tic;
                measure_i = measures{ms_i};
                
                isRFmeasure = any(strcmp(measure_i, RF_Measures));
                if isRFmeasure && (loc_i > 1)
                    continue;  % RF measures don't have a location ==> are always put in loc = 1.
                end

                [vals_mean, vals_median, vals_std, vals_KS, N] = deal(cell(length(pairTypes), 1));
                
%                 [medianProb, meanProb, ksStat, ksProb, ksProb2, vals_Bcc] = deal( 0 );
                                
                Si = S{ loc_i, ms_i };
                                
                tf_Wcc_pair_passedCrit = getPairIdxsThatPassCriteria(criteria, measure_i, loc_i, Wcc_pairIdxs, pairData, allPossPhaseMeasures);                
                tf_Wcc_vals_nonnan = ~any( isnan(Si.val(  Wcc_pairIdxs, 1, :)),3);
                tf_Wcc_idx_use = tf_Wcc_pair_passedCrit & tf_Wcc_vals_nonnan;
                
                Wcc_vals = Si.val(  Wcc_pairIdxs(tf_Wcc_idx_use), 1, : );
                Wcc_vals = Wcc_vals(:);
                  
                                
                
                
                if doRandKStest
                    if useShuffledPhaseControl
                        ctrl_dist = ks_getCtrlDist(Si.ctrlVal, tf_Wcc_idx_use);
                    else                        
                        error('add code to make sure bcc pairs pass criteria, too')                        
                        ctrl_dist = ks_getCtrlDist(single(Si.val(pairTypeIdxs{Bcc_idx})));                                                                                                
                    end
                else
                    ctrl_dist = [];
                end
                

                for pt_i = 1:length(pairTypes)                                    
                    pairIdxs = pairTypeIdxs{pt_i};
                    switch pairTypes{pt_i}
                        case {'Wcc'},
%                             vals = Si.val( pairTypeIdxs{Wcc_idx}(tf_Wcc_idx_use) );
                            vals = Si.val( tf_Wcc_idx_use, 1, : );
                            [vals_mean{pt_i}, vals_median{pt_i}, vals_std{pt_i}, vals_KS{pt_i}, N{pt_i}] = ...
                                getMeanMedianStdKS(vals(:), ctrl_dist);                          
                        case 'Wscc',
                            shuffVals = vertcat( Si.shuffVal {tf_Wcc_idx_use} );
                            [vals_mean{pt_i}, vals_median{pt_i}, vals_std{pt_i}, vals_KS{pt_i}, N{pt_i}] = ...
                                getMeanMedianStdKS(shuffVals, ctrl_dist);                                                                              
                        case 'Wrcc',
                            assert(iscell(pairIdxs));
                            [vals_mean{pt_i}, vals_median{pt_i}, vals_std{pt_i}, vals_KS{pt_i}, N{pt_i}] = ...
                                getMeanMedianStdKS( Si.val, ctrl_dist, pairIdxs );                                                    
                    end
                    
                end
                
%                 one_sided_cmp = switchh(measure_i, {{'cc', 'rho'}, {'dphi', 'dF1'}}, {@ge, @le});
                measure_i_tmp = measure_i;
                if ~isempty(strfind(measure_i_tmp, 'cc'))
                    measure_i_tmp = 'cc';
                end
                one_sided_tail = switchh(measure_i_tmp, {{'cc', 'rho'}, {'dphi', 'dF1', 'dPh_rel'}}, {'right', 'left'});
                null_mean_value = switchh(measure_i_tmp, {{'cc', 'rho'}, {'dphi', 'dF1', 'dPh_rel'}}, [0, 90]);
                
                [medianProb_WTest, meanProb_TTest] = deal(nan);
                % Median Data
                Wcc_median = vals_median{Wcc_idx};
                rand_medians = vals_median{rand_idx};                
                medianProb = getRandomizedProb(Wcc_median, rand_medians, 'both');  
                if length(rand_medians) > 10
                    medianProb_WTest = signrank(Wcc_vals, median(rand_medians));
                end
               
                                                                
                % Mean Data
                Wcc_mean = vals_mean{Wcc_idx};
                rand_means = vals_mean{rand_idx};
                meanProb = getRandomizedProb(Wcc_mean, rand_means, 'both');
                if length(rand_means) > 10
                    [~, meanProb_TTest] = ttest(Wcc_vals, mean(rand_means));
                end

                % STD data
                Wcc_std = vals_std{Wcc_idx};
                rand_stds = vals_std{rand_idx};
                stdProb = getRandomizedProb(Wcc_std, rand_stds, 'both');                
                
                % KS data
                Wcc_KS = vals_KS{Wcc_idx};
                rand_KS = vals_KS{rand_idx};                
                ksProb = getRandomizedProb(Wcc_KS, rand_KS, 'right'); % KS-statistic: alternative hypothesis is that KS-stat is larger.
                
                % chi-square for delta phi:
                if strcmp(measure_i_tmp, 'dphi') && doChiSqrForDphi
                    %%
                    Wcc_vals = Wcc_vals(~isnan(Wcc_vals));
                    binEdges = linspace(0, 180, 13);
                    binEdges(1) = -1; binEdges(end) = 181;
                    if useShuffledPhaseControl
                        ctrl_dist = ks_getCtrlDist(Si.shuffVal(:));
                    else
                        ctrl_dist = ks_getCtrlDist(Si.val(pairTypeIdxs{Bcc_idx}), binEdges);                        
                    end
                    expectedCounts = diff( [0; ctrl_dist.CDF] );
                    expectedCounts = expectedCounts/sum(expectedCounts)*length(Wcc_vals);
                    
                    if isinf(binEdges(1))
                        edges = ctrl_dist.binEdges(2:end);
                    end
                    [~, p_chisqr,chisqr_stat] = chi2gof(Wcc_vals,'edges',binEdges, 'expected',expectedCounts);                
                    
                    3;
                end
                %%
                vals_Wcc  = nonnans( Si.val(  Wcc_pairIdxs  ));
                if useShuffledPhaseControl
                    idx_use = cellfun(@(x) length(x) > 1, Si.shuffVal);
                    vals_Bcc  = vertcat( Si.shuffVal{idx_use}); 
                else
                    vals_Bcc  = nonnans( Si.val(  pairTypeIdxs{Bcc_idx}  )); 
                end
                [ksProb_standard, ksstat_standard] = deal(nan);
                if ~isempty(vals_Bcc) && (numel(vals_Bcc) < 1e6)
                    [h, ksProb_standard, ksstat_standard] = kstest2(vals_Wcc, vals_Bcc(:));
                end
                %%
%                 Wrcc_stds_std = std(Wrcc_stds);
                
%                 medianProb_1s = getRandomizedProb(Wcc_median, rand_medians, one_sided_tail);
%                 meanProb_1s = probFunc(Wcc_mean, rand_means, one_sided_tail);
                                
                showWorking_ext = 0;
                if showWorking_ext

                    if any(strcmp(measure_i, {'dphi', 'dF1'}))
                        ax_cent = 90;
                    else
                        ax_cent = 0;
                    end

    %                 subplot(length(measures), 2, sub2ind([2 length(measures)], 1, ms_i));
    %                 hist(rand_medians, 30); title(sprintf('%s medians', measure_i));
    %                 xlabel(sprintf('P(1) = %.2f, P(2) = %.2f', medianProb_1s, medianProb_2s))
    %                 xlim( max(abs(xlim-ax_cent))*[-1, 1] + ax_cent)
    %                 drawVerticalLine(mean(rand_medians), 'color', 'r');
    %                 drawVerticalLine(Wcc_median, 'color', 'g');

                    subplot(length(measures), 2, sub2ind([2 length(measures)], 1, ms_i));
                    hist(rand_stds, 30); title(sprintf('%s stds', measure_i));
    %                 xlim( max(abs(xlim-ax_cent))*[-1, 1] + ax_cent)
                    drawVerticalLine(mean(rand_stds), 'color', 'r');
                    drawVerticalLine(Wcc_std, 'color', 'g');

                    subplot(length(measures), 2, sub2ind([2 length(measures)], 2, ms_i) );
                    hist(rand_means, 30); title(sprintf('%s means', measure_i));
                    xlabel(sprintf('P(1) = %.2f, P(2) = %.2f', meanProb_1s, meanProb_2s))
                    xlim( max(abs(xlim-ax_cent))*[-1, 1] + ax_cent)
                    drawVerticalLine(mean(rand_means), 'color', 'r');
                    drawVerticalLine(Wcc_mean, 'color', 'g');
                    3;
                end
                
                                
%                 mwwProb = ranksum(vals_Wcc, vals_Bcc);                


%                 fmt_val = ['%11.3f'];
%                 fmt_pval = ['%11.4f'];
%                 fmt_pval_test = ['%11.3g'];
%                 fmt_pval_w = ['%21.3g'];

%%
                value_str   = @(val) sprintf('| %11.3f ', val); % fmt_val ' '];
                pval_str    = @(val) sprintf('| %s %7.4f ', pvalStarStr(val), val);  
                reg_test_str = @(val) sprintf('| %s %7.4g ', pvalStarStr(val), val);  
%                 pval_str_w  = ['| ' fmt_pval_w ' '];                                
%                 value_pval_pvalT_str = [value_str, pval_str];
                
                dphi_extra_str = '';
                if strcmp(measure_i_tmp, 'dphi') && doChiSqrForDphi
                    dphi_extra_str = sprintf('p_chi2 = %.3g (chi2stat = %.4f)', p_chisqr,chisqr_stat.chi2stat);
                end
                
%                 str_template = ['%14s ' repmat( pairTmp_str, 1, 4)  repmat( pairTmp_str_w, 1, 2)  ' | \n'];        
%                 fprintf(str_template, ...
%                     measure_i, medianProb_1s, medianProb_2s, meanProb_1s, meanProb_2s, ksProb, mwwProb);
%                 misc_str = '';
                misc_str = sprintf('stdProb = %.3f', stdProb);
                str_template = ['%14s ' repmat( '%s', 1, 10),  ' %s %s| \n'];        
                %%
                if showWorking
                    %%
                    fprintf(str_template, ...
                        measure_i, value_str( Wcc_median ), pval_str( medianProb), ...
                                   value_str( Wcc_mean   ), pval_str(   meanProb), ...
                                   value_str( Wcc_KS     ), pval_str(    ksProb), ...
                                   reg_test_str( medianProb_WTest), reg_test_str( meanProb_TTest), ...  
                                   reg_test_str( ksProb_standard ), value_str( ksstat_standard ), ...
                                   dphi_extra_str, misc_str);
                end
                
                if saveControlDistribs
%                     S_save.(measure) = struct('vals_Bcc',  vals_Bcc{1} , ...
%                                               'vals_rand_median', vals_median{rand_idx, 1}, ...
%                                               'vals_rand_mean', vals_mean{rand_idx, 1});
                end
                3;
                %%
                measure_fld = measure_i;
                if ~isRFmeasure
                    measure_fld = [locations{loc_i} '_' measure_fld]; %#ok<AGROW>
                end                
%                 measure_fld = [locations{loc_i} '_' measure_i]; 

                pairStats_S.(measure_fld) = struct('median', Wcc_median, 'medianProb', medianProb, 'Wscc_medians', rand_medians, ...
                                                   'mean', Wcc_mean, 'meanProb', meanProb, 'Wscc_means', rand_means, ...
                       'ksStat', Wcc_KS, 'ksProb', ksProb, 'Wscc_KS', rand_KS, ...
                       'medianProb_W_test', medianProb_WTest, 'meanProb_TTest', meanProb_TTest, 'ksProb_standard', ksProb_standard);

                allMeasureNames{all_ms_idx} = measure_fld; all_ms_idx = all_ms_idx+1;
%                 toc;
            end        
            3   ;
        end
        if showWorking
            fprintf('\n\n');
        end
        
        if saveControlDistribs
%             save(statsDatafile, '-struct', 'S_save');
        end
        
        pairStats_S.allMeasureNames_orig = measures;
        pairStats_S.allMeasureNames = allMeasureNames(~cellfun(@isempty, allMeasureNames));
        allColumns = fieldnames(pairStats_S.(measure_fld));
        allColumns( cellfun(@(s) strncmp(s, 'Wscc', 4), allColumns) ) = [];
        pairStats_S.columns = allColumns;
        pairStats_S.minRoe = curMinR_oe;
        save(pairStats_matFile, '-struct', 'pairStats_S');             
        
    end
    

    if plotInColor % || 1
        %%
        wcc_col = 'b';  ctrl_col = 'r';
        wcc_col_line = 'b';
%         bar_norm = 'b'; 
        bar_norm = [.8, .8, .8];
        wcc_sz = 3; ctrl_sz = 3;
        wcc_line_w = 1; ctrl_line_w = 1; ctrl_err_w = 1;
        
    else
        %%
        wcc_col = 'k'; ctrl_col = [0 0 0];
        wcc_col_line = [.6, .6, .6];
        bar_norm = [.8, .8, .8];
        wcc_sz = 6; ctrl_sz = 3;
        wcc_line_w = 2; ctrl_line_w = 1; ctrl_err_w = 2;
    end

    
    
    
    
    if doCCvsTetrodeDepthPlots
        %%
        makeFiguresForPaper = 0 && 1; %~showWorking;
        
%         penDistStats_matFile = getFileName('penDistStats');
        
        
        makeDistFigs = showWorking;
        penDistStats_S.cmpType = 'phase';
        penDistStats_S.gratingType = gratingType;
        penDistStats_S.subtractSpont = curSubtractSpont;
        penDistStats_S.bccType = curBccType;
        
        nResamplesTotal = 200;
%                     nResamplesTotal = 50;
        randInfo.nResamplesTotal = nResamplesTotal;
        randInfo.getBrccPairs = 1;
        randInfo.nUnits = nUnits;
        randInfo.idxMtx = idxMtx;

        
        saveStatsToFile = (nResamplesTotal == 2000);
        
        %             [Wrcc_idxs, nsimp_wrcc,   Brcc_idxs, nsimp_brcc] = getRandomizedPairs(randInfo, loop_i==1, nResamplesPerLoop);
        
%         dF1oDC_cc = abs(diff(pairData_ori.F1oDCs_pref, [], 2));
%         dF1oDC_wcc = dF1oDC_cc(Wcc_oo_pairIdxs);
%         dF1oDC_bcc = dF1oDC_cc(Bcc_oo_pairIdxs);
                

        if makeFiguresForPaper

            fig_idx = 200;
            subSpcM = [0.01 0.015, 0.01];
            subSpcN = [0.02 0.00, 0.02];
        else
            fig_idx = curGratingType*100;
        end
        %             stimTypes = {'spf', 'ori'};
        measureTypes      = {'cc'}; 
        oriMeasureFieldNames = {'cc'};
                
        set(0,'DefaultFigureWindowStyle','docked')

%         locName = 'maxMinFracR';
%         locName = {'p75MinFracR', 'p50MinFracR', 'p25MinFracR', 'minMinFracR'};
%         locName = 'p50MinFracR';
        locName = 'p75MinFracR';
        
        loc_idx = find(strcmp(locName, locations));
        
        cc_idx = find(strcmp(measures, 'cc'), 1);
        cc_all = S{loc_idx, cc_idx}.val;
        cc_bcc = cc_all(Bcc_pairIdxs);
        cc_wcc = cc_all(Wcc_pairIdxs);
        cc_binEdges = S{loc_idx, cc_idx}.binEdges;

        dphi_idx = find(strcmp(measures, 'dphi'), 1);
        dphi_all = S{loc_idx, dphi_idx}.val;
        dphi_bcc = dphi_all(Bcc_pairIdxs);
        dphi_wcc = dphi_all(Wcc_pairIdxs);
        dphi_binEdges = S{loc_idx, dphi_idx}.binEdges;

        measureTypes       = {'CC', 'DPhi'};
        measureFieldNames  = {'cc', 'dphi'};
        fmt = {'%.3f', '%.1f'};
        
        measures_cc  = {cc_all, dphi_all};
        measures_bcc = {cc_bcc, dphi_bcc};
        measures_wcc = {cc_wcc, dphi_wcc};
        ms_binEdges = {cc_binEdges, dphi_binEdges};

        if strcmp(gratingType, 'flashed')
            STAcc_idx = find(strcmp(measures, 'STA_cc'), 1);
            STAcc_all = S{1, STAcc_idx}.val;
            STAcc_bcc = STAcc_all(Bcc_pairIdxs);
            STAcc_wcc = STAcc_all(Wcc_pairIdxs);
            STAcc_binEdges = S{1, STAcc_idx}.binEdges;

            measureTypes = [measureTypes, 'STA_cc'];
            measureFieldNames = [measureFieldNames, 'STA_cc'];
            fmt = [fmt, '%.3f'];
            
            measures_cc  = [measures_cc, STAcc_all];
            measures_bcc = [measures_bcc, STAcc_bcc];
            measures_wcc = [measures_wcc , STAcc_wcc];
            ms_binEdges = [ms_binEdges, STAcc_binEdges];

        end
        
        
        
        
%         if strcmp(stimTypes{ti}, 'ori')

            
            
            doSSpairs = 0;
            [Wrcc_idxs_ori, ~,   Brcc_idxs_ori, ~] = getRandomizedPairs(randInfo, 1, nResamplesTotal);
            [Wrcc_pairIdxs, Brcc_idxs] = deal(Wrcc_idxs_ori, Brcc_idxs_ori);

%         end

        %%
        if doSSpairs
            bcc_SC_type = pairData.SCtype_pref(Bcc_pairIdxs);
            wcc_SC_type = pairData.SCtype_pref(Wcc_pairIdxs);
            idx_use_bcc = find(bcc_SC_type == 2);
            idx_use_wcc = find(wcc_SC_type == 2);

            Bcc_pairIdxs_M = Bcc_pairIdxs_M(idx_use_bcc);
            Wcc_pairIdxs_M = Wcc_pairIdxs_M(idx_use_wcc);
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
        %%
        tf_pairs_cell = pairData.cellIds > 0;
        tf_cc = tf_pairs_cell(:,1) & tf_pairs_cell(:,2);
        tf_cm = xor(tf_pairs_cell(:,1), tf_pairs_cell(:,2));
        tf_mm = ~tf_pairs_cell(:,1) & ~tf_pairs_cell(:,2);
        assert(nnz(tf_cc) + nnz(tf_cm) + nnz(tf_mm) == length(tf_pairs_cell));
        
        tf_sameGrp = pairData.Gids(:,1) == pairData.Gids(:,2);
        
        tf_bcc = tf_cc(B_pairIdxs);
        tf_bcm = tf_cm(B_pairIdxs);
        tf_bmm = tf_mm(B_pairIdxs);

        %%
        

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
            
            title_fsize = 10; % 9
            xlabel_fsize = 10; %9;
            ylabel_fsize = 10; %9;
            legend_fsize = 10;

        end


        nSubplots = showDistributions + showWP + double(showWA);

        for mi = 1:length(measureTypes);

            if makeFiguresForPaper && makeDistFigs
                sub_row_idx = find(strcmp(measureTypes{mi}, column_measures), 1);
                if isempty(sub_row_idx);
                    continue;
                end
            end


            linew = 2;
            x_wcc = measures_wcc{mi};
            x_bcc = measures_bcc{mi};
            if doSSpairs
                x_wcc = x_wcc(idx_use_wcc);
                x_bcc = x_bcc(idx_use_bcc);
            end


            %                     dX_wrcc = cellfun(@(idx) measures_cc{mi}(idx), Wrcc_idxs, 'un', 0);
            %                     dX_brcc = cellfun(@(idx) measures_cc{mi}(idx), Brcc_idxs, 'un', 0);

            X_name = measureTypes{mi};
            if ~makeFiguresForPaper && makeDistFigs
                fig_idx = fig_idx+1;
                figure(fig_idx); clf;

            end

%             tf_b_pairType = tf_bcc;  b_pairType_str = 'Cell-Cell';
%             tf_b_pairType = tf_bcm; b_pairType_str = 'Cell-MU';
            tf_b_pairType = tf_bmm; b_pairType_str = 'MU-MU';
            
            tf_WA          = tf_b_pairType &                 sameAnimal_bcc==1; 
            tf_WP          = tf_b_pairType & samePen_bcc==1;
            tf_BP_WA       = tf_b_pairType & samePen_bcc==0 & sameAnimal_bcc==1;
            tf_BP_BH_WA    = tf_b_pairType & samePen_bcc==0 & sameHemisphere_bcc == 0 & sameAnimal_bcc==1;
            tf_BP          = tf_b_pairType & samePen_bcc==0;
            tf_BA          = tf_b_pairType &                 sameAnimal_bcc==0;
            tf_BP_WA_close = tf_b_pairType & tf_BP_WA & penDist_bcc_nz < .5;
            tf_BP_WA_far   = tf_b_pairType & tf_BP_WA & penDist_bcc_nz > 2;

%             stat_name = 'median'; stat_func = @nanmedian;
            stat_name = 'mean'; stat_func = @nanmean;
            %%% look at all the distributions:
            Dists_C     = {nonnans(x_bcc),             nonnans(x_bcc(tf_BP_WA)),  nonnans(x_bcc(tf_WP)),   nonnans(x_wcc)};
            dists_names = {'BS             ', 'BS-BP-WA ',        'BS-WP      ',   'WS            '};
            i_bs = 1; i_bs_bp_wa = 2; i_bs_wp = 3; i_ws = 4;
            means = cellfun(stat_func, Dists_C);
            bs_mean = stat_func(x_bcc);
            bs_wp_mean = stat_func(x_bcc(tf_WP));
            bs_bp_wa_mean = stat_func(x_bcc(tf_BP_WA));
            ws_mean = stat_func(x_wcc);

            N = cellfun(@length, Dists_C);
            nBS = N(1); nWABP= N(2);  nWP = N(3); nWA = nWABP+nWP;
%                 fprintf('N (BS) = %d. N(WA) = %d. N(WA-BP) = %d. R(BS/WA) = %.1f. R(BS/WA-BP) = %.1f)\n', nBS, nWA, nWABP, nBS/nWA, nBS/nWABP);
            mean_strs = arrayfun(@(m) sprintf([fmt{mi}], m), means, 'un', 0);
            leg_strs = arrayfun(@(nm, n, s) sprintf(['%s (N = %d, mean = %s) .'], nm{1}, n, s{1}), dists_names, N, mean_strs, 'un', 0);

            
            if showDistributions && makeDistFigs
                subplotGap(1,nSubplots,1);
                %                     Dists_C = {dX_bcc, dX_bcc(tf_WA), dX_bcc(tf_WP), dX_bcc(tf_BP_WA), dX_wcc};
                %                     dists_names = {'BS              ',  'BS-WA       ', 'BS-WP       ', 'BS-BP-WA ', 'WS              '};
                %                     Dists_C = {dX_bcc, dX_bcc(tf_BP), dX_bcc(tf_WA), dX_bcc(tf_WP), dX_bcc(tf_BP_WA), dX_wcc};
                %                     dists_names = {'BS              ',  'BS-BP       ', 'BS-WA       ', 'BS-WP       ', 'BS-BP-WA ', 'WS              '};


                h_hists = hist2(Dists_C, ms_binEdges{mi}, 'line', 'norm');
                mks = 'sovs.';
                for i = 1:length(h_hists)
                    set(h_hists(i), 'marker', mks(i), 'linewidth', linew)
                end
                legend(leg_strs, 'location', 'NE', 'fontsize', 8);

                c_order = get(gca, 'colorOrder');
                for i = 1:length(N)
                    drawVerticalLine(means(i), 'color', c_order(i,:), 'linestyle', ':', 'linewidth', 2);
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
            %                     depth_means = zeros(1, nSteps);
            %                     for i = 1:nSteps
            %                         depth_means(i) = nanmean(  dX_bcc(samePen_bcc==1 & depthDist_bcc > dLims(i) ) );
            %
            %                     end




            %% look at within-penetration distributions


            %                     figure(402);
            %                     plot(dLims, depth_means, 'o-');
            %                     opt.nPairsEach = nPairsEach;
            %                     nPairsEach = 250;
            opt.nPairsEach = nPairsEach;
            opt.overlap_f = 5;
            opt.bs_mean = bs_mean;
            opt.sigTestTail = 'both';
            %                     opt.sigTestTail = 'left';


            %                     dX_wrcc = cellfun(@(idx) measures_cc{mi}(idx), Wrcc_idxs, 'un', 0);

            %                     opt.allDistances_sorted = sort(nonnans(depthDist_bcc));
            
            depthDist_bcc_WP = depthDist_bcc(tf_WP);
            dX_bcc_WP = x_bcc(tf_WP);
            dX_brcc_WP = cellfun(@(idx) measures_cc{mi}(idx(tf_WP)), Brcc_idxs, 'un', 0);
            doCheck = 0;
            if doCheck
                %%
                dX_brcc = cellfun(@(idx) measures_cc{mi}(idx), Brcc_idxs, 'un', 0);
                dX_brcc_WP2  = cellfun(@(dX) dX(tf_WP), dX_brcc, 'un', 0);
                assert(isequaln(dX_brcc_WP, dX_brcc_WP2));
            end



            
            penDist_bcc_BP_WA = penDist_bcc_nz(tf_BP_WA);
            dX_bcc_BP_WA = x_bcc(tf_BP_WA);
            dX_brcc_BP_WA = cellfun(@(idx) measures_cc{mi}(idx(tf_BP_WA)), Brcc_idxs, 'un', 0);

            if doCheck
                dX_brcc = cellfun(@(idx) measures_cc{mi}(idx), Brcc_idxs, 'un', 0);
                dX_brcc_BP_WA2 = cellfun(@(dX) dX(tf_BP_WA), dX_brcc, 'un', 0);
                assert(isequaln(dX_brcc_BP_WA, dX_brcc_BP_WA2));
            end
            %%
%             [y_depth_means, y_depth_rand_lo, y_depth_rand_hi] = deal([]);
            %                     if showWP
            opt.overlap_f = 3;
            opt.y_stat_func = stat_func;
            [xs_depth, x_depth_lo, x_depth_hi,   y_depth_means, y_depth_lo, y_depth_hi,  y_depth_pvalues, depth_rand_med, y_depth_rand_lo, y_depth_rand_hi] = ...
                getSlidingWindowMediansVsDist(depthDist_bcc_WP, dX_bcc_WP, dX_brcc_WP, opt);
            %                     end
            %%
            [horiz_means, horiz_rand_lo, horiz_rand_hi] = deal([]);
            if showWA
                opt.overlap_f = 3.5;
                [xs_horiz, x_horiz_lo, x_horiz_hi, horiz_means, horiz_means_lo, horiz_means_hi, horiz_pvalues, horiz_rand_med, horiz_rand_lo, horiz_rand_hi] = ...
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

            idx_notSig = find(y_depth_pvalues > pval_th, 1);
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
                    subplotGap(1,nSubplots, 1+showDistributions);
                end
                %                     h = ploterr(xs_depth, y_depth_means, {x_depth_lo, x_depth_hi}, {depth_median_lo, depth_median_hi}, 'hhxy', .7); hold on;

                if ~makeFiguresForPaper
                    %                             h_rand = plot(xs_depth, depth_rand_med, 'linewidth', 1, 'color', 'b', 'marker', 's', 'markersize', 6); hold on;
                end
                line_w_hiLo = iff(makeFiguresForPaper, 1, 2);
                h_rand_hiLo = plot(xs_depth, [y_depth_rand_lo', y_depth_rand_hi'], 'linewidth', line_w_hiLo, 'color', randHiLo_color, 'linestyle', '-'); hold on;

                %                     h_rand = ploterr(xs_depth, depth_rand_med, [], {depth_rand_lo, depth_rand_hi}, 'hhxy', .7); hold on;
                %                     set(h_rand, 'linewidth', 2, 'color', 'b') ;
                %                     set(h_rand(1), 'marker', 's', 'markersize', 6)
                %                     set(h_rand(2), 'color', [0 0 .7]) ;

                h_sliding = ploterr(xs_depth, y_depth_means, {x_depth_lo, x_depth_hi}, [], 'hhxy', 1); hold on;
                set(h_sliding, 'linewidth', 1.5, 'color', sliding_color) ;
                set(h_sliding(1), 'marker', 'o', 'markersize', 6)
                set(h_sliding(2), 'color', sliding_color_err) ;

                %                     h_hiLo = plot(xs_depth, [depth_median_lo', depth_median_hi'], 'color', [.7 0 0], 'linestyle', ':', 'linewidth', 2);



                %                     , 'o-', 'linewidth', linew, 'color', 'r'); hold on;

                ylims1 = lims([y_depth_means, means, y_depth_rand_lo, y_depth_rand_hi], [.5, .3]);
                ylims2 = lims([horiz_means, means, horiz_rand_lo, horiz_rand_hi], .05);
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

                    depth_medians_sm = gaussSmooth(y_depth_means, 3);
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
                    h_B = drawHorizontalLine(bs_mean, 'linestyle', '--', 'color', 'b', 'linewidth', linew);
                    h_BW = drawHorizontalLine(bs_wp_mean, 'linestyle', '--', 'color', 'r', 'linewidth', linew);
                    h_W = drawHorizontalLine(ws_mean, 'linestyle', '--', 'color', [0, .75, .75], 'linewidth', linew);
                    h_lines_legend = [h_B, h_BW, h_W, h_sliding(1), h_rand_hiLo(1)];
                    legend_str_lines = {sprintf('BS %s (%s)', stat_name, mean_strs{i_bs}), sprintf('BS-WP %s (%s)', stat_name, mean_strs{i_bs_wp}), ...
                        sprintf('WS %s (%s)', stat_name, mean_strs{i_ws}), 'BS-WP(window)', 'BS-WP(window)-rand'};
                else
                    showBS = 1;
                    h_W = drawHorizontalLine(ws_median, 'linestyle', '--', 'color', 'k', 'linewidth', linew);
                    h_lines_legend = [h_W, h_sliding(1), h_rand_hiLo(1)];
                    legend_str_lines = {sprintf('WS %s', stat_name), 'BS-WP (observed)', 'BS-WP (control)'};
                    if showBS
                        h_B = drawHorizontalLine(bs_median, 'linestyle', ':', 'color', .3*[1, 1, 1], 'linewidth', linew);
                        h_lines_legend = [h_B, h_lines_legend];
                        legend_str_lines = [sprintf('BS %s', stat_name), legend_str_lines];
                    end

                end

                show05 = 0;
                show01 = 1;
                alsoShow001 = 0;
                if alsoShow001
                    idx_lt_001 = find( y_depth_pvalues < .001 );
                    idx_lt_01 = find( y_depth_pvalues >= .001 & y_depth_pvalues < .01 );
                else
                    idx_lt_01 = find( y_depth_pvalues < pval_th );
                end
                idx_lt_05 = find( y_depth_pvalues >=  .01 & y_depth_pvalues < .05 );

                %                         idx_lt_001 = find( y_depth_pvalues < .001 );
                %                         idx_lt_01 = find( y_depth_pvalues >= .001 & y_depth_pvalues < .01 );
                %                         if makeFiguresForPaper
                yvals_stars = @(idx) yval_stars * ones(size(idx));
                %                         else
                %                             yvals_stars = @(idx) y_depth_means(idx);
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
                xlabel(sprintf('%s distance between sites (\\mum)%s', stat_name, nPairs_str), 'fontsize', xlabel_fsize);
                X_name_lower = X_name;
                if ~any(strncmp(X_name, {'DSI', 'SF', 'F1'}, 2))
                    X_name_lower(1) = lower(X_name(1));
                end
                ylabel(sprintf('%s diff in %s', stat_name, X_name_lower), 'fontsize', ylabel_fsize);
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

                %                     plot(xs_horiz, horiz_means, 'o-', 'linewidth', linew, 'color', [0 .75 0]); hold on;
                h_sliding = ploterr(xs_horiz, horiz_means, {x_horiz_lo, x_horiz_hi}, [], 'hhxy', 1.2); hold on;
                %                     h = ploterr(xs_horiz, horiz_means, {x_horiz_lo, x_horiz_hi}, {horiz_medians_lo, horiz_medians_hi}, 'hhxy', 1.2); hold on;
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
                    horiz_medians_sm = gaussSmooth(horiz_means, 10);
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

                if ~isempty(idx_lt_001), h_sig_horiz(1) = plot(xs_horiz(idx_lt_001), horiz_means(idx_lt_001), 'k*', 'markersize', 8, 'linewidth', 2); end
                if ~isempty(idx_lt_01), h_sig_horiz(2) = plot(xs_horiz(idx_lt_01), horiz_means(idx_lt_01), 'ks', 'markersize', 8, 'linewidth', 2); end
                if ~isempty(idx_lt_05), h_sig_horiz(3) = plot(xs_horiz(idx_lt_05), horiz_means(idx_lt_05), 'ko', 'markersize', 8, 'linewidth', 2); end

                legend([h_B, h_BP_WA, h_W, h_sliding(1), h_rand], {sprintf('BS %s (%s)', stat_name, mean_strs{i_bs} ), ...
                    sprintf('BS-BP-WA %s (%s)', stat_name, mean_strs{i_bs_bp_wa} ), ...
                    sprintf('WS %s (%s)', stat_name, mean_strs{i_ws} ), ...
                    'BS-BP-WA(window)', 'BS-BP-WA(window)-rand'}, 'location', 'SE', 'fontsize', 8);

                xlabel(sprintf('Median AP-ML distance (mm) of %d pairs in sliding window', nPairsEach), 'fontsize', 9);
                ylabel(sprintf('%s diff in %s', stat_name, X_name));
                title(sprintf('%s (%s gratings) : BS-BP-WA distribution. %s %s', X_name, gratingType, x_half_str, sc_str), 'fontsize', 9);

            end
            3;

        end


        
        
        penDistStats_S.columns = {'dist'};
        penDistStats_S.allMeasureNames = allMeasureFieldNames;
        penDistStats_S.nResamples = nResamplesTotal;
%         if saveStatsToFile
%             save(penDistStats_matFile, '-struct', 'penDistStats_S');
%         end
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
    
    
    if doPlots
        %%
       doFigs = 2;
%        controlVals = 'Wrcc';
%        controlVals = 'Bcc';    
       controlVals = 'Wscc';

%         dphi_normalize = 'wrcc'; % 'wrcc', 'bcc', 'null'
        dphi_normalize = 'null'; % 'wrcc', 'bcc', 'null', 'none'
       
%         useShuffledPhaseCtrl_inPlots = strcmp(controlVals, 'Wscc');

%         addToExistingFigures = ~useShuffledPhaseControl;
        addToExistingFigures = 0;
        addLinesAtMean = 0;
                
        stackMeasuresVertically = 0;
        showMeasuresInSeparateFigures = 0;
        
        pvals_inTitle = {'t-test', 'chi-sqr', 'rand_mean'};
        doAddSubplotLetters = 1;        
            m_spacing = [0 0 0];
            n_spacing = [0 0 0];

        showStatsOfControl = 0;
%         ctrl_col = iff(addToExistingFigures, [0 .85 0], [1 0 0]);
        
        xlabelFontSize = 13;
        ylabelFontSize = 11;

%         if useShuffledPhaseControl
%             newplots = 1;
%             plotWcc = 1;
%             ctrl_color = 'r';
%         else
%             newplots = 1;
%             plotWcc = 0;
%             ctrl_color = 'k';
%         end

%         cellF1oDC = [allCells.F1oDC];
%         Wcc_pairF1oDC = cellF1oDC(Wcc_pairs);       
%         Wcc_nSimp = sum(Wcc_pairF1oDC > 1, 2);
%         idx_ss_wcc = Wcc_nSimp == 2;

       if any(doFigs == [1, 2]) 
           %%

           allMeasuresCanPlot = {'cc',  'dphi', 'dF1', 'MID_cc', 'MID_ovlp_cc', 'MID_fit_cc', 'STA_cc'};
           allMeasureLabels = {'cc',  '\Delta\phi', '\DeltaF1', 'cc_{MID}', 'cc_{MID-overlap}',  'cc_{Gabor}', 'cc_{STA}'};
           
           isRFmeasure = @(ms) any(strcmp(ms, {'MID_cc', 'MID_ovlp_cc', 'MID_fit_cc', 'STA_cc'}));
%            fig1_measures = {'cc', 'dphi', 'dF1'};
           
           fig1_measures = {'cc', 'dphi'};
%            fig2_measures = {'MID_cc', 'MID_ovlp_cc'};
%            fig2_measures = {'MID_cc', 'MID_fit_cc', 'MID_ovlp_cc'};
           fig2_measures = {'MID_cc', 'MID_fit_cc'};
           
           doMIDplots = 1;
           addCtrlValues = 1;
           addCtrlErrorBars = 1;


            colorBySimpleComplexPairs = 0;
            normalizeDistributionsTo1 = 0;
                normalizeDphiDistributionsTo1 = 1;
                rescaleDphiDistributionToExpectedFrac = 1;
            plotCumulativeDistributions = 1;

                
           if strcmp(gratingType, 'flashed') && ~curMatchDB && doMIDplots
%                measuresToPlot_allFigs = {fig1_measures, fig2_measures};
               measuresToPlot_allFigs = {fig1_measures};
%                measuresToPlot_allFigs = {fig2_measures};
           else
               measuresToPlot_allFigs = {fig1_measures};
           end
           fig_ids = [1001, 1002];
           
           for fig_i = 1:length(measuresToPlot_allFigs);
               
               measuresToPlot = measuresToPlot_allFigs{fig_i};
               redoCriteria = 1;
               displayPlot = 1;
               if redoCriteria      
                   
                    criteria = struct;
%                     nPhase_S = 4;                                         

                    % 1. CELL SELECTION CRITERIA
%                     criteria.n_phases = struct('op', @eq, 'value', );
%                     criteria.n_phases = struct('op', @eq, 'value', 60);
%                     GLF_overlap = 20; criteria.GLF_overlap = struct('op', @gt, 'value', GLF_overlap);
%                     criteria.negAmps_overlap = struct('op', @gt, 'value', amps_overlap);
                    minID = 10;   criteria.minID = struct('op', @gt, 'value', minID);

                    % PAIR SELECTION CRITERIA
%                     minFrac = 0.25;     criteria.loc_minFracOfMaxes     = minFrac;                    
%                     minF1oDC_cmp = 1;   criteria.loc_minF1oDC_cmp       = minF1oDC_cmp;

                    % simple/complex pairing
                    
%                     SCpref_type = 'maxR_avP_noSm';
%                     SCpref_type = 'maxR_maxP_noSm';
%                     SCpref_type = 'maxR_avP_sm';
%                     SCpref_type = 'maxR_maxP_sm';
%                     SCpref_type = 'prefR_avP_noSm';
                    SCpref_type = 'prefR_maxP_noSm';
%                     SCpref_type = 'prefR_avP_sm';
%                     SCpref_type = 'prefR_maxP_sm';

                    SC_f1odc_th = 0.25;
                    criteria.(['SCtype_pref_' SCpref_type]) = struct('op', @eq, 'value', 2);
                    criteria.(['F1oDC_' SCpref_type '_maxJackStd']) = struct('op', @lt, 'value', SC_f1odc_th);
                    
                    
                    % spatial frequency tuning similarity 
%                     criteria.D_spf_pref = struct('op', @gt, 'value', 1);
%                     criteria.Dw_spf = struct('op', @gt, 'value', 1);

                    % VALUE SELECTION CRITERIA - odd/even reproducibility
                    if strcmp(gratingType, 'flashed')
                        min_rsqr_oe = 0.25;
                        criteria.MID_cc.min_rsqr_oe = min_rsqr_oe;            
                        criteria.MID_ovlp_cc.min_rsqr_oe = min_rsqr_oe;    
                        
                        criteria.MID_fit_cc.min_rsqr_oe = min_rsqr_oe;    
                        criteria.MID_fit_cc.min_rsqr_fit = min_rsqr_oe;        
                    end
                    
                    if strcmp(curResponseType(''), 'gainCorrected')
                        criteria.minCOV = struct('op', @lt, 'value', .25);
                    end


                    % VALUE SELECTION CRITERIA - jackknive std err limits
                    doJackStdLimits = 1;
                    if doJackStdLimits
                        jackStd_cc_max = .25;
                        jackStd_dphi_max = 25;
                        jackStd_RF_cc = .25;

                        jackStdCrit.cc     = jackStd_cc_max;
                        jackStdCrit.dphi   = jackStd_dphi_max; 
                        jackStdCrit.MID_cc       = jackStd_RF_cc; 
                        jackStdCrit.MID_ovlp_cc  = jackStd_RF_cc; 
                        jackStdCrit.MID_fit_cc   = jackStd_RF_cc; 
                        jackStdCrit.STA_cc       = jackStd_RF_cc; 
                    else
                        jackStdCrit = struct;
                    end
                    
                    
               end
               
%                criteria_str = getCriteriaStr(criteria);
%                criteria = struct;
%                 pairData;

                locationsToPlot = locations;
                if isRFmeasure( measuresToPlot{1} )
                    locationsToPlot = locations(1);
                else
                    
%                     locationsToPlot = {'maxMinFracR3'};
%                     locationsToPlot = {'maxR1xR2'};
                    locationsToPlot = {'maxMinFracR'};
%                     locationsToPlot = {'maxMinFracR2'};
                end


    %             measuresToPlot = {'cc',  'dphi', 'MID_cc', 'MID_fit_cc'};
%                 measuresToPlot = {'cc', 'dphi'}; 
    %             measuresToPlot = {'MID_cc', 'MID_ovlp_cc', 'MID_fit_cc'};


                %             locationsToPlot = {'maxMinFracR', 'maxR1xR2'};
                

                measureLabels = allMeasureLabels( cellfun(@(ms) find(strcmp(ms, allMeasuresCanPlot),1), measuresToPlot) );

                measuresAvailable = cellfun(@(s) any(strcmp(s, measures)), measuresToPlot);
                measuresToPlot = measuresToPlot(measuresAvailable);
                measureLabels = measureLabels(measuresAvailable);

                nMeasuresPlot = length(measuresToPlot);                
                nMeasuresPlotInFig = iff(showMeasuresInSeparateFigures, 1, nMeasuresPlot);

                nLocationsPlot = length(locationsToPlot);
    %             cc_idx = find(strcmp(measures, 'cc'), 1);
    %             dphi_idx = find(strcmp(measures, 'dphi'), 1);


    %             [h_ax7, h_ax7_inset] = deal( zeros(1,nPlots) );
    %             [Wcc_binVals, Ctrl_binVals] = deal(cell(1,nPlots));
                if cmp_opts.usePosNegDphi
                    dPhiBinE = @(nBin) binCent2edge( linspace(-180, 180, nBin) );
                else
                    dPhiBinE = @(nBin) binCent2edge( linspace(0, 180, nBin) );
                end
    %             nLocations = length(locations);
                
                for loc_i = 1:nLocationsPlot
                    figId_base = 20+loc_i*10+(curGratingType-1)*1000 + fig_ids(fig_i) *50;
                    figure(figId_base); 
                    set(figId_base, 'name', locationsToPlot{loc_i})
                    if ~addToExistingFigures
                        clf;
                    end
                    loc_idx = find(strcmp(locationsToPlot{loc_i}, locations), 1);

                    for ms_i = 1:nMeasuresPlot
                        if showMeasuresInSeparateFigures && (ms_i > 1)
                            figure(figId_base+ms_i); clf;
                            set(figId_base+ms_i, 'name', locationsToPlot{loc_i})
                        end
                        
                        measure_i = measuresToPlot{ms_i};
                        ms_idx = find(strcmp(measure_i, measures), 1);    
                        null_mean_value = switchh(measure_i, {'cc', 'dphi', 'dF1'}, [0, 90, 90, 0]);
                        

                        Si = S{loc_idx, ms_idx};                     
                        Si_val_Wcc = Si.val(Wcc_pairIdxs,:,:);
                        n3 = size(Si.val,3);
                        nShuffle = max(cellfun(@(x) size(x,2), Si.shuffVal));

                        [tf_Wcc_pair_passedCrit, general_crit_str, specific_crit_str] = getPairIdxsThatPassCriteria(criteria, measure_i, loc_idx, Wcc_pairIdxs, pairData, allPossPhaseMeasures);
                        tf_Wcc_pair_passedCrit = tf_Wcc_pair_passedCrit(:,:, ones(1, n3));
                        tf_Wcc_vals_nonnan = ~isnan( Si_val_Wcc );
                           
                        tf_Wcc_val_passedCrit = tf_Wcc_pair_passedCrit & tf_Wcc_vals_nonnan;
                        doJackStdLimits_now = doJackStdLimits && isfield(jackStdCrit, measure_i) && isfield(Si, 'jackStd');
                        if doJackStdLimits_now
                            tf_Wcc_vals_jackStdOK = Si.jackStd(Wcc_pairIdxs, :,:)  < jackStdCrit.(measure_i);
                            tf_Wcc_idx_use = tf_Wcc_val_passedCrit & tf_Wcc_vals_jackStdOK;
                            
                            tf_Wcc_shuffVals_jackStdOK = Si.shuffJackStd(Wcc_pairIdxs, :,:)  < jackStdCrit.(measure_i);
                            tf_Wcc_shuffIdx_use = bsxfun(@and, tf_Wcc_val_passedCrit, tf_Wcc_shuffVals_jackStdOK);
                            
                            specific_crit_str = appendToStr(specific_crit_str, ...
                                getPairCriteriaStr('jackStdErr', struct('op', @lt, 'value', jackStdCrit.(measure_i))), '; ');
                        else
                            
                            tf_Wcc_idx_use = tf_Wcc_val_passedCrit;
                            tf_Wcc_shuffIdx_use = tf_Wcc_idx_use(:, ones(1, nShuffle), 1:n3);
                            
                        end
                        
                        nPairs_used = nnz(any(tf_Wcc_idx_use,3));
                        Wcc_vals_orig = Si_val_Wcc ( tf_Wcc_idx_use );
                        Wcc_vals = Wcc_vals_orig(:);

                        Wcc_pairIdxs_rep = Wcc_pairIdxs(:,:,ones(1,n3));
                        Wcc_pairIdxsUse = Wcc_pairIdxs_rep(tf_Wcc_idx_use);
                        %%
%                         Wcc_vals_orig;
%                         [uVals, valsCount] = uniqueCount(Wcc_vals_orig(:));
%                         ii = 1;
%                         nEach = sum(Wcc_vals_orig == uVals(ii), 2);
%                         idx = nEach > 0;
%                         coOccuring = setdiff(Wcc_vals_orig(idx,:), uVals(ii));
%                         figure(55); clf;
%                         hist(coOccuring, 60);
                        
                        %%
                        plotdPhiInIndivBins = 1;

                        binEdges = Si.binEdges;

                        nCCbins_drift = 26;
                        nCCbins_flashed = 18;

                        nDphiBins_drift = 12;
                        if plotdPhiInIndivBins
%                             uPhases = unique(pairData.n_phases(Wcc_pairIdxs(tf_Wcc_idx_use)));
                            uPhases = unique(pairData.n_phases(Wcc_pairIdxsUse));

                            if length(uPhases) == 1
                                nDphiBins_drift = uPhases/2;
                            else
                                nDphiBins_drift = 60;
                            end
                        
%                             nDphiBins_drift = 120;

                        end
                        if displayPlot
                            nDphiBins_drift = 12;
                        end
                        nDphiBins_flashed = 5;
    %                     nCCbins_drift = 40;
    %                     nCCbins_flashed = 22;
                        nMIDbins = 21;
                        nMID_fit_bins = 25;

                        if strncmp(locationsToPlot{loc_i}, 'wgtSum', 3)
                            nDphiBinFactor_drift = 2;
                            nDphiBinFactor_flashed = 5;
                            nCCBinFactor = 2;
                        else
                            nDphiBinFactor_drift = 1;
                            nDphiBinFactor_flashed = 1;
                            nCCBinFactor = 1;                        
                        end                    

                        switch gratingType
                            case 'drifting',
                                switch measure_i
                                    case {'cc', 'rho'},  binEdges = linspace(-1, 1, nCCbins_drift*nCCBinFactor);
                                    case {'dphi', 'dF1'},binEdges = dPhiBinE(nDphiBins_drift*nDphiBinFactor_drift+1);
                                end
                            case 'flashed',
                                switch measure_i
                                    case {'cc', 'rho'}, binEdges = linspace(-1, 1, nCCbins_flashed*nCCBinFactor);
                                    case {'dphi','dF1'}, binEdges = dPhiBinE(nDphiBins_flashed*nDphiBinFactor_flashed);
                                    case {'STA_cc', 'MID_cc', 'MID_ovlp_cc'}, binEdges = linspace(-1, 1, nMIDbins+1);                                
                                    case 'MID_fit_cc', binEdges = linspace(-1, 1, nMID_fit_bins+1);
                                end
                        end

    %                     binEdges = [-1:.05:1]; %S{dsi_si_idx}.binEdges;

                        dphi_measure = binEdges(end) > 50; % ie. is a dphi measure
                        binCents = binEdge2cent(binEdges);
            %         nBins = length(binEdges)-1;



                        [Wcc_cumFracBins, Wcc_binVals] = cumFracBinsFromVals(Wcc_vals, binEdges);
                        if colorBySimpleComplexPairs
                            [uSCtype, scTypeIdxs] = uniqueList( pairData.SCtype_cmp(Wcc_pairIdxsUse) );

                            [~, Wcc_binVals_col_C{1}] = cumFracBinsFromVals(Wcc_vals(scTypeIdxs{3}), binEdges);
                            [~, Wcc_binVals_col_C{2}] = cumFracBinsFromVals(Wcc_vals(scTypeIdxs{1}), binEdges);
                            [~, Wcc_binVals_col_C{3}] = cumFracBinsFromVals(Wcc_vals(scTypeIdxs{2}), binEdges);
                            Wcc_binVals_col = vertcat( Wcc_binVals_col_C{:})';
                        end
                        mean_Wcc_val = mean(Wcc_vals);
                        stderr_Wcc_val = stderr(Wcc_vals);
                        3;

                        allCtrlVals = [];
                        switch controlVals
                            case 'Bcc',                                
                                tf_Bcc_pair_passedCrit = getPairIdxsThatPassCriteria(criteria, measure_i, loc_idx, Bcc_pairIdxs, pairData, measures);
                                tf_Bcc_vals_nonnan = ~any( isnan(Si.val(  Bcc_pairIdxs, 1, :)),3);
                                tf_Bcc_idx_use = tf_Bcc_pair_passedCrit & tf_Bcc_vals_nonnan;
                                
                                Bcc_vals_orig = squeeze( Si.val(  Bcc_pairIdxs(tf_Bcc_idx_use), 1, : ) );
                                Bcc_vals = Bcc_vals_orig(:);
                                
                                [Bcc_cumFracBins, Bcc_binVals] = cumFracBinsFromVals(Bcc_vals, binEdges);
                                [Ctrl_cumFracBins, Ctrl_binVals] = deal(Bcc_cumFracBins, Bcc_binVals);                                
                                allCtrlVals = Bcc_vals;                                                                                                
                            case 'Wrcc', 

                                [Wrcc_cumFracBins, Wrcc_binVals] = cellfun(@(idxs) cumFracBinsFromVals(Si.val(idxs(tf_Wcc_idx_use)), binEdges), Wrcc_pairIdxs, 'un', 0);
                                rand_cumFracBins = cat(1, Wrcc_cumFracBins{:});
                                rand_binVals = cat(1, Wrcc_binVals{:});
                                                                
                                allCtrlVals_C = cellfun(@(idxs) single(Si.val(idxs(tf_Wcc_idx_use))), Wrcc_pairIdxs, 'un', 0);
                                allCtrlVals = cat(1, allCtrlVals_C{:});
%                                 allCtrlVals = Si.val;                                
                                

                            case 'Wscc', 
                                [Wscc_cumFracBins, Wscc_binVals] = deal( cell(1, Npermutes) );

%                                 tf_Wcc_idx_use
%                                 allShuffVals = vertcat(Si.shuffVal{Wcc_pairIdxsUse});
                                allShuffVals_all = vertcat(Si.shuffVal{Wcc_pairIdxs});
                                idx_have_shuffledVals =  ~cellfun(@isempty, Si.shuffVal) ;
                                tf_Wcc_shuffIdx_use_have = tf_Wcc_shuffIdx_use(idx_have_shuffledVals, :,:);
                                
                                if n3 == 1
                                    allShuffVals = [allShuffVals_all(tf_Wcc_shuffIdx_use_have)]; %;  allShuffVals_all(tf_Wcc_idx_use(:,:,2), :, 2)];
                                elseif n3 == 2
                                    allShuffVals = [allShuffVals_all(tf_Wcc_shuffIdx_use_have)]; %;  allShuffVals_all(tf_Wcc_idx_use(:,:,2), :, 2)];
%                                     allShuffVals = [allShuffVals_all(tf_Wcc_shuffIdx_use_have(:,:,1), :, 1);  allShuffVals_all(tf_Wcc_shuffIdx_use_have(:,:,2), :, 2)];
                                end
%                                 allShuffVals = allShuffVals_all
%                                 vertcat(Si.shuffVal{Wcc_pairIdxsUse});
                                if isempty(allShuffVals)
                                    addCtrlValues = 0;
                                    addCtrlErrorBars = 0;
                                else
                                    %%
                                    allCtrlVals_C = cell(1, Npermutes);
%                                     allCtrlVals_M = cell(1, Npermutes);
                                    for i = 1:Npermutes
                                        vals_permute_i = allShuffVals_all( tf_Wcc_shuffIdx_use_have(:,i), i,:);
                                        [Wscc_cumFracBins{i}, Wscc_binVals{i}] = cumFracBinsFromVals( vals_permute_i(:), binEdges);
                                        allCtrlVals_C{i} = vals_permute_i(:);
                                    end
                                    rand_cumFracBins = cat(1, Wscc_cumFracBins{:});
                                    rand_binVals = cat(1, Wscc_binVals{:});
                                    allCtrlVals = allShuffVals(:);

                                end

                        end
                        
                        mean_ctrl_val = nanmean(allCtrlVals);
                        stderr_ctrl_val = nanstderr(allCtrlVals);                        
                                                
                        if exist('rand_binVals', 'var')
                            Ctrl_binVals     = mean(rand_binVals, 1);
                            Ctrl_cumFracBins = mean(rand_cumFracBins, 1);
                            if addCtrlErrorBars
                                Ctrl_binVals_E     = std(rand_binVals, [], 1);
                                Ctrl_cumFracBins_E = std(rand_cumFracBins, [], 1);
                            end
                        else
                            Ctrl_binVals = [];
                            addCtrlErrorBars = 0;
                        end
                        
                        
%                         pvals = zeros(size(Wcc_binVals));
%                         for bi = 1:length(Wcc_binVals)
%                             pvals(bi) = getRandomizedProb(Wcc_binVals(bi), rand_binVals(:,bi), 'both');  
%                         end
%                         fprintf('min pval for %s : %g\n', measure_i, min(pvals))
                        
                        
                        
                        %%%%%% 1. Distribution of differences
                        
                        subM = iff(plotCumulativeDistributions, 2, 1);
                        ms_i_inFig = iff(showMeasuresInSeparateFigures, 1, ms_i);                            
                        
                        if ~stackMeasuresVertically                            
                            h_ax(ms_i, 1) = subplotGap(subM, nMeasuresPlotInFig, 1, ms_i_inFig);
                        else
                            h_ax(ms_i, 1) = subplotGap(nMeasuresPlotInFig, subM, ms_i_inFig, 1);
                        end
                        hold on; box on;                                        

    %                     h_ax(ms_i, 1) = subplot(2, nMeasuresPlotInFig, ms_i);  hold on; box on;                    
    %                     plot(binCents, Wcc_binVals/sum(Wcc_binVals), 'b.-'); hold on;
    %                     plot(binCents, Ctrl_binVals/sum(Ctrl_binVals), 'r.-');
                        if (~dphi_measure && normalizeDistributionsTo1) || ...
                           ( dphi_measure && normalizeDphiDistributionsTo1)
                            Wcc_vals_norm = Wcc_binVals/sum(Wcc_binVals);
                            sumCtrl_vals = sum(Ctrl_binVals);
                            Ctrl_vals_norm = Ctrl_binVals/sumCtrl_vals;
                        else
                            Wcc_vals_norm = Wcc_binVals;
                            sumCtrl_vals = sum(Ctrl_binVals);
                            Ctrl_vals_norm = Ctrl_binVals/sumCtrl_vals*sum(Wcc_vals_norm);
                        end
                        
                        if addCtrlErrorBars
                            Ctrl_vals_norm_E = Ctrl_binVals_E/sumCtrl_vals*sum(Wcc_vals_norm);
                        else
                            Ctrl_vals_norm_E = zeros(size(Ctrl_binVals));
                        end
                        if dphi_measure 
            %%
                            allNPhases = unique(pairData.n_phases);
                            [setNPhases_wcc, setNPhaseCounts_wcc] = uniqueCount([pairData.n_phases(Wcc_pairIdxsUse)]);
                            if size(Wcc_vals_orig, 2) == 2   % strcmp(cmp_opts.phase_oe_action, 'keepBoth')
                                setNPhaseCounts_wcc = setNPhaseCounts_wcc*2;
                            end;
                            [dPhiNull_val_wcc, dPhiNull_Count_wcc] = deltaPhiNull(setNPhases_wcc, setNPhaseCounts_wcc, allNPhases);
                            dPhiNull_count_wcc_binned = rebinnedDPhiCount(dPhiNull_val_wcc, dPhiNull_Count_wcc, binEdges);                            
              
                            if rescaleDphiDistributionToExpectedFrac
                                ylab = 'Rel. proportion of pairs';
                                if ~isempty(dphi_normalize)
                                    switch dphi_normalize
                                        case 'bcc', renormFactor  = Ctrl_vals_norm / sum(Ctrl_vals_norm);
                                        case 'wrcc', renormFactor = Ctrl_vals_norm / sum(Ctrl_vals_norm); 
                                        case 'null', renormFactor = dPhiNull_count_wcc_binned / sum(dPhiNull_count_wcc_binned);
                                        case 'none', renormFactor = 1;                  
                                    end
                                    Wcc_vals_norm = Wcc_vals_norm ./ renormFactor;
                                    Ctrl_vals_norm = Ctrl_vals_norm ./ renormFactor;
                                    if addCtrlErrorBars
        %                                 Wrcc_vals_norm = bsxfun(@rdivide, Wrcc_vals_norm, renormFactor);
                                        Ctrl_vals_norm_E = Ctrl_vals_norm_E ./ renormFactor;
                                    end
                                end
                            end
                        else
                            ylab = 'Number of pairs';
                        end
%                         if ~addToExistingFigures
                        
                        if ~colorBySimpleComplexPairs
                            bar(binCents, Wcc_vals_norm, 1, 'facecolor', bar_norm);
                        else
                            h = bar(binCents, Wcc_binVals_col, 1, 'b', 'stacked');
                            set(h(1), 'facecolor', 'k');
                            set(h(2), 'facecolor', 'w');
                            set(h(3), 'facecolor', [.6 .6 .6]);
                            3;
                        end
%                         end
                        
                        if addCtrlValues
                            stairs2(binEdges(1:end-1), Ctrl_vals_norm, ['-'], 'color', ctrl_col, 'linewidth', 2);
                        end

                        if addCtrlErrorBars                                               
                            h_e1 = errorbar(binCents, Ctrl_vals_norm, Ctrl_vals_norm_E, '.', 'color', ctrl_col);
                            removeFromLegend(h_e1);
                                                        
                        end                       

                        ylims = lims([Wcc_vals_norm, Ctrl_vals_norm + Ctrl_vals_norm_E], .1);
                        
                        addStarsAboveSigBins = 1;
                        
                        if addStarsAboveSigBins
                            idx_pos = find(Wcc_binVals);
                            for ii = idx_pos                                
                                pval = getRandomizedProb(Wcc_binVals(ii), rand_binVals(:,ii), 'both');
                                if pval < 0.05                                    
                                    text(binCents(ii), ylims(2), pvalNumStr(pval), 'horiz', 'cent', 'vert','top')
                                end
                                if ii == 2
%                                     idx_Wcc_idx_use = find(tf_Wcc_idx_use);
%                                     idx_6 = find( any(Wcc_vals_orig == 6, 2) );
%                                     ind_pr_w6 = idx_Wcc_idx_use(idx_6);
%                                     GCs = [pairData.Gids(ind_pr_w6, 1), pairData.cellIds(ind_pr_w6, :)];
                                    
                                    
%                                     Wcc_idxs(tf_Wcc_idx_use)
                                    
                                end
                            end
%                                 Wcc_tst = Wcc_binVals(ii)
%                             Wscc_binVals
                            

%                             Wscc_binCounts =     


                        end

                        if dphi_measure
                            if length(binCents) < 6
                                xticks = binCents;                
                            else
                                xticks = binCents(1):45:binCents(end);
                            end
                            set(gca, 'xtick', xticks);

                        end
                        xlim(binEdges([1, end]))                    
%                         ylim([0 inf]);

                        ylims = lims([Wcc_vals_norm, Ctrl_vals_norm + Ctrl_vals_norm_E], .1);
                        ylim([0 ylims(2)]);
                        ylabel(ylab, 'fontsize', ylabelFontSize);
                        loc_str = iff(any(strcmp(measure_i, RF_Measures)), '', locations{loc_idx});
%                         loc_str = '';
                        
%                         grating_str = iff(~any(strcmp(measure_i, RF_Measures)), sprintf('(%s gratings) ', titleCase(gratingType)), '');
                        grating_str =  sprintf('(%s gratings) ', titleCase(gratingType)) ;
                        nValues_used = nnz(tf_Wcc_idx_use);
                        nPairs_ctrl = length(allCtrlVals);
                        if n3 == 1
                            assert(nValues_used == nPairs_used);
                            nPairs_str = sprintf('(N = %d pairs)', nPairs_used);
                        elseif n3 == 2
                            nPairs_str = sprintf('(N = %d values; N = %d pairs)', nValues_used, nPairs_used);
                        end


                        fudge = 0;
%                         title_str = sprintf('%s Gratings : %s. N = %d', titleCase( gratingType ), loc_str, nValues_used);
                        title_str = sprintf('%s Gratings : %s', titleCase( gratingType ), loc_str);
    %                     title_str = sprintf('%s %s', measureLabels{ms_i}, grating_str);
    
                        title_str_C = {title_str};
                        addNPairsToTitle = true;
                        addExpDetailsToTitle = true;
                        if addNPairsToTitle
                            title_str_C = {title_str_C{:}, nPairs_str};
                        end
                        if addExpDetailsToTitle
                            [~, ~, ~, phase_oe_str] = curPhaseOEmode('');
                            [responseType_str] = titleCase(curResponseType(''));
                            title_str_C = {title_str_C{:}, [responseType_str, '; ', phase_oe_str]};
                        end
                        
                        title(title_str_C); %, 'interpreter', 'none');
%                         title({measureLabels{ms_i}, grating_str, nPairs_str}); %, 'interpreter', 'none');
                        %%
%                         pvals_inTitle = 'U-test';
                        pvals_inTitle = {'t-test', 'chi-sqr', 'rand_mean'};
                        blank_space = '{\fontsize{4} }';
                        pval_str_C = cell(1, length(pvals_inTitle));
                        for p_i = 1:length(pvals_inTitle)
                            switch pvals_inTitle{p_i}
                                case 't-test', [~, pval_t] = ttest(Wcc_vals, null_mean_value); pval_str_i = formattedPvalstr( pval_t/fudge, 'p_T' );
                                case 'U-test', [~, pval_U] = ranksum(Wcc_vals, allCtrlVals); pval_str_i = formattedPvalstr( pval_U, 'p_U' ); 
                                case 'chi-sqr', 
                                    if dphi_measure
                                        %%
                                        binnedChiSqr_DphiTest_drifting = 0;
                                        useBinnedVals = strcmp(gratingType, 'flashed') || binnedChiSqr_DphiTest_drifting;
                                        if useBinnedVals
                                            pval_chi2 = histChiSqrTest(  Wcc_binVals, dPhiNull_count_wcc_binned);
                                        else
                                            pval_chi2 = histChiSqrTest(  Wcc_vals, dPhiNull_val_wcc, dPhiNull_Count_wcc);
                                        end
                                        pval_str_i = formattedPvalstr( pval_chi2, 'p_X' ); 
                                                                                        
                                    else
                                        pval_str_i = {};
                                    end
                                    
                                case 'rand_mean', 
                                    %%
                                    if any(strcmp(controlVals, {'Wrcc', 'Wscc'}))
%                                     if 
                                        allRandMeans = cellfun(@nanmean, allCtrlVals_C);
                                        pval_R = getRandomizedProb(mean_Wcc_val, allRandMeans, 'both');
                                        pval_str_i = formattedPvalstr( pval_R/fudge, 'p_R' );                                            
                                    else
                                        pval_str_i = '';
                                    end
                            end
                            pval_str_C{p_i} = pval_str_i;
                        end
                        pval_str_C = pval_str_C(~cellfun(@isempty, pval_str_C));
                        pval_str = cellstr2csslist(pval_str_C);
                        data_color = iff(showStatsOfControl, 'blue', 'black');
                        
                        fmt = iff(dphi_measure, '%.1f', '%.2f');
                        deg = iff(dphi_measure, '�', '');
                        mean_std_str = sprintf(['\\color{%s} ' fmt deg ' \\pm ' fmt deg '   (N = %d pairs)'], data_color, mean_Wcc_val, stderr_Wcc_val, nValues_used);
                        
                        if showStatsOfControl
                            pvals_inTitle = 't-test';
                            switch pvals_inTitle
                                case 't-test', [~, pval] = ttest(allCtrlVals, null_mean_value); ctrl_pval_str = formattedPvalstr(pval, 'p_T');
%                                 case 'U-test', [~, pval] = ranksum(Wcc_vals, allCtrlVals); ctrl_pval_str = sprintf('p_U = %.2g', pval);
                            end
                            ctrl_std_str_C = {sprintf(['\\color{red} Ctrl: ' fmt deg ' \\pm ' fmt deg '   (N = %d pairs). %s'], mean_ctrl_val, stderr_ctrl_val, nPairs_ctrl, ctrl_pval_str)};
                            title_str_C = {[mean_std_str, '  ' pval_str], ctrl_std_str_C{:} };
                        else                            
                            title_str_C = {mean_std_str, blank_space, pval_str };
                        end
                            
                        
%                         title(title_str_C);
                        if ~plotCumulativeDistributions
                            xlabel([measureLabels{ms_i} ' ' grating_str]);
                        end
                        
                        mean_LineWidth = 2;
                        if addLinesAtMean                            
                            drawVerticalLine(mean_Wcc_val, 'linestyle', '--', 'linewidth', mean_LineWidth, 'color', 'b');
                            3;
                        
                            if showStatsOfControl
                                drawVerticalLine(null_mean_value, 'linestyle', '--', 'linewidth', mean_LineWidth, 'color', 'k');
                            end

                            drawVerticalLine(mean_ctrl_val, 'linestyle', '--', 'linewidth', mean_LineWidth, 'color', 'r');                            
                        end

                        if dphi_measure && strcmp(dphi_normalize, 'null') 
                            if rescaleDphiDistributionToExpectedFrac
                                drawHorizontalLine(1, 'linestyle', ':', 'color', 'k')
                            else
%                                 stairs(dPhiNull_count_wcc_binned)
                                stairs2(binEdges(1:end-1), dPhiNull_count_wcc_binned, [':'], 'color', 'k', 'linewidth', 2);
                            end
                        end
                        3;
                        
                        %%
                        if doAddSubplotLetters
                            addSubplotLetter(subM, nMeasuresPlotInFig, 1, ms_i_inFig, m_spacing, n_spacing, char('A'+ms_i-1));
                        end

                        3;
                        %%%%%% 2. Cumulative Distribution of differences
                        if plotCumulativeDistributions
                            if ~stackMeasuresVertically
                                h_ax(ms_i, 2) = subplotGap(2, nMeasuresPlotInFig, 2, ms_i_inFig);
                            else
                                h_ax(ms_i, 2) = subplotGap(nMeasuresPlotInFig, 2, ms_i_inFig, 2);
                            end
        %                     h_ax(ms_i, 2) = subplot(2, nMeasuresPlotInFig, 2);  
        %                     mWrcc_cumFracBins = mean(Wrcc_cumFracBins, 1);
                            if ~addToExistingFigures                                  
                                plot(binEdges, [0 Wcc_cumFracBins], 'o-',  'color', wcc_col_line, 'linewidth', wcc_line_w, 'markersize', wcc_sz); hold on;
                            end
                            if addCtrlValues
                                plot(binEdges, [0 Ctrl_cumFracBins], ['-'], 'color', ctrl_col, 'markersize', ctrl_sz, 'linewidth', ctrl_line_w);  
                            end
                            if addCtrlErrorBars                         
        %                         plot(binEdges, [0 mean(Wrcc_cumFracs, 1)], 'kd-');
                                h_e2 = errorbar(binEdges, [0 Ctrl_cumFracBins], [0 Ctrl_cumFracBins_E], 'color', ctrl_col, 'linestyle', 'none', 'linewidth', ctrl_err_w);
                                removeFromLegend(h_e2);
                                3;
                            end
                            xlim(binEdges([1, end]));
                            x_lab = measureLabels{ms_i};
                            addCriteriaToXlabel = true;
                            addStatisticsToXlabel = true;
                            xlabel_interp = 'tex';
                            if addCriteriaToXlabel
                                if strcmp(measureLabels{ms_i}, '\Delta\phi')
                                    x_lab = 'DPHI'; %'??';
                                end
                                x_lab = ['** ' upper(x_lab) ' **'];
                                general_crit_str_C = {}; specific_crit_str_C = {};
                                if ~isempty(general_crit_str)
                                    general_crit_str_C = strsplit(general_crit_str, '; ');
                                end
                                if ~isempty(specific_crit_str)
                                    specific_crit_str_C = strsplit(specific_crit_str, '; ');
                                end
                               x_lab = {x_lab, general_crit_str_C{:}, specific_crit_str_C{:}};  %#ok<AGROW>
                               xlabelFontSize = 8;
                               xlabel_interp = 'none';
                            end
                            
                            
                            if addStatisticsToXlabel
                             


            %                             vals = Si.val( pairTypeIdxs{Wcc_idx}(tf_Wcc_idx_use) );
%                                 vals = Si.val( tf_Wcc_idx_use, 1, : );
                                ctrl_dist = ks_getCtrlDist(allCtrlVals);

                                [Wcc_mean, Wcc_median, ~, Wcc_KS] = ...
                                    getMeanMedianStdKS(Wcc_vals, ctrl_dist);                          

%                                 [rand_means, rand_medians, ~, rand_KS] = ...
%                                     getMeanMedianStdKS(allCtrlVals_C, ctrl_dist);                                                                                                           
                                [rand_means, rand_medians, ~, rand_KS] = ...
                                    cellfun(@(vals) getMeanMedianStdKS(vals, ctrl_dist), allCtrlVals_C);
                    
%                                 measure_i_tmp = measure_i;
%                                 if ~isempty(strfind(measure_i_tmp, 'cc'))
%                                     measure_i_tmp = 'cc';
%                                 end
%                                 one_sided_tail = switchh(measure_i_tmp, {{'cc', 'rho'}, {'dphi', 'dF1', 'dPh_rel'}}, {'right', 'left'});
%                                 null_mean_value = switchh(measure_i_tmp, {{'cc', 'rho'}, {'dphi', 'dF1', 'dPh_rel'}}, [0, 90]);
                
                                medianProb = getRandomizedProb(Wcc_median, rand_medians, 'both');  
                                meanProb = getRandomizedProb(Wcc_mean, rand_means, 'both');
                                ksProb = getRandomizedProb(Wcc_KS, rand_KS, 'right'); % KS-statistic: alternative hypothesis is that KS-stat is larger.
  
                                median_stars = strrep(pvalStarStr (medianProb), ' ', '');
                                mean_stars = strrep(pvalStarStr (meanProb), ' ', '');
                                ks_stars = strrep(pvalStarStr (ksProb), ' ', '');
%                                 medianStars =
                                
                               	prob_str = sprintf('Median = %.3f (p = %.4f %s). Mean = %.3f (p = %.4f %s). KS = %.2f (p = %.4f %s)', ...
                                    Wcc_median, medianProb, median_stars, Wcc_mean, meanProb, mean_stars, Wcc_KS, ksProb, ks_stars);
                                
                                x_lab = {x_lab{:}, prob_str};
                            end
                                
                            xlabel(x_lab, 'fontsize', xlabelFontSize, 'interpreter', xlabel_interp);                    
                            ylim([0 1]);
                            if dphi_measure
                                set(gca, 'xtick', xticks);
                            end
        %                     if ms_i == 1
                                ylabel('Cum. frac. of pairs', 'fontsize', ylabelFontSize);
        %                     end

        %                     ht7 = title(sprintf('DSI%s', spont_incl_str));
        %                     set(ht7, 'fontsize', title_fsize)
                            
                            if (ms_i == nMeasuresPlot) && ~addToExistingFigures % && addLegends
                                if useShuffledPhaseControl
                                    leg_strs = {'Original', 'Phase-shuffled'};
                                else
                                    leg_strs = {'Within-site', 'Between-site (ctrl)'};
                                end                            
                                h_leg = legend(leg_strs, 'location', 'NW', 'fontsize', 9);
                            end

                            if addToExistingFigures
                                legend({'Within-site', 'Within-site : randomized', 'Between-site'}, 'location', 'NW', 'fontsize', 10);
                            end
                            if doAddSubplotLetters
                                addSubplotLetter(2, nMeasuresPlotInFig, 2, ms_i_inFig, m_spacing, n_spacing, char('C'+ms_i-1));
                            end

        %                     axis(h_ax7_inset(plot_i), 'tight'); 
        %                     xlim([0 1]);
        %                     ylims = get(h_ax7_inset(plot_i), 'ylim');
        %                     set(h_ax7_inset(plot_i), 'ylim', [0, ylims(2)*1.05]);

        %                     if (plot_i == nPlots) && addLegends
        %                         ax_inset_pos = get(h_ax7_inset(plot_i), 'position');                        
        %                         leg_pos = get(h_leg7, 'position');
        %                         LB = [ax_inset_pos(1), ax_inset_pos(2)+ax_inset_pos(4)+.02];
        %                         set(h_leg7, 'position', [LB, leg_pos(3:4)]);                        
        %                     end

        %                     figure(50+loc_i+(curGratingType-1)*1000); 
        %                     subplot(2,1,1); hold on;
        %                     col = iff(useShuffledPhaseControl, 'r', 'b');
        %                     errorbar(binCents, Ctrl_vals_norm, Ctrl_vals_norm_E, [col '.-']);
        %                     axis tight;
        %                     subplot(2,1,2); hold on;    
        %                     errorbar(binEdges, [0 Ctrl_cumFracBins], [0 Ctrl_cumFracBins_E], [col '.-']);
        %                     axis tight;
                        end
                        3;
                    end % for ms_i = 1:nMeasuresPlot
                    
                    3;
                    
                    
                    
                end % for loc_i = 1:nLocationsPlot

           end % fig_i = 1:length(measuresToPlot_allFigs);
           3;
           
        end
 
        
        
        
    end



end




function removeFromLegend(h)
    set(get(get(h,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); 
end

function [pairs_renumbered, pairsIdxs_selected] = renumberSubsetOfCellPairs(curPairs, curPairIdxs, selectedCells_idx, idxMtx)
    
    if iscell(curPairs)
        [pairs_renumbered, pairsIdxs_selected] = cellfun(@(curP, curPidxs) ...
            renumberSubsetOfCellPairs(curP, curPidxs, selectedCells_idx, idxMtx), curPairs, curPairIdxs, 'un', 0);
        return;
    end

    pairs_renumbered = binarySearch(selectedCells_idx, curPairs, [], 'exact');
    idxPairs_inSubset = find( all(pairs_renumbered, 2) > 0);
    pairs_renumbered = pairs_renumbered(idxPairs_inSubset, :);  % remove entries with '0'
    pairsIdxs_selected = idxMtx(curPairIdxs(idxPairs_inSubset));    
    
end


function [cumFracBins, binVals] = cumFracBinsFromVals(allVals, binEdges)
    
    binVals = histcnt(allVals(:), binEdges(:))';
    binVals = binVals(:)';
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

function [x_mean, x_std, x_median, x_p25, x_p75, n] = getMeanStdStats(x)
     x_mean = nanmean(x);   
     x_std = nanstd(x);
%      x_median = nanmedian(x);
%      x_p25 = prctile(x,25);
%      x_p75 = prctile(x,75);     
     x_prctiles = prctile(x, [25, 50, 75]);
        x_p25 = x_prctiles(1);
        x_median = x_prctiles(2);
        x_p75 = x_prctiles(3);
%         assert(isequal(x_pct, [x_median, x_p25, x_p75]));
     n = nnz(~isnan(x));    
end



function [x_mean, x_std, x_median, x_p25, x_p75, n] = getMeanStdPctile(x)
     x_mean = nanmean(x);   
     x_std = nanstd(x);
     x_median = nanmedian(x);
     x_p25 = prctile(x,25);
     x_p75 = prctile(x,75);
     n = nnz(~isnan(x));    
end

function [x_mean, x_median, x_std, ks_stat, n] = getMeanMedianStdKS(x, ctrl_dist, pairIdxs)
    doKS = (nargin >= 2) && ~isempty(ctrl_dist);
    
    if all(isnan(x(:))) % eg. MID_fit_cc, and didn't calculate
        [x_mean, x_median, x_std, ks_stat, n] = deal(nan);
        return;
    end
    
%     assert(~any(isnan(x(:))));
    %%
    if (nargin < 3) && isvector(x)
         x_mean = nanmean(x);
         x_median = nanmedian(x);
         x_std = nanstd(x);
         n = nnz(~isnan(x));    
         if doKS
            ks_stat = getKSstat(x, ctrl_dist);
         else
            ks_stat = 0;
         end

    else
       
        sizeX = size(x);
        x = single(x);
        if nargin == 3
            pairIdxs_cat = [pairIdxs{:}];            
            X_all = x(pairIdxs_cat);            
        elseif ismatrix(x)
            X_all = x;
        elseif sizeX(3) > 1
            X_all = [x(:,:,1); x(:,:,2)];            
        end
        nPermutes = size(X_all, 2);
        anyNans = any(isnan(X_all(:)));
        
        x_mean = nanmean(X_all, 1);
        x_median = nanmedian(X_all, 1);
        x_std = nanstd(X_all, [], 1);
        n = sum(~isnan(X_all), 1);
        if doKS
            useFastKS = ~anyNans;
            if useFastKS
                ks_stat = getKSstat_Mult(X_all, ctrl_dist);        
            else         
                ks_stat = zeros(1, nPermutes);
                for j = 1:nPermutes         
                    ks_stat(j) = getKSstat(nonnans(X_all(:,j)), ctrl_dist);
                end
            end
        else
            ks_stat = nan;
        end
        
     
     end
    
    
    
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

function pval = getRandomizedProb(val_wcc, val_permute, tail, return0forL0)
                    
    if (nargin < 4)
        return0forL0 = 1;
    end
        
    nPermutes = length(val_permute);
   
    useCorrectProbs = 1;
    if useCorrectProbs
        probFunc = @(L,N) iff(L==0 && return0forL0, 0, (L+1)/(N+2));
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

function pval_str = formattedPvalstr( pval, prefix )
    pval_str = '?';
    if pval > 0.1
        pval_str = sprintf('%s = %.1f', prefix, pval);
    elseif pval > 0.01
        pval_str = sprintf('%s = %.2f', prefix, pval);
    elseif pval > 0.001
        pval_str = sprintf('%s = %.3f', prefix, pval);
    elseif pval > 0
        pval_str = sprintf('%s = 10^{%d}', prefix, round( log10(pval) ));
    elseif pval == 0 
        if strcmp(prefix, 'p_R')
            pval_str = 'p_R < 10^{-4}';
        else
            pval_str = sprintf('%s < 10^{-20}', prefix);
        end
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

function [pairIdxs, pairIdxList, idxMtx] = useSubsetOfPairIdxs(allPT, pairTypes, nUnits)
    
    if isstruct(allPT)
        [    Wcc_pairIdxs,   Wcm_pairIdxs,   Bcc_pairIdxs,   Bcm_pairIdxs,   Bmm_pairIdxs,   Wrcc_pairIdxs,   Wrcm_pairIdxs] = ...
        deal(allPT.Wcc_idxs, allPT.Wcm_idxs, allPT.Bcc_idxs, allPT.Bcm_idxs, allPT.Bmm_idxs, allPT.Wrcc_idxs, allPT.Wrcm_idxs);   
        
        allPairIdxs = {Wcc_pairIdxs, Wrcc_pairIdxs, Bcc_pairIdxs,  Wcm_pairIdxs, Wrcm_pairIdxs, Bcm_pairIdxs, Bmm_pairIdxs,  Wcc_pairIdxs};
    elseif iscell(allPT)
        assert(length(allPT) == 6);
    end
    allPairTypes = {'Wcc', 'Wrcc', 'Bcc',   'Wcm', 'Wrcm', 'Bcm', 'Bmm',   'Wscc'};
    
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

function str = pvalStarStr(pval)
    if isnan(pval) || pval > 0.05
        str = '   ';
    elseif ibetween(pval, .01, .05)
        str = '  *';
    elseif ibetween(pval, .001, .01)
        str = ' **';
    elseif pval < .001
        str = '***';
    end
end

function str = pvalNumStr(pval)
    if isnan(pval) || pval > 0.05
        str = '';
    elseif ibetween(pval, .01, .05)
        str = '1';
    elseif ibetween(pval, .001, .01)
        str = '2';
    elseif pval < .001
        str = '3';
    elseif pval < .0002
        str = '4';
    end
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

    
%{
    -- mean & median ratios
                medianRatio = vals_median{Bcc_idx}/vals_median{Wcc_idx};
                L_median = nnz( vals_median{Wrcc_idx} <= vals_median{Wcc_idx});
                if useCorrectProbs
                    medianProb = L_median/Npermutes;
                else
                    medianProb = (L_median+1)/(Npermutes+2);
                end

                meanRatio = vals_mean{Bcc_idx}/vals_mean{Wcc_idx};                                                
                L_mean = nnz( vals_mean{Wrcc_idx} <= vals_mean{Wcc_idx});                
                if useCorrectProbs
                    meanProb = (L_mean+1)/(Npermutes+2);                    
                else
                    meanProb = L_mean/Npermutes;%(L_mean+1)/(Npermutes+2);                    
                end
                
%}