function printDegreeTuningCellsStats(suppressOutput_flag)

    global allID allLoc

    if exist('suppressOutput_flag', 'var') && isequal(suppressOutput_flag, 'all')
        curGratingType('f');
        curSubtractSpont(0); printDegreeTuningCellsStats(1);
        curSubtractSpont(1); printDegreeTuningCellsStats(1);
        curGratingType('d');
        curSubtractSpont(0); printDegreeTuningCellsStats(1);
        curSubtractSpont(1); printDegreeTuningCellsStats(1);
        return;
    end

    if exist('suppressOutput_flag', 'var') 
        printCellSelectionCriteriaStats = any(suppressOutput_flag == 1);
        printStatsOfGoodCells           = any(suppressOutput_flag == 2);        
        printCellsPerAnimalDistribs     = any(suppressOutput_flag == 3);
    else
        printCellSelectionCriteriaStats = 1;
        printStatsOfGoodCells = 1;
        printCellsPerAnimalDistribs = 1;
    end
    
    showWorking = nargin == 0 || isempty(suppressOutput_flag);
    
    minIsolationDistance = 10;
    useSigBckgResponseForSS = 0;
    
    % before doing this: run generateGratingCellsDatafile with 'all' selected

    curCmpType('degree');    
    gratingType = curGratingType('');  % FLASHED_GRATING = 1; DRIFTING_GRATING = 2;    
    subtractSpont = curSubtractSpont;
    bccType = curBccType();
%     [cmpType, cmpType_s] = curCmpType;        

    opt.applyStdErrorThresholds = true;
    opt.maxOriStdErr_deg = 5;
    opt.maxDSIStdErr = 0.1;
    opt.maxSpfStdErr = 0.5;
    opt.maxF1oDCErr = 0.25;


    ospDatafile = getFileName('osps', '_all');  % use file which includes non-reproducible/selective cells 
                                                % (but still has narrowed down to 1 group per site)         
    
    subtractSpont_str = iff(subtractSpont, '(Spontaneous subtracted)', '(Spontaneous included)');

    if showWorking
        fprintf(' ************************************* \n ****** %s GRATINGS %s ****** \n ***********************************', upper(curGratingType('')), subtractSpont_str);
    end

                                                    
    printMUproperties = 0;
%     printNCellsPerSiteStats = 1;
    
%     fprintf('Loading ... '); tic;    
    S1 = load(ospDatafile);        
    allCells = S1.allCells;
    clear S1;
%     nUnits = length(allCells);    

    switch gratingType, 
        case 'flashed'
            oriUnits = allCells;
            spfUnits = allCells;
            
        case 'drifting',
            oriUnits = allCells( arrayfun(@(s) strncmp(s.stimType, 'Grating:Orientation', 19), allCells ) );
            spfUnits = allCells( arrayfun(@(s) strncmp(s.stimType, 'Grating:Spatial Freq', 19), allCells ) );
    end
    
%     oriUsableCells = o

    allOriGrpCellIds = [oriUnits.Gid] * 10000 + [oriUnits.cellId];
    allOriGrpIds     = [oriUnits.Gid];
    allSpfGrpCellIds = [spfUnits.Gid] * 10000 + [spfUnits.cellId];
    allSpfGrpIds     = [spfUnits.Gid];
        
    if printMUproperties
        oriIdx_MU    = find([oriUnits.cellId] == 0);
    end
    %%
%     [963 + 557 + 470]/[1091 + 650 + 537]
    oriIdx_cells_tf = [oriUnits.cellId] > 0;
    spfIdx_cells_tf = [spfUnits.cellId] > 0;

    oriIdx_cells = find(oriIdx_cells_tf);
    spfIdx_cells = find(spfIdx_cells_tf);

%     allCells2 = allCells([allCells.cellId] > 0);
    oriUnits_ID = nestedFields(oriUnits, 'spkFeatures', 'IsolationDistance');
    spfUnits_ID = nestedFields(spfUnits, 'spkFeatures', 'IsolationDistance');
    
    oriCells_ID = oriUnits_ID(oriIdx_cells);
    spfCells_ID = spfUnits_ID(spfIdx_cells);
    
    oriCellsIsolated_frac = [nnz(oriCells_ID > minIsolationDistance), length(oriIdx_cells)];    
    spfCellsIsolated_frac = [nnz(spfCells_ID > minIsolationDistance), length(spfIdx_cells)];
    
    if showWorking
        fprintf('* Cell Isolation Stats:\n');
        fprintf('   Orientation Units:       %d / %d (%.1f%%) satisfied ID > %d\n', ...
            oriCellsIsolated_frac(1), oriCellsIsolated_frac(2), oriCellsIsolated_frac(1)/oriCellsIsolated_frac(2)*100, minIsolationDistance);
        fprintf('   Spatial Frequency Units: %d / %d (%.1f%%) satisfied ID > %d\n\n', ...
            spfCellsIsolated_frac(1), spfCellsIsolated_frac(2), spfCellsIsolated_frac(1)/spfCellsIsolated_frac(2)*100, minIsolationDistance);
    end
    3;
    
    oriIdx_cells_use = find(oriIdx_cells_tf & oriUnits_ID > minIsolationDistance);
    spfIdx_cells_use = find(spfIdx_cells_tf & spfUnits_ID > minIsolationDistance);
    
%%
    if ~subtractSpont
        oriStatField_use = 'oriStats_si';
        spfStatField_use = 'spfStats_si';
    else
        oriStatField_use = 'oriStats_ss';
        spfStatField_use = 'spfStats_ss';
    end               
    
    useSSforSelectionCriteria = 1;    
    
    if useSSforSelectionCriteria 
        oriStatField_crit = 'oriStats_ss';
        spfStatField_crit = 'spfStats_ss';
    else
        oriStatField_crit = oriStatField_use;
        spfStatField_crit = spfStatField_use;
    end

    oriCells = oriUnits(oriIdx_cells_use);
    oriTuningStats = [oriCells.tuningStats];
    oriCellStats_select = [oriTuningStats.(oriStatField_crit)];
    oriCellStats_use = [oriTuningStats.(oriStatField_use)];
    if printMUproperties
        oriMUStats_select = nestedFields(oriUnits(oriIdx_MU), 'tuningStats', oriStatField_crit);
    end
    
    spfCells = spfUnits(spfIdx_cells_use);
    spfTuningStats = [spfCells.tuningStats];
    spfCellStats_select = [spfTuningStats.(spfStatField_crit)]; 
    spfCellStats_use = [spfTuningStats.(spfStatField_use)]; 
    
                
%     w_ori_global_si = [oriCellStats_si.w_ori_global];
%     w_ori_global_ss = [oriCellStats_ss.w_ori_global];
    

    if printCellSelectionCriteriaStats
        cmpWithAl = 0;
        alpha = .01;
        minRsqr_ori = 0.4;
        minRsqr_spf = 0.4;
        
        
        % Cell selection criteria statistics
        ori_nCells = length(oriIdx_cells_use);
        idx_oriSelective    = [oriCellStats_select.ori_sel_pval] < alpha;
        idx_oriReproducible = [oriCellStats_select.ori_rep_pval] < alpha;
        idx_goodFit         = [oriCellStats_select.rsqr] > minRsqr_ori;

%         if subtractSpont
%             idx_oriSigResponse = [oriCellStats_select.response_size_pval] < alpha;
%         end
        idx_oriSelected = [oriCellStats_select.cellOK];

        frac_pct = @(i,n) sprintf('%3d/%3d (%.1f%%)', i, n, i/n*100);

    %     ori_nAccepted = nnz([oriCellStats_select
        if cmpWithAl
            fcode = '(vs %.1f%%)';
        else
            fcode = '%d';
        end
    
        if showWorking
            %%
            fprintf('\n\nOrientation Batches : Cell Selection criteria\n');
            fprintf([' Ori : OriSelective : %s ' fcode '\n'], frac_pct(  nnz( idx_oriSelective ), ori_nCells), iff(cmpWithAl, 293/360*100, []));
            fprintf([' Ori : Reproducible : %s ' fcode '\n'], frac_pct(  nnz( idx_oriReproducible ), ori_nCells), iff(cmpWithAl, 283/360*100, []));
            
            fprintf([' Ori : OriSelective & Reproducible : %s ' fcode '\n']', frac_pct(  nnz( idx_oriSelective & idx_oriReproducible ), ori_nCells), iff(cmpWithAl, 257/360*100, []));
            fprintf( ' Ori : OriSelective & Reproducible & good fit: %s\n', frac_pct(  nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit ), ori_nCells ) );
            assert(isequal(nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit ), nnz( idx_oriSelected )));
            
            fprintf([' Ori : Total currently used : %s ' fcode '\n']', frac_pct(  nnz( idx_oriSelected ), ori_nCells), iff(cmpWithAl, 257/360*100, []));
            
        end
        
        
        S_criteria_ori = struct('n_Total_ori', ori_nCells, 'n_Selective_ori',  nnz( idx_oriSelective ), 'n_Reproducible_ori', nnz( idx_oriReproducible ), ...
            'n_Sel_Rep_ori', nnz( idx_oriSelective & idx_oriReproducible ), 'n_Sel_Rep_Fit_ori', nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit ), ...
            'n_Used_ori', nnz(idx_oriSelected) ...
            ... 'pct_Used_ori', nnz( idx_oriSelected) / ori_nCells * 100
            );

        data_crit_S.S_criteria_ori = S_criteria_ori;          
    
            
                        
        
        

%         if subtractSpont
%             fprintf([' Ori : Sig Response : %s ' fcode '\n']', frac_pct(  nnz( idx_oriSigResponse ), ori_nCells), iff(cmpWithAl, 259/360*100, []));
%             fprintf([' Ori : Reproducible & Sig Response : %s ' fcode '\n'], frac_pct(  nnz( idx_oriReproducible & idx_oriSigResponse), ori_nCells), iff(cmpWithAl, 229/360*100, []));
%             fprintf([' Ori : Sig Response & OriSelective : %s ' fcode '\n'], frac_pct(  nnz( idx_oriSigResponse & idx_oriSelective), ori_nCells), iff(cmpWithAl, 238/360*100, []));
%             fprintf([' Ori : All 3       : %s ' fcode '\n'], frac_pct(  nnz( idx_oriSelective & idx_oriReproducible & idx_oriSigResponse), ori_nCells), iff(cmpWithAl, 217/360*100, []));
%         end
        
        
%         assert( nnz( idx_oriSelected ) == nnz( idx_oriSelective & idx_oriReproducible & idx_oriSigResponse) );        
        
    %%
%         for z = 1:length(spfCellStats_select), spfCellStats_select(z).stcProps.nFits = nan; end
            idx_withConstrains = find(arrayfun(@(s) isfield(s.stcProps, 'constraints'), spfCellStats_select));
            spfSTCProps = arrayfun(@(s) rmfieldIfExists(s.stcProps, 'constraints'), spfCellStats_select);

        %%
%             spfSTCProps = [spfCellStats_select.stcProps];
            spf_nCells = length(spfIdx_cells_use);
            spf_rep_pval = cat(1, spfCellStats_select.spf_rep_pval); spf_rep_pval = spf_rep_pval(:,1)';
            idx_spfReproducible = spf_rep_pval < alpha;
    %         idx_spfSigResponse = [spfCellStats_select.response_size_pval] < alpha;
            idx_spfGoodFit_rsqr = ([spfCellStats_select.rsqr] > minRsqr_spf);% & ~strncmp({spfCellStats_select.SLNfitResult}, 'Fitting', 5);
            idx_spfPrefInside = ([spfSTCProps.PrefDistFromStim] == 0);
            idx_spfGoodFit_rsqr_inside  = idx_spfGoodFit_rsqr & idx_spfPrefInside;

            idx_spfSelected = [spfCellStats_select.cellOK];        
            if showWorking
                %%
                fprintf('\n\nSpatial Frequency Batches : Cell Selection criteria\n');
                fprintf([' Spf : Reproducible : %s ' fcode '\n'], frac_pct(  nnz( idx_spfReproducible ), spf_nCells), iff(cmpWithAl, 110/193*100, []));        
                fprintf([' Spf : Reproducible & Good Fit (rsqr): %s \n'], frac_pct(  nnz( idx_spfReproducible & idx_spfGoodFit_rsqr), spf_nCells) );
                fprintf([' Spf : Reproducible & Good Fit (rsqr & not outside): %s \n'], frac_pct(  nnz( idx_spfReproducible & idx_spfGoodFit_rsqr_inside), spf_nCells) );
                fprintf([' Spf : Total currently used: %s ' fcode '\n'], frac_pct(  nnz( idx_spfSelected ), spf_nCells), iff(cmpWithAl, 99/193*100, []) );
            end
            assert( nnz( idx_spfSelected ) == nnz( idx_spfReproducible & idx_spfGoodFit_rsqr_inside ));
%                 find(  idx_spfSelected  & ~( idx_spfReproducible & idx_spfGoodFit_rsqr_inside ) )


%         if subtractSpont
%             fprintf([' Spf : Sig Response : %s ' fcode '\n'], frac_pct(  nnz( idx_spfSigResponse ), spf_nCells), iff(cmpWithAl, 144/193*100, []));        
%             fprintf([' Spf : Reproducible & Sig Response) : %s ' fcode '\n'], frac_pct(  nnz( idx_spfReproducible & idx_spfSigResponse), spf_nCells), iff(cmpWithAl, 99/193*100, []));
%         end

        S_criteria_spf = struct('n_total_spf', spf_nCells, 'n_reproducible_spf', nnz( idx_spfReproducible ), ...
            'n_rep_fit_spf', nnz( idx_spfReproducible & idx_spfGoodFit_rsqr_inside ), ...
            'n_used_spf', nnz( idx_spfSelected ) ...
            ... 'pct_used_spf', nnz( idx_spfSelected ) / spf_nCells * 100 ...            
             ...
        );        

        data_crit_S.S_criteria_spf = S_criteria_spf;
        
        data_crit_S.stats = mergeStructs(S_criteria_ori, S_criteria_spf);
        data_crit_S.allMeasureNames = {'stats'};
        data_crit_S.columns = fieldnames(data_crit_S.stats); %fieldnames(data_S.(allMeasureNames{i}));
        data_crit_S.cmpType = curCmpType('');
        data_crit_S.gratingType = curGratingType('');
        data_crit_S.subtractSpont = subtractSpont;
        
        
        criteria_fn = getFileName('cellSelectedStats', [], [], struct('gratingType', gratingType, 'cmpType', 'degree', 'subtractSpont', subtractSpont));
        save(criteria_fn, '-struct', 'data_crit_S');    
       
    end
%%
%     IDs_used_ori = oriUnits_ID(  oriIdx_cells_use (idx_oriSelected ) );
%     IDs_used_spf = spfUnits_ID(  spfIdx_cells_use (idx_spfSelected ) );
    
%     allID.(['IDs_used_ori_' gratingType(1)]) = IDs_used_ori;
%     allID.(['IDs_used_spf_' gratingType(1)]) = IDs_used_spf;
  %%  
    if printStatsOfGoodCells 
        
        doOnlySScells = 1;
        
        stats_matFile = getFileName('cellStats');
        data_S.cmpType = curCmpType('');
        data_S.gratingType = gratingType;
        data_S.subtractSpont = subtractSpont;
        
        if showWorking
            if doOnlySScells
                fprintf(' ******** ONLY SIMPLE CELLS **********\n ');
            end
            fprintf('\n\n  Parameter      |    Mean     |      Std    |    Median   |     P25     |     P75     |   N (N_sing)     |      P5     |    P95      | \n');
        end

    %     measureNames = {'w_ori_global_cell', 'w_ori_local_cell', 'w_ori_global_MU', 'w_ori_local_MU', 'dsi_cell', 'dsi_MU', 'spf_width', 'spf_pref'};
        
        if printMUproperties
            oriMU_measureNames = {'w_ori_global_MU', 'w_ori_local_MU'};
            dirMU_measureNames = {'DSI_global_MU', 'DSI_local_MU'};
        else
            oriMU_measureNames = {};
            dirMU_measureNames = {};
        end
    
        oriMeasureNames = {'w_ori_global', 'w_ori_local', oriMU_measureNames{:}};      %#ok<CCAT>
        dirMeasureNames = {'DSI_global', 'DSI_local', dirMU_measureNames{:}}; %#ok<CCAT>
        spfMeasureNames = {'w_spf', 'f_opt'};
        RnullMeasureNames = {'R_spont_abs', 'R90_total_abs', 'R90_stim_abs', 'R_spont_rel', 'R90_total_rel', 'R90_stim_rel'};
%         oriMeasureErrNames = {'w_ori_global_err', 'w_ori_local_err'};
        oriMeasureErrNames = {'ori_pref_err', 'w_ori_global_err', 'w_ori_local_err'};
        dirMeasureErrNames = {'DSI_global_err'};
        spfMeasureErrNames = {'w_spf_err', 'f_opt_err'};
        simpCompMeasureErrNames = {'oriF1oDC_err', 'spfF1oDC_err'};
        simpCompMeasureNames = {'oriF1oDC', 'spfF1oDC'};
        
        calculateErrorsAsFractions = false;
        
        allMeasureNames_eachStim = [oriMeasureNames, dirMeasureNames, spfMeasureNames, RnullMeasureNames, oriMeasureErrNames, dirMeasureErrNames, spfMeasureErrNames, simpCompMeasureErrNames, simpCompMeasureNames];    
        if ~printMUproperties
            is_MU_measure = cellfun(@(s) ~isempty(strfind(s, '_MU')), allMeasureNames_eachStim);
            allMeasureNames_eachStim = allMeasureNames_eachStim(~is_MU_measure);
        end
        
        %%
%         allMeasureNames = {'w_ori_global', 'w_ori_local', 'DSI_global', 'DSI_local', 'w_spf', 'f_opt', 'R_spont_abs', 'R90_total_abs', 'R90_stim_abs', 'R_spont_rel', 'R90_total_rel', 'R90_stim_rel'}
%         err_names = {     'w_ori_global_err', 'w_ori_global_err_rel', 'w_ori_local_err', 'DSI_global_err'
%     
%     'ori_pref_deg_err'
%             
%         'w_ori_global', 'w_ori_local', 'DSI_global', 'DSI_local', 'w_spf', 'f_opt', 'R_spont_abs', 'R90_total_abs', 'R90_stim_abs', 'R_spont_rel', 'R90_total_rel', 'R90_stim_rel'}
%         
%         h = createLookupTable(
                
        idx_oriSelected = [oriCellStats_select.cellOK];
        ori_cells_ok = [oriCellStats_select.cellOK] & [oriCellStats_use.cellOK];

        if doOnlySScells
            F1oDC_ori = [oriCellStats_use.F1oDC];
            ori_cells_ok = ori_cells_ok & [F1oDC_ori >= 1];
        end
        
        oriCellStats_use_ok = oriCellStats_use(ori_cells_ok);    
        oriCellStats_use_ok_errs = [oriCellStats_use_ok.error_jack];

        if printMUproperties
            ori_mu_ok = [oriMUStats.cellOK];
            oriMUStats_ok = oriMUStats(ori_mu_ok);        
        end
        oriGrpCellIds = allOriGrpCellIds(oriIdx_cells_use(ori_cells_ok));
        oriGrpIds     = allOriGrpIds(oriIdx_cells_use(ori_cells_ok));
        
        spf_cells_ok = [spfCellStats_select.cellOK];
        if doOnlySScells
            F1oDC_spf = [spfCellStats_use.F1oDC];
            spf_cells_ok = spf_cells_ok & [F1oDC_spf >= 1];
        end
        
        spfStats_use_ok = spfCellStats_use(spf_cells_ok);                
        spfStats_use_ok_errs = [spfStats_use_ok.error_jack];

        spfGrpCellIds = allSpfGrpCellIds(spfIdx_cells_use(spf_cells_ok));
        spfGrpIds = allSpfGrpIds(spfIdx_cells_use(spf_cells_ok));

        3;        
        
        for i = 1:length(allMeasureNames_eachStim)
            measureName = allMeasureNames_eachStim{i};
            measureName_display = measureName;

            if strcmp(gratingType, 'flashed') && any(strcmp(measureName, dirMeasureNames))
                continue;
            end
            
            typeNow = 0;
            switch measureName
                case oriMeasureNames{1},  
                    s = '*** ORIENTATION ***'; typeNow = 1;
                case dirMeasureNames{1},  
                    s = '*** DIRECTION *** '; typeNow = 1;
                case spfMeasureNames{1}, 
                    s = '*** SPATIAL FREQUENCY ***'; typeNow = 1;
                case RnullMeasureNames{1},
                    s = '*** RESPONSE AT NULL ORIENTATION ***'; typeNow = 1;
                case oriMeasureErrNames{1},
                    s = '*** Errors ***'; typeNow = 1;
            end
            if showWorking && typeNow
                fprintf('%s\n', s);
            end
            
            isOriOrDirStat = ~any(strcmp(measureName, [spfMeasureNames, spfMeasureErrNames]));            
            isMUmeasure = printMUproperties && ~isempty(strfind(measureName, 'MU'));            
            measureName = strrep(measureName, '_MU', '');            
            isErrMeasure = ~isempty(strfind(measureName, 'err'));
            isF1oDC_measure = ~isempty(strfind(measureName, 'F1oDC'));
            if isErrMeasure
                3;
            end
            
            if isF1oDC_measure 
                isOriOrDirStat = strncmp(measureName, 'ori', 3);
                measureName = measureName(4:end);
            end
%             [vals_mean, vals_std, vals_median, vals_P25, vals_P75, N, N_sing] = deal( 0 )cell(1,1) );
            

            N_sing = 0;
            if isOriOrDirStat                    

                if ~isMUmeasure
                    
                    valName = valMeasureName(measureName);
                    errName = errMeasureName(measureName);
                    
                    haveMeasure = isfield(oriCellStats_use_ok, valName);
                    if haveMeasure
                        if isfield(oriCellStats_use_ok(1).orig, valName)
                            vals_orig_all = nestedFields(oriCellStats_use_ok, 'orig', valName);
                        end
                        vals_orig = [oriCellStats_use_ok.(valName)];
                        vals = vals_orig;
                    end
                    haveErrMeasure = isfield(oriCellStats_use_ok_errs, errName);
                    if haveErrMeasure
                        errs_orig = [oriCellStats_use_ok_errs.(errName)];
                        errs = errs_orig;
                    end
                    
                    oriGrpIds_use = oriGrpIds;
                    if haveErrMeasure && opt.applyStdErrorThresholds
                        3;
                        switch valName
                            case {'ori_pref', 'w_ori_global', 'w_ori_local'}, idx_use = (errs <= opt.maxOriStdErr_deg);
                            case 'DSI_global',                                idx_use = (errs <= opt.maxDSIStdErr);
                            case 'F1oDC',                                     idx_use = (errs <= opt.maxF1oDCErr);
                            otherwise, error('Unknown name');
                        end
                        
                        errs = errs(idx_use);
                        vals = vals(idx_use);
                        assert(~any(isnan(vals)));
                        oriGrpIds_use = oriGrpIds_use(idx_use);
                    end
                    N_sites = length(unique(oriGrpIds_use));
                    
                    
                    if ~isErrMeasure
                        assert(haveMeasure);
                        vals_use = vals; %[oriCellStats_use_ok.(measureName)];
                                                
                                            
                    elseif isErrMeasure
                        assert(haveErrMeasure);

                        vals_use = errs;
                        
                        if calculateErrorsAsFractions && ~strcmp(measureName_use, 'ori_pref')

                            show = 1;
                            if show
                                %%
                                figure(i);
                                plot(errs, vals, '.');
                                xlabel(measureName, 'interp', 'none'); ylabel( measureName_use, 'interp', 'none');
                                [cc, cc_p] = corr(real(errs(:)), vals(:), 'type', 'spearman');
                                title({[titleCase(gratingType) ' gratings : ' measureName_use], sprintf('cc(spearman) = %.2f. p = %.2g', cc, cc_p)}, 'interp', 'none')
                                3;
                                
                            end
                            3;
                            
%                             fprintf('%s\n', measureName)
%                             [cc, pcc] = corr(vals(:), vals_measured(:))
                            vals_use = errs ./ vals;
                            
                        end
                        
                    end
%                     v = [oriCellStats_use_ok.([measureName '_err'])];
                    [uGrp, uGrpCount] = uniqueCount(oriGrpIds);
                    N_sing = nnz(uGrpCount == 1);
                    
                else
                    
                    
                    vals_use = [oriMUStats_ok.(measureName)];
                end

            elseif ~isOriOrDirStat  % is spf stat
                
                valName = valMeasureName(measureName);
                errName = errMeasureName(measureName);
                                                
                
                haveMeasure = isfield(spfStats_use_ok, valName);
                if haveMeasure
                    vals_orig_all = nestedFields(spfStats_use_ok, 'orig', valName);
                    vals_orig = [spfStats_use_ok.(valName)];
                    vals = vals_orig;
                end
                haveErrMeasure = isfield(spfStats_use_ok_errs, errName);
                if haveErrMeasure
                    errs_orig = [spfStats_use_ok_errs.(errName)];
                    errs = errs_orig;
                end
                
                spfGrpIds_use = spfGrpIds;
                if haveErrMeasure && opt.applyStdErrorThresholds
                    3;
                    switch valName
                        case {'w_spf', 'f_opt'}, idx_use = (errs <= opt.maxSpfStdErr);
                        case 'F1oDC',            idx_use = (errs <= opt.maxF1oDCErr);
                        otherwise, error('Unknown name');
                    end
                    
                    errs = errs(idx_use);
                    vals = vals(idx_use);
                    assert(~any(isnan(vals)));
                    spfGrpIds_use = spfGrpIds_use(idx_use);
                end
                N_sites = length(unique(spfGrpIds_use));
                
                
                if ~isErrMeasure
                    vals_use = vals;
                    
%                     if opt.applyStdErrorThresholds 
%                         vals_use = vals;
%                     end
                        
                    
                elseif isErrMeasure
                    vals_use = errs;% = [spfStats_use_ok_errs.(measureName_use)];
                    
                    if calculateErrorsAsFractions &&  ~strcmp(measureName_use, 'f_opt')
%                         vals_measured = [spfStats_use_ok.(measureName_use)];
                        vals_use = errs ./ vals_measured;
                        
                         
                            show = 1;
                            if show
                                %%
                                figure(i);
                                plot(vals, vals_measured, '.');
                                xlabel(measureName, 'interp', 'none'); ylabel( measureName_use, 'interp', 'none');
                                [cc, cc_p] = corr(vals(:), vals_measured(:), 'type', 'spearman');
                                title({[titleCase(gratingType) ' gratings : ' measureName_use], sprintf('cc(spearman) = %.2f. p = %.2g', cc, cc_p)}, 'interp', 'none')
                                3;
                                
                            end
%                         fprintf('%s\n', measureName)
%                         [cc, pcc] = corr(vals(:), vals_measured(:))
                        %                             vals = vals ./ vals_measured;
                    end
                    
                end
                if ~isMUmeasure
                    [uGrp, uGrpCount] = uniqueCount(spfGrpIds);
                    N_sing = nnz(uGrpCount == 1);
                end
            end
            [vals_mean, vals_std, vals_median, vals_P25, vals_P75, N, vals_P5, vals_P95] = getMeanStdPrctileStats(vals_use);  
                
            if isF1oDC_measure 
                3;
            end
            
            w = num2str( 11 );
            if strncmp(measureName, 'w_ori', 3) 
                fcode = ['%' w '.4f'];
            elseif strncmp(measureName, 'DSI', 3) || any(strcmp(measureName, RnullMeasureNames)) || any(strcmp(measureName, spfMeasureNames))
                fcode = ['%' w '.4f'];
            else
               3; 
            end

            fcode_tmp = ['| ' fcode  ' '];
            n_tmp = '%3d (%3d single)';

            str_template = ['%16s ' repmat( fcode_tmp, 1, 5) ' | ' n_tmp  repmat( fcode_tmp, 1, 2)  '\n'];        
            if showWorking
                fprintf(str_template, ...
                    measureName_display, vals_mean, vals_std, vals_median, vals_P25, vals_P75, N, N_sing, vals_P5, vals_P95);
            end
%             medianVar = 
            3;
            data_S.(measureName_display) = struct('mean', vals_mean, 'std', vals_std, 'vals_median', vals_median, 'vals_P25', vals_P25, 'vals_P75', vals_P75, 'N', N, 'N_sites', N_sites, ...
                                                    ... 'vals_P5', vals_P5, 'vals_P95', vals_P95, 
                                                    'n_total', length(vals_orig) ); 
            
        end
        3;
        %%
        all_w_ori_global = [oriCellStats_use_ok.w_ori_global];
        all_circVar = [oriCellStats_use_ok.circVar];
        all_sigma = (180/(2*pi)) * sqrt( 2 * all_circVar );
        idx_use_ori_widths = ~isnan(all_w_ori_global) & ~isnan(all_sigma);
        w_ori_global = all_w_ori_global(idx_use_ori_widths);
        sigma = all_sigma(idx_use_ori_widths);
        
        cc_ori = corr(w_ori_global(:), sigma(:));
        p_fit_ori = polyfit(w_ori_global, sigma, 1);
        miscTuningStats.WOriGlobal_sigma_Rsqr = cc_ori^2;
        miscTuningStats.WOriGlobal_sigma_fit_m = p_fit_ori(1);
        miscTuningStats.WOriGlobal_sigma_fit_c = p_fit_ori(2);
        
        all_w_spf = [spfCellStats_use.w_spf];
%         all_SLN_params = [spfCellStats_use.SLNparams];
        all_w_SLN = nestedFields(spfCellStats_use, 'SLNparams', 'w', 1);
        idx_use_spf_widths = ~isnan(all_w_spf) & ~isnan(all_w_SLN);

        w_spf = all_w_spf(idx_use_spf_widths);
        w_SLN = all_w_SLN(idx_use_spf_widths);
        
        
        cc_spf = corr(w_spf(:), w_SLN(:));
        p_fit_spf = polyfit(w_spf, w_SLN, 1);
        miscTuningStats.WSpf_WSLN_Rsqr = cc_spf^2;
        miscTuningStats.WSpf_WSLN_fit_m = p_fit_spf(1);
        miscTuningStats.WSpf_WSLN_fit_c = p_fit_spf(2);
        
        
        
        
        %%
        
        data_S.allMeasureNames = allMeasureNames_eachStim;
        data_S.oriMeasureNames = oriMeasureNames;
        data_S.dirMeasureNames = dirMeasureNames;
        data_S.spfMeasureNames = spfMeasureNames;
        data_S.RnullMeasureNames = RnullMeasureNames;
        data_S.columns = fieldnames(data_S.(allMeasureNames_eachStim{i}));
        data_S.miscStats.tuningStats = miscTuningStats;
        
        if ~doOnlySScells
            save(stats_matFile, '-struct', 'data_S');     
        end
                
        3;
    end
    
    
    if printCellsPerAnimalDistribs
        cellDistribs_matFile = getFileName('cellDistribs');
        cellDistribs_S.cmpType = 'degree';
        cellDistribs_S.gratingType = gratingType;
        cellDistribs_S.subtractSpont = subtractSpont;
        cellDistribs_S.bccType = bccType;
       3;
       
       
        oriCells_use = oriCells(ori_cells_ok);
        spfCells_use = spfCells(spf_cells_ok);

       
       % cells / site; sites/penetration, penetrations/animal.
       allMeasureNames = {};
       stimTypes = {'OR', 'SF'};
       for si = 1:length(stimTypes)
           stimType = stimTypes{si};
           switch stimType
               case 'OR',
                   allLocData = [oriCells_use.locData];
                   cellSiteIds = [oriCells_use.Gid];
                   
                   allLocData_withUnused = [oriUnits.locData];
               case 'SF',
                   allLocData = [spfCells_use.locData];
                   cellSiteIds = [spfCells_use.Gid];
                   
                   allLocData_withUnused = [spfUnits.locData];
           end
           cellLocIds = [allLocData.LocId];
           cellPenIds = [allLocData.PenId];
           cellCatIds = [allLocData.CatId];
           %%
           
           [uSiteIds, cellSiteIdxs, nCellsPerSite] = uniqueList(cellSiteIds);
           uLocIds = unique(cellLocIds);
           [uPenIds, cellPenIdxs, nCellsPerPen] = uniqueList(cellPenIds);
           [uCatIds, cellCatIdxs, nCellsPerCat] = uniqueList(cellCatIds); %#ok<ASGLU>

           idx_firstCellInSite = cellfun(@(x) x(1), cellSiteIdxs);          
           sitePenIds = cellPenIds(idx_firstCellInSite);
           siteCatIds = cellCatIds(idx_firstCellInSite);
                      
           [uPenIds2, sitePenIdxs, nSitesPerPen] = uniqueList(sitePenIds); %#ok<ASGLU>
           [uCatIds2, siteCatIdxs, nSitesPerCat] = uniqueList(siteCatIds);  %#ok<ASGLU>
           
           idx_firstCellInPen = cellfun(@(x) x(1), cellPenIdxs);          
           penCatIds = cellCatIds(idx_firstCellInPen);
           [uCatIds3, penCatIdxs, nPensPerCat] = uniqueList(penCatIds);  %#ok<ASGLU>
           
           assert(isequal(uPenIds, uPenIds2));
           assert(isequal(uCatIds, uCatIds2));
           assert(isequal(uCatIds, uCatIds3));
           fld_glob = [stimType '_' gratingType];
           
           
           
           prctiles_vals = [5, 50, 95];
%            cellsPerSite_medstats = prctile(nCellsPerSite, prctiles_vals);
%            cellsPerPen_medstats = prctile(nCellsPerPen, prctiles_vals);
%            cellsPerCat_medstats = prctile(nCellsPerCat, prctiles_vals);
%            sitesPerPen_medstats = prctile(nSitesPerPen, prctiles_vals);
%            sitesPerCat_medstats = prctile(nSitesPerCat, prctiles_vals);           
%            pensPerCat_medstats = prctile(nPensPerCat, prctiles_vals);
%       
           allMeasureNames_S = struct('cellsPerSite', nCellsPerSite, ...
                              'cellsPerPen', nCellsPerPen, ...
                              'cellsPerCat', nCellsPerCat, ...
                              'sitesPerPen', nSitesPerPen, ...
                              'sitesPerCat', nSitesPerCat, ...
                              'pensPerCat', nPensPerCat);
           allMeasureNames_eachStim = fieldnames(allMeasureNames_S);
           
           for mi = 1:length(allMeasureNames_eachStim);
               ms_name = allMeasureNames_eachStim{mi};
                x = allMeasureNames_S.(ms_name);
                [vals_p5, vals_med, vals_p95] = dealV(prctile(x, [5, 50, 95]));
                vals_mean = mean(x);
                vals_std = std(x);
                fld_name = [ms_name '_' stimType];
                allMeasureNames = [allMeasureNames, fld_name]; %#ok<AGROW>
                cellDistribs_S.(fld_name) = struct('median', vals_med, 'P5', vals_p5, 'P95', vals_p95, 'mean', vals_mean, 'std', vals_std ); 
                idx = strfind(ms_name, 'Per');
                nm1 = titleCase(ms_name(1:idx-1));
                nm2 = ms_name(idx+3:end);
                fprintf('[%s:%s] N %s / %s : %.2f : %.2f : %.2f [%.2f +- %.2f]\n', stimType, gratingType, nm1, nm2, vals_p5, vals_med, vals_p95, vals_mean, vals_std);
           
             
                miscDistribStats.(['stats_' stimType]) = ...
                    struct('nCells', length(cellSiteIds), ...
                           'nSites', length(uSiteIds), ...
                           'nPens', length(uPenIds), ...
                           'nCats', length(uCatIds));
                       
%                 miscDistribStats.(['nSites_' stimType]) = length(uSiteIds);
%                 miscDistribStats.(['nPens_' stimType]) = length(uPenIds);
%                 miscDistribStats.(['nCats_' stimType]) = length(uCatIds);
                
                     
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
            
         
           allLoc.(fld_glob).locIds_used = uLocIds;
           allLoc.(fld_glob).penIds_used = uPenIds;
           allLoc.(fld_glob).catIds_used = uCatIds;


           allLoc.(fld_glob).locIds_all = unique([allLocData_withUnused.LocId]);
           allLoc.(fld_glob).penIds_all = unique([allLocData_withUnused.PenId]);
           allLoc.(fld_glob).catIds_all = unique([allLocData_withUnused.CatId]);

            
       end
       
       cellDistribs_S.allMeasureNames = allMeasureNames;
       cellDistribs_S.columns = fieldnames(cellDistribs_S.(fld_name));
       cellDistribs_S.miscStats = miscDistribStats;
       save(cellDistribs_matFile, '-struct', 'cellDistribs_S');
       
    end
    
    if 0
        %%
        allLocIds_used = [];
        allPenIds_used = [];
        allCatIds_used = [];

        allLocIds_all = [];
        allPenIds_all = [];
        allCatIds_all = [];
        fn = fieldnames(allLoc);
        %%
        for i = 1:length(fn)
           allLocIds_used = unique([allLocIds_used, allLoc.(fn{i}).locIds_used]);
           allPenIds_used = unique([allPenIds_used, allLoc.(fn{i}).penIds_used]);
           allCatIds_used = unique([allCatIds_used, allLoc.(fn{i}).catIds_used]); 

           allLocIds_all = unique([allLocIds_all, allLoc.(fn{i}).locIds_all]);
           allPenIds_all = unique([allPenIds_all, allLoc.(fn{i}).penIds_all]);
           allCatIds_all = unique([allCatIds_all, allLoc.(fn{i}).catIds_all]); 
           
        end
        
        
        
        
    end
end
    


%     if length(fieldnames(allID)) == 4
%         3;
%         allIDs = [allID.IDs_used_ori_d, allID.IDs_used_ori_f, allID.IDs_used_spf_f];
%         
%         
%     end
%     








function [x_mean, x_std, x_median, x_p25, x_p75, n, x_p5, x_p95] = getMeanStdPrctileStats(x)
     x_mean = nanmean(x);   
     x_std = nanstd(x);
     x_median = nanmedian(x);
     x_p5 = prctile(x,5);
     x_p25 = prctile(x,25);
     x_p75 = prctile(x,75);
     x_p95 = prctile(x,95);

     p = [5, 25, 75, 95];
     x_ps  = prctile(x, p);
     x_p5  = x_ps(1);
     x_p25 = x_ps(2);
     x_p75 = x_ps(3);
     x_p95 = x_ps(4);
     n = nnz(~isnan(x));    
end

function s = rmfieldIfExists(s, fld)
    if isfield(s, fld)
        s = rmfield(s, fld);
    end
end



function valName = valMeasureName(measureName)

%     measureName = strrep(measureName, 'ori_w_', 'w_ori_');

    valName = strrep(measureName, '_err', '');

end

function errName = errMeasureName(measureName)
    errName = measureName;
%     errName = strrep(measureName, 'w_ori_', 'ori_w_');
    errName = strrep(errName, '_err', '');
%     if ~isempty(strfind(errName, '_err'))
%         errName = [measureName, '_err'];
%     end

end
%{
            if ~subtractSpont
                fprintf([' Ori : OriSelective & Reproducible : %s ' fcode '\n']', frac_pct(  nnz( idx_oriSelective & idx_oriReproducible ), ori_nCells), iff(cmpWithAl, 257/360*100, []));
                fprintf( ' Ori : OriSelective & Reproducible & good fit: %s\n', frac_pct(  nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit ), ori_nCells ) );
                assert(isequal(nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit ), nnz( idx_oriSelected )));
            else
                fprintf([' Ori : SigResponse: %s ' fcode '\n'], frac_pct(  nnz( idx_oriSigResponse ), ori_nCells), iff(cmpWithAl, 259/360*100, []));
                fprintf([' Ori : OriSelective & Reproducible : %s ' fcode '\n']', frac_pct(  nnz( idx_oriSelective & idx_oriReproducible ), ori_nCells), iff(cmpWithAl, 257/360*100, []));
                fprintf([' Ori : OriSelective & Reproducible & good fit: %s ' fcode '\n']', frac_pct(  nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit), ori_nCells), iff(cmpWithAl, 257/360*100, []));
                fprintf([' Ori : OriSelective & Reproducible & SigResponse: %s ' fcode '\n']', frac_pct(  nnz( idx_oriSelective & idx_oriReproducible & idx_oriSigResponse), ori_nCells), iff(cmpWithAl, 257/360*100, []));
                fprintf( ' Ori : OriSelective & Reproducible & SigResponse & good fit: %s\n', frac_pct(  nnz( idx_oriSelective & idx_oriReproducible & idx_oriSigResponse & idx_goodFit), ori_nCells ) );                            
                assert(isequal(nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit & idx_oriSigResponse ), nnz( idx_oriSelected )));
            end    




        if ~subtractSpont
                S_criteria_ori = struct('nOriTotal', ori_nCells, 'nOriSelective',  nnz( idx_oriSelective ), 'nOriReproducible', nnz( idx_oriReproducible ), ...
                    'nOriSel_Rep', nnz( idx_oriSelective & idx_oriReproducible ), 'nOriSel_Rep_Fit', nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit ), ...
                    'nUsed', nnz(idx_oriSelected), ...
                    'pctOriSel_Rep_Fit', nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit ) / ori_nCells * 100 );
            else
                S_criteria_ori = struct('nOriTotal', ori_nCells, 'nOriSigResponse', nnz( idx_oriSigResponse ), ...
                            'nOriSel_Rep_Fit', nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit ), ...
                            'nOriSel_Rep_Fit_Sig', nnz( idx_oriSelective & idx_oriReproducible & idx_goodFit & idx_oriSigResponse ));
            end






%}