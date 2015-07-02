function generateGratingCellsDatafile
%     global ori_yes ori_no spf_yes spf_no
    redo_jacks = 0;
    % before running this mfile, first call
%     load('indivCells_movie_fg.mat');
%     load('indivCells_grating.mat');  %% ( if using grating cells )
    global psthStatsSettings;
%     global smoothPhases
    if isempty(psthStatsSettings)
        setPsthGlobalSettings;
    end    
    
%     curTimeWindow([60 90]);
%     curTimeWindow([30 60]);
    curTimeWindow('best');
    
    gratingType = curGratingType('');    
    cmpType = curCmpType('');
    subtractSpont = curSubtractSpont;
    
    saveFileWithAllCells = true;
    saveFileWithOnlyUsableCells = true;
        
    saveSeparateOriSpfFilesForDegreeCmp = 1;
%     includeNonUsableCells = false;
    addR_full = true;
    addR_jackTrials = 1 && strcmp(cmpType, 'phase'); % && strcmp(curResponseType(''), 'raw');
        addRF_jackTrials = 1 && strcmp(cmpType, 'phase') && strcmp(gratingType, 'flashed'); % && strcmp(curResponseType(''), 'raw');
    addR_oeFields = true;
    
%     includeMultiUnits = strcmp(cmpType, 'degree');
    includeMultiUnits = true || strcmp(cmpType, 'degree');
    
    includeFirstCycleOfDriftingGratings = 0;
    
    
    responseType = curResponseType('');
%     useRgainCorrected = strcmp(cmpType, 'phase') && strcmp(responseType, 'gainCorrected');
    
        gainCorrected_p_Rcorr_av_gain_threshold = 0.01;

    if strcmp(cmpType, 'degree')
        multiUnitsSpikes = 'small spikes'; % 'all spikes' or 'small spikes';
    elseif strcmp(cmpType, 'phase')
        multiUnitsSpikes = 'all spikes'; % 'all spikes' or 'small spikes';
    end

    
    [timeWindow, timeWindow_str] = curTimeWindow; %  {'best', [29, 62], [58, 91], 'bestP1', 'bestP2', 'bestM1', 'bestM2', 'stimw'};
        
    include_STAs = strcmp(gratingType, 'flashed') && strcmp(cmpType, 'phase') && ~curMatchDB;

    include_MIDs = strcmp(gratingType, 'flashed') && strcmp(cmpType, 'phase') && ~curMatchDB;    
        MID_jcc_th = 0.4;
%         selectGroupThatMaximizesBestMids = 0;        
%         phase_selectGroupThatMaximizes = 'bestMIDs';
        phase_selectGroupThatMaximizes = 'reproducible';
%         phase_selectGroupThatMaximizes = 'selective'; % to match degree of tuning cells
        
        MID_trialModes = {'all', 'odd', 'even'}; %{'all', 'odd', 'even'};
        
%         RF_timeWindows_str = cellfun(@(w) iff(strcmp(w, 'best'), '', sprintf('__%d_%d', w)), RF_timeWindows, 'un', 0);        

    useR_predFromMID = 0;
        minR_pred_cc = .65;

    include_odd_even_STA_MID = 1;
    include_STAs_MIDs_ofMU = 1;
    
    includeTuningStats = true; %strcmp(cmpType, 'degree') || strcmp(responseType, 'raw');
%         useTuningStats_spont = 'current';
        useTuningStats_spont = 'ss';
    
    
    
    applyCellIsolationCriterion = 1;
    if strcmp(cmpType, 'degree')
        min_ID = 10;
    elseif strcmp(cmpType, 'phase')
        min_ID = 8;
    end

    updateUsedGellsFile = 0;
    addOStoTuningStats = strcmp(gratingType, 'flashed'); 

    addNumSitesAtEachPen = 1;

    windowOffset = curWindowOffset;
    
    if strcmp(cmpType, 'clusters')   % compare clusters
        clustGrp_str = 'clustersPruned';
        
        includeFlashedGratingCells = 1;
        includeDriftingGratingCells = 1;
        selectOneGroupFromEachLocation = false;
        includeMultiUnits = false;
        saveFileWithOnlyUsableCells = false;
    else
        clustGrp_str = 'cells';
        
        includeFlashedGratingCells = strcmp(gratingType, 'flashed');
        includeDriftingGratingCells = strcmp(gratingType, 'drifting');
        selectOneGroupFromEachLocation = true;
        
    end
    curGroupingType(clustGrp_str);
    
        
    excludeSitesWithOneCell = false;
    excludeCellsWithJust1Trial = true; % ie group 4191. only has 1 cell anyway.
    
    restrictToPresOKs = true;

    excludeCellsWithMultipleTempPeriods = true;
    excludeMockExperiments = true;
    excludeSquareGratings = true;
    excludeOrderedStimuli = true && strcmp(responseType, 'gainCorrected');
%     randomizeOrientations = false;
    
    phaseCmp_restrictToMultipleSpf_flashed  = true;
    phaseCmp_restrictToMultipleSpf_drifting = true;   
    
    
    addPhaseTC_CC_fields = strcmp(cmpType, 'phase');
        addPhaseTC_CC_pvals = 0;
        recalculatePhaseTC_CC_usingSmoothedR = 1;
    addDiffFromMUtoStats = true;
    
%     sortByStimType = false;
    sortByGid = true;

    includeCellFeatures = true;
    
    allWaveformFeatures = {'wvfm_scl_mean', 'wvfm_ccw_mean', 'negAmps_ccw_mean', 'negAmps_ccw_cov', 'spikeAmp', 'PCA', 'GLF', 'spkWidthHeight'};
    cellSeparationFeatures = {'IsolationDistance', 'L_ratio'};
    featuresToInclude = [ allWaveformFeatures, cellSeparationFeatures];

    useOspWindowDataFile = true;    
%         featuresToInclude = 'waveformPCA';
%         includeMeanSpikeMN = true;
%         includeMeanSpikeWvfm = true && ~curMatchDB;
    
%     rescaleR = true;

    reapplyStdErrorThresholds = true;
        maxOriStdErr_deg = 5;
        maxDSIStdErr = 0.1;
        maxSpfStdErr = 0.5;
        maxF1oDCErr = 0.25;
    tuningVarTh = struct('w_ori_global', maxOriStdErr_deg, 'w_ori_local', maxOriStdErr_deg, 'ori_pref_deg', maxOriStdErr_deg, ...
                         'dir_pref_deg', maxOriStdErr_deg, 'DSI_global', maxDSIStdErr, ...
                         'w_spf', maxSpfStdErr, 'f_opt', maxSpfStdErr, ...
                         'F1oDC', maxF1oDCErr);
            
    if strcmp(cmpType, 'degree')
%         siteGroupSelection_maximizes = 'reproducible&selective'; % Flashed: 301 usable cells. Drifting: 1164 cells
%         siteGroupSelection_maximizes = 'reproducible';         % Flashed: gives 350 usable cells. Drifting: 1211 cells
        
        cellInclusionCriteron = 'reproducible&selective';
        siteGroupSelection_maximizes = 'reproducible';
    elseif strcmp(cmpType, 'phase')
                
        if strcmp(gratingType, 'flashed') && strcmp(phase_selectGroupThatMaximizes, 'bestMIDs') && include_MIDs 
            cellInclusionCriteron = 'reproducible';
            negLog_cc_p_threshold = 3;
%             cellInclusionCriteron = 'MIDs_ok';            
            siteGroupSelection_maximizes = 'bestMIDs';
%             siteGroupSelection_maximizes = 'reproducible';
        elseif strcmp(phase_selectGroupThatMaximizes, 'selective') 
            cellInclusionCriteron = 'reproducible&selective';
            siteGroupSelection_maximizes = 'reproducible';
            
        elseif strcmp(phase_selectGroupThatMaximizes, 'reproducible') 
            negLog_cc_p_threshold = 3;
            cellInclusionCriteron = 'reproducible';
            siteGroupSelection_maximizes = 'reproducible';
        end
        
    end        
    
    
    
    
    
    
    [smoothPhase_Width, smoothPhase_Method] = curPhaseSmoothing;
    smoothPhases = ~isnan(smoothPhase_Width) && ~(strcmp(smoothPhase_Method, 'Gauss') && (smoothPhase_Width == 0) );    
    
    smoothPhases_now = smoothPhases && strcmp(gratingType, 'drifting') && strcmp(cmpType, 'phase');    
        smoothPhasesFunc = @(R) smoothOSP_Phs(R, smoothPhase_Method, smoothPhase_Width);
%         smoothPhase_Method = 'Gauss';
%         smoothPhase_Width = 1.5;
%         smoothPhase_Width = 3;

%         smoothPhase_Method = 'Fermi';        
%         k_hp = 5;
%         b_hp = .1;
%         smoothPhase_Width = [k_hp b_hp];

%         smoothPhase_Width = 3;

%     smoothPhases = []; % not doing test any more.
%     if ~isempty(smoothPhases)
%         compressPhases = true; %#ok<NASGU>
%         compressAction = smoothPhases{1};
%         compressN = smoothPhases{2};
%         ext = ['_' compressAction '_' num2str(compressN(1))];        
%         dir_ext = 'smooth_test\';
%     else
%         compressPhases = false; %#ok<NASGU>
%         ext = '';
%         dir_ext = '';
%     end
    filename_ext = '';

    if addR_full
%         filename_ext = [filename_ext '_full'];        
    end    
        
%     useEachCellOnlyOnce = false;
%     doRandomResampling = false;

%     addMaxF1oDCfields = false;
    includeF1oDCofAllStim = 1;

    flashed_restrictToNphases = [];  %#ok<NASGU> 
%     drifting_restrictToNphases = [40, 60]; %#ok<NASGU> %[480, 120, 58]; 

    if strcmp(cmpType, 'phase')
        drifting_restrictToNphases = [40, 60]; %#ok<NASGU> %[480, 120, 58]; 
    else
        drifting_restrictToNphases = []; %#ok<NASGU> %[480, 120, 58]; 
    end
                
    pathname = CatV1Path;
    mockExpGids = [3115, 3123, 1557];   %3115, 3123: tests.   1557: rat (?)
    
    % 1. Load group data info
    subtractSpont_str = iff(subtractSpont, ', (Spontaneous SUBTRACTED)','');
    fprintf('\n *** Creating data file for %s comparisons (%s gratings%s)... *** \n', cmpType, gratingType, subtractSpont_str);

    groupData = [];        
    if includeFlashedGratingCells  
        tic; fprintf('Loading flashed grating %s ... ', clustGrp_str);
        
        % group data
        fg_fname = getFileName('Groups', 'movie_fg');
        S_mcg = load(fg_fname);
        groupData = S_mcg.movieGroups_fg;

        % cell calculations        
        S_indivCells = load( getFileName('indiv', 'movie_fg') ); 
        fprintf(' done. '); toc;
        
        
        % include flashed grating "Grating" Experiments? (currently:no)        
    end
        
    if includeDriftingGratingCells
        tic; fprintf('Loading drifting grating %s ... ', clustGrp_str);
        
        % group data                        
        dOri_fname = getFileName('Groups', 'grating_dOr');
        dSpf_fname = getFileName('Groups', 'grating_dSf');
        S_ori = load(dOri_fname);                
        S_spf = load(dSpf_fname);                        
        groupData = [S_ori.gratingGroups_dOr; S_spf.gratingGroups_dSf];
                
        S_ori_cells = load( getFileName('indiv', 'grating_dOr') );
        S_spf_cells = load( getFileName('indiv', 'grating_dSf') );
        S_indivCells = mergeStructs(S_ori_cells, S_spf_cells);
        clear('S_ori_cells', 'S_spf_cells');
        fprintf(' done. '); toc;
        
    end
        
    % remove single spatial frequency experiments for phase tuning curve comparisons
    if strcmp(cmpType, 'phase') && ...
        ((strcmp(gratingType, 'flashed') && phaseCmp_restrictToMultipleSpf_flashed)  || ...
        (strcmp(gratingType, 'drifting') && phaseCmp_restrictToMultipleSpf_drifting))
        idx_multipleSpf = findInStructArray(groupData, 'spPeriod_pix', [], @(x) length(x) > 1);
        groupData = groupData( idx_multipleSpf ); 
    end         
    
    if restrictToPresOKs
        idx_presOK = findInStructArray(groupData, 'siteOK', 'ok');
        if any(~idx_presOK)
            3;
        end
        groupData = groupData( idx_presOK );
    end                          
    
    restrictToNphases = eval([gratingType '_restrictToNphases;']);
    if ~isempty(restrictToNphases)
        okNphases = arrayfun(@(s) any(length(s.spPh_deg) == restrictToNphases), groupData);
        if any(~okNphases)
            3;
        end
        groupData = groupData(okNphases);
    end
        
    if excludeSquareGratings
        isSquare = [groupData.isSquare];
        if any(isSquare)
            3;
        end
        groupData = groupData(~isSquare);
    end

    if excludeMockExperiments
        mockExp = arrayfun(@(s) any(length(s.Gid) == mockExpGids), groupData);
        if any(mockExp)
            3;
        end        
        groupData = groupData(~mockExp);
    end
    
    if excludeOrderedStimuli
        stimOrdered = [groupData.stimOrdered];
        groupData = groupData(~stimOrdered);
    end
        
    if selectOneGroupFromEachLocation
        locData = [groupData.locationData];
        locIds = [locData.LocId];
        
        % don't mix locations from different experiment types (ie. can have one ori batch and one
        % spf batch from the same location:)
        stimTypeCodes = arrayfun(@(s) switchh(s.stimType(1:11), {'Movie:Flash', 'Grating:Ori', 'Grating:Spa'}, [0 .1 .2]), groupData)';
        locIds = locIds + stimTypeCodes;
        
        [uLocIds, locLists] = uniqueList(locIds);

%         nSpikesTot = arrayfun(@(grp) sum(grp.nSpikes), groupData);   % find total # of spikes in each group at a particular location.
%         nCellsTot = arrayfun(@(grp) length(grp.nSpikes), groupData);   
                
        siteGroupSelection_threshold = switchh(siteGroupSelection_maximizes, ...
            {'reproducible', 'reproducible&selective', 'bestMIDs'}, [1, 3, MID_jcc_th]);
    
        cellRep = getCellReproducibilities(groupData, S_indivCells, siteGroupSelection_maximizes);
                
%                 score = nSpikesTot .* nCellsTot;
%         id_loc_max = cellfun(@(lst) indmax(score(lst)), locLists);
%         grps_selected = cellfun(@(lst,idx) lst(idx), locLists, num2cell(id_loc_max));

        id_loc_select = zeros(1, length(locLists));
        for loc_i = 1:length(locLists)
            id_loc_select(loc_i) = selectMostRepGroup(cellRep(locLists{loc_i}), siteGroupSelection_threshold);
        end
        
        3;
        
%         id_loc_max = cellfun(@(lst) indmax(score(lst)), locLists);
        grps_selected = cellfun(@(lst,idx) lst(idx), locLists, num2cell(id_loc_select));
        groupData = groupData(grps_selected);  
    end        

    if excludeSitesWithOneCell
        moreThanOneCell = arrayfun(@(s) length( setdiff(s.cellIds,0) ) > 1, groupData);
        if any(~moreThanOneCell)
            3;
        end    
        groupData = groupData(moreThanOneCell);
    end    
    
    
    if excludeCellsWithJust1Trial
        idx_oneTrial = arrayfun(@(s) s.Gid == 4191, groupData);
        assert(~any(idx_oneTrial))
%         groupData(idx_oneTrial) = [];
    end
    
    if excludeCellsWithMultipleTempPeriods
        idx_singleTempPeriod = arrayfun(@(s) length(s.tempPeriod_sec) == 1, groupData);
        assert(all(idx_singleTempPeriod));
%         groupData = groupData(idx_singleTempPeriod);
    end
    
    

    allCellIds = double([groupData.cellIds]);    
    if ~includeMultiUnits
%         allCellIds(allCellIds == 0) = [];            
        MU_cellId = nan;
    else       
        if curMatchDB
            multiUnitsSpikes = 'small spikes';
        end
        MU_cellId = switchh(multiUnitsSpikes, {'all spikes', 'small spikes'}, [100, 0]);        
    end
    
    nCells = length(allCellIds);
    
    % Create allCells struct array.
%     R_fullFields = iff(addR_fullFields, {'R_full', []}, {});   
%     R_oeFields = iff(addR_oeFields, {'R_oe', []}, {});   
    
    
    allCells = struct([]);
    cells_usable = false(1, nCells);
    cells_usable_ori = false(1, nCells);
    cells_usable_spf = false(1, nCells);
%     allCells = repmat(  struct('Gid', [], 'cellId', [], 'stimType', [], 'ori', [], 'sp', [], 'ph', [], 'ph_max', [], 'tp_sec', [], 'R', [], 'nspikes', [], 
%                         'maxR', [], 'oriSp_maxR_av', [],  'oriSp_maxR_no_av', [], ...
%                         'F1oDC_maxR_av',  [], 'F1oDC_maxR_no_av',  [], ...
%                         'ori_sp_maxMU', [], ...                    
%                        'degPerPix', [], 'stats', [], R_fullFields{:}, 'type', [], 'locData', [], 'PSTH', []),  nCells ,1 );                

    fprintf('\nGathering list of %s: |', clustGrp_str);
    progressBar('init', nCells, 20);
    cell_idx = 1;    
%     meanPhs = zeros(1,nCells);
%     meanPhsAll = zeros(1,nCells);

    
    if include_MIDs
        
        all_mids_filename = ['allMIDs' timeWindow_str '.mat'];

        S_mid = load(all_mids_filename);
        [mid_gids, mid_cellIds] = ...
            deal(S_mid.allGids, S_mid.allCellIds);        
%         [allMIDs, allMIDs_odd, allMIDs_even, mid_gids, mid_cellIds] = ...
%             deal(S_mid.allMIDs, S_mid.allMIDs_odd, S_mid.allMIDs_even, S_mid.allGids, S_mid.allCellIds);        

    end
    
    mean_jack_STA_ccs = nan(1000, 3, 1);
    
    
    if (strcmp(useTuningStats_spont, 'current') && subtractSpont) || strcmp(useTuningStats_spont, 'ss'); 
        oriStatFld_use = 'oriStats_ss';
        spfStatFld_use = 'spfStats_ss';
    elseif (strcmp(useTuningStats_spont, 'current') && ~subtractSpont) || strcmp(useTuningStats_spont, 'si'); 
        oriStatFld_use = 'oriStats_si';
        spfStatFld_use = 'spfStats_si';
    end               
    
    alsoUseSSforSelectionCriteria = 1;    
    
    if alsoUseSSforSelectionCriteria 
        oriStatFld_crit = 'oriStats_ss';
        spfStatFld_crit = 'spfStats_ss';
    else
        oriStatFld_crit = oriStatFld_use;
        spfStatFld_crit = spfStatFld_use;
    end
    tuningStats_date = 735993.628037;

    for Group_i = 1:length(groupData);
        grpData = groupData(Group_i);
        Gid = grpData.Gid;
        cellIds = double(grpData.cellIds);

%         if ~includeMultiUnits 
%             cellIds(cellIds == 0) = [];            
%         end
        oriPref_mu_allSpk = [];
        oriPref_mu_smlSpk = [];

        if strcmp(cmpType, 'clusters')
            cellClusterIds = getClusterCellAssignment(Gid);
        end
        
%         isSmallSpkMUatThisSite = any(cellIds == 0);
        ori_sp_maxMU = [];
        groupLocData = groupData(Group_i).locationData;
        
        if includeCellFeatures
%             groupSpikes = getSpikes(Gid, [], [], 1); % need this either way for cell_ids.
%             groupNegAmps = groupSpikes(:,3:end);            
            haveSomeWaveformFeatures = any(strCcmp(featuresToInclude, allWaveformFeatures));
            haveSomeCellSeparationFeatures = any(strCcmp(featuresToInclude, cellSeparationFeatures));

            if haveSomeWaveformFeatures
                S_wvfm = load(getFileName('mwvfm', Gid));
%                 grpWvfms = S_wvfm.meanWaveforms;                
                grpWvfms = S_wvfm.meanWaveforms;                
                nChannels = 4;
                nT = length(grpWvfms(1).wvfm_raw)/nChannels;
                idx_t0 = find(S_wvfm.t_ms == 0);
                idx_t0s = idx_t0 + [0:nChannels-1]*nT;
                
            end
            
            if haveSomeCellSeparationFeatures
                all_IDs = getIsolationDistance(Gid, cellIds, 'ID_all');
                all_IDs = [all_IDs.ID_allAP];
                all_Lratios = getIsolationDistance(Gid, cellIds, 'L_ratio_allAP');                
            end
               
                        
        end
        
        isDriftingOriBatch = strncmp(grpData.stimType, 'Grating:Orientation', 12);
        isSpatialFreqBatch = strncmp(grpData.stimType, 'Grating:Spatial', 12);
        isFlashedBatch = strncmp(grpData.stimType, 'Movie:Flashed', 12);        
        

        
        [oriPref_mu_smlSpk, dirPref_mu_smlSpk, oriPref_mu_allSpk, dirPref_mu_allSpk] = deal(nan);
        
        if addDiffFromMUtoStats && (isDriftingOriBatch || isFlashedBatch) && strcmp(cmpType, 'degree')

            if any(cellIds == 0)
                varname_small_spk_MU = getName('celldata', Gid, 0);
                v_small_PSTH = S_indivCells.(varname_small_spk_MU).PSTH;
                [L_bin_ssMU, R_bin_ssMU, windowProfile_ssMU] = getLRbins_windowprofile(timeWindow, v_small_PSTH);
                tuningStats_ssMU = getOspDataForPsthWindow(Gid, 0, [], [], L_bin_ssMU, R_bin_ssMU, windowProfile_ssMU, 'tuningStats', tuningStats_date);
                
                if ~allTuningStatsHaveOrig(tuningStats_ssMU)
                    tuningStats_ssMU = getOspDataForPsthWindow(Gid, 0, [], [], L_bin_ssMU, R_bin_ssMU, windowProfile_ssMU, 'tuningStats', 1);
                end

                oriStats_ss_MU = tuningStats_ssMU.(oriStatFld_use);
                ssMU_ok = oriStats_ss_MU.ori_sel_pval < .01 && oriStats_ss_MU.ori_rep_pval < .01;
                if ssMU_ok
                    oriPref_mu_smlSpk = oriStats_ss_MU.ori_pref_deg;
                    dirPref_mu_smlSpk = oriStats_ss_MU.dir_pref_deg;
                end
            end    



            if (any(cellIds == 100) || 1) && ~curMatchDB
                varname_all_spk_MU = getName('celldata', Gid, 100);
                v_all_PSTH = S_indivCells.(varname_all_spk_MU).PSTH;
                [L_bin_asMU, R_bin_asMU, windowProfile_asMU] = getLRbins_windowprofile(timeWindow, v_all_PSTH);
                tuningStats_asMU = getOspDataForPsthWindow(Gid, 100, [], [], L_bin_asMU, R_bin_asMU, windowProfile_asMU, 'tuningStats', tuningStats_date);

                if ~allTuningStatsHaveOrig(tuningStats_asMU)
                    tuningStats_asMU = getOspDataForPsthWindow(Gid, 100, [], [], L_bin_asMU, R_bin_asMU, windowProfile_asMU, 'tuningStats', 1);
                end

                oriStats_as_MU = tuningStats_asMU.(oriStatFld_use);
                asMU_ok = oriStats_as_MU.ori_sel_pval < .01 && oriStats_as_MU.ori_rep_pval < .01;
                if asMU_ok
                    oriPref_mu_allSpk = oriStats_as_MU.ori_pref_deg;
                    dirPref_mu_allSpk = oriStats_as_MU.dir_pref_deg;
                end
            end
            
        end                
        
        3;
%         maxR = max(v_all.OSP.R, [], 3);
%         [~, ori_sp_maxMU] = maxElement(maxR);             

%         if ~includeMultiUnits
%             cellIds = cellIds(cellIds > 0);
%         end
        if isequal(cellIds, 0)  % only multiunit, and no usable cells.
            continue; 
        end
        
        for cell_i = 1:length(cellIds) 
            progressBar;

            cellId = cellIds(cell_i); 
            isMultiUnit = cellId == 0;
            effCellId = iff(isMultiUnit, MU_cellId, cellId); % for retrieving saved data - to distinguish small spike multi-unit (0) from all-spikes multiunit(100)
            if isMultiUnit && ~includeMultiUnits
                continue;
            end
            
            if ~useOspWindowDataFile
                %%
                varname = getName('celldata', Gid, effCellId);                    
                celldata = S_indivCells.(varname);
                 
                PSTH = celldata.PSTH;

                OSP = celldata.OSP;                       
                [oris, sps, phs, degPerPix, tf_Hz] = deal(OSP.ori', OSP.sp, OSP.ph, OSP.degPerPix, OSP.tf_Hz);
            else
                PSTH = getPSTHforCell(Gid, effCellId);
                [oris, sps, phs, degPerPix, tf_Hz] = deal(grpData.ori_deg', grpData.spPeriod_pix, grpData.spPh_deg, grpData.stimulusInfo.degreesPerBlock, 1./grpData.tempPeriod_sec);
            end
           
                
            
%             [nOri, nSp, nPh] = deal(length(oris), length(sps), length(phs)); 
            [R, R_oe, R_full] = deal([]); %#ok<ASGLU>
            
            % Get Correct Time window for flashed gratings.
            [L_bin, R_bin, windowProfile] = getLRbins_windowprofile(timeWindow, PSTH);
            if isnumeric(timeWindow)
                windowProfile = [];
            end
            
            [L_bin_best, R_bin_best, windowProfile_best] = getLRbins_windowprofile('best', PSTH);
            
%              fprintf('Gid = %d, cellId = %d, L = %d, R = %d\n', Gid, effCellId, L_bin, R_bin);
            
            
 
            if strcmp(responseType, 'gainCorrected')
                getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, {'osp_ph'}, 0);
%                 continue;
%                 try
                if ~isMultiUnit
                    getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, {'osp_ph_oe', 'osp_ph_hoe'}, 0);    
                end
%                 catch Merr
%                     fprintf('error with this one');
%                 end
% continue
            end
            
            [windowStats, R] = getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, {'windowStats', 'osp_ph'}, 0);
            
            R_GC_stats_fields = {};
            if strcmp(responseType, 'gainCorrected')
                [R, R_corr_stats] = deal(double(decompress(R.R)), rmfield(R, 'R'));
                R_corr_stats = rmfield(R_corr_stats, 'opt');
                R_corr_stats = rmfield(R_corr_stats, 'stateFit');
                R_GC_stats_fields = {'GC_stats', R_corr_stats};

            end
                

%             tuningStats = [];
            if includeTuningStats
                curTimeWindow('best')
                tuningStats = getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin_best, R_bin_best, windowProfile_best, {'tuningStats'}, tuningStats_date);
                if ~allTuningStatsHaveOrig(tuningStats)
                    tuningStats = getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin_best, R_bin_best, windowProfile_best, {'tuningStats'}, 1);
                end
                curTimeWindow(timeWindow);
%                 tuningStats.( );
            end

            if includeTuningStats && reapplyStdErrorThresholds % && strcmp(cmpType, 'degree') 
                allTuningStatsNames = fieldnames(tuningStats);
                allTuningVarsWithTh = {'w_ori_global', 'w_ori_local', 'ori_pref_deg', 'dir_pref_deg', 'DSI_global', 'w_spf', 'f_opt', 'F1oDC'};
                for i = 1:length(allTuningStatsNames)
%                     isOriStat = ~isempty(strfind(allTuningStatsNames{i}, 'ori'));
                    S = tuningStats.(allTuningStatsNames{i});
                    if ~isfield(S, 'orig')
                        continue;
                    end
                    origVals = S.orig;
                    errVals = S.error_jack;
                    if isfield(errVals, 'ori_pref')
                        [errVals.ori_pref_deg, errVals.dir_pref_deg] = deal(errVals.ori_pref);
                    end
                    %%
                    for j = 1:length(allTuningVarsWithTh)
                        var_name = allTuningVarsWithTh{j};

                        if isfield(S, var_name)
                            val_orig = origVals.(var_name);
                            val_err = errVals.(var_name);
                            th = tuningVarTh.(var_name);

                            if val_err > th
                                S.(var_name) = nan;
                            else
                                S.(var_name) = val_orig;
                            end
                            3;
                        end


                    end

                    tuningStats.(allTuningStatsNames{i}) = S;
                end

            end
         
                
            cellWellIsolated = (~applyCellIsolationCriterion) || isMultiUnit || (all_IDs(cell_i) > min_ID);
            switch cellInclusionCriteron                
                case 'reproducible',
                    cell_ok = cellWellIsolated && windowStats.cc_p > negLog_cc_p_threshold;
%                     cell_ok = OSP.stats.isRep;
                case 'reproducible&selective'
                    % either ok for ori or for spf
                    cell_ok_ori = cellWellIsolated && isfield(tuningStats, oriStatFld_crit) && tuningStats.(oriStatFld_crit).cellOK ... 
                                                   && isfield(tuningStats, oriStatFld_use)  && tuningStats.(oriStatFld_use ).cellOK;
                    cell_ok_spf = cellWellIsolated && isfield(tuningStats, spfStatFld_crit) && tuningStats.(spfStatFld_crit).cellOK ...
                                                   && isfield(tuningStats, spfStatFld_use)  && tuningStats.(spfStatFld_use ).cellOK;                    
                    cell_ok = cell_ok_ori || cell_ok_spf;
%                     cellOKs = structfun(@(s) s.cellOK, tuningStats); 
%                     cell_ok = any(cellOKs);  % for flashed gratings: accept if either ok for ori OR for spf
                case 'MIDs_ok'                    
                    cell_ok = isfield(celldata.STAs, 'jackCC') && celldata.STAs.jackCC > 0;
                    
            end
            
            
                        
%             if isMultiUnit   % if multiunit (and are keeping multiunits), always keep.
%                 cell_ok = includeMultiUnits;
%             end   
            
            
            cells_usable(cell_idx) = cell_ok;
            if strcmp(cmpType, 'degree')
                cells_usable_ori(cell_idx) = cell_ok_ori;
                cells_usable_spf(cell_idx) = cell_ok_spf;
            end
                        
            if ~saveFileWithAllCells && ~cell_ok  % && isMultiUnit
                continue;                
            end
            
            % Get [R]response matrix
%             R = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin, R_bin, windowProfile, {'osp_ph'});
            
            if includeFirstCycleOfDriftingGratings && strcmp(gratingType, 'drifting');
%                 [allBins, psthVals] = dbGetCellSpkStimHists(Gid, effCellId, osp_opt);%                 
%                 [osp, osp_oe] = getOspDataForPsthWindow(Gid, effCellId, [], [], 1, 1, [], {'osp_ph', 'osp_ph_oe'});
                histArgs = {'keepFirstDriftingGratingCycle', 1};
                [R, R_oe] = calcOspForPsthWindow(Gid, effCellId, 1, 1, 0, [], {'osp_ph', 'osp_ph_oe'}, histArgs);
                
            end

            if smoothPhases_now %% && strcmp(gratingType, 'drifting') && strcmp(cmpType, 'phase')
                R = smoothPhasesFunc(R);
            end
             
            
            % Replace the original OSP response with that predicted from the MID? 
            if useR_predFromMID && isfield(OSP, 'R_pred')
                if OSP.R_pred_cc < minR_pred_cc  
                    continue; % skip this cell if MID wasn't a good predictor of cell responsese.  
                end                
                R = OSP.R_pred;                
                R_pred_cc = OSP.R_pred_cc;
            else                
                R_pred_cc = -1;
            end
            
            
            

            % R with odd/even trials
            R_oeFields = {}; %#ok<NASGU>
            if addR_oeFields 
                if ~isMultiUnit && isempty(R_oe) % ie. not multiunit                    
                    R_oe = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin, R_bin, windowProfile, {'osp_ph_oe'});
                    if strcmp(responseType, 'gainCorrected') && isstruct(R_oe)
                        R_oe = R_oe.R;
                        if isstruct(R_oe)
                            R_oe = decompress(R_oe);
                        end
                    end
                    if smoothPhases_now % && strcmp(gratingType, 'drifting') && strcmp(cmpType, 'phase')
                        R_oe = smoothPhasesFunc(R_oe);
                    end
                else
                    R_oe = [];
                end
                R_oeFields = {'R_oe', R_oe};
            end

            % R with ALL trials
            R_fullFields = {};
            R_jackTrialsFields = {};
            if addR_full || addR_jackTrials
                R_full = getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, {'osp_full'}, 0);
                if addR_full
                    R_fullFields = {'R_full', R_full};
                end
                
                if addR_jackTrials    
                    if smoothPhases_now
                        [R_jackTrials_aa, R_jackTrials_oe] = getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, {'osp_ph_jackTrials_sm', 'osp_ph_oe_jackTrials_sm'}, redo_jacks);
                    else
                        [R_jackTrials_aa, R_jackTrials_oe] = getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, {'osp_ph_jackTrials', 'osp_ph_oe_jackTrials'}, redo_jacks);
                    end
                                        
                    R_jackTrialsFields = {'R_jackTrials_aa', R_jackTrials_aa, 'R_jackTrials_oe', R_jackTrials_oe };
                end
            end

               
            
            stimType = [groupData(Group_i).stimType '_' num2str(round(groupData(Group_i).frameLength_ms))];
            if addDiffFromMUtoStats && (isDriftingOriBatch || isFlashedBatch) && strcmp(cmpType, 'degree');

                all_oddEven = {'', '_even', '_odd'};
                for oe_i = 1:length(all_oddEven)
                    oriStatFld_i = [oriStatFld_use all_oddEven{oe_i}];
                    if ~isfield(tuningStats, oriStatFld_i)
                        continue;
                    end
                    tuningStats.(oriStatFld_i).Dori_pref_allSpkMU = circDist(tuningStats.(oriStatFld_i).ori_pref_deg, oriPref_mu_allSpk, 180);
                    tuningStats.(oriStatFld_i).Dori_pref_smlSpkMU = circDist(tuningStats.(oriStatFld_i).ori_pref_deg, oriPref_mu_smlSpk, 180);

                    if isDriftingOriBatch
                        tuningStats.(oriStatFld_i).Ddir_pref_allSpkMU = circDist(tuningStats.(oriStatFld_i).dir_pref_deg, dirPref_mu_allSpk, 360);
                        tuningStats.(oriStatFld_i).Ddir_pref_smlSpkMU = circDist(tuningStats.(oriStatFld_i).dir_pref_deg, dirPref_mu_smlSpk, 360);
                    end
                end
            end
            
            if addOStoTuningStats && (isFlashedBatch) && strcmp(cmpType, 'degree')
                3;
                tuningStats.(oriStatFld_use).OS = mean(R,3);
                if ~isempty(R_oe)
                    tuningStats.([oriStatFld_use '_odd' ]).OS = mean(R_oe(:,:,:,1),3);
                    tuningStats.([oriStatFld_use '_even']).OS = mean(R_oe(:,:,:,2),3);                    
                end
                
                3;
                % add odd/even?
                
            end


            if strcmp(cmpType, 'clusters')
                cellIdForThisClust = find(cellfun(@(ids) any(ids == cellId), cellClusterIds));
                if isempty(cellIdForThisClust)
                    keyboard;
                    beep;
%                         continue; (?)
                end
                cellIdFields = {'clustId', cellId, 'cellId', cellIdForThisClust};                     
            else
                cellIdFields = {'cellId', cellId};
            end
            
            
            
            % Phase Tuning curve - extra data
            phaseTC_CC_fields = {};
            
            if addPhaseTC_CC_fields
%                 [phaseTC_CCs_oe,  phaseTC_CCs_hoe, phaseTC_CCs_fs, phaseTC_CC_ps_oe, phaseTC_CC_ps_hoe, phaseTC_CC_ps_fs] = deal([]); %#ok<ASGLU>
                if recalculatePhaseTC_CC_usingSmoothedR && smoothPhases_now
                    if isempty(R_full)
                        R_full = getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, {'osp_full'}, 0);
                    end
                    R_full = decompress(R_full);
                    R_full_sm = smoothPhasesFunc(R_full);
                    % stimPSTH_vals_allTrials
                    [phaseTC_CCs_oe,  phaseTC_CC_ps_oe ] = getAllPhaseTcCCs_sampled(R_full_sm, 'cc', 'odd_vs_even');
                    [phaseTC_CCs_hoe, phaseTC_CC_ps_hoe] = getAllPhaseTcCCs_sampled(R_full_sm, 'cc', 'half_odd_vs_half_even');
%                     [phaseTC_CCs_fs,  phaseTC_CC_ps_fs ] = getAllPhaseTcCCs_sampled(R_full_sm, 'cc', 'first_vs_second_half');
                else

                    [phaseTC_CCs_oe, phaseTC_CC_ps_oe,   phaseTC_CCs_hoe, phaseTC_CC_ps_hoe] = ...
                        getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, ...
                        {'phaseTC_CCs_oe', 'phaseTC_CC_ps_oe',   'phaseTC_CCs_hoe', 'phaseTC_CC_ps_hoe'});

                end
                
                phaseTC_CC_fields = {'phaseTC_CCs_oe',  single(phaseTC_CCs_oe), ...
                                     'phaseTC_CCs_hoe', single(phaseTC_CCs_hoe), ...  'phaseTC_CCs_fs',  single(phaseTC_CCs_fs)
                                     };          
                if addPhaseTC_CC_pvals || 0
                    phaseTC_CC_fields = [phaseTC_CC_fields, {'phaseTC_CC_ps_oe',  single(phaseTC_CC_ps_oe), ...
                                                             'phaseTC_CC_ps_hoe', single(phaseTC_CC_ps_hoe), ...  'phaseTC_CC_ps_fs',  single(phaseTC_CC_ps_fs)  
                                                             }];    %#ok<AGROW>
                end
            end

            [stimF1oDCs, stimF1oDC_jackStds] = getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, {'stimF1oDCs', 'stimF1oDC_jackStds'});
            if includeF1oDCofAllStim
                f1odc_fields = {'F1oDCs', single(stimF1oDCs), 'F1oDC_jackStds', stimF1oDC_jackStds};
%                     phaseTC_CC_fields = {'phaseTC_cc', single(phaseTC_cc), 'phaseTC_cc_p', single(phaseTC_cc_p)};

%                     [phaseTC_cc, phaseTC_cc_p, phaseTC_dot] = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin, R_bin, windowProfile, {'phaseTC_CCs', 'phaseTC_CC_ps', 'phaseTC_Dots'});
%                     phaseTC_CC_fields = {'phaseTC_cc', phaseTC_cc, 'phaseTC_cc_p', phaseTC_cc_p}; %, 'phaseTC_dot', phaseTC_dot};
            else
                f1odc_fields = {};
            end
            
            
            
            
%             F1oDC_S = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin, R_bin, windowProfile, 'F1oDC');
            idxs = getIdxPrefStim(Gid, R_full);
            pref_types = fieldnames(idxs);
            for i = 1:length(pref_types)
                os_idx = idxs.(pref_types{i});
                F1oDC_S.(['F1oDC_' pref_types{i}])         = stimF1oDCs(        os_idx(1),os_idx(2));
                F1oDC_S.(['F1oDC_jackStd_' pref_types{i}]) = stimF1oDC_jackStds(os_idx(1),os_idx(2));
            end
            F1oDC_S.types = pref_types;
                                    
            include_STAs_thisUnit = include_STAs && (include_STAs_MIDs_ofMU || ~isMultiUnit);
            switchDim3ToCell = @(X) reshape( mat2cell(X, size(X,1), size(X,2), ones(1, size(X,3))), [], size(X,3));
                
            if include_STAs
                STA_S = struct('L', nan);
                if include_STAs_thisUnit %&& isfield(celldata, 'STAs') && ~isempty(celldata.STAs) ;

                     [STA_aa, STA_oe] = getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, {'STA', 'STA_oe'});
                     STA_odd = double(STA_oe(:,:,1)); STA_even = double(STA_oe(:,:,2));
                     STA_oe_av = (STA_odd + STA_even)/2;
                     if strcmp(responseType, 'raw')
                         %%
                         STA_aa_d = double(STA_aa); 
                         maxDiff = max(abs(STA_aa_d(:)-STA_oe_av(:) ) );
                         cc = corr(STA_aa_d(:), STA_oe_av(:) );
                         
                         assert(maxDiff < 1e-5 || 1-cc < 1e-2);
                     end
                                          
                     STA_aa_S = struct('STA', single(STA_aa));
                     STA_odd_S = struct('STA', single(STA_odd));
                     STA_even_S = struct('STA', single(STA_even));
                     if addRF_jackTrials
                         
                         [STA_all_jacks, STA_odd_jacks, STA_even_jacks] = getOspDataForPsthWindow(Gid, effCellId, [], [], L_bin, R_bin, windowProfile, {'STA_jack', 'STA_odd_jack', 'STA_even_jack'}, redo_jacks);
                         STA_aa_S.STA_jacks = switchDim3ToCell( single(STA_all_jacks) );
                         STA_aa_S.jackCC = getRF_jackCC(STA_all_jacks);
                        
                         STA_odd_S.STA_jacks = switchDim3ToCell( single(STA_odd_jacks) );
                         STA_odd_S.jackCC = getRF_jackCC(STA_odd_jacks);
                         
                         STA_even_S.STA_jacks = switchDim3ToCell( single(STA_even_jacks) );
                         STA_even_S.jackCC = getRF_jackCC(STA_even_jacks);                         
                     end
                     
                     STA_S.STA = STA_aa_S;
                     STA_S.STA_odd = STA_odd_S;
                     STA_S.STA_even = STA_even_S;
                     [r_STA_oe, p_STA_oe] = pearsonR(STA_odd, STA_even);
                     STA_S.rsqr_oe = r_STA_oe^2;
                     STA_S.pval_oe = p_STA_oe;
                     STA_S.L       = length(STA_oe_av);
                     

                end
                STA_fields = {'STAdata', STA_S };                        
            else
                STA_fields = {};
            end

            

            if include_MIDs 
                MID_data = struct('L', nan);
                mid_idx = find(mid_gids == Gid & mid_cellIds == cellId, 1);
                if ~isempty(mid_idx)                                            

                    haveTrialModes = false(1, length(MID_trialModes));
                    for tm_i = 1:length(MID_trialModes)
                        trialMode_str = iff(strcmp(MID_trialModes{tm_i}, 'all'), '', ['_' MID_trialModes{tm_i}]);

                        fld_name_src = ['allMIDs' timeWindow_str trialMode_str];
                        fld_name_dst = ['MID' trialMode_str];
%                             trialMode_str = iff(strcmp(MID_trialModes{ti}, 'all'), '', ['_' MID_trialModes{ti}]);
                        if isfield(S_mid, fld_name_src)
                            Si = S_mid.(fld_name_src)(mid_idx);
                            if ~isempty(Si.MID)
                                haveTrialModes(tm_i) = true;
                                if isempty(Si.MID_fit) && ~isempty(Si.gparams)
                                    Si.MID_fit = single(getMIDfitToGabor(Gid, Si.gparams));
                                end                                        
                                MID_jack_fields = {}; 
                                if addRF_jackTrials
                                    
                                    if isempty(Si.MID_fit_jacks{1}) && ~isempty(Si.gparams_jacks{1})
                                        %%
                                        for jack_i = 1:length(Si.MID_fit_jacks)
                                            Si.MID_fit_jacks{jack_i} = single(getMIDfitToGabor(Gid, Si.gparams_jacks{jack_i}));
                                        end
                                    end
                                    MID_jack_fields = {'MID_jacks', {switchDim3ToCell( Si.MID_jacks )}, 'MID_fit_jacks', {Si.MID_fit_jacks}, ...
                                        'p_jacks', {Si.gparams_jacks} };
                                    
                                end

                                MID_data.(fld_name_dst) = struct(...
                                    'MID', single(Si.MID), 'jackCC', Si.jackMeanCC, 'MID_select', Si.above2std_gabor, 'MID_fit', Si.MID_fit, 'rsqr_fit', Si.rsqr, 'p', Si.gparams, ...
                                    MID_jack_fields{:});
                                MID_data.L = length(Si.MID);
                            end
                        end
                    end
                    idx_odd = strcmp(MID_trialModes, 'odd');
                    idx_even = strcmp(MID_trialModes, 'even');
                    if haveTrialModes(idx_odd) && haveTrialModes(idx_even)
                        M_odd = MID_data.MID_odd;
                        M_even = MID_data.MID_even;
                        idx_use = M_odd.MID_select & M_even.MID_select;
                        if isempty(idx_use) || (nnz(idx_use) / numel(idx_use) < .01)
                            [rsqr_oe_MID, pval_oe_MID] = deal(nan);
                        else
                            [r_oe_MID, pval_oe_MID] = pearsonR( M_odd.MID(idx_use), M_even.MID(idx_use) );
                            rsqr_oe_MID = r_oe_MID.^2;
                        end
                        MID_data.rsqr_oe = rsqr_oe_MID;
                        MID_data.pval_oe = pval_oe_MID;
                    end

%                         end                        
                end
                MID_fields = {'MIDdata', MID_data };                
                if isfield(MID_data, 'rsqr_oe') && ~isfield(MID_data, 'pval_oe')
                    3;
                end
            else
                MID_fields = {};
            end


            if includeCellFeatures                    
                spkFeatures = struct;

                if haveSomeWaveformFeatures
                    wvfm_idx = find([grpWvfms.cellId] == effCellId);

                    if any(strcmp(featuresToInclude, 'wvfm_scl_mean'))
                        spkFeatures.wvfm_scl_mean = grpWvfms(wvfm_idx).wvfm_scl;
                    end
                    if any(strcmp(featuresToInclude, 'wvfm_ccw_mean'))
                        spkFeatures.wvfm_ccw_mean = grpWvfms(wvfm_idx).wvfm_ccw;
                    end
                    if any(strcmp(featuresToInclude, 'negAmps_ccw_mean'))
                        spkFeatures.negAmps_ccw_mean = grpWvfms(wvfm_idx).wvfm_ccw(idx_t0s);
                    end                    
                    if any(strcmp(featuresToInclude, 'negAmps_ccw_cov'))
                        spkFeatures.negAmps_ccw_cov = grpWvfms(wvfm_idx).covNegAmp_ccw;
                    end                
                    if any(strcmp(featuresToInclude, 'PCA')) && isfield(grpWvfms(wvfm_idx), 'PCAcuw_mean')
                        spkFeatures.PCAcuw_mean = grpWvfms(wvfm_idx).PCAcuw_mean;
                        spkFeatures.PCAcuw_cov  = grpWvfms(wvfm_idx).PCAcuw_cov;
                    end      
                    if any(strcmp(featuresToInclude, 'GLF')) && isfield(grpWvfms(wvfm_idx), 'GLFcuw_mean')
                        spkFeatures.GLFcuw_mean = grpWvfms(wvfm_idx).GLFcuw_mean;
                        spkFeatures.GLFcuw_cov  = grpWvfms(wvfm_idx).GLFcuw_cov;
                    end

                    voltageFactor = 1/4;
                    rawNegAmps = grpWvfms(wvfm_idx).wvfm_raw(idx_t0s)*voltageFactor;
                    idx_bestChannel = indmax( abs(rawNegAmps));
                    if any(strcmp(featuresToInclude, 'spikeAmp'))
                        spkFeatures.spikeAmp = rawNegAmps( idx_bestChannel );
                    end

                    if any(strcmp(featuresToInclude, 'spkWidthHeight'))                            
                        spkWH = grpWvfms(wvfm_idx).spikeWidthHeight(idx_bestChannel);
                        spkFeatures.spkWidthHeight = spkWH;
                    end                        

                end

                if haveSomeCellSeparationFeatures

                    if any(strcmp(featuresToInclude, 'IsolationDistance'))
                        spkFeatures.IsolationDistance = all_IDs(cell_i);
                    end

                    if any(strcmp(featuresToInclude, 'L_ratio'))
                        spkFeatures.L_ratio = all_Lratios(cell_i);
                    end
                end                            

                spkFeatureFields = {'spkFeatures', spkFeatures};

            else
                spkFeatureFields = {};
            end

%                 oris = unique(oris)-min(oris);
            oris = mod(oris, 360);

            cellStruct = struct('Gid', Gid, cellIdFields{:}, 'stimType', stimType, ...
                'ori', oris, 'sp', sps, 'ph', phs, 'tf_Hz', tf_Hz, 'R', R, 'nspikes', grpData.nSpikes(cell_i), ... 
                ...'oriSp_maxR_avP', oriSp_maxR_avP, 'F1oDC_maxR_avP',  F1oDC_maxR_avP, 
                ...'F1oDC_maxR_avP_sm', cellF1oDC_opt_sm, 'F1oDC_jackStd', F1oDC_jackStd, ...                
                'F1oDC', F1oDC_S, ...
                 ...'ori_sp_maxMU', ori_sp_maxMU, 
                'degPerPix', degPerPix, 'windowStats', windowStats, 'tuningStats', tuningStats', ... 'stats', stats, ...
                'type', gratingType, 'locData', groupLocData, 'PSTH', PSTH, 'R_pred_cc', R_pred_cc, ...
                R_GC_stats_fields{:}, phaseTC_CC_fields{:}, f1odc_fields{:}, R_oeFields{:}, R_fullFields{:}, ...
                R_jackTrialsFields{:}, spkFeatureFields{:}, STA_fields{:}, MID_fields{:}, 'cell_ok', cell_ok );

            if isempty(allCells)
                clear('allCells');
                allCells(nCells,1) = blankStruct(cellStruct);
            end

            allCells(cell_idx) = cellStruct;                    
            

            cell_idx = cell_idx + 1;

        end
    end

%     progressBar('done');
    
    getOspDataForPsthWindow('save');
    
    allCells(cell_idx:end) = [];
    cells_usable(cell_idx:end) = [];
    cells_usable_ori(cell_idx:end) = [];
    cells_usable_spf(cell_idx:end) = [];
    
    
%     cellF1oDCs_opt = [allCells.F1oDC_maxR_avP];
%     cellF1oDCs_opt_sm = [allCells.F1oDC_maxR_avP_sm];    
    if isempty(allCells)
        fprintf('No Cells! exiting...');
        return;
    end

    if sortByGid
        
        [allCells, idx_newOrder] = sortByFields(allCells, {'Gid', 'cellId', 'clustId'});

    elseif sortByStimType
        nosp = cell2mat(arrayfun(@(s) [size(s.R)], allCells, 'un', 0));        
        [tmp, idx_newOrder] = sortrows(nosp);
        allCells = allCells(idx_newOrder);       
    end
    cells_usable = cells_usable(idx_newOrder);
    cells_usable_ori = cells_usable_ori(idx_newOrder);
    cells_usable_spf = cells_usable_spf(idx_newOrder);
    
    [pathname, filename_base] = getFileName('osps', filename_ext);
    filename_base = strrep(filename_base, '.mat', '');
    if saveFileWithAllCells 
        
        S_all.allCells = allCells;
        [pathname, filename_all] = getFileName('osps', [filename_ext '_all']);        
        stat_str = getStatStr(S_all, clustGrp_str);
        
        save([pathname filename_all], '-struct', 'S_all', '-v6');
        fprintf('\nSaved data for ALL %s grating %s (%s)\n      to file %s\n', ...
            gratingType, clustGrp_str, stat_str, filename_all);
        
    end
    if saveFileWithOnlyUsableCells
        
                
        if strcmp(cmpType, 'degree') && saveSeparateOriSpfFilesForDegreeCmp
            
%             ori_ok = nestedFields(allCells, 'stats', 'tuningStats', 'oriStats_si', 'cellOK');
%             spf_ok = nestedFields(allCells, 'stats', 'tuningStats', 'spfStats_si', 'cellOK');
% 
%             switch gratingType,
%                 case 'flashed',
%                     is_ori = nestedFields(allCells, 'stats', 'tuningStats', oriStatFld, 'cellOK');
%                     is_spf = nestedFields(allCells, 'stats', 'tuningStats', oriStatFld, 'cellOK');
%                 case 'drifting',
%                     is_ori =  strncmp({allCells.stimType}, 'Grating:Orientation', 12);
%                     is_spf =  strncmp({allCells.stimType}, 'Grating:Spatial Frequency', 12);
%             end                                                
            
            S_usable_ori.allCells = allCells(cells_usable_ori);       
            if addNumSitesAtEachPen
                S_usable_ori.allCells = addSiteStatisticsData(S_usable_ori.allCells);
            end
            stat_str = getStatStr(S_usable_ori, clustGrp_str);            
            [pathname, filename_usable_ori] = getFileName('osps', [filename_ext '_ori']);
            save([pathname filename_usable_ori], '-struct', 'S_usable_ori', '-v6');
            fprintf('\nSaved data for USABLE *ORI* %s grating %s (%s)\n       to file %s\n', ...
                gratingType, clustGrp_str, stat_str, filename_usable_ori);

            S_usable_spf.allCells = allCells(cells_usable_spf);            
            if addNumSitesAtEachPen
                 S_usable_spf.allCells = addSiteStatisticsData(S_usable_spf.allCells);
            end
            stat_str = getStatStr(S_usable_spf, clustGrp_str);            
            [pathname, filename_usable_spf] = getFileName('osps', [filename_ext '_spf']);
            save([pathname filename_usable_spf], '-struct', 'S_usable_spf', '-v6');
            fprintf('\nSaved data for USABLE *SPF* %s grating %s (%s)\n       to file %s\n', ...
                gratingType, clustGrp_str, stat_str, filename_usable_spf);            
            3;
            
        else % save all usable cells into a single file

            S_usable.allCells = allCells(cells_usable);
            stat_str = getStatStr(S_usable, clustGrp_str);
            filename_usable = [filename_base '.mat'];
            
            save([pathname filename_usable], '-struct', 'S_usable', '-v6');
            fprintf('\nSaved data for USABLE %s grating %s (%s)\n       to file %s\n', ...
                gratingType, clustGrp_str, stat_str, filename_usable);
            
        end
     
        
        if updateUsedGellsFile && strcmp(timeWindow, 'best') && strcmp(gratingType, 'flashed') && strcmp(cmpType, 'phase')
            
            %%
            usedFilename = [pathname 'usedCells.mat'];
            idx_usable_cells = [S_usable.allCells.cellId] > 0;
            S_used.usedGids = [S_usable.allCells(idx_usable_cells).Gid];
            S_used.usedCellIds = [S_usable.allCells(idx_usable_cells).cellId];
            
            idx_all_cells = [S_all.allCells.cellId] > 0;
            S_used.allGids = [S_all.allCells(idx_all_cells).Gid];
            S_used.allCellIds = [S_all.allCells(idx_all_cells).cellId];
            %%
            save(usedFilename, '-struct', 'S_used');
        end
        
    end
    
  
    toc;


end



function cellRep = getCellReproducibilities(groupData, S_indivCells, inclusionCriterion)
    nGrps = length(groupData);
    cellRep = cell(1,nGrps);    
    
    for grp_i = 1:nGrps
        cellIds = groupData(grp_i).cellIds;
        cellIds(cellIds == 0) = [];        
        
        cellRep{grp_i} = zeros(1,length(cellIds));        
        for cell_i = 1:length(cellIds)
            nm = getName('celldata', groupData(grp_i).Gid, cellIds(cell_i));            
            switch inclusionCriterion
                case 'reproducible', 
                    stats = S_indivCells.(nm).OSP.stats;
                    val = stats.allWindowStats.cc_p;
                case 'reproducible&selective'
                    stats = S_indivCells.(nm).OSP.stats;
                    if isfield(stats.tuningStats, 'oriStats_si')
                        val = stats.tuningStats.oriStats_si.cellOK;
                    elseif isfield(stats.tuningStats, 'spfStats_si')
                        val = stats.tuningStats.spfStats_si.cellOK;
                    end
                case 'bestMIDs'
                    STAs_fld = S_indivCells.(nm).STAs;
                    if isempty(STAs_fld) || ~isfield(STAs_fld, 'jackCC');
                        val = -1;
                    else                        
                        val = STAs_fld.jackCC;                    
                    end
            end            
            cellRep{grp_i}(cell_i) = val;
        end
    end
end
                
function grp_idx = selectMostRepGroup(cellReps, th)
    nGrps = length(cellReps);
    if nGrps == 1
        grp_idx = 1; 
    else
        nCellsAboveTh = cellfun(@(x) nnz(x > th), cellReps);
        [nMax, idx_max] = max(nCellsAboveTh);
        idx_GroupsWithMaxNRepCells = find(nCellsAboveTh == nMax);
        if length(idx_GroupsWithMaxNRepCells) == 1
            grp_idx = idx_max;
        else % multiple groups at this site tied for the highest # of cells above threshold. 
            % pick the one with the highest mean value.
            meanRep = cellfun(@(x) mean(x(x > th)), cellReps(idx_GroupsWithMaxNRepCells) );
            grp_idx = idx_GroupsWithMaxNRepCells ( indmax(meanRep) );
        end
    end
end


function [S2, idx_newOrder] = sortByFields(S, flds)
    Nrecs = length(S);
    Nflds = length(flds);
    vals = zeros(Nrecs, Nflds);        
    for i = 1:Nflds
        if ~isfield(S(1), flds{i})
            continue;
        end
        field_vals = [S(:).(flds{i})];
        vals(:,i) = field_vals;
    end
    [~, idx_newOrder] = sortrows(vals);
    S2 = S(idx_newOrder);
end

function stat_str = getStatStr(S, clustGrp_str)
        
        allGroupIds = [S.allCells.Gid];
        nSites = length(unique(allGroupIds));
        allCellIds = [S.allCells.cellId];
        nCells = nnz(allCellIds > 0);
        nMUs  = nnz(allCellIds == 0);
        
        stat_str = sprintf('Total of %d units (%d %s, %d multi-units) from %d sites', ...
            length(allCellIds), nCells, clustGrp_str, nMUs, nSites);
end




function allCells = addSiteStatisticsData(allCells)
%%
    locData = [allCells.locData];
    CatIds = [locData.CatId]';
    PenIds = [locData.PenId]';
    LocIds = [locData.LocId]';
    Gids = [allCells.Gid]';
        
    assert(length(unique(Gids)) == length(unique(LocIds)));
    cellIds = [allCells.cellId]';    
    
    idx_notMU = cellIds > 0;
    
    %%
    cell_nSitesAtPen = zeros(size(allCells));
    cell_BccFracAtPen = zeros(size(allCells));
    
    [uPenIds, penIdList] = uniqueList(PenIds);
    for pen_i = 1:length(uPenIds)
        idxUnits_thisPen = penIdList{pen_i};
        idxCells_thisPen = idxUnits_thisPen(idx_notMU(idxUnits_thisPen));
        
        Gids_thisPen = Gids(idxCells_thisPen);
        [uGids, gidCount] = uniqueCount(Gids_thisPen);
        
        BccFrac = getFracBetweenSitePairs(gidCount);        
        nSites_thisPen = length( uGids );
        
        cell_nSitesAtPen(idxUnits_thisPen) = nSites_thisPen;
        cell_BccFracAtPen(idxUnits_thisPen) = BccFrac;
        
    end

    
    cell_nSitesInAnimal = zeros(size(allCells));
    cell_BccFracInAnimal = zeros(size(allCells));
    [uCatIds, catIdList] = uniqueList(CatIds);
    for cat_i = 1:length(uCatIds)
        idxUnits_thisCat = catIdList{cat_i};
        idxCells_thisCat = idxUnits_thisCat(idx_notMU(idxUnits_thisCat));
        
        Gids_thisCat = Gids(idxCells_thisCat);
        [uGids, gidCount] = uniqueCount(Gids_thisCat);
        
        BccFrac = getFracBetweenSitePairs(gidCount);        
        nSites_thisCat = length( uGids );
        
        cell_nSitesInAnimal(idxUnits_thisCat) = nSites_thisCat;
        cell_BccFracInAnimal(idxUnits_thisCat) = BccFrac;
        
    end

    
    for cell_i = 1:length(allCells)
        allCells(cell_i).locData.nSitesAtPen = cell_nSitesAtPen(cell_i);
        allCells(cell_i).locData.BccFracAtPen = cell_BccFracAtPen(cell_i);

        allCells(cell_i).locData.nSitesInAnimal = cell_nSitesInAnimal(cell_i);
        allCells(cell_i).locData.BccFracInAnimal = cell_BccFracInAnimal(cell_i);
    end
    
    %%

end



function meanJackCC = getRF_jackCC(RF_jacks)
     [nx, ny, nJacks] = size(RF_jacks);
     jack_RF_vec = reshape(RF_jacks, [nx*ny, nJacks]);

     [~, jack_ccs] = pearsonRm(jack_RF_vec);
     meanJackCC = mean(jack_ccs);
end


function tf = allTuningStatsHaveOrig(tuningStats)

    fn = fieldnames(tuningStats);
    tf = true;
    for j = 1:length(fn)
        if ~isfield( tuningStats.(fn{j}), 'orig')
            tf = false;
        end
    end
end

%     function x = wo0(x)
%         
%     end


%     % before running this mfile, first call
%     load('indivCells_movie_fg.mat');
%     load('indivCells_grating.mat');  %% ( if using grating cells )
%     save results_GratingCells.mat celldata*
%     S = load(results_GratingCells.mat)
% 
%     % then call this function with:
%     generateGratingDatafiles(S)




%{
S_fg = load('flashedGratingCelldata.mat');
fgCells = S_fg.allCells;
fgstats = [fgCells.stats];
fg_ori_sel_p = [fgstats.orientationSelectivePval];
fg_ori_rep_p = [fgstats.orientationReproduciblePval];
fg_spf_rep_p = [fgstats.spatFreqReproduciblePval];

S_dg = load('driftingGratingCelldata.mat');
dgCells = S_dg.allCells;
dgstats = [dgCells.stats];
dg_ori_sel_p = [dgstats.orientationSelectivePval];
dg_ori_rep_p = [dgstats.orientationReproduciblePval];
dg_spf_rep_p = [dgstats.spatFreqReproduciblePval];


S_fg = load('flashedGratingCelldata.mat');
[allCells, Wcc_origIdxs, Wcm_origIdxs, Bcc_origIdxs, Bcm_origIdxs] = elements(S_fg);

f1m2 = [Wcc_origIdxs.frac1AtMax2];
f2m1 = [Wcc_origIdxs.frac2AtMax1];
fm12 = [Wcc_origIdxs.fracMax12];


%}


%     if have 2 cells - 1 pair - but 1 independent value, should contribute 1 - each pair contributes 1 
%     if have 3 cells - 3 pairs - but 2 independent values, should contribute 2 total - each pair contributes 2/3 
%     if have 4 cells - 6 pairs - but 3 independent values, should contribute 3 total - each pair contributes 3/6 = 2/4 
%     if have 5 cells - 10 pairs - but 4 independent values, should contribute 4 total - each pair contributes 4/10 

%                 n = length(setdiff(cellIds, 0)); % number of non-MU cells in this group.
%                 if any(n == [0, 1])
%                     contrib = 0;
%                 else
%                     contrib = 2/n; %= (n-1)/(.5*(n^2-n)); % # independent variables / number of within-site cell-cell pairs.
%                 end


%         % Pick out which cells between sites can be compared.
%         tic;
%         cmpOriPh = @(x, y)  length(x) == length(y);
%         cmpSpf   = @(x, y)  (length(x) == length(y));
%         okPairs = false( nGroups );
%         for i = 1:nGroups
%            for j = i+1:nGroups
%                 okPairs(i,j) = cmpOriPh(oris{i}, oris{j}) && cmpSpf(sps_deg{i}, sps_deg{j}) && cmpOriPh(phs{i}, phs{j});
%                 % (bottom-left half of matrix will be zero.)
%            end
%         end        
% 
%         [Wcc_origIdxs, Wcm_origIdxs, Bcc_origIdxs, Bcm_origIdxs] ...
%            = getCellGroupPairs(cellGroupIds, cellIds, okPairs);  
%        toc; 






% Select flashed/drifting groups
%         if strcmp(gratingType, 'flashed')
%             idx_selected = findInStructArray(gratingGroups, 'stimType', [], @(s)   strncmp(s, 'Movie:Flash', 11) || strncmp(s, 'Grating:Flash', 13)  );
%         elseif strcmp(gratingType, 'drifting')
%             idx_selected = findInStructArray(gratingGroups, 'stimType', [], @(s) ~(strncmp(s, 'Movie:Flash', 11) || strncmp(s, 'Grating:Flash', 13)) );            
%         end             
%         gratingGroups = gratingGroups(idx_selected);
%         
%         % Always remove single grating experiments (just 1 ori / sp)
%         idx_notSingle = findInStructArray(gratingGroups, 'stimType', [], @(s) ~strncmp(s, 'Grating:Single Grating', 14) );
%         gratingGroups = gratingGroups( idx_notSingle ); 
%         
%         if excludeCellsWithMultipleTempPeriods
%             idx_singleTempPeriod = findInStructArray(gratingGroups, 'tempPeriod_sec', [], @(x) length(x) == 1);
%             gratingGroups = gratingGroups( idx_singleTempPeriod ); % remove multiple temporal frequency experiments
%         end
%                 
%         % Combine with movie groups (flashed grating) if are considering movie groups as well.
%         if ~isempty(groupData)
%             groupData = combineStructArrays(groupData, gratingGroups);
%         else
%             groupData = gratingGroups;
%         end
%         if exist('S_indivCells', 'var') %|| ~isempty(S_indivCells)            
%             S_indivCells = mergeStructs(S_indivCells, S_grating_data);
%         else
%             S_indivCells = S_grating_data;
%         end




%     useRanksToExcludeBadCells = false;
%     rankThreshold = 3;
%     addRankFields = true;


%     if useRanksToExcludeBadCells || addRankFields
%         S_rnk = load(['cellRanks_' gratingType '.mat']);
%         [rankGids, rankCellIds, rankVals] = deal([S_rnk.ranks.Gid], [S_rnk.ranks.cellId], [S_rnk.ranks.rank]);        
%     end
%     rankFields = iff(addRankFields, {'manualRank', []}, {});

%             if useRanksToExcludeBadCells
%                 ind = find( (rankGids == Gid) & ( rankCellIds == cellId));
%                 if isempty(ind) || rankVals(ind) < rankThreshold
%                     continue; %ie. skip the cell.
%                 end
%             end

%                 if addRankFields
%                     ind = find(Gid == rankGids & cellId == rankCellIds);
%                     if ~isempty(ind)
%                         rnk = rankVals(ind);
%                     else
%                         rnk = -1;
%                     end
%                     rankFields = {'manualRank', rnk};
%                 end


%         switch selectGroupThatMaximizes
%             case 'numCells', score = nCellsTot;
%             case 'numSpikes', score = nSpikesTot ;
%             case 'numCellsXnumSpikes', score = nSpikesTot .* nCellsTot;
%             case 'mostReproducibleCells', 


%             if any(strcmp(featuresToInclude, 'meanMD'))
%                 S_prop = load(getFileName('properties', Gid));                
%             end
%                 featuresToInclude = {'meanNegAmps', 'meanWvfm', 'meanMD', 'spikeAmp'};
                
                
%             switch featuresToInclude
%                 case 'negAmps',     features = groupNegAmps;
%                 case 'waveformPCA', features = getGroupWaveformCoefficients('PCA', Gid, 2, 'separate', 'unnorm', 'ccw', 'match');
%             end
% %             nFeatures = size(features, 2);
%             spkCellIds = groupSpikes(:,2);
%             
%             if includeMeanSpikeWvfm                
%                 allWvfms = getSpikeWaveforms(Gid, [], 'unnorm', 'raw', 'nomatch');
% %                 [channelMeans, channelCov] = getChannelMeansAndCovariance(Gid);            
%             end



%             if ~iscell(R)
%               ....
%             else
%                 allCells(cell_idx) = struct('Gid', Gid, 'cellId', cellId, 'stimType', [], 'ori', 0, 'sp', 0, 'ph', 0, 'ph_max', 0, 'tp_sec', 0, 'R', 0, 'nspikes', [], ...
%                     'maxR', 0, 'oriSp_maxR_av', 0,  'oriSp_maxR_no_av', 0, ...
%                     'F1oDC_maxR_av',  0, 'F1oDC_maxR_no_av',  0, ...
%                     'ori_sp_maxMU', 0, ...
%                     'degPerPix', 0, 'stats', {stats}, 'type', gratingType, 'locData', groupLocData, 'PSTH', PSTH);
%             end



%                     idx_thisCell = (spkCellIds == cellId);
%                     cellFeatures = features(idx_thisCell,:);
%                     featureData = struct('mean', mean(cellFeatures,1), 'allFeatures', cellFeatures);
% %                     featureFields = {'elecrodeData', elecData};                    
%                     
%                     if includeMeanSpikeMN 
%                         cellNegAmps = groupNegAmps(idx_thisCell,:)';
%                         [M,C] = getChannelMeansAndCovariance(Gid);
%                         mndists = sqrt( sum( cellNegAmps .* (C \ cellNegAmps), 1) );
%                         meanMNdist = mean(mndists);
%                         featureData.meanMNdist = meanMNdist;
%                     end
%                     if includeMeanSpikeWvfm                    
%                         meanWvfm = mean( allWvfms(:,:, idx_thisCell), 3);
%                         featureData.meanWaveform = meanWvfm(:);
%                     end


%             if addMaxF1oDCfields
%                 R2 = reshape(R, [nOri*nSp, nPh]); 
%                 allF1oDCs = reshape(getF1oDC(R2), [nOri, nSp]);
%                 [maxF1oDC, oriSp_maxF1oDC] = maxElement(allF1oDCs);
%             end


%     switch cmpType
%         case 'degree', 
%             inclusionCriteron = 'reproducible&selective';
%             inclusionThreshold = 1;  % 1 --> 'true';
%         case 'phase',  % if strcmp(gratingType, 'flashed') || strcmp(cmpType, 'phase');
%             inclusionCriteron = 'reproducible';        
%             inclusionThreshold = 3;  % 3 --> pval < 1e-3            
% %         case 'clusters'
% %             inclusionCriteron = 'multiple clusters',
% %             inclusionThreshold = 3;
%             
%     end



%                 if plotF1sAndDcs
%                     figure(8); clf;
%                     [allF1s, allDCs] = getF1oDC(R2);
%                     for i = 1:2
%                         subplot(1,2,i);
%                         if i == 1
%                            r = mean(R,3);
%                         elseif i == 2
%                            r = max(R, [], 3);
%                         end
%                         plot(r(:), allF1oDCs(:), 'bo', 'markerfacecolor', 'b', 'markersize', 6); hold on;
%                         plot(r(:), allF1s(:),    'go', 'markerfacecolor', 'g', 'markersize', 5);
%                         plot(r(:), allDCs(:),    'ro');
%                         indMaxF1oDC = indmax(allF1oDCs(:)); plot(r(indMaxF1oDC), allF1oDCs(indMaxF1oDC), 'bs', 'markersize', 10);
%                         indMaxF1    = indmax(allF1s(:));    plot(r(indMaxF1 ),   allF1s(indMaxF1),    'gs', 'markersize', 10);
% 
%                         plot([min(allDCs(:)), max(allDCs(:))], [min(allDCs(:)), max(allDCs(:))], 'r:')
%                         plot(0,0);
%                         if i == 1
%                             xlabel('mean response at ori/sp');  
%                             legend('F1/DC', 'F1', 'DC', 'location', 'bestoutside');
%                         elseif i == 2
%                            xlabel('max response at ori/sp');  
%                            legend('F1/DC', 'F1', 'DC', 'location', 'bestoutside');
%                         end
% 
%                     end                       
%                     3;
%                 end

%     if showF1oDCdistribution
%         F1oDCs = [allCells.F1oDC_maxR];
%         hist(F1oDCs);
%         title(['Distribution of F1 / DC of phase tuning curve (at pref. ori/sf)   (N = ' num2str(length(F1oDCs)) ')']);
%         xlabel('F1 / DC');
%         xlim([0 2]);
%         return;
%     end


%{
                cellStruct = struct('Gid', Gid, cellIdFields{:}, 'stimType', stimType, ...
                    'ori', oris, 'sp', sps, 'ph', phs, 'ph_max', ph_max, 'tf_Hz', tf_Hz, 'R', R, 'nspikes', grpData.nSpikes(cell_i), 'contrib', contrib, ...
                    'maxR', maxR_no_av, 'oriSp_maxR_av', oriSp_maxR_av,  'oriSp_maxR_no_av', oriSp_maxR_no_av, ...
                    'F1oDC_maxR_av',  F1oDC_maxR_av, 'F1oDC_maxR_no_av',  F1oDC_maxR_no_av, ...
                    ... 'oriSp_maxF1oDC', oriSp_maxF1oDC, 'maxF1oDC', maxF1oDC, ...
                     'ori_sp_maxMU', ori_sp_maxMU, 'degPerPix', degPerPix, 'stats', stats, ...
                    'type', gratingType, 'locData', groupLocData, 'PSTH', PSTH, ...
                    phaseTC_CC_fields{:},  R_fullFields{:}, spkFeatureFields{:} );
%}


  %{      
      filename_base = strrep(filename_base, '_S.mat', '.mat'); % remove _S or _C from filename for all files
    filename_base = strrep(filename_base, '_C.mat', '.mat');
      if saveSimpleComplexCellOnlyFiles
            idx_simple  = cellF1oDCs_opt_sm >= 1;        
            S_usable_simpleCells.allCells = allCells(cells_usable & idx_simple & (allCellIds > 0));
            stat_str = getStatStr(S_usable_simpleCells, clustGrp_str);
            filename_simple = [filename_base '_S.mat'];
%             save([pathname filename_simple], '-struct', 'S_usable_simpleCells', '-v6');
            fprintf('\nSaved data for USABLE %s grating SIMPLE %s (%s)\n       to file %s\n', ...
                gratingType, clustGrp_str, stat_str, filename_simple);

            idx_complex = cellF1oDCs_opt_sm <  1;                        
            S_usable_complexCells.allCells = allCells(cells_usable & idx_complex & (allCellIds > 0));
            stat_str = getStatStr(S_usable_complexCells, clustGrp_str);
            filename_complex = [filename_base '_C.mat'];
%             save([pathname filename_complex], '-struct', 'S_usable_complexCells', '-v6');
            fprintf('\nSaved data for USABLE %s grating COMPLEX %s (%s)\n       to file %s\n', ...
                gratingType, clustGrp_str, stat_str, filename_complex);            
        end
%}


%{


            if ~strcmp(timeWindow, 'best') && isfield(celldata.OSP, 'fixedWindowOSPs')
                %%
                windows_ms_avail = floor(celldata.OSP.fixedWindowOSPs.windows_ms);
                wind_idx = findRows(timeWindow, windows_ms_avail);                
                if isempty(wind_idx)                
                    error('Window not found');
                end
                R_oe = celldata.OSP.fixedWindowOSPs.OSPs{wind_idx};
                R = mean(R_oe, 4);
                
            else
                R_full = OSP.R_full;
                R_oe = [];
                R = OSP.R;
            end
%}

%{

            if ~isempty(tuningStats)
                cellF1oDC_opt_sm_all = structfun(@(s) s.F1oDC, tuningStats);
%                 assert(length(unique(cellF1oDC_opt_sm_all))==1);
                             
            else
                R_os = mean(R, 3);
                [~, i_max] = maxElement(R_os);
                ptc = squeeze( R(i_max(1), i_max(2), :));
                f1odc = getF1oDC(ptc);
                cellF1oDC_opt_sm = f1odc;
            end
            %}

%{
%             if strcmp(gratingType, 'flashed');
%                 [L_bin, R_bin] = dealV(PSTH.timeWindow_bins);
%                 windowProfile = PSTH.windowProfile;
%             else
%                 [L_bin, R_bin] = deal(1);
%                 windowProfile = [];
%             end
%}

%{

            
%                 if isSpatialFreqBatch || isFlashedBatch
%                     tuningStats.spfStats_si.F1oDC = F1oDC_maxR_avP;                    
%                 end

%                 [maxR_no_av, oriSp_maxR_no_av] = maxElement(R);
%                 [ori_i_no_av, sp_i_no_av] = elements(oriSp_maxR_no_av);
%                 phaseTC_atMax_no_av = squeeze(R(ori_i_no_av, sp_i_no_av, :));
%                 F1oDC_maxR_no_av = getF1oDC(phs, phaseTC_atMax_no_av, 360 ); 

%                 meanPhs(cell_idx) = mean(phaseTC_atMax_av);
%                 meanPhsAll(cell_idx) = mean(R(:));

%                 ph_max_ind = indmax(phaseTC_atMax_av);
%                 ph_max = round(phs(ph_max_ind));



            %}

%{
%             meanR = mean(R, 3);                                                
%             [maxR_av, oriSp_maxR_avP] = maxElement(meanR);
%             [ori_i_av, sp_i_av]      = elements(oriSp_maxR_avP);                
%             phaseTC_atMax_avP = squeeze(R(ori_i_av, sp_i_av, :));                    
%             F1oDC_maxR_avP = getF1oDC(phs, phaseTC_atMax_avP, 360 ); % phase_max is always 360 for phases.                
            
%}


%{


%                      s.STA_odd_jack{w_idx} = getJackknifedSTA(r_full, Gid, 'odd', nJackSegments_STA, j_method);
%                      s.STA_even_jack{w_idx} = getJackknifedSTA(r_full, Gid, 'odd', nJackSegments_STA, j_method);
%                      end
                     
                     
%                     if strcmp(timeWindow, 'best')
% %                         if timeWindow_str
% %                             [L_bin, R_bin] = dealV(PSTH.timeWindow_bins);
%                         STA_oe = celldata.STAs.STA;
%                         haveSTA = 1;
%                     elseif isfield(celldata.STAs, 'fixedWindowSTAs')
% 
%                         windows_ms_avail = floor(celldata.STAs.fixedWindowSTAs.windows_ms);
%                         wind_idx = findRows(timeWindow, windows_ms_avail);
%                         haveSTA = ~isempty(wind_idx);
%                         if haveSTA
%                             STA_oe = celldata.STAs.fixedWindowSTAs.STAs{wind_idx};
%                         end                            
%                     end

%}

%{

                    [phaseTC_CCs_oe, phaseTC_CC_ps_oe, phaseTC_CCs_hoe, phaseTC_CC_ps_hoe, phaseTC_CCs_fs, phaseTC_CC_ps_fs] = ...
                        getOspDataForPsthWindow(Gid, cellId, [], [], L_bin, R_bin, windowProfile, ...
                        {'phaseTC_CCs_oe', 'phaseTC_CC_ps_oe', 'phaseTC_CCs_hoe', 'phaseTC_CC_ps_hoe', 'phaseTC_CCs_fs', 'phaseTC_CC_ps_fs'});

%}

%{

%                             timeWindow_str_dst = iff(strcmp(RF_timeWindow_str, 'best'), '', RF_timeWindow_str);
%                             timeWindow_str_src = iff(strcmp(RF_timeWindow_str, 'best'), '', ['_' RF_timeWindow_str]);
%                             timeWindow_str_src = iff(strcmp(RF_timeWindow_str {wi}, 'best'), '', sprintf('__%d_%d', RF_timeWindows{wi}));
%                             timeWindow_str_dst = iff(strcmp(RF_timeWindows{wi}, 'best'), '_best', sprintf('_%d_%d', RF_timeWindows{wi}));

%                         RF_timeWindows_str{wi}
%                             MID_window_fldname = ['MID' timeWindow_str_dst];

%}