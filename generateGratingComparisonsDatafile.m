function generateGratingComparisonsDatafile(cmpOriSpfType)

    global suspectGCs suspectGids suspectCellIds suspectProfileCC
    
    global maxMinFracR_cmp
    maxMinFracR_cmp = [];
    
%     global smoothPhases;
    % checklist before running this algorithm:
    %   - complete all calculations for all flashed/drifting-grating cells using 'doCalculationsForAllCells' 
    %   - compile the flashed grating pair data using 'generateGratingPairDatafiles' 
    %    
    
    %%% odd /even trial options
    odd_even_trials_useAllTrialsCrit = true;
%     opt.getPhaseCmpJackStdErr = 1;
    opt.getPhaseCmpJackStdErr_ptc = 1;
    opt.getPhaseCmpJackStdErr_rf = 1;

    doRandMIDs  = 1;
    doRandMID_fits = 1; % bottleneck when we do it, so only do it when necessary
    opt.doRandMID_fits_jacks = 0; % this takes even longer (~20-40 minutes).
    
    flashed_doRF_measures = true;
    [opt.phase_oeState, opt.phase_oe_mode, opt.phase_oe_action] = curPhaseOEmode;
%     opt.phase_oeState = 'oe';    
%     opt.phase_oe_mode = 'diff';
%     opt.phase_oe_action = 'keepBoth'; % 'keepBoth' from odd/even, or 'average' 


    [opt.degree_oeState, opt.degree_oe_mode, opt.degree_oe_action] = curDegreeOEmode;
%     opt.degree_oeState = 'aa';
%     opt.degree_oe_mode = 'diff';
%     opt.degree_oe_action = 'keepBoth'; 

    useExternalMinReqOE = 1;
    
    if useExternalMinReqOE        
%         curMinR_oe(0.25, 'hoe') 
        
        [opt.PTC_minReq_r_oe, opt.PTC_oe_corrMode] = curMinR_oe;
        opt.PTC_require_oe_corr = ~isnan(opt.PTC_minReq_r_oe);
        
    else
        curMinR_oe(nan);
        opt.PTC_require_oe_corr = 1; %~isnan(curMinR_oe);
        opt.PTC_oe_corrMode = 'hoe'; %~isnan(curMinR_oe);
        opt.PTC_minReq_r_oe = 0.5; %curMinR_oe;
    %     opt.PTC_minReq_r_oe = 0.25; %curMinR_oe;
    end
    
    opt.PTC_require_oe_corr_p = 0; %~isnan(curMaxP_oe);
%     opt.PTC_oe_p_corrMode = 'hoe';
%     opt.PTC_maxReq_p_oe = 1; %curMaxP_oe;
    

    opt.makeDphiNullWhenAmbiguous = 1;
    opt.usePosNegDphi = 0;

    %%% MID criteria options
    opt.MID_minReq_jackCC = -1; %0.2;
    opt.MID_minReq_rsqr_oe = .25;
    opt.MID_minReq_rsqr_fit = .25; %0.2
    opt.MID_ovlp_minFrac = .01;


    opt.applyStdErrorThresholds = true;
    opt.maxOriStdErr_deg = 5;
    opt.maxDSIStdErr = 0.1;
    opt.maxSpfStdErr = 0.5;
    opt.maxF1oDCErr = 0.25;

    
    %%%% misc options
    opt.F1oDC_field = 'F1oDC_maxR_avP_sm';
    


    opt.getDegreeDataForPhaseCmp = 1;
    firstDegreeStat = true;

%     useGainCorrectedResponses = strcmp(curResponseType(''), 'gainCorrected');
    
    
    GC = [];
    GC1 = [];
    
    allBCCs = [];
    allBCC_idx = 1;

    allWCCs = [];
    allWCC_idx = 1;
    recordAllCCs = 0;



        opt.doRandMID_fitsWithoutEnvelope = 0;
    
    doR90_measures = 1;
    
    %%% TASKS
%     compare_dOriSp = false;    
        calcPvalues = false;
        shufflePhases = false;
        shuffleStimuli = false;
%         shiftStimuliOri = false;
        shiftStimuliOri = curShiftOri;
        shuffleCells = false;
        flipMIDs = 0;

    randomlySwapPairs = false;
    includePhaseTC = 0;
    includeOriSp_cmp = 0;
    checkForImaginaryVals = false;

    subtractSpont = curSubtractSpont;
    responseType = curResponseType('');
    
    %%% 1. Which type of gratings to do
    gratingType = curGratingType('');    
    cmpType = curCmpType('');    
    if strcmp(cmpType, 'psth') && strcmp(gratingType, 'drifting')
        error('Can''t do psth comparisons with drifting grating data');
    end

    filename_ext = '';
    if strcmp(cmpType, 'degree') 
        if nargin == 0
            generateGratingComparisonsDatafile('ori');               
            if ~curPreserveAligned
                generateGratingComparisonsDatafile('spf');  % for aligned/anti-aligned, only direction (ori) is relevant
            end
            return;
        end
    else 
        cmpOriSpfType = '';
    end        
    
    cmpOriSpfType_str     = switchh(cmpOriSpfType, {'ori', 'spf', ''}, {'(ORIENTATION ONLY)', '(SPATIAL FREQ. ONLY)', ''});
    cmpOriSpfType_fileext = switchh(cmpOriSpfType, {'ori', 'spf', ''}, {'_ori', '_spf', ''});
    filename_ext = [filename_ext, cmpOriSpfType_fileext];
    
    
    weightingBy = 'product';
    selectingBy = 'minNormRate';
%     ptc_criteria = 'all_phases';   %'Only 1 phase', 'all_phases' 'Any phases > 0',
%     ptc_criteria = 'Any phases > 0';   %'Only 1 phase', 'all_phases' 'Any phases > 0',
%     ptc_criteria = 'Only 1 phase';   %'Only 1 phase', 'all_phases' 'Any phases > 0',
    ptc_criteria = '';
        
    cmpType_str = switchh(cmpType, {'phase', 'degree', 'psth'}, {'Phase Tuning', 'Degree of Tuning', 'Psth Similarity'});
    subtractSpont_str = iff(subtractSpont, '(SPONT Subtracted)', '');
    shuff_txt = iff(shufflePhases || shuffleCells || shiftStimuliOri, '(SHUFFLED)', '');  
    bccType_txt = curBccType;
    with_str = switchh(cmpType, {'phase', 'degree'}, {'', 'with \n'});        
    
    [phaseSmoothW, phaseSmoothMethod] = curPhaseSmoothing(''); phaseSmooth_str = curPhaseSmoothing('');
    %%
    smooth_str = iff(strcmp(gratingType, 'drifting') && strcmp(cmpType, 'phase') && ~isnan(phaseSmoothW) && ~isempty(phaseSmooth_str), ...
        sprintf('  --> Using SMOOTHED phase tuning curves (%s; %.1f)\\n', phaseSmoothMethod, phaseSmoothW), '');
    fprintf('\n *** Creating %s Comparison data file for %s grating stimuli %s %s %s %s %s ... ***\n', cmpType_str, gratingType, with_str, cmpOriSpfType_str, subtractSpont_str, shuff_txt, bccType_txt);
    fprintf(smooth_str);
    
    %%
    cells_fileName = getFileName('osps', filename_ext);
    fprintf('Loading cells from file : %s\n', basename( cells_fileName ));
    tic; S1 = load(cells_fileName); toc;
    allCells = S1.allCells;        
    clear('S1');
    nUnits = length(allCells);
    for ci = 1:nUnits
        allCells(ci).R = double(allCells(ci).R); % double precision.
    end
    
    if shuffleCells
        allCells = allCells(randperm(nUnits));
    end
    if shuffleStimuli
        rand('state', 0);
        tic; fprintf('Shuffling stimuli...');
        for ci = 1:nUnits
            R = double(allCells(ci).R);
            R_perm = nan(size(R));
            [nOris, nSpfs, nPhs] = size(R);            
            ri = randi(nOris*nSpfs);
%             idx = reshape( [ri, setdiff(1:nOri*nSpf, ri)], nOri, nSpf);
            idx = reshape( [ri: nOris*nSpfs, 1:ri-1], nOris, nSpfs);
            for ori_k = 1:nOris
                for sp_k = 1:nSpfs
                    [ori_k1, sp_k1] = ind2sub([nOris, nSpfs], idx(ori_k, sp_k));
%                     assert(all([ori_k1 sp_k1] == [ori_k sp_k]))
                    R_perm(ori_k, sp_k,:) = R(ori_k1, sp_k1, :);
                end
            end            
            assert(sum(R(:)) == sum(R_perm(:)));
            
            allCells(ci).R = R_perm;
        end
        fprintf('done.'); toc;
    end

    if shiftStimuliOri
        
        rand('state', 0);
        tic; fprintf('Shuffling stimuli...');
        for ci = 1:nUnits
            [nOris, nSpfs, nPhs] = size(allCells(ci).R);
%             ori_shift_idx = randi(nOris);
%             spf_shift_idx = randi(nSpfs);
%             new_ori_idxs = [ori_shift_idx: nOris, 1:ori_shift_idx-1];
%             new_spf_idxs = [spf_shift_idx: nSpfs, 1:spf_shift_idx-1];
%             new_spf_idxs = 1:nSpfs;
            new_ori_idxs = randperm(nOris);
            new_spf_idxs = randperm(nSpfs);
            
            allCells(ci).R = allCells(ci).R(new_ori_idxs,new_spf_idxs,:);            
            
            allCells(ci).tuningStats.oriStats_ss.OS  = allCells(ci).tuningStats.oriStats_ss.OS(new_ori_idxs,new_spf_idxs,:);
            allCells(ci).tuningStats.oriStats_ss.r_k_dir  = allCells(ci).tuningStats.oriStats_ss.r_k_dir(new_ori_idxs);            
        end
        fprintf('done.'); toc;
    end
    
    if flipMIDs
       
        rand('state', 0);
        tic; fprintf('Flipping & Rotating STAs & MIDs...');
        for ci = 1:nUnits
            
            rot_i = randi(4)-1; %0, 90, 180, 270,
            flip_i = randi(2)-1; % flip - yes/no

            allCells(ci).STA = flip_rot(allCells(ci).STA, flip_i, rot_i);
            
            data = allCells(ci).MIDdata;
            if isempty(fieldnames(data))
                continue;
            end
            
            allCells(ci).MIDdata.MID = flip_rot(data.MID, flip_i, rot_i);
            allCells(ci).MIDdata.MID_fit = flip_rot(data.MID_fit, flip_i, rot_i);
            allCells(ci).MIDdata.MID_select = flip_rot(data.MID_select, flip_i, rot_i);
                                            
        end
        fprintf('done.'); toc;
        
        
        
    end
    
    if opt.doRandMID_fitsWithoutEnvelope
        
        tic; fprintf('Recalculating MID_fits without envelopes');
        idx_sigmas = [4, 5]; % A, mu_x, mu_y, sig_x, sig_y, 
        for ci = 1:nUnits
            
            Gid = allCells(ci).Gid;
            MID_types = {'MID', 'MID_odd', 'MID_even'};
            for mi = 1:length(MID_types)
            
                if isfield(allCells(ci).MIDdata, MID_types{mi});
                    p = allCells(ci).MIDdata.MID_types{mi}.p;
                    MID_fit_cur = getMIDfitToGabor(Gid, p);
                    assert(isequal(MID_fit_cur, allCells(ci).MIDdata.(MID_types{mi}).MID_fit));
                    p(idx_sigmas) = 1e10;
                    
                    MID_fit_new = getMIDfitToGabor(Gid, p);
                    allCells(ci).MIDdata.(MID_types{mi}).MID_fit = MID_fit_new;
                    
                end
                
            end
            
            
            if isempty(fieldnames(data))
                continue;
            end
            
            allCells(ci).MIDdata.MID = flip_rot(data.MID, flip_i, rot_i);
            allCells(ci).MIDdata.MID_fit = flip_rot(data.MID_fit, flip_i, rot_i);
            allCells(ci).MIDdata.MID_select = flip_rot(data.MID_select, flip_i, rot_i);
                                            
        end
        fprintf('done.'); toc;
        
    end
    
    if odd_even_trials_useAllTrialsCrit
        
       for ci = 1:nUnits
           ts = allCells(ci).tuningStats;
           
          
           if isfield(ts, 'oriStats_ss') && isfield(ts, 'oriStats_ss_even') && isfield(ts, 'oriStats_ss_odd')
               ori_fields = {'w_ori_global', 'w_ori_local', 'ori_pref_deg', 'DSI_global'};

               for field_i = 1:length(ori_fields)
                   
                   if ~isnan(ts.oriStats_ss.(ori_fields{field_i})) ... 
                           ... && ( isnan( ts.oriStats_ss_odd.(ori_fields{field_i}) ) || isnan( ts.oriStats_ss_even.(ori_fields{field_i}) ) )

                       ts.oriStats_ss_odd.( ori_fields{field_i}) = ts.oriStats_ss_odd.orig.( ori_fields{field_i});
                       ts.oriStats_ss_even.(ori_fields{field_i}) = ts.oriStats_ss_even.orig.(ori_fields{field_i});
                       
                   elseif isnan(ts.oriStats_ss.(ori_fields{field_i}))
                       ts.oriStats_ss_odd.(  ori_fields{field_i}) = nan;
                       ts.oriStats_ss_even.( ori_fields{field_i}) = nan;
                       
                   end
                   
                   
               end

           
           end
           
%            nPassed_ss = 0;
%            nPassed_ss_even = 0;
%            nPassed_ss_odd = 0;
           if isfield(ts, 'spfStats_ss') && isfield(ts, 'spfStats_ss_even') && isfield(ts, 'spfStats_ss_odd')
               spf_fields = {'w_spf', 'f_opt'};
               for field_i = 1:length(spf_fields)
                   
                    nPassed_ss.(spf_fields{field_i})(ci) = ~isnan(ts.spfStats_ss.(spf_fields{field_i}));
                    nPassed_ss_odd.(spf_fields{field_i})(ci) = ~isnan(ts.spfStats_ss_odd.(  spf_fields{field_i}));
                    nPassed_ss_even.(spf_fields{field_i})(ci) = ~isnan(ts.spfStats_ss_even.(  spf_fields{field_i}));
                                      
                   
                   if ~isnan(ts.spfStats_ss.(spf_fields{field_i})) ... 
                           ... && ( isnan( ts.spfStats_ss_odd.(spf_fields{field_i}) ) || isnan( ts.spfStats_ss_even.(spf_fields{field_i}) ) )

                       ts.spfStats_ss_odd.( spf_fields{field_i}) = ts.spfStats_ss_odd.orig.( spf_fields{field_i});
                       ts.spfStats_ss_even.(spf_fields{field_i}) = ts.spfStats_ss_even.orig.(spf_fields{field_i});
                       
                   elseif isnan(ts.spfStats_ss.(spf_fields{field_i}))
                       ts.spfStats_ss_odd.(  spf_fields{field_i}) = nan;
                       ts.spfStats_ss_even.( spf_fields{field_i}) = nan;
                   end
                   
                   
               end
           end
           
           allCells(ci).tuningStats = ts;
       end
        
        
        
    end
    
    tic;
%     if strcmp(cmpType, 'phase') && useGainCorrectedResponses
%         opt.getPhaseCmpJackStdErr = false;
%         for ci = 1:nUnits
%             useThisCell = allCells(ci).R_corr.p_Rcorr_av_gain > .01;
%             if useThisCell
%                 allCells(ci).R = allCells(ci).R_corr.R;
%                 allCells(ci).R_oe = allCells(ci).R_corr_oe.R;
%                 allCells(ci).R_hoe = allCells(ci).R_corr_hoe.R;
%             else
%                 allCells(ci).R = [];
%                 allCells(ci).R_oe = [];
%                 allCells(ci).R_hoe = [];
%             end
%             
%         end
        
%     end
    
    
    if strcmp(cmpType, 'phase') && opt.getPhaseCmpJackStdErr_ptc  % calculate all jackknife-average phase tuning curves now
        3;
%         fprintf('Calculating jackknife of phase tuning curves ... \n')
        fprintf('Decompressing phase tuning curve jackknives ... \n'); tic;
        cell_jack_field = ['R_jackTrials_' opt.phase_oeState];
        for ci = 1:nUnits
%             if isfield(allCells(ci), cell_jack_field) && ~isempty(allCells(ci).(cell_jack_field));
%                 allCells(ci).R_stim_jackTrials = allCells(ci).(cell_jack_field);
%             else
%                 R_full_jackknifeTrials = getPhaseTuningJackknifedTrials(allCells(ci).R_full, opt.phase_oeState);
%                 allCells(ci).R_stim_jackTrials = R_full_jackknifeTrials;
%             end
            allCells(ci).R_stim_jackTrials = decompress( allCells(ci).(cell_jack_field) );
        end            
        toc;
    end

    F1oDC_types = allCells(1).F1oDC.types;
%     F1oDCs = [allCells.F1oDC_maxR_avP_sm];
%     idxCell = [allCells.cellId] > 0;
%     idx_simpCell = F1oDCs > 1 & idxCell;
%     idx_use1 = idxCell;
%     idx_use1 = idx_simpCell;
%     nP_sf_pref = nPassed_ss.f_opt(idx_use1);
%     nP_sf_w = nPassed_ss.f_opt(idx_use1);
% 
%     nP_sf_pref_even = nPassed_ss_even.f_opt(idx_use1);
%     nP_sf_w_even = nPassed_ss_even.f_opt(idx_use1);
% 
%     nP_sf_pref_odd = nPassed_ss_odd.f_opt(idx_use1);
%     nP_sf_w_odd = nPassed_ss_odd.f_opt(idx_use1);

    3;
    
%     if ~isempty(smoothPhases)
%         ospDatafile_noSmooth = [gratingType 'GratingCells.mat'];            
%         S1_noSmooth = load([pathname ospDatafile_noSmooth]);    
%         allCells_noSmooth = S1_noSmooth.allCells;        
%     end

    pairs_fileName = getFileName('pairs', filename_ext);
    fprintf('Loading pairs from file : %s\n', basename( pairs_fileName ));
    tic; S2 = load(pairs_fileName); toc;
    [Wcc_pairIdxs, Wcm_pairIdxs, Bcc_pairIdxs, Bcm_pairIdxs, Bmm_pairIdxs,  Wrcc_pairIdxs, Wrcm_pairIdxs ] = ...
        deal(S2.Wcc_idxs, S2.Wcm_idxs, S2.Bcc_idxs, S2.Bcm_idxs, S2.Bmm_idxs, S2.Wrcc_idxs, S2.Wrcm_idxs);
    clear('S2');
    %%% 1. Which pairtypes to do:
    
%     pairTypes = {'Wcc', 'Wrcc', 'Bcc', 'Wcm', 'Wrcm', 'Bcm'}; 
%     pairTypes = {'Bcc'}; 
%     pairTypes = {'Wcc', 'Wcm'}; 
%     pairTypes = {'Wcc', 'Wrcc', 'Wcm'}; 

            
    if strcmp(cmpType, 'degree')
        curPairTypes('Wcc', 'Wrcc', 'Bcc');
    elseif strcmp(cmpType, 'phase')
        curPairTypes('Wcc', 'Wscc');
%         curPairTypes('Wcc', 'Wcm', 'Bcc', 'Bcm', 'Bmm', 'Wscc')
%         curPairTypes('Wcc', 'Wrcc', 'Bcc');
    end
        

    [pt_ids, pairTypes] = curPairTypes;
    
%     pairTypes_str = [pairTypes{:}];
%     pairTypes = {'Wcc'};     
    if strcmp(curCmpType(''), 'degree')
        pairTypes(strcmp(pairTypes, 'Wscc')) = [];
    end

    doWcc   = any(strcmp(pairTypes, 'Wcc'));
    doWscc  = any(strcmp(pairTypes, 'Wscc'));
    doWcm   = any(strcmp(pairTypes, 'Wcm'));
    doBcc   = any(strcmp(pairTypes, 'Bcc'));
    doBcm   = any(strcmp(pairTypes, 'Bcm'));
    doBmm   = any(strcmp(pairTypes, 'Bmm'));
    
    doWrcc  = any(strcmp(pairTypes, 'Wrcc'));
    doWrcm = any(strcmp(pairTypes, 'Wrcm'));    

    Wrcc_pairIdxs = uniqueInts(Wrcc_pairIdxs);
    Wrcm_pairIdxs = uniqueInts(Wrcm_pairIdxs);

    allPairIdxs = [iff(doWcc || doWscc, Wcc_pairIdxs(:), []); 
                   iff(doWcm, Wcm_pairIdxs(:), []);
                   iff(doBcc, Bcc_pairIdxs(:), []);
                   iff(doBcm, Bcm_pairIdxs(:), []);
                   iff(doBmm, Bmm_pairIdxs(:), []);
                   iff(doWrcc, Wrcc_pairIdxs, []);
                   iff(doWrcm, Wrcm_pairIdxs, []) ];
    allPairIdxs = unique(allPairIdxs);
    nPairs = length(allPairIdxs);        
    pairsToSkip = zeros(nPairs,1);
    pairsToSkip_i = 1;
    
    
    idxMtx = zeros(nUnits, nUnits, 'uint32');
    idxMtx(sort(allPairIdxs)) = 1:length(allPairIdxs);    

    nStimMax = iff(strcmp(gratingType, 'flashed'), 360, 72);
    
    if recordAllCCs
        if doWcc
            allWCCs = zeros(length(Wcc_pairIdxs), nStimMax);
            allWCC_idx = 1;
        end
        if doBcc
            allBCCs = zeros(length(Wrcc_pairIdxs) + length(Bcc_pairIdxs), nStimMax);
            allBCC_idx = 1;
        end
    end
    
    
    switch cmpType
        case 'phase',
    %%% 2. Which locations to measure at:
    % ( all locations: 'maxR1',  'maxR2',  'mean12',   'maxR1xR2',  'maxMU',  'maxMinFracR', 'wgted sum', 'mtx x mtx'  )  
%     locations = {'maxMinFracR', 'wgtedSum_all', 'wgtedSum_top10', 'wgtedSum_top20', 'wgtedSum_above90', 'wgtedSum_above50'}; 
%     locations = {'wgtedSum_all', 'wgtedSum_1phase', 'wgtedSum_2phase', 'wgtedSum_4phase', 'wgtedSum_allphase'}; 
%     locations = {'maxMinFracR', 'p75MinFracR', 'p50MinFracR', 'p25MinFracR'}; 
%     locations = {'maxMinFracR', 'maxR1xR2', 'wgtedSum_all'}; 
%     locations = {'maxMinFracR', 'maxR1xR2', 'p50MinFracR'}; 
%     locations = {'maxMinFracR', 'p75MinFracR', 'p50MinFracR', 'p25MinFracR', 'wgtedSum_all', 'wgtedSum_all_oe'}; 
%         locations = {'wgtedSum_all'}; 
    
%     locations = {'maxMinFracR'}; 
    locations = {'maxMinFracR', 'maxMinFracR2',  'maxMinFracR3', 'maxR1xR2'}; 
%     locations = {'maxR1',  'maxR2', 'maxR1xR2', 'maxMinFracR', 'maxMU'}; 
%     locations = {'maxMinFracR', 'maxR1xR2', 'minMinFracR', 'p75MinFracR', 'p50MinFracR', 'p25MinFracR'}; % all locations: 'maxR1',  'maxR2',  'mean12',   'maxR1xR2',  'maxMU',  'maxMinFracR', 'wgted sum', 'mtx x mtx'    
    %     locations = {'maxR1xR2', 'maxMinFracR', 'maxMU'};  
    %     locations = {'maxR1',  'maxR2',  'mean12',   'maxR1xR2',  'maxMU', 'maxMinFracR'};
    
    
        case 'degree',
            locations = {'NA'};
    end
    
    nLocations = length(locations);    

    [doMaxR1, doMaxR2, doMean12, doMaxR1xR2, doMaxMU, doMaxMinFracR, doMinMinFracR] = ...
        dealV( cellfun(@(loc) any(strcmp(locations, loc)), {'maxR1', 'maxR2', 'mean12', 'maxR1xR2', 'maxMU', 'maxMinFracR', 'minMinFracR'}) );
    
    if recordAllCCs 
        loc_idx_first_wgtSum = find(strncmp(locations, 'wgtedSum', 6), 1);
        if isempty(loc_idx_first_wgtSum)
            error('To record all ccs, at least one location must be a weighted sum location');
        end
        location_record = locations{loc_idx_first_wgtSum};              
    end

    if strcmp(cmpType, 'phase')
        oe_keepBoth = strcmp(opt.phase_oe_action, 'keepBoth');
        opt.oe_keep1 = strcmp(opt.phase_oe_action, 'keep1') || oe_keepBoth;
        opt.oe_keep2 = strcmp(opt.phase_oe_action, 'keep2') || oe_keepBoth;
        opt.oe_average = strcmp(opt.phase_oe_action, 'average');
    elseif strcmp(cmpType, 'degree')
        oe_keepBoth = strcmp(opt.degree_oeState, 'oe') && strcmp(opt.degree_oe_action, 'keepBoth');
    end
                    
    opt.n3 = iff(oe_keepBoth, 2, 1);
    opt.oe_keepBoth = oe_keepBoth;
% 
    %%% 3. Which statisitical measures to perform:
%     measures = {'cc'};                        % all measures = {'cc', 'rho', 'tau', 'dphi', 'dF1'};        
%     measures = {'cc', 'rho', 'dphi'};                        % all measures = {'cc', 'rho', 'tau', 'dphi', 'dF1'};        
%     measures = {'cc', 'dphi' };                        % all measures = {'cc', 'rho', 'tau', 'dphi', 'dF1'};        

%     measures = {'dphi'};        

    cellClosenessMeasures = {'negAmps_dist', 'negAmps_overlap', 'fullWaveform_cc', 'fullWaveform_ed', 'channelWvfm_meanCC', ...
        'PCA_overlap', 'GLF_overlap', 'minID', 'maxL_ratio', 'diffFWHM', 'diffPtPWidth'};    
%     cellClosenessMeasures = {'negAmps_dist', 'negAmps_cc', 'fullWaveform_cc', 'channelWvfm_meanCC', 'PCA_overlap', 'GLF_overlap'};    
    
    
    % bins:
    % 1. for cc, rho
    xBins0_1   = linspace( 0, 1, 18);
    xBins1_1   = linspace(-1, 1, 18);
    xBins1_1_36   = linspace(-1, 1, 36);
    [xBins_dphi, xBins_dF1, Dw_ori_bins, Dw_ori_bins, D_dir_bins, Dw_spf_bins, D_spf_bins] = deal([]);
    
    % 2. for dphi & dF1, 
    switch cmpType,
        case 'phase'            
            
            PhaseTCMeasures = {'cc', 'rho', 'dphi', 'dF1'};
            RF_Measures = {'STA_cc', 'MID_cc', 'MID_ovlp_cc', 'MID_fit_cc', 'dPh_rel'};
            allMeasures = PhaseTCMeasures;
            if flashed_doRF_measures && strcmp(gratingType, 'flashed')
                allMeasures = [allMeasures, RF_Measures];
            end
            
            cc_bins = xBins1_1;
            rho_bins = xBins1_1;
            D_rel_ph_bins = [0:5:180];
            
            if opt.usePosNegDphi
                dphiBinE = @(nBin) binCent2edge( linspace(-180, 180, nBin) );
            else
                dphiBinE = @(nBin) binCent2edge( linspace(0, 180, nBin) );
            end
            
            switch gratingType
                case 'flashed', 
                    n_dPhiBins = switchh(opt.phase_oeState, {'aa', 'oe'}, [5, 9]);
                    n_dF1Bins = n_dPhiBins;
                    xBins_dphi = dphiBinE(n_dPhiBins);
                    xBins_dF1  = dphiBinE(n_dF1Bins);                                
                case 'drifting', 
                    n_dPhiBins = 11; xBins_dphi = dphiBinE(n_dPhiBins);
                    n_dF1Bins  = 11;  xBins_dF1  = dphiBinE(n_dF1Bins); % = linspace(-1e-5, 180+1e-5, n_dF1Bins+1);                

            end
            doOriSpfCmps = false;      
            PhaseTC_bins = {cc_bins, rho_bins, xBins_dphi, xBins_dF1};
            RF_bins = {xBins1_1_36, xBins1_1_36, xBins1_1_36, xBins1_1_36, D_rel_ph_bins};
            allMeasureBins = [PhaseTC_bins];
            if flashed_doRF_measures && strcmp(gratingType, 'flashed')
                allMeasureBins = [allMeasureBins, RF_bins];
            end
%             measures = allMeasures;
%             measures = {'cc', 'dphi'};
            measures = {'cc', 'dphi'};
            RF_ms = {'MID_cc', 'MID_fit_cc', 'MID_ovlp_cc', 'STA_cc', 'dPh_rel'};
            if flashed_doRF_measures && strcmp(gratingType, 'flashed')
                measures = [measures, RF_ms];
            end
%             measures = {'cc', 'rho', 'dphi'};
%             measures = {'STA_cc', 'MID_cc', 'MID_fit_cc', 'MID_ovlp_cc'};
%             measures = {'MID_cc', 'MID_fit_cc'};            
            
            if strcmp(gratingType, 'drifting') || (curMatchDB == 1)
                measures = setdiff(measures, RF_Measures);
                if isempty(measures)
                    error('Can''t do RF measures with drifting gratings / with DB calculations');
                end
            end                
            doingRFmeasures = any(strCcmp(measures, RF_Measures));
            doOSmeasures = false;
        case 'degree',
%             measures = {'Dw_ori_glob', 'Dw_ori_loc', 'D_dsi', 'D_ori_pref', 'D_dir_pref', 'Dw_spf', 'D_spf_pref'};

            oriMeasures = {'Dw_ori_glob', 'Dw_ori_loc', 'D_ori_pref', 'D_F1oDC_ori', 'D_F1pair_ori'};
            D_ori_bins = [0:5:90];
            D_F1oDC_bins = [0:.05:2];
            D_F1pair_bins = [0:.25:2];
            oriMeasure_bins = {D_ori_bins, D_ori_bins, D_ori_bins,   D_F1oDC_bins, D_F1pair_bins};
            
            dirMeasures = {'D_dsi_glob', 'D_dsi_glob_unrec', 'D_dsi_loc', 'D_dsi_loc_unrec', 'D_dir_pref', 'D_aligned_pair'};
            dsi_bins = [0:.05:1];
            dsi_unrec_bins = [0:.05:10];
            D_dir_bins = [0:5:180];
            D_aligned_bins = [1:1:222];
            dirMeasure_bins = {dsi_bins, dsi_unrec_bins, dsi_bins, dsi_unrec_bins, D_dir_bins, D_aligned_bins};
            
            spfMeasures = {'Dw_spf'    'D_w_SLN', 'D_Bn', 'D_skew',  'D_spf_pref',  'D_F1oDC_spf', 'D_F1pair_spf'};
            Dw_spf_bins = [0:.5:15]; % max: ~14.57
            D_spf_bins = [0:.25:6.25]; %max: 5.98
            spfMeasure_bins = {Dw_spf_bins, Dw_spf_bins, Dw_spf_bins, Dw_spf_bins, D_spf_bins, D_F1oDC_bins, D_F1pair_bins};

            OriSpMeasures_ori = {'cc_OSP', 'cc_ori'};
            OriSpMeasures_spf = {'cc_spf'};
            cc_bins = [-1:.05:1];
            OriSpMeasure_bins_ori = {cc_bins, cc_bins};
            OriSpMeasure_bins_spf = {cc_bins};
            
            R90_measures = {'dR_spont_abs', 'dR90_total_abs', 'dR90_stim_abs', 'dR_spont_rel', 'dR90_total_rel', 'dR90_stim_rel'};                        
            DR_abs_bins = [0:.1:80]; % max: ~75
            DR_rel_bins = [0:.05:3.5]; % max: ~3.12
            R90_measure_bins = {DR_abs_bins, DR_abs_bins, DR_abs_bins, DR_rel_bins, DR_rel_bins, DR_rel_bins};
            
            doDirMeasures = strcmp(gratingType, 'drifting');
            if ~doDirMeasures
                dirMeasures = {};
                dirMeasure_bins = {};
            end
            
            doOSmeasures = true; %strcmp(gratingType, 'flashed');
            if ~strcmp(cmpOriSpfType, 'ori') || ~doOSmeasures
                OriSpMeasures_ori = {};
                OriSpMeasure_bins_ori = {};
            end
            if ~strcmp(gratingType, 'flashed');
                idx_rm = strcmp(OriSpMeasures_ori, 'cc_OSP');
                OriSpMeasures_ori(idx_rm) = [];
                OriSpMeasure_bins_ori(idx_rm) = [];
            end
            if ~strcmp(cmpOriSpfType, 'spf') || ~doOSmeasures
                OriSpMeasures_spf = {};
                OriSpMeasure_bins_spf = {};
            end
%             OriSpMeasures      = [OriSpMeasures_ori, OriSpMeasures_spf];
%             OriSpMeasure_bins = [OriSpMeasure_bins_ori, OriSpMeasure_bins_spf];
                        
            
            if subtractSpont || curPreserveSimpleComplex || curPreserveAligned || ~strcmp(curBccType, 'full')
                R90_measures = {};
                R90_measure_bins = {};
            end
            
            
%             allMeasures    = [oriMeasures,     dirMeasures,     spfMeasures,     R90_measures];
%             allMeasureBins = [oriMeasure_bins, dirMeasure_bins, spfMeasure_bins, R90_measure_bins];
            switch cmpOriSpfType
                case 'ori'
                    allMeasures    = [oriMeasures,     dirMeasures,     OriSpMeasures_ori,     R90_measures];
                    allMeasureBins = [oriMeasure_bins, dirMeasure_bins, OriSpMeasure_bins_ori, R90_measure_bins];
                case 'spf'
                    allMeasures    = [spfMeasures     OriSpMeasures_spf];
                    allMeasureBins = [spfMeasure_bins OriSpMeasure_bins_spf];
            end


            
            measures = allMeasures;
            
            doOriSpfCmps = true;
            

        case 'clusters'
            
            measures = {};
            allGids = unique([allCells.Gid]);
            nGids = length(allGids);
            
            S_refr = getRefrClusterData(allGids);
            
%             allGroupClustIds = cell(1,nGids);
%             allGroupClustNSpk = cell(1,nGids);
%             allGroupCrossPruningMtxs = cell(1, nGids);
%             allGroupNRefrSpikes = cell(1, nGids);
%             for gi = 1:nGids
%                 clusterData = getClusterData(allGids(gi), 2);
%                 if isempty(fieldnames(clusterData))
%                     continue;
%                 end
%                 allGroupClustIds{gi} = clusterData.clustIds;
%                 allGroupClustNSpk{gi} = clusterData.clustNSpikes;
%                 allGroupNRefrSpikes{gi} = clusterData.isiData.pairNRefrSpikes;                
%                 allGroupCrossPruningMtxs{gi} = clusterData.crossRefrData.crossPruningData.nSpksRemoved;                                
%             end
            3;
            
            
        case 'psth', 
            
            n_dPhiBins = 20;
            xBins_dphi = linspace(0, 130, n_dPhiBins+1);     % max: 113

            n_dF1Bins = 30;
            xBins_dF1 = linspace(0, 40, n_dF1Bins+1);  % max: 39
            
            doOriSpfCmps = false;
    end
    
    
    nMeasures = length(measures);
    opt.measures = measures;

    for ms_i = 1:nMeasures
        ind = find(strcmp(measures{ms_i}, allMeasures));
        assert(length(ind) == 1);
        measureBins{ms_i} = allMeasureBins{ind}; %#ok<AGROW>
%         assert(~isempty(measureBins{ms_i}))
    end                                        
    [doDot, doCC, doRho, doTau, doDPhi, doDF1, doSTAcc, doMIDcc, doMID_fit_cc, doMID_ovlp_cc, doDRelPhase] = ...
        dealV( cellfun(@(ms) any(strcmp(measures, ms)), {'dot', 'cc', 'rho', 'tau', 'dphi', 'dF1', 'STA_cc', 'MID_cc', 'MID_fit_cc', 'MID_ovlp_cc', 'dPh_rel'}) );

%     [RF_timeWindow, RF_timeWindow_str] = curTimeWindow;
%     STA_fld = ['STA' RF_timeWindow_str];  STA_odd_fld = [STA_fld '_odd']; STA_odd_fld = [STA_fld '_even'];
%     MID_fld = ['MID' RF_timeWindow_str];  MID_odd_fld = 'MID_odd'; MID_odd_fld = [MID_fld '_even'];

    switch opt.phase_oeState, 
        case 'aa', hasMID = @(CellData)  isfield(CellData, 'MIDdata') && isfield(CellData.MIDdata, 'MID');
                   hasSTA = @(CellData)  isfield(CellData, 'STAdata') && isfield(CellData.STAdata, 'STA');
        case 'oe', hasMID = @(CellData)  isfield(CellData, 'MIDdata') && ( isfield(CellData.MIDdata, 'MID_odd') && isfield(CellData.MIDdata, 'MID_even') );
                   hasSTA = @(CellData)  isfield(CellData, 'STAdata') && ( isfield(CellData.STAdata, 'STA_odd') && isfield(CellData.STAdata, 'STA_even') );
            
    end
    
    dDir_MU_useAllSpkMU = 0;
    dDirMU_field = iff(dDir_MU_useAllSpkMU, 'Ddir_pref_allSpkMU', 'Ddir_pref_smlSpkMU');                
    
    anyWgtedSum = any(strncmp(locations, 'wgtedSum', 6), 1);
    if anyWgtedSum
        fprintf('Multiple PTC: weight by %s. Selection by %s. Phase-tuning selection: %s\n', weightingBy, selectingBy, ptc_criteria);    
    end

    fprintf('-Using these pair types : %s\n', cellstr2csslist(pairTypes));
    fprintf('-Computing these measures : %s\n', cellstr2csslist(measures));
    if strcmp(cmpType, 'phase')
        fprintf('-at these locations : %s\n', cellstr2csslist(locations));
    end
        
    
    if strcmp(cmpType, 'phase')
        switch opt.phase_oeState
            case 'aa', trial_str = 'All-All';
                       keep_str = '';
            case 'oe', trial_str = switchh(opt.phase_oe_mode, {'same', 'diff'}, {'Odd-Odd', 'Odd-Even'});
                       keep_str = switchh(opt.phase_oe_action, {'keep1', 'keep2', 'keepBoth', 'average'}, ...
                           {'(Keeping 1st value)', '(Keeping 2nd value)', '(Keeping both values)', '(Averaging values)'});
                
        end
        fprintf('Phase Tuning Curves / MIDs : comparing %s trials %s\n', trial_str, keep_str);
    elseif strcmp(cmpType, 'degree')
        switch opt.degree_oeState
            case 'aa', trial_str = 'All-All';
                keep_str = '';
            case 'oe', trial_str = switchh(opt.degree_oe_mode, {'same', 'diff'}, {'Odd-Odd', 'Odd-Even'});
                keep_str = switchh(opt.degree_oe_action, {'keepBoth', 'average'}, {'(Keeping both values)', '(Averaging values)'});
        end
        fprintf('Ori/Spf Tuning comparison: use %s trials %s\n', trial_str, keep_str);                
    end
    3;

    if strcmp(cmpType, 'phase') && opt.PTC_require_oe_corr
        oe_str = switchh(opt.PTC_oe_corrMode, {'oe', 'hoe', 'fs'}, {'odd/even', 'half-odd/half-even', 'first/second half'});
        fprintf('* Phase Tuning curves required to have %s cc of at least %.2f\n', oe_str, opt.PTC_minReq_r_oe);
    end
    
    subtractSpont = curSubtractSpont;   
    alwaysUseTuningWithSS = true;
    spont_suffix = iff(subtractSpont || alwaysUseTuningWithSS, '_ss', '_si');
    
    oriStats_fld = ['oriStats' spont_suffix];
    spfStats_fld = ['spfStats' spont_suffix];
    
    
    if doWscc        
        rand('state', 100);
        opt.nPhaseShuffles = 5000;
        keepPhaseOrder = 1;
        
        nShuffleGroups = 10; % to save memory space - do shuffling in chunks
        nPerShuffleGrp = opt.nPhaseShuffles / nShuffleGroups;        
        assert(nPerShuffleGrp == round(nPerShuffleGrp))
        shuffleGrpIdxs = arrayfun(@(i) [1:nPerShuffleGrp] + (i-1)*nPerShuffleGrp, 1:nShuffleGroups, 'un', 0);
        
        assert(~any(strncmp(locations, 'wgtedSum', 7)));
%         assert(~doBcc);
        
        all_nPhases = unique(arrayfun(@(s) length(s.ph), allCells));        
        lcm_nPhases = lcm(all_nPhases(1), all_nPhases(2)); assert( length(all_nPhases) == 2);
%         maxNPh = max( all_nPhases );        
        
        cellGids = arrayfun(@(s) s.Gid, allCells);
        cellCellIds = arrayfun(@(s) s.cellId, allCells);
        [~, groupCellCount] = uniqueCount(cellGids);
        maxNCell = max([groupCellCount; cellCellIds])+1; % 1 for multunits: add 1 to all indices
        cellPairs = nchoosek(1:maxNCell,2);
        nGroupPermMax = size(cellPairs, 1);
        cellPairIdx = zeros(maxNCell, maxNCell);        
        cellPairIdx(sub2indV([maxNCell, maxNCell], cellPairs)        ) = 1:nGroupPermMax;
        cellPairIdx(sub2indV([maxNCell, maxNCell], fliplr(cellPairs))) = 1:nGroupPermMax;        
        
        if (doCC || doRho || doDPhi || doDF1)  % calculate randperm indices only if needed for phase tuning curves
        
            Wscc_randperm_idxs = cell(length(all_nPhases), nGroupPermMax, 2);   
            Wscc_randperm_idxFirst = cell(length(all_nPhases), nGroupPermMax, 1);   
            listAllOrders = zeros(1, length(all_nPhases));

            for nph_i = 1:length(all_nPhases)                        
                nph = all_nPhases(nph_i);                        

                nOrdersPossible = iff(keepPhaseOrder, nph, factorial(nph));
                listAllOrders(nph_i) = opt.nPhaseShuffles > nOrdersPossible;


                if listAllOrders(nph_i)
                    i1 = [1:nph]';
                    idx1 = uint8( i1(:,ones(1,nOrdersPossible)) );
                    if keepPhaseOrder                    
                        idx2 = mod(bsxfun(@plus, i1-1, i1'-1), nph)+1;
                    else
                        idx2 = perms(1:nph)';
                    end
                    idx2 = uint8(idx2);
                    assert(max(all_nPhases < 2^8));
                end                


                for perm_i = 1:nGroupPermMax
                    if listAllOrders(nph_i)
                        whichPerm = randi(nOrdersPossible, 1, opt.nPhaseShuffles);
                        Wscc_randperm_idxs{nph_i, perm_i, 1} = {idx1, idx2, whichPerm};
                        idxFirst = arrayfun(@(i) find(i == whichPerm, 1), 1:nOrdersPossible, 'un', 0);
                        Wscc_randperm_idxFirst{nph_i, perm_i} = [idxFirst{:}];

                    else
                        Wscc_randperm_idxs{nph_i, perm_i, 1} = zeros(nph, opt.nPhaseShuffles, 'uint8');
                        Wscc_randperm_idxs{nph_i, perm_i, 2} = zeros(nph, opt.nPhaseShuffles, 'uint8');
                        for grp_i = 1:nShuffleGroups
                            idx1 = uint8( ord(rand(nph, nPerShuffleGrp), 1) );
                            idx2 = uint8( ord(rand(nph, nPerShuffleGrp), 1) );
                            Wscc_randperm_idxs{nph_i, perm_i, 1}(:, shuffleGrpIdxs{grp_i}) = idx1;
                            Wscc_randperm_idxs{nph_i, perm_i, 2}(:, shuffleGrpIdxs{grp_i}) = idx2;
                        end
                    end
                    
                    
                    
                end

            end
        end
        3;
        
        %{
        Wscc_randperm_idxs = cell(length(all_nPhases), nGroupPermMax, nShuffleGroups, 2);                
        for nph_i = 1:length(all_nPhases)                        
            nph = all_nPhases(nph_i);            
            
%             nPermPossible = factorial(nph);
%             if nPhaseShuffles > nPermPossible
%             else            
            for perm_i = 1:nGroupPermMax                
                for grp_i = 1:nShuffleGroups
                    i1 = uint8( ord(rand(nph, nPerShuffleGrp), 1) );
                    i2 = uint8( ord(rand(nph, nPerShuffleGrp), 1) );
                    Wscc_randperm_idxs{nph_i, perm_i, grp_i, 1} = i1;
                    Wscc_randperm_idxs{nph_i, perm_i, grp_i, 2} = i2;
                end
            end
        end
        %}
        3;
        
        
        
    end

%     nInterp = iff(gratingType == 1, 5, 2);
    nInterp_deltaPhi = 1;
        
    doF1oDC = true;
        doF1oDCforMultipleStim = false;
    calcLowerSumPhs = false;

    % Miscellaneous variables    
%     doWeightedSums = any(strcmp(locations, 'wgtdSum'));  
%     doMtxXMtx      = any(strcmp(locations, 'R1xR2'));  
           
    

    [pathname, cmpDatafile] = getFileName('comparisons', filename_ext);    
            
    pairData = []; %(nPairs,1) = struc;
    pairData2 = [];
    allStatsC = cell(nLocations, nMeasures);                

%     emptyStatStruct = struct('dot', nan, 'cc', nan, 'cc_p', nan, 'rho', nan, 'rho_p', nan, 'tau', nan, 'tau_p', nan, 'dphi', nan, 'dF1', nan, ...
%         'F1oDC1', nan, 'F1oDC2', nan, 'frac1ofMax', nan, 'frac2ofMax', nan, 'sumPhs1', nan, 'sumPhs2', nan, 'cc_atMaxDphi', nan);

    switch cmpType
        case 'phase',  
            dPhi_method = 'cross-correlation';
            dF1_method = 'F1 phase';
        case 'psth', 
            dPhi_method = 'dist between maxes';
            dF1_method =  'dist between COMs';
            doF1oDC = false;
            calcLowerSumPhs = true;
    end            
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    showWorking = false;
    if showWorking
        nWorkCols = 5;
        nWorkRows = 6;
        binIds = round(linspace(1, length(xBins1_1)-1, nWorkCols));
        upTo = ones(nWorkCols, 1);
        binsL = xBins1_1(binIds);
        binsR = xBins1_1(binIds+1);
        figure(900); clf; hold on;
    end
        
    
    function phaseTuningStats = getPairPhaseTuningStats(phs, ph1, ph2, ph1_jackknives, ph2_jackknives, phaseTuningStats, pairPermIdx)
%         [dt, cc, rho, tau, dphi, dF1, cc_p, rho_p, tau_p, F1oDC1, F1oDC2, lowerSumPhs, cc_atMaxDphi]
        
        oneIsZeroOrUniform = any( [all(ph1==0), all(ph2==0), all(diff(ph1)==0), all(diff(ph2)==0)]);
        if oneIsZeroOrUniform
            return;
        end
        

        [nPh, nSamples] = size(ph1);                
        
        if shufflePhases %% if want to shuffle phases for entire file:
            if nSamples == 1
                ph1 = ph1(randperm(nPh));
                ph2 = ph2(randperm(nPh));
            else
                
% %                 idx1 = reshape(randperm(numel(ph1)), size(ph1));
% %                 idx2 = reshape(randperm(numel(ph2)), size(ph2));
% 
%                 idx1 = randperm(nPh);
%                 idx2 = randperm(nPh);
%                 
%                 ph1 = ph1(idx1,:);
%                 ph2 = ph2(idx2,:);                
            end
        end        
        
        
        if exist('pairPermIdx', 'var') && ~isempty(pairPermIdx)
            
            nph_idx = find(all_nPhases == nPh, 1);
            
            if listAllOrders(nph_idx)
                [perm_idx1, perm_idx2, result_perm_idx] = deal(Wscc_randperm_idxs{nph_idx, pairPermIdx, 1}{:});
                nSamples = size(perm_idx1,2);
            else
                perm_idx1 = Wscc_randperm_idxs{nph_idx, pairPermIdx, 1};
                perm_idx2 = Wscc_randperm_idxs{nph_idx, pairPermIdx, 2};
                nSamples = opt.nPhaseShuffles;
            end
            ph1 = ph1(perm_idx1);
            ph2 = ph2(perm_idx2);        
            if ~isempty(ph1_jackknives)
                %%
                ph1_jackknives_orig = ph1_jackknives(:);
                ph2_jackknives_orig = ph2_jackknives(:);
                %%
                ph1_jackknives = cellfun(@(ph1) ph1(perm_idx1), ph1_jackknives_orig, 'un', 0);  % ph1_jackknives(perm_idx1);
                ph2_jackknives = cellfun(@(ph2) ph2(perm_idx2), ph2_jackknives_orig, 'un', 0);  % ph1_jackknives(perm_idx1);
                
                % every column consists of the nTrials jackknifed phase tuning curves. 
                
            end
                

            fn = fieldnames(phaseTuningStats);
            if (nSamples == opt.nPhaseShuffles)
                        
                nShuffleGroups = 20;
                nPerShuffleGrp = opt.nPhaseShuffles / nShuffleGroups;        
                shuffleGrpIdxs = arrayfun(@(i) [1:nPerShuffleGrp] + (i-1)*nPerShuffleGrp, 1:nShuffleGroups, 'un', 0);

                phaseTuningStats_grp = phaseTuningStats;                
                for fld_i = 1:length(fn)
                    phaseTuningStats_grp.(fn{fld_i}) = zeros(1, nSamples);
                end

                for samp_grp_i = 1:nShuffleGroups
                    grp_idx = shuffleGrpIdxs{samp_grp_i};
                    phaseTuningStats_i = getPairPhaseTuningStats(phs, ph1(:,grp_idx), ph2(:,grp_idx), phaseTuningStats); % already shuffled above, so don't need to shuffle again

                    for fld_i = 1:length(fn)
                        phaseTuningStats_grp.(fn{fld_i})(grp_idx) = phaseTuningStats_i.(fn{fld_i});
                    end

                end
                phaseTuningStats = phaseTuningStats_grp;
            else
                phaseTuningStats = getPairPhaseTuningStats(phs, ph1, ph2, ph1_jackknives, ph2_jackknives,  phaseTuningStats); % already shuffled above, so don't need to shuffle again

                for fld_i = 1:length(fn)
                    field_vals = phaseTuningStats.(fn{fld_i});
                    if length(field_vals) > 1
                        phaseTuningStats.(fn{fld_i}) = field_vals(result_perm_idx);
                    end
                end
                
                
            end
            
            
            return;                                    
        end
        
        
        
%         ph1 = ph1(:);
%         ph2 = ph2(:);
        
%         s0 = 0;
%         dt = s0; cc = s0; rho = s0; tau = s0; dphi = s0; dF1 = s0; cc_p = s0; rho_p = s0; tau_p = s0; lowerSumPhs = 0; F1oDC1 = 0; F1oDC2 = 0; 
%         cc_atMaxDphi = s0;
        if nSamples == 1
            pearsonFunc = @pearsonR;
            spearmanFunc = @spearmanRho;
            tauFunc = @kendallTau;
        else
            pearsonFunc = @(x,y) pearsonR(x,y,1);
            spearmanFunc = @spearmanRho_v;
            tauFunc = @kendallTau_v;
            3;
        end
            
        if doDot,  phaseTuningStats.dot = normDotProd(ph1, ph2);      end
        % cc, rho, tau
        if calcPvalues; % want p-values as well            
            if doCC,  [phaseTuningStats.cc, phaseTuningStats.cc_p] = pearsonFunc(ph1, ph2); end
            if doRho, [phaseTuningStats.rho, phaseTuningStats.rho_p] = spearmanFunc(ph1, ph2); end
            if doTau, [phaseTuningStats.tau, phaseTuningStats.tau_p] = tauFunc(ph1, ph2);  end            
        else  
            if doCC,  phaseTuningStats.cc = pearsonFunc(ph1, ph2);      end
            if doRho, phaseTuningStats.rho = spearmanFunc(ph1, ph2);  end
            if doTau, phaseTuningStats.tau = tauFunc(ph1, ph2);   end      
            
            
            if opt.getPhaseCmpJackStdErr_ptc && doCC && ~isempty(ph1_jackknives) && ~isempty(ph2_jackknives)
                if nSamples == 1
                    cc_jackknives = cellfun(pearsonFunc, ph1_jackknives, ph2_jackknives);
                    phaseTuningStats.cc_jackStd = jackknifeStdErr(cc_jackknives, phaseTuningStats.cc);
                else
                    cc_jackknives = cellfun(pearsonFunc, ph1_jackknives, ph2_jackknives, 'un', 0);
                    cc_jackknives_allShifts = [cat(1, cc_jackknives{:})];
                    cc_jackStd = arrayfun(@(i) jackknifeStdErr(cc_jackknives_allShifts(:,i), phaseTuningStats.cc(i)), 1:nSamples );
                    phaseTuningStats.cc_jackStd = cc_jackStd;
                end
                
            end
            
            
        end   
        if doCC 
            assert( all( phaseTuningStats.cc >= -1 & phaseTuningStats.cc <= 1) );
%             phaseTuningStats.cc(phaseTuningStats.cc < -1) = -1;
        end                                    
            
        
        % F1/DC
%         if (doF1oDC || doDPhi) && (nSamples == 1) % if nargout > 7  % also provide F1/DC
%             F1oDC1 = getF1oDC( single(phs), single( ph1'), 360 )';
%             F1oDC2 = getF1oDC( single(phs), single( ph2'), 360 )';
%             phaseTuningStats.F1oDC1 = F1oDC1;
%             phaseTuningStats.F1oDC2 = F1oDC2;
%         end
               
        
%         if (doDPhi || doDF1) && (nSamples > 1)
%             idx_use = find( (F1oDC1 > 1e-10) & (F1oDC2 > 1e-10)); % either all zeros, or [0 eps 0 eps].
%             nUse = length(idx_use);
%         end
        
        if doDPhi
            assert(nInterp_deltaPhi == 1);
            
            if length(phs) == size(ph1, 1)                
%                 dphi = bestCircShift_Matlab(phs, ph1, ph2);  
%                 dphi2 = bestCircShift(phs, double(ph1), double(ph2));  
%                 if any(isnan(dphi)) && ~opt.makeDphiNullWhenAmbiguous
%                     [dphi, all_ccShifts] = bestCircShift_Matlab(phs, double(ph1), double(ph2), [], opt.usePosNegDphi);  
                    dPhiFunc = @(p1, p2) bestCircShift_Matlab(phs, p1, p2, [], opt.usePosNegDphi);  
                    
%                     dphi = bestCircShift_Matlab(phs, ph1, ph2, [], opt.usePosNegDphi);  
                    phaseTuningStats.dphi = dPhiFunc(ph1, ph2);  
                    
                    if opt.getPhaseCmpJackStdErr_ptc  && ~isempty( ph1_jackknives ) && ~isempty(ph2_jackknives)
                        if nSamples == 1
                            dphi_jackknives = cellfun(dPhiFunc, ph1_jackknives, ph2_jackknives);
                            phaseTuningStats.dphi_jackStd = jackknifeStdErr(dphi_jackknives, phaseTuningStats.dphi);
                        else
                            dphi_jackknives = cellfun(dPhiFunc, ph1_jackknives, ph2_jackknives, 'un', 0);
                            dphi_jackknives_allShifts = [cat(1, dphi_jackknives{:})];
                            dphi_jackStd = arrayfun(@(i) jackknifeStdErr(nonnans(dphi_jackknives_allShifts(:,i)), phaseTuningStats.dphi(i)), 1:nSamples );
                            phaseTuningStats.dphi_jackStd = dphi_jackStd;
%                             assert (all( dphi_jackStd == dphi_jackStd(1)) || isnan(dphi_jackStd(1)) )
                        end
                        
                        
                    end

                    
%                     assert(isequaln(dphi2, dphi))
                if size(phs, 1) == size(phs, 2)
                    dphi3 = min(phs, 360-phs);
                    assert(isequal(sort(dphi3), sort(dphi)))
                end
                    
                    
%                     dphi = bestCircShift_Matlab(phases, double(ph_tc1), double(ph_tc2), 0, opt.usePosNegDphi)
%                 end                
%                 dphi = bestCircShift_Matlab(phs, double(ph1), double(ph2));  

            else  % for mtx x mtx case     
%                 dphi = 0;            
                phaseTuningStats.dphi = 0;
            end
            
%             phaseTuningStats.cc_atMaxDphi = cc_atMaxDphi;
        end             
        
        if doDF1 && length(phs) == size(ph1, 1)
            
            if nSamples == 1       
                if any([phaseTuningStats.F1oDC1, phaseTuningStats.F1oDC2] < 1e-5)
                    dF1 = nan;
                else                    
                    ph1_F1 = firstNHarmonics(ph1, 1); % keep only first harmonic
                    ph2_F1 = firstNHarmonics(ph2, 1);
                    dF1 = bestCircShift_Matlab(phs, ph1_F1, ph2_F1);
%                     dF1 = deltaPhi( phs, ph1_F1, ph2_F1, dPhi_method);
                end
            elseif nSamples > 1
%%
                idx_use = 1:nSamples; nUse = length(idx_use);
                ph1_F1 = firstNHarmonics(ph1(:,idx_use), 1); % keep only first harmonic
                ph2_F1 = firstNHarmonics(ph2(:,idx_use), 1);
                
                dF1 = nan(1,nSamples);  
                if nUse > 0
%                     dF1(idx_use)  = deltaPhi( phs, ph1_F1, ph2_F1, dPhi_method); 
                    dF1(idx_use) = bestCircShift_Matlab(phs, ph1_F1, ph2_F1);
                end
                
                chk = 0;
                if chk
                    dF1b = nan(1,nSamples);                    
                    for i = 1:nUse;
                        ph1_F1 = firstNHarmonics(ph1(:,idx_use(i)), 1); % keep only first harmonic
                        ph2_F1 = firstNHarmonics(ph2(:,idx_use(i)), 1);
                        dF1b(idx_use(i)) = bestCircShift_Matlab( phs, ph1_F1, ph2_F1);                        
%                         dF1b(idx_use(i)) = deltaPhi( phs, ph1_F1, ph2_F1, dPhi_method);
                    end                                        
                    assert(isequaln(dF1, dF1b));                    
                end
                
            end
                                            
            if doDPhi
%                 dF1(isnan(dphi)) = nan;            
            end
            phaseTuningStats.dF1 = dF1;
        end

                
        if calcLowerSumPhs  % also return smaller sum of phases.
            sum1 = sum(ph1);  sum2 = sum(ph2);
            if sum1 <= sum2,
                lowerSumPhs = sum1;
            else
                lowerSumPhs = sum2;
            end
            phaseTuningStats.lowerSumPhs = lowerSumPhs;
        end

        if 0&&(length(phs) > 5) && all([F1oDC1 F1oDC2] > .5) && (rand < .1)

            [phi, f_cos1, t1] = getF1phase(deg2rad(phs), ph1, 2*pi);
            [phi, f_cos2, t2] = getF1phase(deg2rad(phs), ph2, 2*pi);
            [phi, f_cos1n, t1n] = getF1phase(deg2rad(phs), ph1/max(ph1), 2*pi);
            [phi, f_cos2n, t2n] = getF1phase(deg2rad(phs), ph2/max(ph2), 2*pi);
            
            figure(19); clf;
            plot(phs, ph1, 'bo-', phs, ph2, 'gs-'); hold on;
            plot(rad2deg(t1), f_cos1, 'b:'); plot(rad2deg(t2), f_cos2, 'g:');
            figure(20); clf;
            plot(phs, ph1/(max(ph1)), 'bo-', phs, ph2/(max(ph2)), 'gs-'); hold on;
            plot(rad2deg(t1n), f_cos1n, 'b:'); plot(rad2deg(t2n), f_cos2n, 'g:');
            xlim([0 360]);
            3;
        end
        
%         prod_cc_rho = phaseTuningStats.cc .* phaseTuningStats.rho;
%         if 0 && prod_cc_rho < .1
%             figure(19); clf;
%             plot(phs, ph1, 'bo-', phs, ph2, 'gs-'); hold on;
% %             figure(20); clf;            
%         end
        
%         if 0&& (dphi == 0) && (length(phs) == 60)
%             figure(123); clf;            
%             plot(phs, ph1, 'bo-', phs, ph2, 'gs-');
%             3;
%         end        
        
    end

    function pairPhaseStruct_av = getPairPhaseTuningStats_OE(phases, ph_tc1_a,    ph_tc1_b,    ph_tc2_a,    ph_tc2_b,  ...
                                                                     ph_jacks1_a, ph_jacks1_b, ph_jacks2_a, ph_jacks2_b, blankPhaseStruct, varargin)
        %%
        pairPhaseStruct_a = getPairPhaseTuningStats(phases, ph_tc1_a, ph_tc2_a,  ph_jacks1_a,  ph_jacks2_a,  blankPhaseStruct, varargin{:});
        pairPhaseStruct_b = getPairPhaseTuningStats(phases, ph_tc1_b, ph_tc2_b,  ph_jacks1_b,  ph_jacks2_b,  blankPhaseStruct, varargin{:});
        pairPhaseStruct_av = blankPhaseStruct;
        
        fn = fieldnames(pairPhaseStruct_av);        
        
        for fld_i = 1:length(fn)
            fld = fn{fld_i};
            pairPhaseStruct_av.(fld) = combineValues(pairPhaseStruct_a.(fld), pairPhaseStruct_b.(fld), fld, opt);                            
        end
    end


    function oriPairStats = getOriDegreeTuningStats_for2Cells(cell1_stats, cell2_stats, oriPairStats, opt)        
                
        switch opt.degree_oeState
            case 'aa', 
                oriPairStats = getOriDegreeTuningStats_for2Stats(cell1_stats.(oriStats_fld), cell2_stats.(oriStats_fld), oriPairStats);
            case 'oe',
                fn_odd = [oriStats_fld '_odd'];
                fn_even = [oriStats_fld '_even'];                           
                
%                 case 'same',  [idx1a, idx2a,   idx1b, idx2b] = deal(1,1, 2,2); 
%                 case 'diff',  [idx1a, idx2a,   idx1b, idx2b] = deal(1,2, 2,1); 
                switch opt.degree_oe_mode
                    case 'same',  [fn_1a, fn_2a,   fn_1b, fn_2b] = deal(fn_odd, fn_odd,   fn_even, fn_even); % 1,2 : cell 1/cell 2
                    case 'diff',  [fn_1a, fn_2a,   fn_1b, fn_2b] = deal(fn_odd, fn_even,   fn_even, fn_odd); % a,b : first/second comparison
                end
                if (isnan(cell1_stats.(fn_odd).w_ori_global) || isnan(cell1_stats.(fn_even).w_ori_global))  && ...
                        ~isnan(  cell1_stats.(oriStats_fld).w_ori_global );
                    3;
                end
                
                oriPairStats_a = getOriDegreeTuningStats_for2Stats(cell1_stats.(fn_1a), cell2_stats.(fn_2a), oriPairStats);
                oriPairStats_b = getOriDegreeTuningStats_for2Stats(cell1_stats.(fn_1b), cell2_stats.(fn_2b), oriPairStats);
                %%                
                fn = fieldnames(oriPairStats_a);
                for fi = 1:length(fn)
                    fld = fn{fi};
                    oriPairStats.(fld) = combineValues(oriPairStats_a.(fld), oriPairStats_b.(fld), fld, opt);
                end
                
%                 st1 = struct2array(oriPairStats_a);
%                 st2 = struct2array(oriPairStats_b);
%                 st_av = (st1 + st2)/2
                
                3;
                
        end        
        
    end

    function oriPairStats = getOriDegreeTuningStats_for2Stats(s1, s2, oriPairStats)        

        
        doD_ori_pref_thisPair = 1;
        doD_w_ori_global_thisPair = 1;
        doD_w_ori_local_thisPair = 1;
        doD_dsi_thisPair = 1;
        doD_F1oDC_thisPair = 1;

        if opt.applyStdErrorThresholds && strcmp(opt.degree_oeState, 'aa') 
            err1 = s1.error_jack;
            err2 = s2.error_jack;
            if any([err1.ori_pref, err2.ori_pref] > opt.maxOriStdErr_deg)
                doD_ori_pref_thisPair = 0;
            end
            if any([err1.w_ori_global, err2.w_ori_global] > opt.maxOriStdErr_deg)
                doD_w_ori_global_thisPair = 0;
            end
            if any([err1.w_ori_local, err2.w_ori_local] > opt.maxOriStdErr_deg)
                doD_w_ori_local_thisPair = 0;
            end
            if any([err1.DSI_global, err2.DSI_global] > opt.maxOriStdErr_deg)
                doD_dsi_thisPair = 0;
            end
            if any([err1.F1oDC, err2.F1oDC] > opt.maxF1oDCErr)
                doD_F1oDC_thisPair = 0;
            end

        end
            
                
%                     opt.maxOriStdErr_deg = 5;
%     opt.maxDSIStdErr = 0.1;
%     opt.maxSpfStdErr = 0.5;
        
        
        oriMax = 180;
        dirMax = 360;
        
        if doD_w_ori_global_thisPair
            oriPairStats.Dw_ori_glob = circDist( s1.w_ori_global,   s2.w_ori_global,   oriMax);
        end
        if doD_w_ori_local_thisPair
            oriPairStats.Dw_ori_loc  = circDist( s1.w_ori_local(1), s2.w_ori_local(1), oriMax);
        end
        if doD_ori_pref_thisPair
            oriPairStats.D_ori_pref  = circDist( s1.ori_pref_deg,   s2.ori_pref_deg,   oriMax);
        end
        if doDirMeasures
            if doD_dsi_thisPair
                oriPairStats.D_dsi_glob = abs(s1.DSI_global - s2.DSI_global );
                if isfield(s1, 'DSI_global_unrec')
                    oriPairStats.D_dsi_glob_unrec = abs(s1.DSI_global_unrec - s2.DSI_global_unrec );
                end
                oriPairStats.D_dsi_loc  = abs(s1.DSI_local  - s2.DSI_local  );
                if isfield(s1, 'DSI_local_unrec')
                    oriPairStats.D_dsi_loc_unrec  = abs(s1.DSI_local_unrec  - s2.DSI_local_unrec  );
                end
            end            
            
            if doD_ori_pref_thisPair
                oriPairStats.D_dir_pref    = circDist( s1.dir_pref_deg, s2.dir_pref_deg, dirMax);            

                cell1_aligned_factor = getCellAligned_factor(s1.(dDirMU_field)); % aligned = 1. anti-aligned: 10, unaligned: 100
                cell2_aligned_factor = getCellAligned_factor(s2.(dDirMU_field));                        
                oriPairStats.D_aligned_pair = cell1_aligned_factor + cell2_aligned_factor; % a/a : 2. a/aa : 11. aa/aa: 20
            
            end
        end
        
        if doOSmeasures
            if strcmp(gratingType, 'flashed')
                oriPairStats.cc_OSP = pearsonR(s1.OS, s2.OS);
            end
            if length(s1.r_k_dir) == length(s2.r_k_dir)
                oriPairStats.cc_ori = pearsonR(s1.r_k_dir, s2.r_k_dir);            
            end
        end
        
        if doR90_measures
            oriPairStats.dR_spont_abs   = abs( s1.R_spont_abs   - s2.R_spont_abs   );
            oriPairStats.dR90_total_abs = abs( s1.R90_total_abs - s2.R90_total_abs );
            oriPairStats.dR90_stim_abs  = abs( s1.R90_stim_abs  - s2.R90_stim_abs  );
            oriPairStats.dR_spont_rel   = abs( s1.R_spont_rel   - s2.R_spont_rel   );
            oriPairStats.dR90_total_rel = abs( s1.R90_total_rel - s2.R90_total_rel );
            oriPairStats.dR90_stim_rel  = abs( s1.R90_stim_rel  - s2.R90_stim_rel  ); 
        end
        
        if doD_F1oDC_thisPair
            oriPairStats.D_F1oDC_ori  = abs( s1.F1oDC - s2.F1oDC);
            F1_type1 = s1.F1oDC > 1;
            F1_type2 = s2.F1oDC > 1;
            oriPairStats.D_F1pair_ori   = F1_type1 + F1_type2; % 0 = Complex/Complex;  1 = Simple/Complex;  2 = Simple/Simple        
            3;
        end    

        
    end
    

    function spfPairStats = getSpfDegreeTuningStats_for2Cells(cell1_stats, cell2_stats, spfPairStats)
 
        switch opt.degree_oeState
            case 'aa', 
                spfPairStats = getSpfDegreeTuningStats_for2Stats(cell1_stats.(spfStats_fld), cell2_stats.(spfStats_fld), spfPairStats);
            case 'oe',
                fn_odd = [spfStats_fld '_odd'];
                fn_even = [spfStats_fld '_even'];                           
                
%                 case 'same',  [idx1a, idx2a,   idx1b, idx2b] = deal(1,1, 2,2); 
%                 case 'diff',  [idx1a, idx2a,   idx1b, idx2b] = deal(1,2, 2,1); 
                switch opt.degree_oe_mode
                    case 'same',  [fn_1a, fn_2a,   fn_1b, fn_2b] = deal(fn_odd, fn_odd,   fn_even, fn_even); % 1,2 : cell 1/cell 2
                    case 'diff',  [fn_1a, fn_2a,   fn_1b, fn_2b] = deal(fn_odd, fn_even,   fn_even, fn_odd); % a,b : first/second comparison
                end
                           
                spfPairStats_a = getSpfDegreeTuningStats_for2Stats(cell1_stats.(fn_1a), cell2_stats.(fn_2a), spfPairStats);
                spfPairStats_b = getSpfDegreeTuningStats_for2Stats(cell1_stats.(fn_1b), cell2_stats.(fn_2b), spfPairStats);
                  
                fn = fieldnames(spfPairStats_a);
                for fi = 1:length(fn)
                    if oe_keepBoth
                        spfPairStats.(fn{fi}) = cat(3, spfPairStats_a.(fn{fi}), spfPairStats_b.(fn{fi}));
                    else
                        spfPairStats.(fn{fi}) = (spfPairStats_a.(fn{fi}) + spfPairStats_b.(fn{fi}))/2;
                    end
                end
                
%                 spfPairStats = structfun(@(a,b) (a+b)/2, spfPairStats_a, spfPairStats_b, 'un', 0);
        end        

    end

    function spfPairStats = getSpfDegreeTuningStats_for2Stats(s1, s2, spfPairStats)
%         if ~isfield(cell1_stats, 'spfStats_si')
%             return;
%         end
%         s1 = cell1_stats.spfStats_si;
%         s2 = cell2_stats.spfStats_si;

        doD_spf_pref_thisPair = 1;
        doD_w_spf_thisPair = 1;
        doD_F1oDC_thisPair = 1;

        if opt.applyStdErrorThresholds && strcmp(opt.degree_oeState, 'aa')
            err1 = s1.error_jack;
            err2 = s2.error_jack;
            if any([err1.w_spf, err2.w_spf] > opt.maxSpfStdErr)
                doD_w_spf_thisPair = 0;
            end
            if any([err1.f_opt, err2.f_opt] > opt.maxSpfStdErr)
                doD_spf_pref_thisPair = 0;
            end
            if any([err1.F1oDC, err2.F1oDC] > opt.maxF1oDCErr)
                doD_F1oDC_thisPair = 0;
            end
        end

        if doD_w_spf_thisPair
            spfPairStats.Dw_spf     = abs( s1.w_spf       - s2.w_spf);
            spfPairStats.D_w_SLN    = abs( s1.SLNparams.w - s2.SLNparams.w);
            spfPairStats.D_Bn       = abs( s1.Bn          - s2.Bn);
            spfPairStats.D_skew     = abs( s1.SLNparams.s - s2.SLNparams.s);
        end
        
        if doD_spf_pref_thisPair
            spfPairStats.D_spf_pref = abs( log2( s1.f_opt / s2.f_opt) );  % == abs( s1.SLNparams.f_opt_log - s2.SLNparams.f_opt_log );
        end
        
        if isnan( spfPairStats.D_spf_pref )
            3;
        end
        
        if doD_F1oDC_thisPair
            spfPairStats.D_F1oDC_spf  = abs( s1.F1oDC - s2.F1oDC);
            F1_type1 = s1.F1oDC > 1;
            F1_type2 = s2.F1oDC > 1;
            spfPairStats.D_F1pair_spf   = F1_type1 + F1_type2; % 0 = Complex/Complex;  1 = Simple/Complex;  2 = Simple/Simple        
        end    
        
        if doOSmeasures
            if length(s1.spf_tc) == length(s2.spf_tc)
                spfPairStats.cc_spf = pearsonR(s1.spf_tc, s2.spf_tc);
            end
        end

        
    end


   

    function psthStats = getStatsFor2PSTHs(psth1, psth2)
        
        [bins, bins2, vals1, vals2] = deal(psth1.bins, psth2.bins, psth1.vals, psth2.vals);
        assert(isequal(bins, bins2));
        
        restrictToWindowUnion = false;
        showWorking = false;        
        
        if restrictToWindowUnion
            idxs1 = [find(bins' > psth1.timeWindow(1), 1, 'first')-1 : find(bins' < psth1.timeWindow(2), 1, 'last')+3];
            idxs2 = [find(bins' > psth2.timeWindow(1), 1, 'first')-1 : find(bins' < psth2.timeWindow(2), 1, 'last')+3];
            
            idx_start = min([idxs1(:); idxs2(:)]);
            idx_end   = max([idxs1(:); idxs2(:)]);
            uidx = idx_start:idx_end;
            
            [bins, vals1, vals2]  = deal(bins(uidx), vals1(uidx), vals2(uidx));            
        end
        
        [dt, cc, rho, tau, dphi, dF1] = getPairPhaseTuningStats(bins, vals1, vals2);        
        psthStats = struct('dot', dt, 'cc', cc, 'rho', rho, 'tau', tau, 'dphi', dphi, 'dF1', dF1, ...
            'F1oDC1', 0, 'F1oDC2', 0,  'frac1ofMax', 0, 'frac2ofMax', 0, 'sumPhs1', 0, 'sumPhs2', 0, 'cc_atMaxDphi', 0);        
        
        if showWorking
            figure(154); clf
            plot2PSTHs(psth1.bins, psth1.vals, psth2.vals);
            drawVerticalLine(bins2(idxs1([1, end])), 'color', 'b');
            drawVerticalLine(bins2(idxs2([1, end])), 'color', 'r');
            drawVerticalLine(bins2(uidx([1, end])), 'color', 'g');
            title(sprintf('dot = %.2f, cc = %.2f, \rho = %.2f, \Delta \phi = %.2f, \DeltaF1 = %.2f', dt, cc, rho, dphi, dF1))
            3;            
        end        
        
    end


    function [curPairStats, oriSpCmpInds, phaseTCs, shuffPairStats] = getPhaseTuningStatsFor2OSPs(...
            phases, R1, R2, pairInfo, blankPhaseStruct, blankShuffPairStats, R1_oe, R2_oe)                                         

        curPairStats(1:nLocations) = blankPhaseStruct;
        oriSpCmpInds = zeros(2,nLocations, 'uint8');
        phaseTCs = zeros(length(phases),nLocations, 2, 'single');
        
        pairIdx = pairInfo.pairIdx;
        shuffPairStats(1:nLocations) = blankShuffPairStats;
        
        haveR_oe = exist('R1_oe', 'var') && exist('R2_oe', 'var');
%         if (nargin < 8) 
            [R1_forIdx, R2_forIdx] = deal(R1, R2);
%         end        
        
        compressMethod = 'mean'; %options: 'max', 'mean'
        
        [nOri, nSpf, nPh] = size(R1);           
        R1_stim_ph = reshape(R1, [nOri*nSpf, nPh])';
        R2_stim_ph = reshape(R2, [nOri*nSpf, nPh])';
        
        R1_forIdx_stim_ph = reshape(R1_forIdx, [nOri*nSpf, nPh])';
        R2_forIdx_stim_ph = reshape(R2_forIdx, [nOri*nSpf, nPh])';
        
        if strcmp(compressMethod, 'max')
            mR1 = max(R1_forIdx_stim_ph, [], 1);  
            mR2 = max(R2_forIdx_stim_ph, [], 1);  
        elseif strcmp(compressMethod, 'mean')
            mR1 = mean(R1_forIdx_stim_ph, 1);  
            mR2 = mean(R2_forIdx_stim_ph, 1);  
        end        
        fR1 = mR1/max(mR1(:));
        fR2 = mR2/max(mR2(:));                                  
        minFracR = min(fR1, fR2);        
                    
        sumPhs1 = sum(R1_stim_ph,1); normSumRates1 = sumPhs1/max(sumPhs1);
        sumPhs2 = sum(R2_stim_ph,1); normSumRates2 = sumPhs2/max(sumPhs2);

%         nRespPh1 = sum(R1_stim_ph>0,1);
%         nRespPh2 = sum(R2_stim_ph>0,1);                                                    

        if opt.PTC_require_oe_corr
            R1_phaseTC_oe = pairInfo.phaseTC_cc1(:);
            R2_phaseTC_oe = pairInfo.phaseTC_cc2(:);
            min_phaseTC_oe = min(R1_phaseTC_oe, R2_phaseTC_oe);            
            neitherIsNan = ~isnan(R1_phaseTC_oe) & ~isnan(R2_phaseTC_oe);
            
            phaseTC_r_oe_ok = (min_phaseTC_oe > opt.PTC_minReq_r_oe) & neitherIsNan;

            fR1(~phaseTC_r_oe_ok) = nan;
            fR2(~phaseTC_r_oe_ok) = nan;
            maxMinFracR_noOE = max(minFracR);
            minFracR(~phaseTC_r_oe_ok) = nan;
            maxMinFracR = max(minFracR);
            minMinFracR = min(minFracR);
            if maxMinFracR_noOE > .8 && maxMinFracR < .5
                3;
            end
            maxMinFracR_cmp(end+1,:) = [maxMinFracR_noOE, maxMinFracR];
            
            if (maxMinFracR_noOE > .98 & maxMinFracR < .1)
                3;
                
            end
            
        end


        if opt.PTC_require_oe_corr_p
            R1_phaseTC_p_oe = pairInfo.phaseTC_cc_p1(:);
            R2_phaseTC_p_oe = pairInfo.phaseTC_cc_p2(:);
            max_phaseTC_p_oe = max(R1_phaseTC_p_oe, R2_phaseTC_p_oe)';
            phaseTC_p_oe_ok = max_phaseTC_p_oe < opt.PTC_maxReq_p_oe;

            fR1(~phaseTC_p_oe_ok) = nan;
            fR2(~phaseTC_p_oe_ok) = nan;
            minFracR(~phaseTC_p_oe_ok) = nan;
        end
        
        
        for loc_i = 1:length(locations)                        
            
            locationName = locations{loc_i};            
            explicitOE_measure = ~isempty(strfind(locations{loc_i}, '_oe'));
            explicitAA_measure = ~isempty(strfind(locations{loc_i}, '_aa'));            
                            
            if explicitOE_measure
                locationName = strrep(locationName, '_oe', '');
                isOEmeasure = 1;
            elseif explicitAA_measure 
                locationName = strrep(locationName, '_aa', '');
                isOEmeasure = 0;
            elseif ~explicitAA_measure && ~explicitOE_measure
                isOEmeasure = strcmp(opt.phase_oeState, 'oe');
            end
                            
            if strncmp(locationName, 'wgtedSum', 8);
                locationName = 'wgtedSum';
            end      

            
            if isOEmeasure && ~exist('R1_stim_ph_oe', 'var')
                %%
                R1_stim_ph_oe = permute(reshape(R1_oe, [nOri*nSpf, nPh, 2]), [2 1 3]);
                R2_stim_ph_oe = permute(reshape(R2_oe, [nOri*nSpf, nPh, 2]), [2 1 3]);
                
                sumPhs1_oe = sum(R1_stim_ph_oe,1); normSumRates1_oe = sumPhs1_oe/max(sumPhs1_oe(:));
                sumPhs2_oe = sum(R2_stim_ph_oe,1); normSumRates2_oe = sumPhs2_oe/max(sumPhs2_oe(:));
                nRespPh1_oe = min(sum(R1_stim_ph_oe>0,1), [], 3);
                nRespPh2_oe = min(sum(R2_stim_ph_oe>0,1), [], 3);                
            end
            
            
            switch locationName
                case {'maxR1', 'maxR2', 'mean12', 'maxR1xR2', 'maxMU', 'maxMinFracR', 'maxMinFracR2', 'maxMinFracR3', 'minMinFracR', 'p75MinFracR', 'p50MinFracR', 'p25MinFracR'}  % 'max F1oDC 1x2', 'max F1oDC MU'
                    
                    if strncmp(locations{loc_i}, 'maxMU', 5) && (~isfield(pairInfo, 'ori_sp_maxMU') || isempty(pairInfo.ori_sp_maxMU)), continue; end % for maxMU, some sites don't have multi-units.
                    
                    [fracR_sorted, ordFracR] = sort(minFracR, 'descend');
                    
                    switch locationName
                        case 'maxR1',       [fracR_val, stim_idx] = max(fR1);
                        case 'maxR2',       [fracR_val, stim_idx] = max(fR2);
                        case 'maxR1xR2',    [fracR_val, stim_idx] = max(fR1 .* fR2);
                        case 'maxMU',       stim_idx = sub2indV([nOri, nSpf], pairInfo.ori_sp_maxMU); fracR_val = minFracR(stim_idx);
                        case 'maxMinFracR', [fracR_val, stim_idx] = max(minFracR); assert(stim_idx == ordFracR(1));
                        case 'maxMinFracR2', [stim_idx] = ordFracR(2);
                        case 'maxMinFracR3', [stim_idx] = ordFracR(3);
                        case 'p75MinFracR', [fracR_val, stim_idx] = min(abs( minFracR - prctile(minFracR, 75)));
                        case 'p50MinFracR', [fracR_val, stim_idx] = min(abs( minFracR - prctile(minFracR, 50)));
                        case 'p25MinFracR', [fracR_val, stim_idx] = min(abs( minFracR - prctile(minFracR, 25)));
                        case 'minMinFracR', [fracR_val, stim_idx] = min(minFracR);
                    end                    
%                     oriSpInds.(locationName) = stim_idx;
                    
                    ori_sp_idx = ind2subV([nOri, nSpf], stim_idx);
                    oriSpCmpInds(:,loc_i) = ori_sp_idx(:);     
                    cell1_oe_corr = pairInfo.phaseTC_cc1(stim_idx);
                    cell2_oe_corr = pairInfo.phaseTC_cc2(stim_idx);
                    min_ptcOEcorr = min(cell1_oe_corr, cell2_oe_corr);                    
%                     curPairStats(loc_i).ptcOEcorr = cat(3, cell1_oe_corr, cell2_oe_corr);                    
%                     ori_i = ori_sp(1); sp_i = ori_sp(2);
                    blankPhaseStruct_i = blankPhaseStruct;                    
                    blankPhaseStruct_i.F1oDC1 = pairInfo.F1oDCs1(stim_idx);
                    blankPhaseStruct_i.F1oDC2 = pairInfo.F1oDCs2(stim_idx);
                    blankPhaseStruct_i.frac1ofMax = fR1(stim_idx);
                    blankPhaseStruct_i.frac2ofMax = fR2(stim_idx);

               
                    
                    if ~isOEmeasure
                        
                        ph_tc1 = R1_stim_ph( :, stim_idx);  
                        ph_tc2 = R2_stim_ph( :, stim_idx);

                        phaseTCs(:, loc_i, 1) = ph_tc1;
                        phaseTCs(:, loc_i, 2) = ph_tc2;

                        
                        if pairInfo.sameGroup && opt.getPhaseCmpJackStdErr_ptc
                            ph1_jacks = pairInfo.Cell1.R_stim_jackTrials(stim_idx,:);
                            ph2_jacks = pairInfo.Cell2.R_stim_jackTrials(stim_idx,:);
                        else
                            ph1_jacks = [];
                            ph2_jacks = [];
                        end
                        
%                         oneIsZeroOrUniform = any( [all(ph_tc1<=0), all(ph_tc2<=0), all(diff(ph_tc1)==0), all(diff(ph_tc2)==0)]);
%                         invalid_frac = (fracR_val == 0) || isnan(fracR_val);
%                         if ~oneIsZeroOrUniform && ~invalid_frac
                            pairPhaseStats = getPairPhaseTuningStats(phases, ph_tc1(:), ph_tc2(:),  ph1_jacks, ph2_jacks,  blankPhaseStruct_i);

                            curPairStats(loc_i) = pairPhaseStats;
                            if doWscc && pairInfo.sameGroup
%                                 shuffPairStats(loc_i) = getPairPhaseTuningStats(phases, ph_tc1(:), ph_tc2(:),  [], [],  blankShuffPairStats, pairInfo.pairIdx);
                                shuffPairStats(loc_i) = getPairPhaseTuningStats(phases, ph_tc1(:), ph_tc2(:),  ph1_jacks, ph2_jacks,  blankShuffPairStats, pairInfo.pairIdx);
                            end
                        
%                         end
                        
%                         if showExampleOfPhaseTuningCurves ...
                        
                    elseif isOEmeasure
                        
%%
                        ph_tc1_odd  = R1_stim_ph_oe(:,stim_idx, 1);
                        ph_tc1_even = R1_stim_ph_oe(:,stim_idx, 2);
                        
                        ph_tc2_odd  = R2_stim_ph_oe(:,stim_idx, 1);
                        ph_tc2_even = R2_stim_ph_oe(:,stim_idx, 2);  
                        
                        
                        if pairInfo.sameGroup && opt.getPhaseCmpJackStdErr_ptc
                            ph_jacks1_odd  = pairInfo.Cell1.R_stim_jackTrials(stim_idx,:,1);
                            ph_jacks1_even = pairInfo.Cell1.R_stim_jackTrials(stim_idx,:,2);

                            ph_jacks2_odd  = pairInfo.Cell2.R_stim_jackTrials(stim_idx,:,1);
                            ph_jacks2_even = pairInfo.Cell2.R_stim_jackTrials(stim_idx,:,2);                            
                        else
                            ph_jacks1_odd = [];
                            ph_jacks1_even = [];

                            ph_jacks2_odd = [];
                            ph_jacks2_even = [];
                        end
                        
                        %%
                        switch opt.phase_oe_mode
                            case 'same', 
                                [ph_tc1_a,    ph_tc2_a,     ph_tc1_b,    ph_tc2_b]    = deal(ph_tc1_odd,    ph_tc2_odd,      ph_tc1_even,    ph_tc2_even);
                                [ph_jacks1_a, ph_jacks2_a,  ph_jacks1_b, ph_jacks2_b] = deal(ph_jacks1_odd, ph_jacks2_odd,   ph_jacks1_even, ph_jacks2_even);
                            case 'diff', 
                                [ph_tc1_a,    ph_tc2_a,     ph_tc1_b,    ph_tc2_b]    = deal(ph_tc1_odd,    ph_tc2_even,     ph_tc1_even,    ph_tc2_odd);
                                [ph_jacks1_a, ph_jacks2_a,  ph_jacks1_b, ph_jacks2_b] = deal(ph_jacks1_odd, ph_jacks2_even,  ph_jacks1_even, ph_jacks2_odd);
                        end                                        
                        
%                         oneIsZeroOrUniform = any( [all(ph_tc1_a<=0), all(ph_tc2_a<=0), all(diff(ph_tc1_a)==0), all(diff(ph_tc2_a)==0), ...
%                                                    all(ph_tc1_b<=0), all(ph_tc2_b<=0), all(diff(ph_tc1_b)==0), all(diff(ph_tc2_b)==0) ]);
                                               
%                         invalid_frac = (fracR_val == 0) || isnan(fracR_val);
%                         if ~oneIsZeroOrUniform && ~invalid_frac                                               
                            
                            pairPhaseStats = getPairPhaseTuningStats_OE(phases, ph_tc1_a,    ph_tc1_b,     ph_tc2_a,    ph_tc2_b,  ...
                                                                                ph_jacks1_a, ph_jacks1_b,  ph_jacks2_a, ph_jacks2_b, blankPhaseStruct_i);
                           
                            curPairStats(loc_i) = pairPhaseStats;
                            if doWscc && pairInfo.sameGroup                           
                                shuffPairStats(loc_i) = getPairPhaseTuningStats_OE(phases, ph_tc1_a,    ph_tc1_b,     ph_tc2_a,    ph_tc2_b,  ...
                                                                                           ph_jacks1_a, ph_jacks1_b,  ph_jacks2_a, ph_jacks2_b,  blankShuffPairStats, pairInfo.pairIdx);
                            end                            
%                         end
                                               
                    end
                    
                    3;
                    curPairStats(loc_i).ptcOEcorr = min_ptcOEcorr;
                    
                    cc_vals = squeeze(curPairStats(loc_i).cc);
                    cc_val = mean(cc_vals);
                    dphi_vals = squeeze(curPairStats(loc_i).dphi);
                    dphi_val = mean(dphi_vals);
                    GCs = [pairInfo.Cell1.Gid, pairInfo.Cell1.cellId, pairInfo.Cell2.cellId];
                    if ~isempty(suspectGCs)
                        idx_suspect = findRows(GCs, suspectGCs);
                        isSuspectPair = any(dphi_vals == 6) && ~isempty(idx_suspect);
                    end
                        maxMinFracR_val = max(minFracR);
                    
%                         if  strcmp(locationName, 'maxMinFracR') && min_ptcOEcorr > 0.7 && length(phases) == 8 && fracR_val > 0.6
%                         if  strcmp(locationName, 'maxMinFracR') && min_ptcOEcorr < -0.5 && length(phases) == 8
%                     if  strcmp(locationName, 'maxMinFracR') && any(dphi_vals == 6) && length(phases) == 60 && isSuspectPair
%                     if  strcmp(locationName, 'maxMinFracR') && length(phases) == 60 && any(isnan(dphi_vals)) && ~isnan(cc_val)
                    if  strcmp(locationName, 'maxMinFracR') && length(phases) == 8 && any(dphi_vals==180) &&  (pairInfo.SCtype_pref == 2) && 0
%                     if  strcmp(locationName, 'maxMinFracR') && length(phases) == 4 && any(dphi_vals==0) && maxMinFracR_val > .5
%                             ibetween(dphi_val, [30, 40]) % ibetween(any(dphi_vals == 6) && 0% (cc_val > .8) && 0;
                        
                        %%
                        if 1
                            
                            showProfiles = 1;
                            if showProfiles
                                ori = pairInfo.oris;
                                spf = pairInfo.spfs;
                                
                                %                                 figure(1);clf; imagesc(ori, spf, mean(R1, 3)); hold on;
                                figure(1); clf; imageOSP(pairInfo.Cell1, 'mean:ph', 'OSP', 'nolabels');
                                colorbar;
                                title(sprintf('FracR = %.2f', pairPhaseStats.frac1ofMax))
                                %                                         title(' ');
                                hold on;
                                plot(ori_sp_idx(2), ori(ori_sp_idx(1)), 'wo', 'markersize', 10, 'linewidth', 3); colorbar;

                                figure(2); clf; imageOSP(pairInfo.Cell2, 'mean:ph', 'OSP', 'nolabels');
                                colorbar;
                                hold on;
                                %                                 figure(2);clf; imagesc(ori, spf, mean(R2, 3)); hold on;
                                plot(ori_sp_idx(2), ori(ori_sp_idx(1)), 'wo', 'markersize', 10, 'linewidth', 3); colorbar;
                                title(sprintf('FracR = %.2f', pairPhaseStats.frac2ofMax))
                                %                                 title(' ');
%                                 suspectProfileCC(idx_suspect) = corr(pairInfo.Cell1.R(:), pairInfo.Cell2.R(:));
                            end
                            nrm = @(x) x/max(x);
                            plotComparedTuningCurves = 1;
                            if isOEmeasure
                                %                                 if plotComparedTuningCurves
                                %%
                                ptc_1 = {ph_tc1_odd,  ph_tc1_even};
                                ptc_2 = {ph_tc2_even, ph_tc2_odd };
                                cc1 = corr(ph_tc1_odd(:), ph_tc1_even(:)); 
                                cc2 = corr(ph_tc2_odd(:), ph_tc2_even(:));
                                ss = sprintf('cc1 = %.2f. cc2 = %.2f', cc1, cc2);
                                h_axx = zeros(2,2);
                                sortedCCs= cell(1,2);
                                for ii = 1:2
                                    figure(30+ii); clf;
                                    p1 = ptc_1{ii}; p2 = ptc_2{ii};
                                    h_axx(ii,1) = subplot(2,1,1); h_p1 = plot(phases, nrm(p1), 'bo-', phases, nrm(p2), 'gs-');
                                    ccs_fromShifts = getCCsFromShifts(p1, p2);
                                    title(sprintf('dphi = %.1f', dphi_vals(ii)));                                    

                                    [mx_cc, ind_max] = max(ccs_fromShifts);
                                    h_axx(ii,2)= subplot(2,1,2); plot(phases, ccs_fromShifts, 'k.-'); hold on;
                                    plot(phases(ind_max), mx_cc, 'ro');
                                    if isnan(dphi_vals(ii))
                                        assert(nnz(ccs_fromShifts > mx_cc *.9999) > 1);
                                    end
                                    xlabel(ss);
                                    sortedCCs{ii} = sort(ccs_fromShifts, 'descend');                                    
                                    set(h_axx(ii,:), 'xtick', [0:90:360], 'xlim', [0 360]);
                                end
                                
                                3;
                            else
                                %%
                                figure(5); clf;
                                h_axx(1) = subplot(2,1,1); h_p1 = plot(phases, nrm(ph_tc1(:)), 'bo-', phases, nrm(ph_tc2(:)), 'gs-');
                                ccs_fromShifts = getCCsFromShifts(ph_tc1, ph_tc2);
                                [mx_cc, ind_max] = max(ccs_fromShifts); 
                                h_axx(2)= subplot(2,1,2); plot(phases, ccs_fromShifts, 'k.-'); hold on;
                                plot(phases(ind_max), mx_cc, 'ro');
                                assert(nnz(ccs_fromShifts > mx_cc *.99) == 1);
                                set(h_axx, 'xtick', [0:90:360], 'xlim', [0 360]);
                                
                                3;
                            end
                            
                            
                            
                            
                            3;
                            %                                 else %
                            %                                     figure(3); plot(phases, ph_tc1_odd, 'bo-', phases, ph_tc1_even, 'gs-');
                            %                                     title(sprintf('cell 1 (odd/even)', dphi_vals(1)))
                            %                                     figure(4); plot(phases, ph_tc2_odd, 'bo-', phases, ph_tc2_even, 'gs-');
                            %                                     title(sprintf('ptc 2 (odd/even)', dphi_vals(2)))
                            %                                 end
                            %                                 if plotComparedTuningCurves
                            %                                 3;
                            %                                 fig_ids = [3,4];
                            %                                 ph_cells_odd = {ph_tc1_odd, ph_tc2_odd};
                            %                                 ph_cells_even = {ph_tc1_even, ph_tc2_even};
                            %                                 cols = {'b', [0 .7 0]};
                            %                                 for i = 1:2
                            %                                     wrp = @(x) [x; x(1)];
                            %                                     phases_ext = [phases, 360];
                            %                                     figure(fig_ids(i)); clf;
                            %                                     ph_odd = ph_cells_odd{i};
                            %                                     ph_even = ph_cells_even{i};
                            %                                     M_spc = [.0 .05 0.15 ];
                            %                                     h_ax(1) = subplotGap(2,1,1,1, M_spc); plot(phases_ext, wrp(ph_odd), 'o-', 'color', cols{i}, 'linewidth', 2);
                            %                                     cc_oe1 = corr(ph_odd(:), ph_even(:));
                            %     %                                 title(sprintf('cc (odd vs even) = %.1f', cc_oe1))
                            %                                     title('Odd Trials')
                            %                                     h_ax(2) = subplotGap(2,1,2,1, M_spc); plot(phases_ext, wrp(ph_even), 'o-', 'color', cols{i}, 'linewidth', 2);
                            %                                     title('Even Trials')
                            %                                     set(h_ax, 'xtick', [0:90:360], 'xlim', [0 360]);
                            %
                            %                                     hh(i) = annotation('textbox', [0 0, 1, .1], 'string', ...
                            %                                         sprintf('cc (odd vs even) = %.2f', cc_oe1), ...
                            %                                         'horiz', 'center', 'verticalAlignment', 'top', 'linestyle', 'none');
                            %                                 end
                            %
                            
                            
                            3;
                        end
                        %%
                        imCC = 0;
                        if imCC
                            
                            figure(140);
                            clf;
                            [h_axx, h_imm] = imageOSP(pairInfo.Cell2, 'mean:ph', 'OSP', 'nolabels');
                            colorbar;
                            set(h_imm, 'cdata', pairInfo.phaseTC_cc2);
                            set(h_axx, 'clim', [-1 1]);
                            
                            %%
                            indx = find(pairInfo.phaseTC_cc2 > .9);
                            indx = indx(3);
                            
                            ph_odd  = R2_stim_ph_oe(:,indx, 1);
                            ph_even = R2_stim_ph_oe(:,indx, 2);
                            
                            wrp = @(x) [x; x(1)];
                            phases_ext = [phases, 360];
                            figure(115); clf;
                            M_spc = [.0 .05 0.15 ];
                            h_ax(1) = subplotGap(2,1,1,1, M_spc); plot(phases_ext, wrp(ph_odd), 'o-', 'color', cols{1}, 'linewidth', 2);
                            cc_oe1 = corr(ph_odd(:), ph_even(:));
                            %                                 title(sprintf('cc (odd vs even) = %.1f', cc_oe1))
                            title('Odd Trials')
                            h_ax(2) = subplotGap(2,1,2,1, M_spc); plot(phases_ext, wrp(ph_even), 'o-', 'color', cols{1}, 'linewidth', 2);
                            title('Even Trials')
                            set(h_ax, 'xtick', [0:90:360], 'xlim', [0 360]);
                            
                            hh(i) = annotation('textbox', [0 0, 1, .1], 'string', ...
                                sprintf('cc (odd vs even) = %.2f', cc_oe1), ...
                                'horiz', 'center', 'verticalAlignment', 'top', 'linestyle', 'none');
                            
                        end
                        
                        %                             count6 = count6 +1;
                        3;
                    end
                         3;                    
                    
                    show = 0;
                    if show
                        if strcmp(locationName, 'maxMinFracR')
                            figure(50); clf;
                            subplotGap(1,3,1, [], [.05 0, .05]);
                            plot(ext(phases), wrp(ph_tc1/max(ph_tc1)), 'bo-', ...
                                 ext(phases), wrp(ph_tc2/max(ph_tc2)), 'gs-' );
                             set(gca, 'xtick', ext(phases), 'xlim', lims(ext(phases)))
                            title(sprintf('%s : cc = %.2f', locationName, pairPhaseStats.cc));

                        end

                        if strcmp(locationName, 'p75MinFracR')
                            figure(50); 
                            subplotGap(1,3,2, [], [.05 0, .05]);
                            plot(ext(phases), wrp(ph_tc1/max(ph_tc1)), 'bo-', ...
                                 ext(phases), wrp(ph_tc2/max(ph_tc2)), 'gs-' );
                             set(gca, 'xtick', ext(phases), 'xlim', lims(ext(phases)))
                            title(sprintf('%s : cc = %.2f', locationName, pairPhaseStats.cc));
                            3;
                        end

                        if strcmp(locationName, 'p50MinFracR')
                            figure(50); 
                            if all(ph_tc2 == 0)
                                3;
                            end
                            subplotGap(1,3,3, [], [.05 0, .05]);
                            plot(ext(phases), wrp(ph_tc1/max(ph_tc1)), 'bo-', ...
                                 ext(phases), wrp(ph_tc2/max(ph_tc2)), 'gs-' );
                             set(gca, 'xtick', ext(phases), 'xlim', lims(ext(phases)))
                            title(sprintf('%s : cc = %.2f', locationName, pairPhaseStats.cc));
                            3;
                        end
                    end                                        
                    
                                        
                case 'wgtedSum', % Weighted sum (summing over all ori,sp, and weighting by firing rates)                    
                    
                    
                    if ~exist('R1_stim_ph', 'var') || ~exist('sumPhs1_av', 'var')                     
%                         R1_stim_ph = reshape(R1, [nOri*nSpf, nPh])';
%                         R2_stim_ph = reshape(R2, [nOri*nSpf, nPh])';
                        sumPhs1_av = sum(R1_stim_ph,1); normSumRates1 = sumPhs1_av/max(sumPhs1_av);
                        sumPhs2_av = sum(R2_stim_ph,1); normSumRates2 = sumPhs2_av/max(sumPhs2_av);
                        
                        nRespPh1_av = sum(R1_stim_ph>0,1);
                        nRespPh2_av = sum(R2_stim_ph>0,1);                        
                    end
                    
                    if isOEmeasure
                        nRespPh1_calc = nRespPh1_oe;
                        nRespPh2_calc = nRespPh2_oe;                        
                    else
                        nRespPh1_calc = nRespPh1_av;
                        nRespPh2_calc = nRespPh2_av;                        
                    end
                    
%                         'wgtedSum_all_oe'
                    switch locations{loc_i}
                        case {'wgtedSum_all',  'wgtedSum_all_oe'}
                                                respAllPhs1 = nRespPh1_calc>0;
                                                respAllPhs2 = nRespPh2_calc>0;
                        case 'wgtedSum_1phase', respAllPhs1 = nRespPh1_calc == 1;
                                                respAllPhs2 = nRespPh2_calc == 1;
                        case 'wgtedSum_2phase', respAllPhs1 = nRespPh1_calc == 2;
                                                respAllPhs2 = nRespPh2_calc == 2;
                        case 'wgtedSum_4phase', respAllPhs1 = nRespPh1_calc == 4;
                                                respAllPhs2 = nRespPh2_calc == 4;
                        case 'wgtedSum_allphase',respAllPhs1 = nRespPh1_calc == nPh;
                                                 respAllPhs2 = nRespPh2_calc == nPh;
                    end

                    idx_bothResponded = find( respAllPhs1 & respAllPhs2 );
                    nStimBothResponded = length(idx_bothResponded);

                    minNormSumRates = min(normSumRates1(idx_bothResponded), normSumRates2(idx_bothResponded));                        
                    prodRates = sumPhs1(idx_bothResponded) .* sumPhs2(idx_bothResponded);

                    switch weightingBy
                        case 'product', wgts = prodRates;
                        case 'minNormRate',  wgts = minNormSumRates;
                    end
                    switch selectingBy
                        case 'product', responseMeasure = prodRates;
                        case 'minNormRate',  responseMeasure = minNormSumRates;
                    end
                                                
                        
%                         [bestPairResponse, idx_bestPairResponses] = sort(responseMeasure, 'descend');
%                         bestPairResponse1 = max(responseMeasure);
%                         [bestProdResponse, idx_bestProdResponses] = sort(prodRates, 'descend');
                    if nStimBothResponded > 0

                        if isOEmeasure
                            R_cmp1_odd  = R1_stim_ph_oe(:,idx_bothResponded, 1);
                            R_cmp1_even = R1_stim_ph_oe(:,idx_bothResponded, 2);

                            R_cmp2_odd  = R2_stim_ph_oe(:,idx_bothResponded, 1);
                            R_cmp2_even = R2_stim_ph_oe(:,idx_bothResponded, 2);                                                                

                            switch opt.phase_oe_mode
                                case 'same', [R_cmp1_a, R_cmp2_a,  R_cmp1_b, R_cmp2_b] = deal(R_cmp1_odd, R_cmp2_odd,   R_cmp1_even, R_cmp2_even);
                                case 'diff', [R_cmp1_a, R_cmp2_a,  R_cmp1_b, R_cmp2_b] = deal(R_cmp1_odd, R_cmp2_even,  R_cmp1_even, R_cmp2_odd );
                            end
                            
                            pairPhaseStruct_all_a = getPairPhaseTuningStats(phases, R_cmp1_a, R_cmp2_a, blankPhaseStruct);
                            pairPhaseStruct_all_b = getPairPhaseTuningStats(phases, R_cmp1_b, R_cmp2_b, blankPhaseStruct);
                            pairPhaseStruct_all = blankPhaseStruct;
                            fn = fieldnames(pairPhaseStruct_all);
                            for fld_i = 1:length(fn)
                                pairPhaseStruct_all.(fn{fld_i}) = (pairPhaseStruct_all_a.(fn{fld_i}) + pairPhaseStruct_all_b.(fn{fld_i}))/2;
                            end
                        else
                            R_cmp1 = R1_stim_ph(:,idx_bothResponded);
                            R_cmp2 = R2_stim_ph(:,idx_bothResponded);
                            pairPhaseStruct_all = getPairPhaseTuningStats(phases, R_cmp1, R_cmp2, blankPhaseStruct);
                        end                                                                                        


                    else
                        pairPhaseStruct_all = blankPhaseStruct;
                    end

                    idx_use = 1:nStimBothResponded;


                    if recordAllCCs && strcmp(locations{loc_i}, location_record) % only do once per location.                            
                        [~, idx_sort] = sort(prodRates, 'descend');
                        ccs = nan(1,nStimMax);
                        n = min(nStimMax, length(idx_sort));
                        ccs(1:n) = pairPhaseStruct_all.cc(idx_sort(1:n));                        

                        if pairInfo.sameLocation
                            assert(size(allWCCs, 1) >= allWCC_idx);
                            allWCCs(allWCC_idx,:) = ccs;
                            allWCC_idx = allWCC_idx+1;                                                                
                        else
                            assert(size(allBCCs, 1) >= allBCC_idx);
                            allBCCs(allBCC_idx,:) = ccs;
                            allBCC_idx = allBCC_idx+1;
                        end
                    end
                                                    
%                     switch locations{loc_i}
%                         case 'wgtedSum_all',     idx_use = 1:nStimBothResponded;
%                         case 'wgtedSum_top10',   idx_use = idx_bestPairResponses(1:min(10, nStimBothResponded) );
%                         case 'wgtedSum_top20',   idx_use = idx_bestPairResponses(1:min(20, nStimBothResponded) );
%                         case 'wgtedSum_above90', idx_use = find(responseMeasure > bestPairResponse1*.9);
%                         case 'wgtedSum_above50', idx_use = find(responseMeasure > bestPairResponse1*.5);                    
                    
                    pairPhaseStruct = blankPhaseStruct;                        
    %                     f = @(x) wgt_nanmean(x, wgts_use);
                    useWeightByFiringRate = 0;
                    if useWeightByFiringRate
                        wgts_use = wgts(idx_use);
                        if doCC,   pairPhaseStruct.cc    = wgt_nanmean(pairPhaseStruct_all.cc(idx_use), wgts_use);  end
                        if doRho,  pairPhaseStruct.rho   = wgt_nanmean(pairPhaseStruct_all.rho(idx_use), wgts_use);  end
                        if doTau,  pairPhaseStruct.tau   = wgt_nanmean(pairPhaseStruct_all.tau(idx_use), wgts_use);  end
                        if doDPhi, pairPhaseStruct.dphi  = wgt_nanmean(pairPhaseStruct_all.dphi(idx_use), wgts_use);  end
                        if doDF1,  pairPhaseStruct.dF1   = wgt_nanmean(pairPhaseStruct_all.dF1(idx_use), wgts_use);  end
                        if doF1oDC && doF1oDCforMultipleStim
                            pairPhaseStruct.F1oDC1= wgt_nanmean(pairPhaseStruct_all.F1oDC1(idx_use), wgts_use);
                            pairPhaseStruct.F1oDC2= wgt_nanmean(pairPhaseStruct_all.F1oDC2(idx_use), wgts_use);
                        end
                    else
                        if doCC,   pairPhaseStruct.cc    = nanmean(pairPhaseStruct_all.cc(idx_use));  end
                        if doRho,  pairPhaseStruct.rho   = nanmean(pairPhaseStruct_all.rho(idx_use));  end
                        if doTau,  pairPhaseStruct.tau   = nanmean(pairPhaseStruct_all.tau(idx_use));  end
                        if doDPhi, pairPhaseStruct.dphi  = nanmean(pairPhaseStruct_all.dphi(idx_use));  end
                        if doDF1,  pairPhaseStruct.dF1   = nanmean(pairPhaseStruct_all.dF1(idx_use));  end
                        if doF1oDC && doF1oDCforMultipleStim 
                            pairPhaseStruct.F1oDC1= nanmean(pairPhaseStruct_all.F1oDC1(idx_use));  
                            pairPhaseStruct.F1oDC2= nanmean(pairPhaseStruct_all.F1oDC2(idx_use));  
                        end
                    end
                    if doDPhi && (pairPhaseStruct.dphi > 180 )   % due to numerical rounding error
                        pairPhaseStruct.dphi = 180;
                    end                    
                    if doDF1 && (pairPhaseStruct.dF1 > 180 )
                        pairPhaseStruct.dF1 = 180;
                    end
                    pairPhaseStruct.lowerSumPhs = nStimBothResponded;

                    if ~pairInfo.sameGroup && strcmp(locations{loc_i}, 'wgtedSum_all')
%                        pairPhaseStruct.cc_atMaxDphi = pearsonR(mean(R1_stim_ph, 1)', mean(R2_stim_ph, 1)');                        
                    end                    
%                     if strcmp(locations{loc_i}, 'wgtedSum_all') && isfield(opts, 'STA_cc') && (abs(STA_cc) < .18) && (pairPhaseStruct.cc > .27)
%                         3;
%                     end
                    
                    
                    curPairStats(loc_i)  = pairPhaseStruct;                     
                             
                case 'mtx x mtx'
                    f = 0.1;
                    OSP_sigInds = (R1(:) > f * mean(R1(:)) )  |  (R2(:) > f * mean(R2(:)) );
                    [dt, cc, rho, tau, tmpxx, tmpyy, cc_p, rho_p, tau_p ] = getPhaseTuningStats(phases, R1(OSP_sigInds), R2(OSP_sigInds));
                    meanPhs1 = squeeze( mean(mean(R1, 1), 2) );
                    meanPhs2 = squeeze( mean(mean(R2, 1), 2) );
                    dphi = deltaPhi( phases, meanPhs1, meanPhs2, 'max' );
                    dF1  = deltaPhi( phases, meanPhs1, meanPhs2, 'angle' );
                    
                    curPairStats(loc_i) = struct('dot', dt, 'cc', cc, 'cc_p', cc_p, 'rho', rho, 'rho_p', rho_p, 'tau', tau, 'tau_p', tau_p, 'dphi', dphi, 'dF1', dF1, 'F1oDC1', nan, 'F1oDC2', nan, 'lowerSumPhs', nan);
                    
                otherwise
                    error('Unknown "location" parameter');
            end
        end
                           
                                      
            if showWorking
                figure(900);
                cc = curPairStats(1).cc;
                binId = find( binsL < cc & cc < binsR );
                if ~isempty(binId) && upTo(binId) < nWorkRows;
                    subplot(nWorkRows, nWorkCols, (upTo(binId)-1)*nWorkCols + binId);
                    
                    [ori_i, sp_i] = elements( oriSpInds.maxR1xR2 );            
                    ph_tc1 = squeeze( R1( ori_i, sp_i, : ) );
                    ph_tc2 = squeeze( R2( ori_i, sp_i, : ) );                    
                    
                    plot(phases, ph_tc1, 'ro-', phases, ph_tc2, 'go-');
                    text(0,0, num2str(cc), 'verticalalignment', 'bottom');
                    xlim([0 phases(end)])
                    set(gca, 'xtick', []);
                    upTo(binId) = upTo(binId)+1;
                end
            end                
            
            %{
%             if all(upTo == nWorkRows);
%                 3;
%             end
%                 ind_i = mxPhsOutputId;
% %                 [ori_i, sp_i] = elements( outputOriSpInds(ind_i,:) );
%                 ph_tc1 = squeeze( R1( ori_i, sp_i, : ) );
%                 ph_tc2 = squeeze( R2( ori_i, sp_i, : ) );
%                 [dphi, idx_maxPhi1, idx_maxPhi2] = deltaPhi(phs, ph_tc1, ph_tc2);
%                 if ~isnan(idx_maxPhi1)
%                     thesePhs = phases([idx_maxPhi1, idx_maxPhi2]);
%                     mxPhs(mxPhsI,:) = thesePhs;
%                     mxPhsI = mxPhsI + 1;
% 
%                 end
%             end
        %}

        4;
    end


    function [curPairRFStats, shuffPairRFStats, stats] = getReceptiveFieldCorrelationsFor2Cells(Cell1, Cell2, curPairRFStats, shuffPairRFStats, blankSingleRFStat)
        
        sameGroup = Cell1.Gid == Cell2.Gid;
%         shuffPairStats = struct;

        min_rqr_oe_STA = nan;
        
        if doSTAcc 
            
%             [nSTA1, nSTA2] = deal(nan);
            cell1_hasSTA = hasSTA(Cell1);
            cell2_hasSTA = hasSTA(Cell2);
            
            
            if (cell1_hasSTA && cell2_hasSTA) 
                min_rqr_oe_STA = min(Cell1.STAdata.rsqr_oe, Cell2.STAdata.rsqr_oe);                        

                cmp_stas = (length(Cell1.STAdata.L) == length(Cell2.STAdata.L)) ...   (nSTA1 == nSTA2)
                    && (min_rqr_oe_STA > opt.MID_minReq_rsqr_oe);

                if cmp_stas
                    [curPairRFStats.STA_cc, curPairRFStats.STA_cc_jackStd] = rfCC(Cell1.STAdata, Cell2.STAdata, 'STA', 'STA', blankSingleRFStat, opt);        
                end


                if doWscc && doRandMIDs && sameGroup 
                    if cmp_stas       
                        STA_cc_rand = randFlipSignHalf(curPairRFStats.STA_cc, opt.nPhaseShuffles);
                        shuffPairRFStats.STA_cc = STA_cc_rand;
                    end
                    
                end
            end
            shuffPairRFStats.STA_cc_jackStd = curPairRFStats.STA_cc_jackStd ( ones(1, opt.nPhaseShuffles, opt.n3) );            
        end

        
        [nMID1, nMID2] = deal(nan);
        cell1_hasMID = hasMID(Cell1);
        cell2_hasMID = hasMID(Cell2);

        
        if cell1_hasMID                        
            MIDdata1 = Cell1.MIDdata;
            nMID1 = MIDdata1.L;                        
        end
        if cell2_hasMID                        
            MIDdata2 = Cell2.MIDdata;
            nMID2 = MIDdata2.L;
        end

        if (cell1_hasMID && cell2_hasMID) && (nMID1 == nMID2)                                 
            
            jackCC_1 = min(MIDdata1.MID_odd.jackCC, MIDdata1.MID_even.jackCC);                        
            rsqr_fit_1 = min(MIDdata1.MID_odd.rsqr_fit, MIDdata1.MID_even.rsqr_fit);
            rqr_oe_1 = MIDdata1.rsqr_oe;

            jackCC_2 = min(MIDdata2.MID_odd.jackCC, MIDdata2.MID_even.jackCC);
            rsqr_fit_2 = min(MIDdata2.MID_odd.rsqr_fit, MIDdata2.MID_even.rsqr_fit);
            rqr_oe_2 = MIDdata2.rsqr_oe;

            min_jackCC = min(jackCC_1, jackCC_2);
            max_jackCC = max(jackCC_1, jackCC_2);

            min_rsqr_fit = min(rsqr_fit_1, rsqr_fit_2);
            max_rsqr_fit = max(rsqr_fit_1, rsqr_fit_2);

            min_rsqr_oe = min(rqr_oe_1, rqr_oe_2);
            max_rsqr_oe = max(rqr_oe_1, rqr_oe_2);

        else
            [min_jackCC, max_jackCC, min_rsqr_fit, max_rsqr_fit, min_rsqr_oe, max_rsqr_oe] = deal(nan);
        end
        
        
        if doMIDcc
            cmp_mids = (nMID1 == nMID2) && (min_rsqr_oe > opt.MID_minReq_rsqr_oe); % && (rsqr_oe < opt.MID_minReq_rsqr_oe);

            if cmp_mids % (nMID1 == nMID2) && (min_jackCC > opt.MID_minReq_jackCC)                            
                [curPairRFStats.MID_cc, curPairRFStats.MID_cc_jackStd] = rfCC(MIDdata1, MIDdata2, 'MID', 'MID', blankSingleRFStat, opt); 

%                 [MID_cc1, MID_cc2] = deal(MID_cc_S.cc_a, MID_cc_S.cc_b);
%                 pearsonR(MIDdata1.MID, MIDdata2.MID);
% 
%                 maxF1oDC_here = max(Cell1.F1oDC_maxR_avP, Cell1.F1oDC_maxR_avP);
%                 if (MID_cc > -0.1) && (maxF1oDC_here < 1) && (min_rsqr_fit > .5)
                if 0 && (min_jackCC > 0.4) && (max_rsqr_fit < .2)
                    %%
                    mid1 = MIDdata1.MID; mid2 = MIDdata2.MID;
                    osp1 = Cell1.R; osp2 = Cell2.R;
                    L1 = max(abs(mid1(:)));
                    L2 = max(abs(mid2(:)));
                    figure(91); clf;
                    subplotGap(2,2, 1,1); imagesc(mid1); axis square; set(gca, 'xtick', [], 'ytick', []); caxis([-L1, L1]); title(sprintf('%d : %d', Gid1, cellId1));
                    subplotGap(2,2, 1,2); imagesc(mid2); axis square; set(gca, 'xtick', [], 'ytick', []); caxis([-L2, L2]); title(sprintf('%d : %d', Gid2, cellId2));
                    subplotGap(2,2, 2,1); imageOSP(osp1, 'mean:ph', 'OSP', 'noLabels'); set(gca, 'xtick', [], 'ytick', []);  colorbar;
                    subplotGap(2,2, 2,2); imageOSP(osp2, 'mean:ph', 'OSP', 'noLabels'); set(gca, 'xtick', [], 'ytick', []);  colorbar;
%                                 subplotGap(1,2, 1,1); imagesc(mid1); axis square; set(gca, 'xtick', [], 'ytick', []); caxis([-L1, L1]); title(sprintf('%d : %d', Gid1, cellId1));
%                                 subplotGap(1,2, 1,2); imagesc(mid2); axis square; set(gca, 'xtick', [], 'ytick', []); caxis([-L2, L2]); title(sprintf('%d : %d', Gid2, cellId2));
                    xlabel(sprintf('cc = %.2f', MID_cc))
                    3;
                end
                
            end


            if doWscc && doRandMIDs && sameGroup 
                if cmp_mids
                    MID_cc_rand = randFlipSignHalf(curPairRFStats.MID_cc, opt.nPhaseShuffles);
                    shuffPairRFStats.MID_cc = MID_cc_rand;
                end
                
            end
            
            shuffPairRFStats.MID_cc_jackStd = curPairRFStats.MID_cc_jackStd ( ones(1, opt.nPhaseShuffles, opt.n3) );

        end

        [min_frac_overlap, max_frac_overlap, min_ovlp] = deal(nan);      
%         MID_ovlp_cc_rand = []; 
        if doMID_ovlp_cc
            cmp_mids = (nMID1 == nMID2) && (min_rsqr_oe > opt.MID_minReq_rsqr_oe) ...
                 && (min_rsqr_fit > opt.MID_minReq_rsqr_fit); %(min_jackCC > opt.MID_minReq_jackCC); % && (rsqr_oe < opt.MID_minReq_rsqr_oe);
            if cmp_mids
%                 mid_overlap = ~( MIDdata1.MID_select & MIDdata2.MID_select );
%                 if nnz(mid_overlap) >= numel(MIDdata1.MID)/200  % overlap > 0.5% of total # pixels
%                     c1 = ;
%                     c2 = pearsonR(MIDdata1.MID, MIDdata2.MID);
%                     fprintf('%.2f, %.2f\n', c1, c2);
%                     curPairStats.MID_ovlp_cc = pearsonR(MIDdata1.MID(mid_overlap), MIDdata2.MID(mid_overlap));
                    [curPairRFStats.MID_ovlp_cc, curPairRFStats.MID_ovlp_cc_jackStd, MID_ovlp_cc_S] = rfCC(MIDdata1, MIDdata2, 'MID', 'MID_overlap', blankSingleRFStat, opt);
                    [mid_overlap1, mid_overlap2, min_ovlp1, min_ovlp2] = deal(...
                        MID_ovlp_cc_S.mid_overlap1, MID_ovlp_cc_S.mid_overlap2, MID_ovlp_cc_S.min_ovlp1, MID_ovlp_cc_S.min_ovlp2);

                  frac_overlap1 = nnz(mid_overlap1)/numel(mid_overlap1);
                  frac_overlap2 = nnz(mid_overlap2)/numel(mid_overlap2);
                  min_frac_overlap = min(frac_overlap1, frac_overlap2);
                  max_frac_overlap = max(frac_overlap1, frac_overlap2);           
                  min_ovlp = min(min_ovlp1, min_ovlp2);
            end


            if doWscc && sameGroup && cmp_mids
                MID_ovlp_cc_rand = randFlipSignHalf(curPairRFStats.MID_ovlp_cc, opt.nPhaseShuffles); 
                shuffPairRFStats.MID_ovlp_cc = MID_ovlp_cc_rand;
            end
            
            shuffPairRFStats.MID_ovlp_cc_jackStd = curPairRFStats.MID_ovlp_cc_jackStd ( ones(1, opt.nPhaseShuffles, opt.n3) );
        end                    
        stats = struct('min_rqr_oe_STA', min_rqr_oe_STA, ...
                       'min_jackCC', min_jackCC, 'max_jackCC', max_jackCC, ...
                       'min_rsqr_fit', min_rsqr_fit, 'max_rsqr_fit', max_rsqr_fit, ...
                       'min_rsqr_oe', min_rsqr_oe, 'max_rsqr_oe', max_rsqr_oe, ...
                       'min_frac_overlap', min_frac_overlap, 'max_frac_overlap', max_frac_overlap, ...
                       'min_ovlp', min_ovlp);

        
        % MID - GABOR fits:           
        if doMID_fit_cc 
            cmp_mid_fits = (nMID1 == nMID2) && ...
                    (min_rsqr_oe > opt.MID_minReq_rsqr_oe) && (min_rsqr_fit > opt.MID_minReq_rsqr_fit);
            if cmp_mid_fits
                [curPairRFStats.MID_fit_cc, curPairRFStats.MID_fit_cc_jackStd] = rfCC(MIDdata1, MIDdata2, 'MID', 'MID_fit', blankSingleRFStat, opt);
%                 [MID_fit_cc1, MID_fit_cc2] = deal(MID_fit_cc_S.cc_a, MID_fit_cc_S.cc_b);
                
                if mean(abs(curPairRFStats.MID_fit_cc)) > .5 && 0
                    %%
                    figure(53);  clf;
                    subplot(2,2,1); imagesc(MIDdata1.MID_odd.MID); axis square; set(gca, 'xtick', [], 'ytick', []);  ylabel(sprintf('[%d, %d] (odd)', Gid1, cellId1))
                    subplot(2,2,2); imagesc(MIDdata1.MID_even.MID); axis square; set(gca, 'xtick', [], 'ytick', []); 
                    subplot(2,2,3); imagesc(MIDdata2.MID_even.MID); axis square; set(gca, 'xtick', [], 'ytick', []); ylabel(sprintf('[%d, %d] (even)', Gid2, cellId2))
                    xlabel(sprintf('cc1: %.2f', MID_cc1));
                    subplot(2,2,4); imagesc(MIDdata2.MID_odd.MID); axis square; set(gca, 'xtick', [], 'ytick', []);                                
                    xlabel({sprintf('cc2: %.2f', MID_cc2), sprintf('mean: %.2f', curPairRFStats.MID_cc) });
                    3;

                    figure(54);  clf;
                    MID1_odd  = MIDdata1.MID_odd.MID;  MID1_odd( ~mid_overlap1) = 0;
                    MID2_even = MIDdata2.MID_even.MID; MID2_even(~mid_overlap1) = 0;
                    MID1_even = MIDdata1.MID_even.MID; MID1_even(~mid_overlap2) = 0;
                    MID2_odd  = MIDdata2.MID_odd.MID;  MID2_odd( ~mid_overlap2) = 0;

                    subplot(2,2,1); imagesc(MID1_odd); axis square; set(gca, 'xtick', [], 'ytick', []);  ylabel(sprintf('[%d, %d] (odd)', Gid1, cellId1))
                    subplot(2,2,2); imagesc(MID1_even); axis square; set(gca, 'xtick', [], 'ytick', []); 
                    subplot(2,2,3); imagesc(MID2_even); axis square; set(gca, 'xtick', [], 'ytick', []); ylabel(sprintf('[%d, %d] (even)', Gid2, cellId2))
                    xlabel(sprintf('cc1: %.2f', MID_ovlp_cc1));
                    subplot(2,2,4); imagesc(MID2_odd); axis square; set(gca, 'xtick', [], 'ytick', []);                                
                    xlabel({sprintf('cc2: %.2f', MID_ovlp_cc2), sprintf('mean: %.2f', MID_ovlp_cc) });
                    3;


                    figure(55);  clf;
                    subplot(2,2,1); imagesc(MIDdata1.MID_odd.MID_fit); axis square; set(gca, 'xtick', [], 'ytick', []);  ylabel(sprintf('[%d, %d] (odd)', Gid1, cellId1))
                    subplot(2,2,2); imagesc(MIDdata1.MID_even.MID_fit); axis square; set(gca, 'xtick', [], 'ytick', []); 
                    subplot(2,2,3); imagesc(MIDdata2.MID_even.MID_fit); axis square; set(gca, 'xtick', [], 'ytick', []); ylabel(sprintf('[%d, %d] (even)', Gid2, cellId2))
                    xlabel(sprintf('cc1: %.2f', MID_fit_cc1));
                    subplot(2,2,4); imagesc(MIDdata2.MID_odd.MID_fit); axis square; set(gca, 'xtick', [], 'ytick', []);                                
                    xlabel({sprintf('cc2: %.2f', MID_fit_cc2), sprintf('mean: %.2f', curPairRFStats.MID_fit_cc) });                                
                    3;

                end

                
            end
            
            if doWscc && sameGroup
                if cmp_mid_fits && doRandMID_fits
                    MID_fit_cc_rand = MID_fits_randPhase_for2Cells(MIDdata1, MIDdata2, Cell1.Gid, opt);
                    shuffPairRFStats.MID_fit_cc = MID_fit_cc_rand;
                end
                
            end            
            shuffPairRFStats.MID_fit_cc_jackStd = curPairRFStats.MID_fit_cc_jackStd( ones(1, opt.nPhaseShuffles, opt.n3) );

        end                    



        if doDRelPhase
            bothMID_fits_ok = cell1_hasMID && cell2_hasMID && ...
                    (min_rsqr_oe > opt.MID_minReq_rsqr_oe) && (min_rsqr_fit > opt.MID_minReq_rsqr_fit);

            if bothMID_fits_ok         
                curPairRFStats.dPh_rel = dRelativePhase(MIDdata1, MIDdata2, opt);
%                             rel_ph1 = rad2deg(MIDdata1.MID.p(8));
%                             rel_ph2 = rad2deg(MIDdata2.MID.p(8));
%                             dPh_rel = circDist(rel_ph1, rel_ph2, 360);                
            end
%             curPairStats.dPh_rel = dPh_rel;

            if doWscc && sameGroup 
                if bothMID_fits_ok && ~any(isnan(curPairRFStats.dPh_rel))
                    dPh_rel_rand = dRelativePhase_rand(opt.nPhaseShuffles, opt);
                    shuffPairRFStats.dPh_rel = dPh_rel_rand;
                    
                end
            end                        
           

        end
    end



    function doCalculationsForPairs(pairIdxs, pairLabel)

        nCurPairs = length(pairIdxs);                
        pairs = ind2subV([nUnits, nUnits], pairIdxs);
        fprintf(['\nCalculating ' pairLabel ' data: (total of ' num2str(nCurPairs) ' pairs)']);

        progressBar('init-', nCurPairs, 60);
        switch cmpType, 
            case 'phase'                                                
                measure_flds = measures;
                assert(~calcPvalues);
                
                F1oDC_flds = iff(doF1oDC, {'F1oDC1', 'F1oDC2'}, {});
                ccMaxDphi_flds = iff(doDPhi || 1, {'cc_atMaxDphi'}, {});
                lowerSumPhs_flds = iff(calcLowerSumPhs, {'lowerSumPhs'}, {});
                frac1ofMax_flds = {'frac1ofMax', 'frac2ofMax'};
                ptcOEcorr_flds = {'ptcOEcorr'};
                jackStd_flds = {};
                if opt.getPhaseCmpJackStdErr_ptc
                    jackStd_flds = cellfun(@(s) [s '_jackStd'], measure_flds, 'un', 0);
                end
                
                allFields = [measure_flds, F1oDC_flds, frac1ofMax_flds, ccMaxDphi_flds, lowerSumPhs_flds, ptcOEcorr_flds, jackStd_flds];
                blankPhaseTuningStats = structFromFieldNames(allFields, nan);                     
                for mj = 1:length(measures)
                    blankPhaseTuningShuff.(measures{mj}) = nan(1, opt.nPhaseShuffles, 1 ); % not opt.n3
                    if opt.getPhaseCmpJackStdErr_ptc
                        blankPhaseTuningShuff.([measures{mj} '_jackStd']) = nan(1, opt.nPhaseShuffles, 1 ); % not opt.n3
                    end
                end
                

                blankRF_fields = RF_Measures;
                if opt.getPhaseCmpJackStdErr_rf
                    blankRF_fields = [blankRF_fields,   cellfun(@(s) [s '_jackStd'], RF_Measures, 'un', 0)];
                end
                blankRFStats       = structFromFieldNames(blankRF_fields, nan);
                blankRFShuffStruct = structFromFieldNames(blankRF_fields, nan(1, opt.nPhaseShuffles, opt.n3 ));
                
                blankSingleRFStat = struct('cc', nan, 'jackStdErr', nan);

            case 'degree'
                doDegree_Ori = strcmp(cmpOriSpfType, 'ori');
                doDegree_Spf = strcmp(cmpOriSpfType, 'spf');
                if doDirMeasures
                    dirMeasureFields = {'D_dir_pref', nan, 'D_dsi_glob', nan, 'D_dsi_loc', nan, 'D_dsi_glob_unrec', nan, 'D_dsi_loc_unrec', nan, 'D_aligned_pair', nan};
                else
                    dirMeasureFields = {};
                end
                if doOSmeasures
                    OSmeasureFields_ori = {'cc_OSP', nan, 'cc_ori', nan};
                    OSmeasureFields_spf = {'cc_spf', nan};
                else
                    OSmeasureFields_ori = {};
                    OSmeasureFields_spf = {};
                end
                
                dF1oDC_ori_fields = {'D_F1oDC_ori', nan, 'D_F1pair_ori', nan};
                dF1oDC_spf_fields = {'D_F1oDC_spf', nan, 'D_F1pair_spf', nan};
                oriFields = {'D_ori_pref', nan, 'Dw_ori_glob', nan, 'Dw_ori_loc', nan};
                spfFields = {'Dw_spf', nan, 'D_w_SLN', nan, 'D_Bn', nan, 'D_skew', nan, 'D_spf_pref', nan};
                if doR90_measures
                    R90Fields = {'dR_spont_abs', nan, 'dR90_total_abs', nan, 'dR90_stim_abs', nan, 'dR_spont_rel', nan, 'dR90_total_rel', nan, 'dR90_stim_rel', nan};
                else
                    R90Fields = {};
                end

                if doDegree_Ori
                    blankDegreeTuningStatStruct = struct(oriFields{:}, dirMeasureFields{:}, dF1oDC_ori_fields{:}, OSmeasureFields_ori{:}, R90Fields{:});
                elseif doDegree_Spf
                    blankDegreeTuningStatStruct = struct(spfFields{:}, OSmeasureFields_spf{:}, dF1oDC_spf_fields{:});
                end
        end        
        featureDistsFields_blank = [cellClosenessMeasures; num2cell(nan(1, length(cellClosenessMeasures)))];
%         featureDistsStruct_blank = struct(featureDistsFields_blank{:});
               
        pairInfo.pairIdx = [];
        
        pairData_method = 2; % 1 = a little faster, simpler, more memory intensive. 2 = a little slower, more memory efficient
        
        [SCtype_pref_fields, F1oDC_maxStdJack_fields, minF1oDC_fields, maxF1oDC_fields] = deal( cell(2, length(F1oDC_types)) );
        [SCtype_pref_vals,   F1oDC_maxJackStd_vals,   minF1oDC_vals,   maxF1oDC_vals] = deal( cell(1, length(F1oDC_types)) );
        
        SCtype_pref_fields(1,:) = cellfun(@(s) ['SCtype_pref_' s ], F1oDC_types, 'un', 0);
        F1oDC_maxStdJack_fields(1,:) = cellfun(@(s) ['F1oDC_' s '_maxJackStd'], F1oDC_types, 'un', 0);
        minF1oDC_fields(1,:) = cellfun(@(s) ['minF1oDC_' s], F1oDC_types, 'un', 0);
        maxF1oDC_fields(1,:) = cellfun(@(s) ['maxF1oDC_' s], F1oDC_types, 'un', 0);
        
        for pair_ind = 1:nCurPairs
            progressBar;
%             if pair_ind == 15, return; end;
            [ind1, ind2]  = deal(pairs(pair_ind,1), pairs(pair_ind,2));
            [ind1_id, ind2_id] = deal(ind1, ind2);
            if randomlySwapPairs && rand < .5
                [ind1_id, ind2_id] = deal(ind2_id, ind1_id);
            end

%             [frac1AtMax2, frac2AtMax1] = deal(allPairsToDo(pair_ind).frac1AtMax2, allPairsToDo(pair_ind).frac2AtMax1);            
%             if (frac1AtMax2 > frac2AtMax1)
%                 [ind1, ind2] = deal(ind2, ind1); % make ind1 the one with the higher fraction.
%             end            
%             flds = {'Gid', 'cellId', 'contrib', 'maxR', 'stats', 'ph'};            
            
            Cell1 = allCells(ind1_id);  
            Cell2 = allCells(ind2_id);    
            
%             [Gid1, cellId1, maxR1, stats1, phs1, ori_sp_maxMU1, locData1, phaseTC1_cc, phaseTC1_cc_p, phaseTC1_dot] = deal( ...
%                 Cell1.Gid, Cell1.cellId, Cell1.maxR, Cell1.stats, Cell1.ph, Cell1.ori_sp_maxMU, Cell1.locData, Cell1.phaseTC_cc, Cell1.phaseTC_cc_p, Cell1.phaseTC_dot);    
%             [Gid2, cellId2, maxR2, stats2, phs2, ori_sp_maxMU2, locData2, phaseTC2_cc, phaseTC2_cc_p, phaseTC2_dot] = deal(...
%                 Cell2.Gid, Cell2.cellId, Cell2.maxR, Cell2.stats, Cell2.ph, Cell2.ori_sp_maxMU, Cell2.locData, Cell2.phaseTC_cc, Cell2.phaseTC_cc_p, Cell2.phaseTC_dot);            

            [Gid1, cellId1, windowStats1, locData1] = deal( ...   ori_sp_maxMU1, 
                Cell1.Gid, Cell1.cellId, Cell1.windowStats, Cell1.locData);   % Cell1.ori_sp_maxMU, 
            [Gid2, cellId2, windowStats2, locData2] = deal(...   ori_sp_maxMU1, 
                Cell2.Gid, Cell2.cellId, Cell2.windowStats, Cell2.locData);  % Cell2.ori_sp_maxMU,           
            
            nPh = length(Cell1.ph);
%             assert(all( [Cell1.PSTH.bins(1), Cell2.PSTH.bins(1)] < -295));
            
            GC = [Gid1, cellId1, Gid2, cellId2];

%             if isequal(GC, [4790, 1, 4790, 2])
%                 3;
%             end
            
%             if ~isempty(smoothPhases)
%                 R1_noSmooth = allCells_noSmooth(ind1_id).R;  
%                 R2_noSmooth = allCells_noSmooth(ind2_id).R;  
%                 Rs_noSmooth = {R1_noSmooth, R2_noSmooth};
%             else
                Rs_noSmooth = {};
%             end
            sameGroup = Gid1 == Gid2;
            if sameGroup
                [sameAnimal, sameHemisphere, samePenetration, sameLocation] = deal(1);
                [penetrationDistance, dAP, dML, depthDistance] = deal(0);
            else
                sameAnimal      = locData1.CatId == locData2.CatId;
                sameHemisphere  = locData1.CatId == locData2.CatId && locData1.HemiId == locData2.HemiId;
                samePenetration = locData1.PenId == locData2.PenId;
                sameLocation    = locData1.LocId == locData2.LocId;
                
                if sameLocation
                    penetrationDistance = 0;
                    dAP = 0; dML = 0;
                    depthDistance = 0;
                elseif samePenetration
                    penetrationDistance = 0;
                    dAP = 0; dML = 0;
                    depthDistance = abs(locData1.depth - locData2.depth); 
                elseif sameAnimal 
                    
                    if any([locData1.AP, locData2.AP, locData1.ML, locData2.ML]) == 0
                        dAP = nan;
                        dML = nan;
                        penetrationDistance = nan;
                    else
                        dAP = abs(locData1.AP - locData2.AP); 
                        dML = abs(locData1.ML - locData2.ML);
                        penetrationDistance = hypot( dAP, dML );
                    end
                    
                    depthDistance = nan;
                else
                    penetrationDistance = nan;
                    dAP = nan; dML = nan;
                    depthDistance = nan;
                end
            
            end                

            
%             F1oDC_types = fieldnames(Cell1.F1oDC.types);
            F1oDCs1 = Cell1.F1oDC;
            F1oDCs2 = Cell2.F1oDC;
            
            for ti = 1:length(F1oDC_types)
                type_i = F1oDC_types{ti};
                
                cells_f1odc         = [F1oDCs1.(['F1oDC_' type_i]),         F1oDCs2.(['F1oDC_' type_i])        ];
                cells_f1odc_jackStd = [F1oDCs1.(['F1oDC_jackStd_' type_i]), F1oDCs2.(['F1oDC_jackStd_' type_i])];

                F1oDC_pref_maxJackStd = max(cells_f1odc_jackStd);
                multiunits = [cellId1, cellId2] == 0;
                simplecells = cells_f1odc > 1;

                if any(multiunits)
                    % [0, 1, 1.5, 2] = [mC-cC,  mC-cS,  mS-cC,  mS-cS]
                    SCtype_pref = sum( simplecells );  
                    if (SCtype_pref == 1) && any(multiunits & simplecells) 
                        SCtype_pref = 1.5; % mS-cC --> 1.5   (mC-cS --remains--> 1.)  
                    end
                else
                    if any(isnan(cells_f1odc))
                        SCtype_pref = nan;
                    else
                        SCtype_pref = sum(simplecells); % [0, 1, 2] = [CC, SC/CS, SS];
                    end
                end            

                SCtype_pref_vals{ti} = SCtype_pref;
                F1oDC_maxJackStd_vals{ti} = F1oDC_pref_maxJackStd;
                minF1oDC_vals{ti} = min(cells_f1odc);
                maxF1oDC_vals{ti} = max(cells_f1odc);
                
            end

            SCtype_pref_fields(2,:) = SCtype_pref_vals;
            F1oDC_maxStdJack_fields(2,:) = F1oDC_maxJackStd_vals;
            minF1oDC_fields(2,:) = minF1oDC_vals;
            maxF1oDC_fields(2,:) = maxF1oDC_vals;

%             doPhaseTuningComparisons = strcmp(cmpType, 'phase');
%             doDegreeOfTuningComparisons = strcmp(cmpType, 'degree');
            switch cmpType, 
                case 'phase'
%                     [curPairStats, shuffPairStats] = 
                    pairInfo.sameGroup = Gid1==Gid2;
                    pairInfo.sameLocation = sameLocation;
%                     if isfield(Cell1, 'STAdata')
%                         opts.STA_cc = pearsonR(Cell1.STA(:), Cell2.STA(:));
%                     end
                    if isfield(Cell1, 'R_oe') && isfield(Cell2, 'R_oe')
                        R_oe = {Cell1.R_oe, Cell2.R_oe};
                    else
                        R_oe = {};
                    end
                    
%                     pairInfo.ori_sp_maxMU = iff(~isempty(ori_sp_maxMU1), ori_sp_maxMU1, ori_sp_maxMU2); % maybe decide which one to use (1 vs 2)?
                    if doWscc
                        pairInfo.pairIdx = cellPairIdx(cellId1+1, cellId2+1);
                    end
                    pairInfo.SCtype_pref = SCtype_pref;
                    
                    pairInfo.phaseTC_cc1 = Cell1.phaseTC_CCs_oe;
                    pairInfo.phaseTC_cc2 = Cell2.phaseTC_CCs_oe;                    
                    if opt.PTC_require_oe_corr
                        if strcmp(opt.PTC_oe_corrMode, 'oe')
                            pairInfo.phaseTC_cc1 = Cell1.phaseTC_CCs_oe;
                            pairInfo.phaseTC_cc2 = Cell2.phaseTC_CCs_oe;
                        elseif strcmp(opt.PTC_oe_corrMode, 'hoe')
                            pairInfo.phaseTC_cc1 = Cell1.phaseTC_CCs_hoe;
                            pairInfo.phaseTC_cc2 = Cell2.phaseTC_CCs_hoe;                            
                        elseif strcmp(opt.PTC_oe_corrMode, 'fs')
                            pairInfo.phaseTC_cc1 = Cell1.phaseTC_CCs_fs;
                            pairInfo.phaseTC_cc2 = Cell2.phaseTC_CCs_fs;                            
                        end
                    end
                    
                    
                    if opt.PTC_require_oe_corr_p
                        if strcmp(opt.PTC_oe_p_corrMode, 'oe')                        
                            pairInfo.phaseTC_cc_p1 = Cell1.phaseTC_CC_ps_oe;
                            pairInfo.phaseTC_cc_p2 = Cell2.phaseTC_CC_ps_oe;
                        elseif strcmp(opt.PTC_oe_p_corrMode, 'hoe')
                            pairInfo.phaseTC_cc_p1 = Cell1.phaseTC_CC_ps_hoe;
                            pairInfo.phaseTC_cc_p2 = Cell2.phaseTC_CC_ps_hoe;                            
                        elseif strcmp(opt.PTC_oe_p_corrMode, 'fs')
                            pairInfo.phaseTC_cc_p1 = Cell1.phaseTC_CC_ps_fs;
                            pairInfo.phaseTC_cc_p2 = Cell2.phaseTC_CC_ps_fs;                            
                        end
                    end
                    
                    if doF1oDC && isfield(Cell1, 'F1oDCs')
                        pairInfo.F1oDCs1 = Cell1.F1oDCs;
                        pairInfo.F1oDCs2 = Cell2.F1oDCs;                        
                    end
                    pairInfo.oris = Cell1.ori;
                    pairInfo.spfs = Cell1.sp;
                    pairInfo.phases = Cell1.ph;
                    pairInfo.Cell1 = Cell1;
                    pairInfo.Cell2 = Cell2;
%                     if isnan(Cell1.F1oDC_maxR_avP_sm) || isnan(Cell2.F1oDC_maxR_avP_sm)
%                         3;
%                     end
%                     if 
                                        
                    [curPairStats, oriSpCmpInds, phaseTCs, shuffPairStats] = getPhaseTuningStatsFor2OSPs(...
                        Cell1.ph, Cell1.R, Cell2.R, pairInfo, blankPhaseTuningStats, blankPhaseTuningShuff, R_oe{:});
%                         blankPhaseTuningStatStruct, opts, R_oe{:});
%                     assert(isnan(curPairStats(2).cc_atMaxDphi));                    
%                     curPairStats(2).cc_atMaxDphi = STA_cc;
%                     if isnan(curPairStats.cc)
%                         3;
%                     end

                    [nOri, nSpf, nPh] = size(Cell1.R);
%                     oriSpCmpIdxs = sub2indV([nOri, nSpf], oriSpCmpInds'); 
%                     oriSpCmpIdxs = oriSpCmpIdxs(:)';
                    doChecks = 0;
                    if doChecks
                        stat_nums = cellfun(@(fld) curPairStats.(fld), measures, 'un', 0);
                        stat_nums = [stat_nums{:}];

                        if any(isnan(stat_nums)) && ~all(isnan(stat_nums))
                            3;  
                        end
                        if all(isnan(stat_nums))
                            pairsToSkip(pairsToSkip_i) = pair_ind;
                            pairsToSkip_i = pairsToSkip_i+1;
                        end
                        if doDPhi && doDF1
                            dphiF1 = [curPairStats.dphi, curPairStats.dF1];
                            assert( all( dphiF1 <= 180 | isnan(dphiF1)) );
                        end
                    end
                    
                    if strcmp(gratingType, 'flashed') && doingRFmeasures
                        %%
                        [RF_curPairStats, RF_shuffPairStats, RF_stats] = getReceptiveFieldCorrelationsFor2Cells(Cell1, Cell2, blankRFStats, blankRFShuffStruct, blankSingleRFStat);
                        %%
                        RF_flds = fieldnames(RF_curPairStats);
                        for fld_i = 1:length(RF_flds)
                            curPairStats(1).(RF_flds{fld_i}) = RF_curPairStats.(RF_flds{fld_i});
                        end
                        
                        RF_shuff_flds = fieldnames(RF_shuffPairStats);
                        for fld_i = 1:length(RF_shuff_flds)
                            shuffPairStats(1).(RF_shuff_flds{fld_i}) = RF_shuffPairStats.(RF_shuff_flds{fld_i});
                        end                                        
                    else
                        RF_stats = struct;
                    end
                    
                    if opt.getDegreeDataForPhaseCmp
                        if firstDegreeStat
                            blankDegreeTuningStatStruct = struct;
                        end                    
                        curDegreePairStats = getSpfDegreeTuningStats_for2Cells(Cell1.tuningStats, Cell2.tuningStats, blankDegreeTuningStatStruct);
                       
                        if firstDegreeStat
                            firstDegreeStat = false;
                            blankDegreeTuningStatStruct = blankStruct( curDegreePairStats, nan );
                        end
                    end                    
                    
                    
                case 'degree',                    
                    curPairStats = blankDegreeTuningStatStruct;
                    if doDegree_Ori
                        curPairStats = getOriDegreeTuningStats_for2Cells(Cell1.tuningStats, Cell2.tuningStats, curPairStats);
                    end
                    if doDegree_Spf
                        curPairStats = getSpfDegreeTuningStats_for2Cells(Cell1.tuningStats, Cell2.tuningStats, curPairStats);
                    end
%                     curPairStats = mergeStructs(oriPairStats, spfPairStats);
                          
                case 'clusters',
                    curPairStats = struct;                    
                    [Clust1, Clust2] = deal(Cell1, Cell2);
                    if (Gid1 == Gid2)
                        clustId1 = Clust1.clustId;
                        clustId2 = Clust2.clustId;
                        
                        Gid_idx = find(S_refr.gids == Gid1);
                        clustIdx1 = find(S_refr.groupClustIds{Gid_idx} == clustId1, 1);
                        clustIdx2 = find(S_refr.groupClustIds{Gid_idx} == clustId2, 1);                        
                        
                        if clustIdx1 > clustIdx2
                            [clustIdx1, clustIdx2] = deal(clustIdx2, clustIdx1);
                        end
                        
                        
                        clust1Nspk = S_refr.groupClustNSpk{Gid_idx}(clustIdx1);
                        clust2Nspk = S_refr.groupClustNSpk{Gid_idx}(clustIdx2);                        
                        
                        pairRefrPct = S_refr.groupNRefrSpikes{Gid_idx}(clustIdx1, clustIdx2);% / (clust1Nspk + clust2Nspk);
                        pairPrunePct = S_refr.groupCrossPruningMtxs{Gid_idx}(clustIdx1, clustIdx2);% / (clust1Nspk + clust2Nspk);
                        pairMinIsis = S_refr.groupMinIsis{Gid_idx}(clustIdx1, clustIdx2);
                        pairRefrPeriod_ms = S_refr.groupRefrPeriods{Gid_idx}(clustIdx1, clustIdx2);
                        pairCCG16ms = S_refr.groupCcg16ms{Gid_idx}(clustIdx2, clustIdx1); % must have idx1 > idx2
                        assert(~isempty(pairCCG16ms));                        
                        
                        pairCcgCenterRatio16 = S_refr.groupCcgCenterRatio16{Gid_idx}(clustIdx1, clustIdx2);
                        pairCcgLrRatio16 = S_refr.groupCcgLrRatio16{Gid_idx}(clustIdx1, clustIdx2);

                        pairCcgCenterRatio4 = S_refr.groupCcgCenterRatio4{Gid_idx}(clustIdx1, clustIdx2);
                        pairCcgLrRatio4 = S_refr.groupCcgLrRatio4{Gid_idx}(clustIdx1, clustIdx2);
                        
                    else
                        3;
                    end
                        
                                                                                        
                    
                case 'psth'                
                    curPairStats = getStatsFor2PSTHs(Cell1.PSTH, Cell2.PSTH);                    

            end
%             if any([curPairStats(2).dF1] == 0)
%                 3;
%             end
           
                        
            
            if (Gid1 == Gid2)
                fet1 = Cell1.spkFeatures;
                fet2 = Cell2.spkFeatures;
                
                if any(strcmp(cellClosenessMeasures, {'negAmps_dist'}))
                    negAmps1 = fet1.negAmps_ccw_mean;
                    negAmps2 = fet2.negAmps_ccw_mean;                    
                    featureDists.negAmps_dist = norm(negAmps1-negAmps2);
                end

                if any(strcmp(cellClosenessMeasures, {'negAmps_cc'}))
                    negAmps1 = fet1.negAmps_ccw_mean;
                    negAmps2 = fet2.negAmps_ccw_mean;
                    featureDists.negAmps_cc = 1-pearsonR(negAmps1, negAmps2);
                end

                if any(strcmp(cellClosenessMeasures, {'negAmps_overlap'}))
                    M1 = fet1.negAmps_ccw_mean;  C1 = fet1.negAmps_ccw_cov;
                    M2 = fet2.negAmps_ccw_mean;  C2 = fet2.negAmps_ccw_cov;                
                    featureDists.negAmps_overlap = - (quadProdGaussians(M1, C1, M2, C2, 'log')...
                                                     -quadProdGaussians(M1, C1, M1, C2, 'log')); % normalize by product if had the same mean
                end

                if any(strcmp(cellClosenessMeasures, {'fullWaveform_cc'}))
                    wvfm_ccw1 = fet1.wvfm_ccw_mean;
                    wvfm_ccw2 = fet2.wvfm_ccw_mean;
                    featureDists.fullWaveform_cc = 1-pearsonR(wvfm_ccw1, wvfm_ccw2);
                end

                if any(strcmp(cellClosenessMeasures, {'fullWaveform_ed'}))
                    wvfm_ccw1 = fet1.wvfm_ccw_mean;
                    wvfm_ccw2 = fet2.wvfm_ccw_mean;
                    featureDists.fullWaveform_ed = norm(wvfm_ccw1-wvfm_ccw2);
                end            

                if any(strcmp(cellClosenessMeasures, {'channelWvfm_meanCC'}))
                    [nT, nChannels] = deal(43, 4);                                
                    wvfm_scl1 =  reshape(fet1.wvfm_scl_mean(:), [nT, nChannels]) ;
                    wvfm_scl2 =  reshape(fet2.wvfm_scl_mean(:), [nT, nChannels]) ;
                    channelCCs = pearsonR(wvfm_scl1, wvfm_scl2, 1);
                    featureDists.channelWvfm_meanCC = 1-mean(channelCCs);                
                end

                if any(strcmp(cellClosenessMeasures, {'PCA_overlap'}))
                    if isfield(fet1, 'PCAcuw_mean')
                        M1 = fet1.PCAcuw_mean;  C1 = fet1.PCAcuw_cov;
                        M2 = fet2.PCAcuw_mean;  C2 = fet2.PCAcuw_cov;                
                        featureDists.PCA_overlap = -( quadProdGaussians(M1, C1, M2, C2, 'log')-quadProdGaussians(M1, C1, M1, C2, 'log')); % normalize by integral if had the same mean
                    else
                        featureDists.PCA_overlap = nan;
                    end                    
                end
                

                if any(strcmp(cellClosenessMeasures, {'GLF_overlap'})) 
                    if isfield(fet1, 'GLFcuw_mean')
                        M1 = fet1.GLFcuw_mean;  C1 = fet1.GLFcuw_cov;
                        M2 = fet2.GLFcuw_mean;  C2 = fet2.GLFcuw_cov;                
                        featureDists.GLF_overlap = -( quadProdGaussians(M1, C1, M2, C2, 'log')-quadProdGaussians(M1, C1, M1, C2, 'log')); % normalize by product if had the same mean
                    else
                        featureDists.GLF_overlap = nan;
                    end                                            
                end            
                
                if any(strcmp(cellClosenessMeasures, {'minID'}))
                    ID1 = fet1.IsolationDistance; ID2 = fet2.IsolationDistance;
                    if any(isnan([ID1, ID2]))
                        minID = nan;
                    else
                        minID = min(ID1, ID2);
                    end
                    featureDists.minID = minID;
                end
                
                if any(strcmp(cellClosenessMeasures, {'maxL_ratio'}))
                    featureDists.maxL_ratio = min(fet1.L_ratio, fet2.L_ratio);
                end

                if any(strcmp(cellClosenessMeasures, {'diffFWHM'}))
                    FWHM1 = fet1.spkWidthHeight.FWHM;
                    FWHM2 = fet2.spkWidthHeight.FWHM;
                    featureDists.diffFWHM = abs(FWHM1 - FWHM2);
                end

                if any(strcmp(cellClosenessMeasures, {'diffPtPWidth'}))
                    PtP_width1 = fet1.spkWidthHeight.PtP_width;
                    PtP_width2 = fet2.spkWidthHeight.PtP_width;
                    featureDists.diffPtPWidth = abs(PtP_width1 - PtP_width2);
                end
                
                featureDistsFields = [fieldnames(featureDists), struct2cell(featureDists)]';
            else                
                featureDistsFields = featureDistsFields_blank;
%                 struct('negAmps_dist', nan, 'negAmps_cc', nan, 'negAmps_overlap', nan, 'fullWaveform_cc', nan,...
%                     'fullWaveform_ed', nan, 'channelWvfm_meanCC', nan, 'PCA_overlap', nan, 'GLF_overlap', nan);
            end
            
            f = @(x) x; % rep p-values are already in -log10 form
            
            clustIdFields = {};
            phaseTCfields = {};
            degreeCmpFields = {};
            clustersCmpFields = {};
            STA_fields = {};
            MID_fields = {};
            SpfCmp_fields = {};
            minCOV_fields = {};
            
            if strcmp(responseType, 'gainCorrected')
                cov_field = iff(strcmp(gratingType, 'flashed'), 'coeff_of_var', 'coeff_of_var_cycAv');
                covs = [Cell1.GC_stats.(cov_field), Cell2.GC_stats.(cov_field)];
                minCOV_fields = {'minCOV', min(covs)};
            end
            
            if strcmp(cmpType, 'phase')
                if any(isnan([curPairStats.F1oDC1]))
                    3;
                end
                SCtype_cmp = sum( [curPairStats.F1oDC1; curPairStats.F1oDC2] > 1 );
                TC_fields = iff(includePhaseTC, {'phaseTCs', phaseTCs}, {});                
                if strcmp(gratingType, 'flashed') && doingRFmeasures
                    if doSTAcc
                        STA_fields = {'min_rqr_oe_STA', RF_stats.min_rqr_oe_STA};
                    end
                    if doMIDcc || doMID_fit_cc || doMID_ovlp_cc
                        MID_fields = {'min_rsqr_fit', RF_stats.min_rsqr_fit, 'max_rsqr_fit', RF_stats.max_rsqr_fit, ...
                            'min_rsqr_oe', RF_stats.min_rsqr_oe, 'max_rsqr_oe', RF_stats.max_rsqr_oe, ...     
                            'min_frac_ovlp', RF_stats.min_frac_overlap, 'max_frac_ovlp', RF_stats.max_frac_overlap, ...  
                            'min_ovlp', RF_stats.min_ovlp, ...
                            ...'min_OSI', RF_stats.min_OSI, 'max_OSI', RF_stats.max_OSI
                            };                
                    end
                end
                oriSpCmp_fields = iff(includeOriSp_cmp, {'oriSp_cmp', oriSpCmpInds}, {});
                    
                if opt.getDegreeDataForPhaseCmp
                     SpfCmp_fields = {'Dw_spf', curDegreePairStats.Dw_spf, 'D_spf_pref', curDegreePairStats.D_spf_pref};
                end
                if calcLowerSumPhs
                    lowerSumPhs_flds = {'loc_minSumPhs', [curPairStats.lowerSumPhs]}; 
                end
                
                phaseTCfields = {... 'loc_F1oDCs_cmp',  [curPairStats.F1oDC1; curPairStats.F1oDC2], ...   all the "loc_*" entries each have one column per location.
                                 'loc_minF1oDC_cmp',  min([curPairStats.F1oDC1; curPairStats.F1oDC2], [], 1), ...                                 
                                 'loc_maxF1oDC_cmp',  max([curPairStats.F1oDC1; curPairStats.F1oDC2], [], 1), ...                                 
                                 'SCtype_cmp', SCtype_cmp, ...    
                                 'loc_cc_atMaxDphi', [curPairStats.cc_atMaxDphi], ...  % fix if do more than one location.
                                 'loc_minFracOfMaxes', min([curPairStats.frac1ofMax; curPairStats.frac2ofMax], [], 1), ...  
                                 ...'loc_minPhaseTC_cc',  min([phaseTC1_cc(oriSpCmpIdxs); phaseTC2_cc(oriSpCmpIdxs)], [], 1), ...  
                                 ...'loc_minPhaseTC_cc_p',  min([phaseTC1_cc_p(oriSpCmpIdxs); phaseTC2_cc_p(oriSpCmpIdxs)], [], 1), ...  
                                 ...'loc_minPhaseTC_dot',  min([phaseTC1_dot(oriSpCmpIdxs); phaseTC2_dot(oriSpCmpIdxs)], [], 1), ...  
                                 ...'loc_minSumPhs',      min([curPairStats.lower sumPhs1; curPairStats.sumPhs2],[], 1), ...         
                                 lowerSumPhs_flds{:}, ...  
                                 'loc_ptcOEcorr',    [curPairStats.ptcOEcorr], ...
                                 oriSpCmp_fields{:}, ... 
                                 TC_fields{:}, ...
                                 SpfCmp_fields{:} ...
                                }; 
            end

            if strcmp(cmpType, 'clusters')
                clustersCmpFields = {'pairRefrPct', pairRefrPct, 'pairPrunePct', pairPrunePct, 'pairMinIsis', pairMinIsis, 'pairCCG16ms', pairCCG16ms, ... 
                    'pairRefrPeriod_ms', pairRefrPeriod_ms,  'pairCcgCenterRatio16', pairCcgCenterRatio16, 'pairCcgLrRatio16', pairCcgLrRatio16, ...
                    'pairCcgCenterRatio4', pairCcgCenterRatio4, 'pairCcgLrRatio4', pairCcgLrRatio4};
                clustIdFields = {'clustIds', [clustId1 clustId2], 'nSpikes', [clust1Nspk, clust2Nspk]};
            end
            
%             if strcmp(cmpType, 'degree')
%                 degreeCmpFields = {'oriPairStats', oriPairStats, 'spfPairStats', spfPairStats};
%             else
%                 degreeCmpFields = {};
%             end
            
            statsForThisPair = struct(...               
                'Gids', [Gid1, Gid2], ...
                'cellIds', [cellId1, cellId2], ...   
                clustIdFields{:}, ...
                ...'F1oDCs_pref', [Cell1.(opt.F1oDC_field), Cell2.(opt.F1oDC_field)],...
                ...'minF1oDC_pref', [min([Cell1.(opt.F1oDC_field), Cell2.(opt.F1oDC_field)])], ...
                ...'maxF1oDC_pref', [max([Cell1.(opt.F1oDC_field), Cell2.(opt.F1oDC_field)])], ...
                ...'SCtype_pref', SCtype_pref, 'F1oDC_pref_maxJackStd', F1oDC_pref_maxJackStd, ...    
                SCtype_pref_fields{:}, F1oDC_maxStdJack_fields{:}, ...
                minF1oDC_fields{:}, maxF1oDC_fields{:}, ...
                'animalCmp', sameAnimal, 'hemiCmp', sameHemisphere, 'penetrCmp', samePenetration, 'locCmp', sameLocation, ...
                'penDist', penetrationDistance, 'depthDist', depthDistance, 'dAP', dAP, 'dML', dML, ...
                'min_rep_cc_p', min( f(windowStats1.cc_p), f(windowStats2.cc_p) ), ...
                ... 'min_rep_rho_p', min( f(stats1.allWindowStats.rho_p), f(stats2.allWindowStats.cc_p) ), ...
                'n_phases', nPh, ...                                
                phaseTCfields{:}, ...
                STA_fields{:}, ...
                MID_fields{:}, ...
                degreeCmpFields{:}, ...
                clustersCmpFields{:}, ...
                featureDistsFields{:}, ...
                minCOV_fields{:} ...
                );

            statsForThisPair = structfun(@single, statsForThisPair, 'un', 0);
            
            if checkForImaginaryVals 
                vals = struct2cell(statsForThisPair); vals = cellfun(@(x) x(:)', vals, 'un', 0); vals = vals(cellfun(@isnumeric, vals));
                assert(~ any(imag([vals{:}]) > 0));
            end
            
            if isempty(pairData)                
                pairData_fields = fieldnames(statsForThisPair);
                if pairData_method == 1
                    clear('pairData');
                    pairData(nPairs,1) = blankStruct(statsForThisPair);
                else
                    % pairData = blankStruct(statsForThisPair);                
                    size1 = structfun(@(x) size(x,1), statsForThisPair);
                    if any(size1 > 1)
                        fieldnms = fieldnames(statsForThisPair);
                        error('Field "%s" has more than 1 row"', fieldnms{find(size1>1,1)});
                    end
                    pairData = structfun(@(x) zeros(nPairs, size(x,2), size(x,3) ), statsForThisPair, 'un', 0);
                    
                end
            end
                        
            pair_idx = idxMtx( ind1, ind2);
            if pairData_method == 1
                c = struct2cell(pairData( pair_idx ));
                assert(isempty([c{:}]));

                pairData( pair_idx ) = statsForThisPair;            
            else
                for fld_i = 1:length(pairData_fields)
                    pairData.(pairData_fields{fld_i})(pair_idx,:,:) = statsForThisPair.(pairData_fields{fld_i});                
                end
            end
            
            for loc_i = 1:nLocations
                pairStats_loc_i = curPairStats(loc_i);
                if doWscc
                    shuffStats_loc_i = shuffPairStats(loc_i);
                else
                    shuffStats_loc_i = struct;
                end
                
                for m_i = 1:nMeasures 
                    if strcmp(cmpType, 'phase') && any(strcmp(measures{m_i}, RF_Measures)) && (loc_i > 1)
                        continue;
                    end
                    
                    is_pval = calcPvalues && any(strcmp(measures{m_i}, {'cc', 'rho', 'tau'}));
                    jackStdErr_field = [measures{m_i} '_jackStd'];
                    haveJackStdErr = isfield(pairStats_loc_i, jackStdErr_field);
                    haveShuffJackStdErr = isfield(shuffStats_loc_i, jackStdErr_field);
                    if pair_ind == 1
                        allStatsC{loc_i,m_i}.measure = measures{m_i};
                        allStatsC{loc_i,m_i}.location = locations{loc_i};
                        if ~isfield(allStatsC{loc_i,m_i}, 'val')
                            allStatsC{loc_i,m_i}.val = zeros(nPairs, 1, opt.n3); 
                        end
                        if doWscc && ~isfield(allStatsC{loc_i,m_i}, 'shuffVal')
                            allStatsC{loc_i,m_i}.shuffVal = cell(nPairs, 1); %zeros(nPairs, opt.nPhaseShuffles, opt.n3, 'single'); 
                            allStatsC{loc_i,m_i}.ctrlVal = cell(nPairs, 2); 
                        end
                        if is_pval && ~isfield(allStatsC{loc_i,m_i}, 'pval')
                            allStatsC{loc_i,m_i}.pval = zeros(nPairs, 1); 
                        end
                        if haveJackStdErr % opt.getPhaseCmpJackStdErr_ptc &&
                            allStatsC{loc_i,m_i}.jackStd = zeros(nPairs, 1, opt.n3); 
                        end
                        if haveShuffJackStdErr % opt.getPhaseCmpJackStdErr_ptc &&
                            allStatsC{loc_i,m_i}.shuffJackStd = zeros(nPairs, opt.nPhaseShuffles, opt.n3); 
                        end
                        
                    end
                    curPairIdx = idxMtx(ind1, ind2);
                    val = single( pairStats_loc_i.(measures{m_i}) );
                    allStatsC{loc_i,m_i}.val( curPairIdx, 1, : ) = val;
                    
                    if doWscc 
                        shuffVals = single( shuffStats_loc_i.(measures{m_i}) );
                        allStatsC{loc_i,m_i}.shuffVal{ curPairIdx } = shuffVals;
                                                
                        if ~any(isnan(shuffVals(:))) && ~isempty(shuffVals)

                            switch measures{m_i}
                                case {'cc', 'dphi', 'dF1'},
                                    nPhases = statsForThisPair.n_phases;
                                    nPhases_idx = find(all_nPhases == nPhases, 1);                        
                                    idx_first = Wscc_randperm_idxFirst{nPhases_idx, pairInfo.pairIdx};

                                    ctrlVals = shuffVals( 1, idx_first, : );
                                    val_weight = lcm_nPhases / nPhases;
                                    assert(length(unique(ctrlVals)) == length(unique(shuffVals)));                                    
                                    
                                case {'MID_cc', 'MID_ovlp_cc', 'STA_cc'},                                    
                                    ctrlVals = [val, -1*val];
                                    val_weight = 1;                                    
                                    assert(length(unique(shuffVals)) == 2*opt.n3 );
                                case {'MID_fit_cc', 'dPh_rel'},                                    
                                    ctrlVals = shuffVals;
                                    val_weight = 1;
%                                     assert(length(unique(shuffVals)) >= (opt.nPhaseShuffles *.99) );
                                otherwise
                                    error('Unhandled case');
                            end
                        
                            allStatsC{loc_i,m_i}.ctrlVal{ curPairIdx, 1 } = val_weight;
                            allStatsC{loc_i,m_i}.ctrlVal{ curPairIdx, 2 } = ctrlVals;
                        else
                            3;
                        end
                            
                            
                    end                    
                                        
                    if is_pval
                        pval = pairStats_loc_i.([measures{m_i} '_p']);
                        allStatsC{loc_i,m_i}.pval( curPairIdx ) = pval;
                    end
                    
                    if haveJackStdErr
                        jackStdErr = pairStats_loc_i.(jackStdErr_field);
                        allStatsC{loc_i,m_i}.jackStd( curPairIdx, :,: ) = jackStdErr;
                    end
                    
                    if haveShuffJackStdErr
                        shuffJackStdErr = shuffStats_loc_i.(jackStdErr_field);
                        allStatsC{loc_i,m_i}.shuffJackStd( curPairIdx, :,: ) = shuffJackStdErr;
                    end
                    
                    
                end
            end
                    
        end
                
        
        % convert pairData : array of struct --> struct of arrays (much more space-efficient)
        %%
        if pairData_method == 1
            for fld_i = 1:length(pairData_fields)
                pairData_array.(pairData_fields{fld_i}) = cat(1, pairData.(pairData_fields{fld_i}));
            end
            pairData = pairData_array;
        end
        
%         tf = isequalwithequalnans(pairData2, pairData);
%         profile viewer;
        3;
        
        %%
        % add bin data to each statistic
        for loc_i = 1:nLocations            
            for m_i = 1:nMeasures       
                if strcmp(cmpType, 'phase') && any(strcmp(measures{m_i}, RF_Measures)) && loc_i > 1
                    continue;
                end
                
                
                binEdges = measureBins{m_i};
                
                if strncmp(locations{loc_i}, 'wgtedSum', 6)
                    binEdges = linspace(binEdges(1), binEdges(end), 35);
                end                
                nBins = length(binEdges)-1; % **

                if ~isfield(allStatsC{loc_i,m_i}, 'bins')
                    allStatsC{loc_i,m_i}.bins = zeros(nPairs, 1, 'uint8');
%                     allStatsC{loc_i,m_i}.binMask = false(nBins, nPairs);
                    allStatsC{loc_i,m_i}.binEdges = binEdges;
                end
                
%                 bins = zeros(nBins, nCurPairs);
                vals = [allStatsC{loc_i,m_i}.val( idxMtx(pairIdxs) )];                
                valLims = lims(vals(~isnan(vals)));
                if isempty(valLims)
                    valLims = [0 1];
                end 
                if (valLims(1) < binEdges(1)) || (valLims(2) > binEdges(end))
                                        
                    warning('Extending binSize for measure : %s', measures{m_i});
                    dBin = binEdges(2)-binEdges(1);
                    bE1 = binEdges(1);
                    bEn = binEdges(end);
                    nbins = length(binEdges);
                    keyboard;
                    while (valLims(1) < bE1)                         
                        bE1 = bE1-dBin;
                        nbins = nbins +1;
                    end
                    while (valLims(2) > bEn)
                        bEn = bEn+dBin;
                        nbins = nbins +1;
                    end
                    binEdges = linspace(bE1, bEn, nbins);
                    
                    assert( (valLims(1) >= binEdges(1)) && (valLims(2) <= binEdges(end)) )
                end                    
                
                [n, whichBins] = histcnt(vals, binEdges);
%                 for bi = 1:nBins
%                     bins(bi, (whichBins == bi)) = true;                    
%                 end
%                 
                allStatsC{loc_i,m_i}.bins( idxMtx(pairIdxs) ) = whichBins;
                assert(all(whichBins(~isnan(vals)) > 0));
%                 allStatsC{loc_i,m_i}.binMask(:, idxMtx(pairIdxs) ) = bins;                
            end
        end

%         progressBar('done');     
        
        
        
    end

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%----------    A. CALCULATIONS    ------- %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
%     if ~exist('compareSpatialData.mat', 'file') || redo

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Compute the phase statistics for all cell pairs %%%%%%

        showWorking = false;
        doCalculationsForPairs(allPairIdxs, 'All cell pairs');
        
        pairsToSkip = pairsToSkip(1:pairsToSkip_i-1);
%         fprintf('%d ', pairsToSkip);  % these are pairs where all the
%         measures are "NaN"s. could save space and remove these pairs,
%         but would be too much trouble - and are very few of these pairs
%         anyway.
%         fprintf('\nN = %d\n', length(pairsToSkip));
%         profile viewer
        varsToSave = {'pairData', 'allStatsC',   'pairTypes', 'measures', 'locations', 'opt'};
        if recordAllCCs
            
            %%
            allWCCs(allWCC_idx:end,:) = [];
            oe_mode_str = iff(strcmp(opt.phase_oeState, 'oe'), ['_' opt.phase_oe_mode], '');
            s_cc.([gratingType '_' opt.phase_oeState oe_mode_str]) = allWCCs; 
            allCCs_filename = [CatV1Path 'allCCs_vs_rank' curMatchDB('') '.mat'];
            append_str = iff(exist(allCCs_filename, 'file'), {'-append'}, {});
            save(allCCs_filename, '-struct', 's_cc', append_str{:})                            
            3;
                        
        end
        
        tic;
        fprintf('\nSaving data to %s ... ', cmpDatafile );
        save([pathname cmpDatafile], varsToSave{:}, '-v6'); % '-v6' to skip file compression.
        fprintf(' done. ');
%         count6
        toc;
                
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%----------      B. PLOTTING      ------- %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% 1. COMPARE ORI & SPATIAL FREQ.
    %{
    function normalizedhist(data, controldata, nBins, colr)
        [nData, binCenters] = hist(data, nBins);
        nControlData = hist(controldata, binCenters);
        h = bar(binCenters, nData ./ nControlData);
        if exist('colr', 'var')
            set(h, 'color', colr);
        end
    end
    
    function plotOriSpComprisons(allOriMaxes, allOrisShown, dOri, alldOris, allSpMaxes, allSpsShown, dSp, alldSps, figbase, colr, labl) %#ok<INUSL>
        do2Dplot = false;
        nCols = iff(do2Dplot, 3, 2);
        figure(figbase); clf;
        
        % 1a) Plot Individual Orientations:
        subplot(2,nCols,1);
        normalizedhist(allOriMaxes, allOrisShown, 10, colr);
        title('Orientation: Individual'); 
        
%                 set(findobj(gca,'Type','patch'),'FaceColor',colr);
        subplot(2,nCols,nCols+1); hist(dOri, 19); title('Orientation: Pair differences');
%         xlim([0 90]);
                set(findobj(gca,'Type','patch'),'FaceColor',colr);
        
        % 2a) Plot Individual Spatial Periods
        subplot(2,nCols,1);
        normalizedhist(log2(allSpMaxes), log2(allOrisShown), 9, colr);
        title('log_2 ( Spatial period: Indiv) ');
        
        
        % 2b) Plot Spatial Period octaves:                
        subplot(2,nCols,nCols+2); hist(dSp, 9); title('Spatial period: Octave differences');
                set(findobj(gca,'Type','patch'),'FaceColor',colr);
        
        if do2Dplot
            subplot(2,nCols,3);       hist2d([ allOris, log2(allSps)], 10, 10, 'plot'); xlabel('ori'); ylabel('sp')
            subplot(2,nCols,nCols+3); hist2d([ dOri,    dSp         ], 10, 10, 'plot'); xlabel('ori'); ylabel('sp');
        end
        suptitle_2([labl, ' : N = ' num2str(length(allOris)/2) ]);
    end

    %}

%     function ccs = RF_randInvertPhase(RFdata1, RFdata2, rfType, nRand)        
%         
%     %     cc_pos = pearsonR(MID1, MID2);
%     %     cc_inv = pearsonR(MID1, -MID2);
% 
%         cc_pos = rfCC(RFdata1, RFdata2, rfType);
%         cc_inv = -cc_pos;
% 
%         idx = randi(2, 1, nRand);
%         ccs = zeros(1, nRand);
% 
%         ccs(idx == 1) = cc_pos;
%         ccs(idx == 2) = cc_inv;
%     
%     end

    

end


%{
Considerations for choosing SC type

When plotting, would like to have S-S at the bottom,
so S-S should have smaller corresponding value(?)

same colors:
	SS & cS-mS - 0
	SC & cC-mS - 1
 	     cS-mC - 1.5	
	CC & cC-mC - 2

vs: want to have CC first, so that matches minF1/DC series better


how many simple cells
	CC & cC-mC - 0
 	SC & cS-mC - 1
	     cC-mS - 1.5
	SS & cS-mS - 2

%}
    function c = combineValues(a, b, measure_name, opt)
   
        doAverage = opt.oe_average || (~any(strcmp(measure_name, opt.measures)) && isempty(strfind(measure_name, 'jackStd')) );
        if doAverage
            c = (a + b)/2;
        elseif opt.oe_keepBoth
            c = cat(3, a, b);  % keep both values, concatenated in 3rd dimension.
        elseif opt.oe_keep1
            c = a;
        elseif opt.oe_keep2
            c = b;
        end            
        
    end



    function [oriMax, spMax, oris, sps] = findMaxOfOSP(ind)
        osp = allCells(ind);
        oris = osp.oris;
        sps = osp.sps;
        mOSP = mean(osp.R, 3);
        
        [tmp, max1inds] = maxElement(mOSP);
        if (tmp == 0)
            oriMax = NaN;
            spMax  = NaN;
        else    
            oriMax = oris(max1inds(1));
            spMax  = sps( max1inds(2));
        end
    end
    
    function [allOriMaxes, allSpMaxes, allOrisShown, allSpsShown, dOri, dSp] = getdOridSpDistribution(allPairs) %#ok<DEFNU>
        if strcmp(gratingType, 'drifting')
            maxOri = 360;
        elseif strcmp(gratingType, 'flashed')
            maxOri = 180;
        end
        nPairsToDo = length(allPairs);
        allOrisShown = cell(nPairsToDo,1);
        allSpsShown  = cell(nPairsToDo,1);
        
        ori1max = zeros(nPairsToDo, 1); sp1max = zeros(nPairsToDo, 1);
        ori2max = zeros(nPairsToDo, 1); sp2max = zeros(nPairsToDo, 1);
        
                
        dOri = zeros(nPairsToDo, 1);  dSp = zeros(nPairsToDo, 1);
        for i = 1:nPairsToDo
            [ind1, ind2] = elements(allPairs(i,:));
            
            [ori1max(i), sp1max(i), ospOris1, ospSps1] = findMaxOfOSP(ind1);
            [ori2max(i), sp2max(i), ospOris2, ospSps2] = findMaxOfOSP(ind2);
            allOrisShown{i} = [ospOris1(:); ospOris2(:)];
            allSpsShown{i}  = [ospSps1(:); ospSps2(:)];
                        
            dOri(i) = circDist( ori1max(i), ori2max(i), maxOri);
            dSp(i)  = abs(  log2 (sp1max(i)/sp2max(i)) );
        end
        
        allOriMaxes = [ori1max; ori2max];
        allSpMaxes =  [sp1max;  sp2max];
        allOrisShown = [allOrisShown{:}];
        allSpsShown  = [allSpsShown{:}];        
    end

    
%     function y = singleIfNumeric(x)
%         if isnumeric(x)
%             y = single(x);
%         else 
%             y = x;
%         end
%     end
    
    function S = getRefrClusterData(Gids)
        redo = 0;
        [filepath, filename1] = getFileName('clusterData2', Gids(1));
        filename = strrep(filename1, ['Group_' num2str(Gids(1))], 'allGroups');        
        if ~exist([filepath, filename], 'file') || redo
                
            allGids = sort(getAllGids);
            nGids = length(allGids);
            S.gids = allGids;
            S.groupClustIds = cell(1,nGids);
            S.groupClustNSpk = cell(1,nGids);
            S.groupCrossPruningMtxs = cell(1, nGids);
            S.groupNRefrSpikes = cell(1, nGids);
            S.groupRefrPeriods = cell(1, nGids);
            S.groupMinIsis = cell(1, nGids);
            S.groupCcgCenterRatio = cell(1, nGids);
            S.groupCcgLrRatio     = cell(1, nGids);            
            S.groupCcg16ms        = cell(1, nGids);            
            
                range_ms = 16; nbins = 80;
                [allRanges_ms, allNbins] = getGlobals('isi_allRanges_ms', 'isi_allNbins');        
                range_idx = find(allRanges_ms==range_ms, 1);
                nbins_idx = find(nbins == allNbins, 1);            

            fprintf('Generating Pairwise-Refractory data file...\n');
            progressBar('init-', nGids);
            for gi = 1:nGids
                progressBar;
                clusterData = getClusterData(allGids(gi), 2);
                if isempty(fieldnames(clusterData))
                    continue;
                end
                S.groupClustIds{gi} = clusterData.clustIds;
                S.groupClustNSpk{gi} = clusterData.clustNSpikes;
                S.groupNRefrSpikes{gi} = clusterData.isiData.pairNRefrSpikes;                
                S.groupCrossPruningMtxs{gi} = clusterData.crossRefrData.crossPruningData.nSpksRemoved;                                
                S.groupRefrPeriods{gi} = clusterData.crossRefrData.crossPruningData.refrPeriod_ms;
%                 pairRefracPeriod_ms = S_refr.groupMinIsis{Gid_idx}(clustIdx1, clustIdx2);
                
                allMin_isis_C = arrayfun(@(s) min(abs(s.isi_ms)), clusterData.crossRefrData.refrSpikesInfo, 'un', 0);                                
                allMin_isis_C(cellfun(@isempty, allMin_isis_C)) = {0};
                allMin_isis = cell2mat(allMin_isis_C);
                allMin_isis = allMin_isis+allMin_isis';                
                S.groupMinIsis{gi} = allMin_isis;
                                
                masterBinIdxs = clusterData.ccgData.mbIdxs{range_idx, nbins_idx};                                
                groupCcg16ms = masterBin2Bin(clusterData.ccgData.ccgs_mb, masterBinIdxs);
                S.groupCcg16ms{gi} = groupCcg16ms;
                
                [centerRatio16, LRratio16] = getCenterLRratios(clusterData.ccgData, 16);
                [centerRatio4, LRratio4] = getCenterLRratios(clusterData.ccgData, 4);                
                S.groupCcgCenterRatio16{gi} = centerRatio16;
                S.groupCcgLrRatio16{gi}     = LRratio16;
                S.groupCcgCenterRatio4{gi} = centerRatio4;
                S.groupCcgLrRatio4{gi}     = LRratio4;
                
                
            end
            save([filepath, filename], '-struct', 'S')
            
        else
            S = load([filepath, filename]);            
        end
        
        idxs_retrieve = binarySearch(sort(S.gids), Gids, [], 'exact');
        S.gids = S.gids(idxs_retrieve);
        S.groupClustIds = S.groupClustIds(idxs_retrieve);
        S.groupClustNSpk = S.groupClustNSpk(idxs_retrieve);
        S.groupNRefrSpikes = S.groupNRefrSpikes(idxs_retrieve);
        S.groupCrossPruningMtxs = S.groupCrossPruningMtxs(idxs_retrieve);            
        S.groupMinIsis = S.groupMinIsis(idxs_retrieve);            
        S.groupRefrPeriods = S.groupRefrPeriods(idxs_retrieve);
        S.groupCcg16ms = S.groupCcg16ms(idxs_retrieve);
        S.groupCcgCenterRatio16 = S.groupCcgCenterRatio16(idxs_retrieve);
        S.groupCcgLrRatio16 = S.groupCcgLrRatio16(idxs_retrieve);        
        S.groupCcgCenterRatio4 = S.groupCcgCenterRatio4(idxs_retrieve);
        S.groupCcgLrRatio4 = S.groupCcgLrRatio4(idxs_retrieve);
    
    end
    
    
    function [centerRatio, LR_ratio, ccgBinVals] = getCenterLRratios(ccgData, windowSize_ms)
                        
        lastBin_ms = ccgData.mbEdges(end);
        idxNegMs = indmin(abs(ccgData.mbEdges-(-windowSize_ms)));
        idxPosMs = indmin(abs(ccgData.mbEdges-( windowSize_ms)));
        idx0Bin  = indmin(abs(ccgData.mbEdges-(0)));                                          
        idxNotEmpty = ~cellfun(@isempty, ccgData.ccgs_mb);
        [rs,cs] = find(idxNotEmpty);
        [centerRatio, LR_ratio] = deal( zeros(size( ccgData.ccgs_mb)) );                
        for i = 1:length(rs)                    
            r = rs(i); c = cs(i);
            binVals = ccgData.ccgs_mb{r,c};
            NegRate = sum(binVals(idxNegMs:idx0Bin-1))/windowSize_ms;
            PosRate = sum(binVals(idx0Bin:idxPosMs))/windowSize_ms;
            bckgRate = sum(binVals)/(lastBin_ms*2);
            centerRatio(r,c) = (PosRate+NegRate)-(bckgRate);         centerRatio(c,r) = centerRatio(r,c);
            LR_ratio(r,c) = log10((PosRate+1e-10)/(NegRate+1e-10));  LR_ratio(c,r)    = LR_ratio(r,c);
        end                
        
    end
    
    function X = flip_rot(X, flip_i, rot_i)
        if flip_i
            X = fliplr(X);            
        end
        if rot_i > 0
            X = rot90(X, rot_i);
        end
    end

    
    function RF_stats = getPairReceptiveFieldStats(RF1_S, RF2_S, rf_subType, RF_stats, opt)
                       %%
        rf_subType_jack = rf_subType;
        if any(strcmp(rf_subType, {'STA', 'MID', 'MID_fit'}))
            RF1 = RF1_S.(rf_subType);
            RF2 = RF2_S.(rf_subType);
        elseif strcmp(rf_subType, {'MID_overlap'})                        
            L = size(RF1_S.MID,1);
            mid_overlap = ( RF1_S.MID_select & RF2_S.MID_select );
            frac = nnz(mid_overlap) / L^2;

            frac_overlap_with_1 = nnz(mid_overlap & RF1_S.MID_select)/nnz(RF1_S.MID_select);
            frac_overlap_with_2 = nnz(mid_overlap & RF2_S.MID_select)/nnz(RF2_S.MID_select);                        
            min_ovlp = min(frac_overlap_with_1, frac_overlap_with_2);

            if (frac < opt.MID_ovlp_minFrac)
                return;
            end      
            RF1 = RF1_S.MID( mid_overlap);
            RF2 = RF2_S.MID( mid_overlap);

            [RF_stats.mid_overlap1, RF_stats.mid_overlap2] = deal(mid_overlap);
            [RF_stats.min_ovlp1, RF_stats.min_ovlp2] = deal(min_ovlp);
            rf_subType_jack = 'MID';
        end

        if ~isempty(RF1) && ~isempty(RF2)
            cc = pearsonR(RF1, RF2);
            RF_stats.cc = cc;

            if opt.getPhaseCmpJackStdErr_rf
                jack_fld = [rf_subType_jack '_jacks'];
                if isfield(RF1_S, jack_fld) && isfield(RF2_S, jack_fld)
                    RF1_jacks = RF1_S.(jack_fld);
                    RF2_jacks = RF2_S.(jack_fld);
                    
                    if strcmp(rf_subType, 'MID_overlap')
                        RF1_jacks = cellfun(@(X) X(mid_overlap), RF1_jacks, 'un', 0);
                        RF2_jacks = cellfun(@(X) X(mid_overlap), RF2_jacks, 'un', 0);
                    end
                    
                    cc_jackknives = cellfun(@pearsonR, RF1_jacks, RF2_jacks);
                    RF_stats.jackStdErr = jackknifeStdErr(cc_jackknives, cc);
                end
            end
        end

    
    end

    
%         function pairPhaseStruct_av = getPairPhaseTuningStats_OE(phases, ph_tc1_a,    ph_tc1_b,    ph_tc2_a,    ph_tc2_b,  ...
%                                                                      ph_jacks1_a, ph_jacks1_b, ph_jacks2_a, ph_jacks2_b, blankPhaseStruct, varargin)

    function pairRFStruct_av = getPairReceptiveFieldStats_OE(RF1_a_S, RF1_b_S,  RF2_a_S, RF2_b_S, rf_subType, blankRF_struct, opt)
        %%
        RF_stats_a = getPairReceptiveFieldStats(RF1_a_S, RF2_a_S, rf_subType, blankRF_struct, opt);
        RF_stats_b = getPairReceptiveFieldStats(RF1_b_S, RF2_b_S, rf_subType, blankRF_struct, opt);
        
        pairRFStruct_av = blankRF_struct;
        
        fn = fieldnames(RF_stats_a);        
        for fld_i = 1:length(fn)
            fld = fn{fld_i};
            pairRFStruct_av.(fld) = combineValues(RF_stats_a.(fld), RF_stats_b.(fld), fld, opt);                            
        end

    end
    
    
    function [cc, jackStdErr, RF_stats] = rfCC(Cell1, Cell2, rf_type, rf_subType, blankRF_struct, opt)
%     [cc, cc_a, cc_b, mid_overlap_a, mid_overlap_b, min_ovlp_a, min_ovlp_b] 

%         fn_rf_odd = [rf_type '_odd'];
%         fn_rf_even = [rf_type '_even'];

%         haveALL =  isfield(Cell1, rf_type) && isfield(Cell2, rf_type);
%         haveOE = isfield(Cell1, fn_rf_odd) && isfield(Cell1, fn_rf_even) && ...
%                         isfield(Cell2, fn_rf_odd) && isfield(Cell2, fn_rf_even);
%         if haveALL && ~haveOE
%             3;
%         end
        cc = nan;
        jackStdErr = nan;
        RF_stats = struct;
        
%         S = struct('cc_a', nan, 'cc_b', nan, 'jackStdErr', nan);
%         [cc, cc_a, cc_b] = deal(  nan(1,1,opt.n3)  );
%         [mid_overlap_a, mid_overlap_b, min_ovlp_a, min_ovlp_b] = deal(nan);
        
        
        L = Cell1.L; L2 = Cell2.L;        
        if (L2 ~= L)
            return;
        end
%         Lsqr = L^2;
        
        if nargin < 4
            rf_subType = '';
        end
        
        switch opt.phase_oeState
            case 'aa'
                if isfield(Cell1, rf_type) && isfield(Cell2, rf_type)
                    RF1_S = Cell1.(rf_type);
                    RF2_S = Cell2.(rf_type);
                    
                    RF_stats = getPairReceptiveFieldStats(RF1_S, RF2_S, rf_subType, blankRF_struct, opt);                       
                    
                end
                                    
            case 'oe'
                fn_rf_odd  = [rf_type '_odd'];
                fn_rf_even = [rf_type '_even'];
                if isfield(Cell1, fn_rf_odd) && isfield(Cell1, fn_rf_even) && ...
                   isfield(Cell2, fn_rf_odd) && isfield(Cell2, fn_rf_even)                         
                
                    switch opt.phase_oe_mode
                        case 'diff',
                            RF1_a_S = Cell1.(fn_rf_odd);
                            RF1_b_S = Cell1.(fn_rf_even);
                            RF2_b_S = Cell2.(fn_rf_odd);
                            RF2_a_S = Cell2.(fn_rf_even);

                        case 'same',
                            RF1_a_S = Cell1.(fn_rf_odd);
                            RF1_b_S = Cell1.(fn_rf_even);
                            RF2_a_S = Cell2.(fn_rf_odd);
                            RF2_b_S = Cell2.(fn_rf_even);
                    end                    
                    
                    RF_stats = getPairReceptiveFieldStats_OE(RF1_a_S, RF1_b_S,  RF2_a_S, RF2_b_S, rf_subType, blankRF_struct, opt);
                    
                end
        end
        cc = RF_stats.cc;
        jackStdErr = RF_stats.jackStdErr;
        
    end

    
    
    
    function [ccs, cc_jackStds] = MID_fits_randPhase_for2Cells(MID1data, MID2data, Gid, opt)
        setCto0 = 1;
        dSamp = 1;
    
        [xs, ys] = getStimulusXY(Gid);
        if dSamp > 1
            xs = xs(1:dSamp:end);
            ys = ys(1:dSamp:end);
        end
        
        [xs_grid, ys_grid] = meshgrid(xs, ys);
        XY = [xs_grid(:), ys_grid(:)];
        
        phis1_rand_all = rand(1, opt.nPhaseShuffles)*2*pi;
        phis2_rand_all = rand(1, opt.nPhaseShuffles)*2*pi;
        
        
       
        switch opt.phase_oeState
            case 'aa', 
%                 p1 = MID1data.MID.p;
%                 p2 = MID2data.MID.p;
%                 ccs = MID_fits_randPhase_for2MIDs(XY, p1, p2, phis1_rand_all, phis2_rand_all, useEnvelope);

                MID1 = MID1data.MID;
                MID2 = MID2data.MID;
                [ccs, cc_jackStds] = MID_fits_randPhase_for2MIDs(XY, MID1, MID2, phis1_rand_all, phis2_rand_all, opt);
                
                
            case 'oe', 
%                 p1_odd = MID1data.MID_odd.p;
%                 p1_even = MID1data.MID_even.p;
%                 p2_odd = MID2data.MID_odd.p;
%                 p2_even = MID2data.MID_even.p;
% 
%                                 
%                 switch opt.phase_oe_mode 
%                     case 'same', [p1a, p1b,  p2a, p2b] = deal(p1_odd, p2_odd,   p1_even, p2_even);
%                     case 'diff', [p1a, p1b,  p2a, p2b] = deal(p1_odd, p2_even,  p1_even, p2_odd );
%                 end              
%                 
%                 ccs_1 = MID_fits_randPhase_for2MIDs(XY, p1a, p1b, phis1_rand_all, phis2_rand_all, useEnvelope);
%                 ccs_2 = MID_fits_randPhase_for2MIDs(XY, p2a, p2b, phis1_rand_all, phis2_rand_all, useEnvelope);
                
                MID1_odd  = MID1data.MID_odd;
                MID1_even = MID1data.MID_even;
                MID2_odd  = MID2data.MID_odd;
                MID2_even = MID2data.MID_even;

                                
                switch opt.phase_oe_mode 
                    case 'same', [MID1a, MID1b,  MID2a, MID2b] = deal(MID1_odd, MID2_odd,   MID1_even, MID2_even);
                    case 'diff', [MID1a, MID1b,  MID2a, MID2b] = deal(MID1_odd, MID2_even,  MID1_even, MID2_odd );
                end              
                
                [ccs_1, cc_jackStds_1] = MID_fits_randPhase_for2MIDs(XY, MID1a, MID1b, phis1_rand_all, phis2_rand_all, opt);
                [ccs_2, cc_jackStds_2] = MID_fits_randPhase_for2MIDs(XY, MID2a, MID2b, phis1_rand_all, phis2_rand_all, opt);
                
                
                ccs         = combineValues(ccs_1, ccs_2, 'MID_fit_cc', opt);
                cc_jackStds = combineValues(cc_jackStds_1, cc_jackStds_2, 'MID_fit_cc', opt);
                
        end        
    end
    
    function [ccs, cc_jackStds] = MID_fits_randPhase_for2MIDs(XY, MID1, MID2, phis1_rand_all, phis2_rand_all, opt)
    
        
        useEnvelope = ~opt.doRandMID_fitsWithoutEnvelope;
         
        p1 = MID1.p;
        p2 = MID2.p;
        ccs = MID_fits_randPhase_for2Gabors(XY, p1, p2, phis1_rand_all, phis2_rand_all, useEnvelope);
        
        if isfield(MID1, 'p_jacks') && ~isempty(MID1.p_jacks) ...
                && isfield(MID2, 'p_jacks') && ~isempty(MID2.p_jacks) && opt.doRandMID_fits_jacks
            p1_jacks = MID1.p_jacks;
            p2_jacks = MID2.p_jacks;
            
            cc_jackknives = cellfun(@(p1, p2) MID_fits_randPhase_for2Gabors(XY, p1, p2, phis1_rand_all, phis2_rand_all, useEnvelope), p1_jacks, p2_jacks, 'un', 0);
            cc_jackknives_allShifts = [cat(1, cc_jackknives{:})];
            cc_jackStds = arrayfun(@(i) jackknifeStdErr(cc_jackknives_allShifts(:,i), ccs(i)), 1:opt.nPhaseShuffles );
        else
            cc_jackStds = [];
        end
                
    end

    function ccs = MID_fits_randPhase_for2Gabors(XY, p1, p2, phis1_rand_all, phis2_rand_all, useEnvelope)
%         ccs = MID_fits_randPhase(MID1data, MID2data, Gid, nRand)
                        
        [A1, mu_x1, mu_y1, sig_x1, sig_y1, theta1, k1, phi1, C1] = dealV(p1);
        [A2, mu_x2, mu_y2, sig_x2, sig_y2, theta2, k2, phi2, C2] = dealV(p2);
        
        nRand = length(phis1_rand_all);
%         cc_actual = pearsonR(MID1data.MID_fit, MID2data.MID_fit);
        
        if ~useEnvelope
            [sig_x1, sig_y1, sig_x2, sig_y2] = deal(1e10);
        end

        ccs = zeros(1,nRand);         
        doPlot = 0;
        
        

%                 figure(101);
%                 subplot(1,2,1); imagesc(MID1data.MID_fit)
%                 subplot(1,2,2); imagesc(MID2data.MID_fit);
%         3;
%         nMemMax_MB = 100; 
%         nMaxAtATime = (nMemMax_MB*1024^2) / ( numel(xs_grid)*4*2 );
%         nMaxAtATime = 500;
        nPerGroup = 50;
        nPerGroup = min(nPerGroup, nRand);
        nGroups = nRand/nPerGroup;
                
        grpIdxs = arrayfun(@(i) [1:nPerGroup] + (i-1)*nPerGroup, 1:nGroups, 'un', 0);               

%         ph = linspace(0, 2*pi, sqrt(nRand));
        
%         [phis1_rand_all, phis2_rand_all] = meshgrid(ph, ph);
%         phis1_rand_all = phis1_rand_all(:);
%         phis2_rand_all = phis2_rand_all(:);
        
        for grp_i = 1:nGroups
            
%             phis1_rand = rand(1,nPerGroup)*2*pi;
%             phis2_rand = rand(1,nPerGroup)*2*pi;
            phis1_rand = phis1_rand_all( grpIdxs{grp_i} );
            phis2_rand = phis2_rand_all( grpIdxs{grp_i} );
                        
            Z1 = gabor(A1, mu_x1, mu_y1, sig_x1, sig_y1, theta1, k1, phis1_rand, C1, XY);
            Z2 = gabor(A2, mu_x2, mu_y2, sig_x2, sig_y2, theta2, k2, phis2_rand, C2, XY);
            
            ccs( grpIdxs{grp_i} ) = pearsonR(Z1, Z2, 1);            
        end

        
        %%
%         figure(4);
%         ccs_mtx = reshape(ccs, 100, 100);
%         imagesc(ph, ph, ccs_mtx); axis xy; colorbar;
%         3;
        
        
    
    end

    function ccs = randFlipSignHalf(cc, nRand)
        %%
        nCCs = size(cc,3);
        ccs = zeros(1, nRand, nCCs);
        idx = randi(2, 1, nRand, nCCs);        
    
        for i = 1:nCCs
            ccs(1,idx(:,:,i) == 1,i)  = cc(1,1,i);
            ccs(1,idx(:,:,i) == 2,i) = -cc(1,1,i);
        end
        
    end    
    
    function dPh_rel = dRelativePhase(MIDdata1, MIDdata2, opt)

        dPh_rel = nan;
        
    
        switch opt.phase_oeState
            case 'aa', 
                
                if ~isfield(MIDdata1, 'MID') || ~isfield(MIDdata2, 'MID') ...
                    || isempty(MIDdata1.MID.p) || isempty(MIDdata2.MID.p)
                    return;
                end                
                
                rel_ph1 = rad2deg(MIDdata1.MID.p(8));
                rel_ph2 = rad2deg(MIDdata2.MID.p(8));
                dPh_rel = circDist(rel_ph1, rel_ph2, 360);
                    
            case 'oe'
                
                if ~isfield(MIDdata1, 'MID_odd') || ~isfield(MIDdata1, 'MID_odd') || ...
                   ~isfield(MIDdata2, 'MID_even') || ~isfield(MIDdata2, 'MID_even') || ...
                    isempty(MIDdata1.MID_odd.p) || isempty(MIDdata1.MID_even.p) || ...
                    isempty(MIDdata2.MID_odd.p) || isempty(MIDdata2.MID_even.p)
                    return;
                end
                
                rel_ph1_odd = rad2deg(MIDdata1.MID_odd.p(8));
                rel_ph1_even = rad2deg(MIDdata1.MID_even.p(8));
                rel_ph2_odd = rad2deg(MIDdata2.MID_odd.p(8));
                rel_ph2_even = rad2deg(MIDdata2.MID_even.p(8));
                
                switch opt.phase_oe_mode
                    case 'same',
                        dPh_rel1 = circDist(rel_ph1_odd, rel_ph2_odd, 360);
                        dPh_rel2 = circDist(rel_ph1_even, rel_ph2_even, 360);
                    case 'diff',
                        dPh_rel1 = circDist(rel_ph1_odd, rel_ph2_even, 360);
                        dPh_rel2 = circDist(rel_ph1_even, rel_ph2_odd, 360);
                end
                                
                dPh_rel = combineValues(dPh_rel1, dPh_rel2, 'dPh_rel', opt);
                
        end    
        
    end
                            
    function dPh_rel_rand = dRelativePhase_rand(nPhaseShuffles, opt)
    
        switch opt.phase_oeState
            case 'aa',                 
                dPh_rel_rand = circDist(rand(1, nPhaseShuffles)*360, rand(1, nPhaseShuffles)*360, 360);  
                                    
            case 'oe'
                
                dPh_rel_rand1 = circDist(rand(1, nPhaseShuffles)*360, rand(1, nPhaseShuffles)*360, 360);
                dPh_rel_rand2 = circDist(rand(1, nPhaseShuffles)*360, rand(1, nPhaseShuffles)*360, 360);
                
                dPh_rel_rand = combineValues(dPh_rel_rand1, dPh_rel_rand2, 'dPh_rel', opt);
                                 
        end        
    end
    
    
    function f = getCellAligned_factor(dori_mu)
        if isnan(dori_mu)
            f = nan;
        elseif dori_mu < 45
            f = 1;
        elseif dori_mu > 135
            f = 10;
        elseif (dori_mu >= 45 && dori_mu <= 135)
            f = 100;
        end
    end
  
    function ccs = getCCsFromShifts(x1, x2)
        x1 = x1(:); x2 = x2(:);
        %%
        ccs = zeros(1, length(x1));
        for i = 1:length(x1)
            x2_i = circshift(x2, i-1);
            ccs(i) = pearsonR(x1, x2_i);      
        end
    end
    
%     function r = flicked(r)
%         if ~all(r==0) && all(diff(r)==0)
%             i = randi(length(r));
%             r(i) = r(i) + 1e-5;
%         end
%     end


%                 switch distMeth
%                     case 'dist'
%                         featureDist = norm(electM1-electM2);
%                         
%                     case '1DgaussOverlap'
%                         amps1 = Cell1.spkFeatures.allFeatures  ;	amps2 = Cell2.spkFeatures.allFeatures;
%                         nspk1 = size(amps1,1);         nspk2 = size(amps2,1);
%                         ds1 = zeros(1,nspk1);           ds2 = zeros(1,nspk2);
%                         for i1 = 1:nspk1
%                             ds1(i1) = norm(electM1 - projectionFromPointToLine(electM1, electM2, amps1(i1,:)));                    
%                         end
%                         for i2 = 1:nspk2
%                             ds2(i2) = norm(electM1 - projectionFromPointToLine(electM1, electM2, amps2(i2,:)));                    
%                         end
%                         m1 = mean(ds1); m2 = mean(ds2);
%                         s1 = std(ds1);  s2 = std(ds2);                                        
%                         dbug = false;
%                         if dbug
%                             featureDist = gaussiansOverlap(m1, s1, m2, s2, 1);
%                             displayNHist({ds1, ds2});
%                             title(sprintf('s = %.2g', featureDist));
%                             3;
%                         else
%                             featureDist = gaussiansOverlap(m1, s1, m2, s2);
%                         end                        
% 
%                     case '4DgaussOverlap'
%                         amps1 = Cell1.spkFeatures.allFeatures  ;	amps2 = Cell2.spkFeatures.allFeatures;                        
%                         M1 = double( mean(amps1,1) );
%                         M2 = double( mean(amps2,1) );
%                         C1 = double( cov( amps1)    );
%                         C2 = double( cov( amps2)    );
%                         
% %                         getGroupWaveformCoefficients('PCA', Gid, 2, 'separate', [], [], matchSpiker);
%                         3;
%                         
% %                         featureDist(1,1) = norm(M1-M2);                        
%                         featureDist = quadProdGaussians(M1, C1, M2, C2)/quadProdGaussians(M1, C1, M1, C2); % normalize by product if had the same mean
% %                         a1 = quadProdGaussians(M1, C1, M1, C2);
% %                         a2 = quadProdGaussians(M2, C1, M2, C2);
% %                         assert(abs(a1-a2)/a1 < 1e-10);
% %                         3;
%                     case 'waveformCC',
%                         wvfm1 = Cell1.spkFeatures.meanWaveform(:);
%                         wvfm2 = Cell2.spkFeatures.meanWaveform(:);
%                         featureDist = pearsonR(wvfm1, wvfm2);
% 
%                     case 'waveformED';
%                         wvfm1 = Cell1.spkFeatures.meanWaveform(:);
%                         wvfm2 = Cell2.spkFeatures.meanWaveform(:);
%                         featureDist = norm(wvfm1-wvfm2);
%                         
%                     case 'perceptron'
%                         amps1 = Cell1.spkFeatures.allFeatures  ;	amps2 = Cell2.spkFeatures.allFeatures;
%                         allFeatures = [amps1; amps2];
%                         allIds  = [ones(size(amps1,1),1); 2*ones(size(amps2,1),1)];
%                         [Wt, Etot] = perceptron(allFeatures', allIds);
% %                         fprintf('%5.2g \n ', Etot);
%                         featureDist = Etot;
%                     
%                 end


%{
    interpolation for delta -phi
    if nInterp > 1
        [r1_tmp, phs_itp_tmp, ph1_itp] = fourierInterp(ph1, [], nInterp);
        [r2_tmp, phs_itp_tmp, ph2_itp] = fourierInterp(ph2, [], nInterp);                                        
        phs_itp = phs_itp_tmp/length(phs)*360;
    else 
        phs_itp = phs;
        ph1_itp = ph1; ph2_itp = ph2;
    end
%}

%{
        % 1-4.(a)  Get ori/sp indices for 4 possible basic outputs
        if doMaxR1 || doMean12,  [~, oriSpInds.maxR1] = maxElement(mR1);   end
        if doMaxR2 || doMean12,  [~, oriSpInds.maxR2] = maxElement(mR2);   end
        if doMaxR1xR2,           [~, oriSpInds.maxR1xR2] = maxElement(fR1 .* fR2); end
        if doMaxMinFracR,        [~, oriSpInds.maxMinFracR] = maxElement(min(fR1, fR2)); end
        if doMaxMU,                  oriSpInds.maxMU = pairInfo.ori_sp_maxMU; end
%}

%                     wrp = @(x) [x(:); x(1)];
%                     ext = @(x) [x(:); x(end)+diff(x(1:2))];

%{
        if doDPhi
            assert(nInterp_deltaPhi == 1);
            if length(phs) == size(ph1, 1)
                if nSamples == 1                
                    if any([F1oDC1 F1oDC2] < 1e-10)
                        dphi = nan;
                        cc_atMaxDphi = nan;
                    else                  
                        [dphi, cc_atMaxDphi] = deltaPhi( phs, ph1, ph2, dPhi_method);                         
                    end
                    idx_use = 1;
                else                              
                    dphi = nan(1,nSamples);
                    cc_atMaxDphi = nan(1,nSamples);
                    if nUse > 0
                        [dphi(idx_use), cc_atMaxDphi(idx_use)]  = deltaPhi( phs, ph1(:,idx_use), ph2(:,idx_use), dPhi_method);                     
                    end
                    
                    chk = 0;
                    if chk
                        dphi2 = nan(1,nSamples);
                        cc_atMaxDphi2 = nan(1,nSamples);                                                            
                        for i = 1:nUse;
                            [dphi2(idx_use(i)), cc_atMaxDphi2(idx_use(i))] = deltaPhi( phs, ph1(:,idx_use(i)), ph2(:,idx_use(i)), dPhi_method); 
                        end                                        
                        assert(isequalwithequalnans(dphi, dphi2));
                        assert(isequalwithequalnans(cc_atMaxDphi, cc_atMaxDphi2));
                    end
                end
                dphi2 = nan(1,nSamples);
                dphi2b = nan(1,nSamples);
                dphi2(idx_use) = bestCircShift(phs, ph1, ph2);  
                dphi2b(idx_use) = bestCircShift_Matlab(phs, ph1, ph2);
                assert(isequalwithequalnans(dphi, dphi2));
                assert(isequalwithequalnans(dphi, dphi2b));
                
            else
                dphi = 0; % for mtx x mtx case
                cc_atMaxDphi = nan;
            end
            phaseTuningStats.dphi = dphi;
            phaseTuningStats.cc_atMaxDphi = cc_atMaxDphi;
        end             
%}

%{

                        if isfield(Cell1, 'STA') && isfield(Cell2, 'STA') && ...
                                (length(Cell1.STA) == length(Cell2.STA))
                            STA_cc = pearsonR(Cell1.STA, Cell2.STA); 
                        else
                            STA_cc = nan;
                        end


%}

%{
                            if showExampleOfPhaseTuningCurves 
                                L_f1 = lims([pairPhaseStats.F1oDC1, pairPhaseStats.F1oDC2]);
                                dphi_val = pairPhaseStats.dphi;
                                cc_val = pairPhaseStats.cc;
                                ph_tc2_orig = ph_tc2;
                                ph_tc1_orig = ph_tc1;
                                r = max(ph_tc1)/max(ph_tc2); if r < 1, r = 1/r; end
                                if ibetween(dphi_val, [30, 40]) && (L_f1(1) > 1.1 && L_f1(2) < 1.7)% && (r < 1.8)
                                    %%
                                    phases_ext = [phases, 360];
                                    w = @(x) [x; x(1)];                             
                                    ph_tc1 = ph_tc1_orig*1.2;

                                    figure(30); clf; 
                                    col2 = [0 .7 0];
                                    h1 = plot(phases_ext, w(ph_tc1), 'b.-', phases_ext, w(ph_tc2), 'g.-');
                                    set(h1(2), 'color', col2);
                                    title('Phase Tuning curves');
                                    set(gca, 'xtick', 0:45:360); xlim([0 360]);
                                    xlabel('Phase'); ylabel('Firing Rate (Hz)');
                                    legend('Cell 1', 'Cell 2');

                                    ph_tc1_norm = ph_tc1-mean(ph_tc1); ph_tc1_norm = ph_tc1_norm/norm(ph_tc1_norm);
                                    ph_tc2_norm = ph_tc2-mean(ph_tc2); ph_tc2_norm = ph_tc2_norm/norm(ph_tc2_norm);


                                    n_shift = dphi_val / diff(phases(1:2));
                                    ph_tc1_norm_l = [ph_tc1_norm(n_shift+1:end); ph_tc1_norm(1:n_shift)];
                                    ph_tc1_norm_r = [ph_tc1_norm(end-n_shift+1:end); ph_tc1_norm(1:end-n_shift)];
                                    cc_l = corr(ph_tc1_norm_l, ph_tc2_norm);
                                    cc_r = corr(ph_tc1_norm_r, ph_tc2_norm);
                                    if cc_l > cc_r
                                        ph_tc1_norm_shft = ph_tc1_norm_l;
                                    else
                                        ph_tc1_norm_shft = ph_tc1_norm_r;
                                    end


                                    figure(31); clf;
                                    h2 = plot(phases_ext, w(ph_tc1_norm), 'b.-', phases_ext, w(ph_tc2_norm), 'g.-');
                                    set(h2(2), 'color', col2)
                                    title(sprintf('cc = %.2f', cc_val));
                                    set(gca, 'xtick', 0:45:360); xlim([0 360]);
                                    xlabel('Phase'); ylabel('Normalized Rate');
                                    drawHorizontalLine(0, 'linestyle', ':', 'color', 'k');

                                    figure(32); clf;
                                    h3 = plot(phases_ext, w(ph_tc1_norm), 'b.-', phases_ext, w(ph_tc2_norm), 'g.-', phases_ext, w(ph_tc1_norm_shft), 'b:');
                                    set(h3(3), 'linewidth', 2);
                                    set(h3(2), 'color', col2);
                                    title(sprintf('\\Delta\\phi = %.0f\\circ', dphi_val));
                                    set(gca, 'xtick', 0:45:360); xlim([0 360]);
                                    xlabel('Phase'); ylabel('Normalized Rate');
                                    drawHorizontalLine(0, 'linestyle', ':', 'color', 'k');
                                    3;
                                end
                            end
%}


%{


            multiunits = [cellId1, cellId2] == 0;
            cells_f1odc = [Cell1.(opt.F1oDC_field), Cell2.(opt.F1oDC_field)];
            simplecells = cells_f1odc > 1;
            if any(multiunits)
                % [0, 1, 1.5, 2] = [mC-cC,  mC-cS,  mS-cC,  mS-cS]
                SCtype_pref = sum( simplecells );  
                if (SCtype_pref == 1) && any(multiunits & simplecells) 
                    SCtype_pref = 1.5; % mS-cC --> 1.5   (mC-cS --remains--> 1.)  
                end
            else
                if any(isnan(cells_f1odc))
                    SCtype_pref = nan;
                else
                    SCtype_pref = sum(simplecells); % [0, 1, 2] = [CC, SC/CS, SS];
                end
            end

%}