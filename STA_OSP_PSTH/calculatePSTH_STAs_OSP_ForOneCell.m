function cellData = calculatePSTH_STAs_OSP_ForOneCell(Gid, cellId, oldCellData)

    global nCellsCompleted;
    global GC
    cellId = double(cellId);
    GC = [Gid, cellId];
    
    showResults = (nargout == 0);
%     showResults = 1;

%     fixedTimeWindows = [];
    fixedTimeWindows = [30 60; 60 90];       % actual windows: [29.2 - 62.5]  &  [58.3 - 91.6]
        nFixedTimeWindows = size(fixedTimeWindows, 1);            
    
    % Build on / replace previous data for this cell
    if exist('oldCellData', 'var') && ~isempty(oldCellData) && isstruct(oldCellData) && ~isfield(oldCellData, 'id')
        cellData = oldCellData;
        if (oldCellData.Gid ~= Gid) || (oldCellData.cellId ~= cellId)
            keyboard;
        end        
    else
        cellData = struct('Gid', Gid, 'cellId', cellId, 'PSTH', [], 'OSP', [], 'STAs', [], 'MIDs', []);
    end
    
    % Information about this experiment
    siteData = siteDataFor('Gid', Gid);
%     Did = siteData.Did;
    stimType_full = siteData.stimType;
    frameLength_ms = siteData.frameLength_ms;
    degreesPerPixel = siteData.stimulusInfo.degreesPerBlock;
    [stimType, stimSubType1] = strtok(stimType_full, ':');
    [stimSubType, stimTypeDetails] = strtok(stimSubType1, ':');
    
%     Did = dbLookup('Did',  'Gid', Gid);    
    
%     tableName = getDatabaseTableForDid(Did, stimType);
%     frameLength_ms = getFrameLength('Did', Did);
%     isFlashedGratingStimulus = (strcmp(stimType, 'Movie')   && strcmp(stimSubType, 'Flashed_Gratings') ) ...
%                           || (strcmp(stimType, 'Grating') && strcmp(stimSubType, 'Flashed Grating Batch') );
    isFlashedGratingMovie = (strcmp(stimType, 'Movie')   && strcmp(stimSubType, 'Flashed_Gratings') );
    isGratingStimulus =  isFlashedGratingMovie || ...
                         (strcmp(stimType, 'Grating') && ~strcmp(stimSubType, 'Single Grating') );

    gratingType = flashedOrDrifting(Gid, 's');
%     if isFlashedGratingMovie || (strcmp(stimType, 'Grating') && strcmp(stimSubType, 'Flashed Grating Batch') );
%         gratingType = 'flashed';
%     elseif (strcmp(stimType, 'Grating')) && ~strcmp(stimSubType, 'Single Grating');
%         gratingType = 'drifting';
%     end

%     osp_method = 'stimulus';  % options: 'stimulus' / 'spike'.
    % Decide which calculations to perform;
    calculatePSTH = (isFlashedGratingMovie); % && (isempty(cellData.PSTH) || ~isfield(cellData.PSTH, 'bckgRate') || (length(cellData.PSTH.bckgRate) < 2));
%     calculatePSTH = (frameLength_ms >= 100) && (~iscell(cellData.PSTH) || isempty(cellData.PSTH) || (length(cellData.PSTH) < 5));
    
    calculateOSP = true;%isGratingStimulus && (isempty(cellData.OSP)) || ~isfield(cellData.OSP, 'stats');
    matchDB = curMatchDB;
    calculateSTAs = ~matchDB && strcmp(gratingType, 'flashed'); %any(strcmp(stimType, {'Movie', 'Noise', 'Mseq'})) && isempty(cellData.STAs);
    
    % Display which calculations are to be performed:
    grpType = curGroupingType('');  
    grpType_sing = strrep(grpType, 'clusters', 'cluster'); grpType_sing = strrep(grpType_sing, 'cells', 'cell');
    fprintf('Group: %d, %s: %s ...', Gid, grpType_sing, iff(cellId == 100, '*', num2str(cellId)));
    if any([ calculateOSP, calculateSTAs, calculatePSTH ])
        fprintf([' Doing calculations: (' iff(calculatePSTH, ' PSTH;'), iff(calculateOSP, ' OSP;'), iff(calculateSTAs, ' STA;'),  ' ) ']);
        nCellsCompleted = nCellsCompleted+ 1;
    else
        fprintf([' Cell is already completed ' ]);        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%      1. CALCULATE PSTH  (so that can find best time window) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PSTHdata = [];
    
    if calculatePSTH
        % spkTsRelToFrame with   delay = 0 ms   for the PSTH.        
        [PSTHdata, stats] = getPSTHforCell(Gid, cellId);
%         figure(4);
%         plotThisPSTH(PSTHdata);
        cellData.PSTH = PSTHdata;
        
        timeWindow_ms = PSTHdata.timeWindow_ms;
        windowProfile = PSTHdata.windowProfile;
        bckgRate      = PSTHdata.bckgRate;
    else
        windowProfile = [];
        if isGratingStimulus  
            timeWindow_ms  = [0 0]; % just see how many spikes at each phase of the drifting grating stimulus.            
        else
            timeWindow_ms  = [30 70];
        end
        if strcmp(gratingType, 'drifting')
            PSTHdata.timeWindow_bins = [1, 1];
        end
%         timeWindow  = [25 65; 40 80; 55 95];
    end

    if calculateSTAs || calculateOSP
%         timeWindow_ms = [30 50];
%         windowProfile = [130 500 160 70];
%         if strncmp(osp_method, 'spike', 3)
%             relContrOfFrameToSpike = getParsedSpikes( 'frame', Gid, cellId, timeWindow_ms(1,:), windowProfile );
%         end        
%         degreesPerPixel = dbGetDegreesPerPixel(Did);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%      2. GENERATE O/S/P MAP   (if grating stimulus) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    if calculateOSP        
        addStandardOSP = strcmp(gratingType, 'flashed') && strcmp(grpType, 'clustersPruned');
                % for merging - is nice to have a OSP with a fixed/standardized time window, 
                % so the OSPs for different clusters can be added together.                         
        addFullOSP = true;
        addBckgSamples = true;
        keepAllTrials_flag = strcmp(grpType, 'cells'); %strcmp(gratingType, 'drifting');
        calcDegreeOfTuning = 1; %~matchDB ; %&& strcmp(grpType, 'cells');
        
        [bckgSamplesFields, OSP_full_fields, OSP_standard_fields] = deal({});        
        
        [R, R_full, oris, sps, phs, tf_Hz, meanRate, bckgSamples] = getOriSpfPhaseProfile_simple(Gid, cellId, PSTHdata, keepAllTrials_flag);
        
        if addStandardOSP
            PSTHdata_standard = struct('bins', PSTHdata.bins, 'meanRate', PSTHdata.meanRate);
            [R2, R_full2] = getOriSpfPhaseProfile_simple(Gid, cellId, PSTHdata_standard, []);
            OSP_standard_fields = {'R2', {R2}, 'R_full2', {R_full2}};
        end
        %%
        if (nFixedTimeWindows > 0) && strcmp(gratingType, 'flashed');
%%
            binC = PSTHdata.bins';
%             binE = binCent2edge(binC)';
            binWidth = diff(binC(1:2));
            allTimeWindows_bins = binarySearch(binC, fixedTimeWindows);
            allTimeWindows_ms = bsxfun(@plus, binC(allTimeWindows_bins),  (binWidth/2)*[-1, 1]);
  %%          
            PSTH_data_extra = PSTHdata;
            PSTH_data_extra.windowProfile = [];            
            
            R_full_extra = cell(1, nFixedTimeWindows);
            keepAllTrials_flag_extra = 0; % just want odd/even trials;
            for extra_i = 1:nFixedTimeWindows
%                 fixedTimeWindows_C{extra_i} = fixedTimeWindows(extra_i, :);
                PSTH_data_extra.timeWindow_bins = allTimeWindows_bins(extra_i,:);
                [~, R_full_extra{extra_i}] = getOriSpfPhaseProfile_simple(Gid, cellId, PSTH_data_extra, keepAllTrials_flag_extra);
            end
            3;
            extraOSP_fields = {'fixedWindowOSPs', struct('windows_bins', allTimeWindows_bins, 'windows_ms', allTimeWindows_ms, 'OSPs', {R_full_extra})};
        else
            extraOSP_fields = {};
        end
        
        %%
        if strcmp(gratingType, 'drifting') % if flashed grating, have already calculated the stats when calculating the PSTH. Drifting gratings: calculate now.
            stats.allWindowStats = getPSTHwindowData(Gid, cellId);
            stats.isRep = stats.allWindowStats.cc_p > 3;
        end
        doExtraWindows = 1;
        if calcDegreeOfTuning
            LR_bins = PSTHdata.timeWindow_bins;
            
            redo_flag = 735174;
%             stats.tuningStats = getOspDataForPsthWindow(Gid, cellId, [], [], LR_bins(1), LR_bins(2), windowProfile, {'tuningStats'}, redo_flag);
%             flds = fieldnames(stats.tuningStats);
%             if ~any( cellfun(@(s) ~isempty(strfind(s, '_odd')), flds) )
%                 stats.tuningStats = getOspDataForPsthWindow(Gid, cellId, [], [], LR_bins(1), LR_bins(2), windowProfile, {'tuningStats'}, redo_flag);
%             end            
            
%             tuningStats = calcDegreeOfTuningStats(R_full, bckgSamples, Gid, cellId);
            tuningStats = getOspDataForPsthWindow(Gid, cellId, [], [], LR_bins(1), LR_bins(2), windowProfile, {'tuningStats'});
            stats.tuningStats = tuningStats;
            
            if doExtraWindows
               3; 
                
            end
            
            3;
        end
        stats.meanRate = meanRate;
        
        if addFullOSP
            OSP_full_fields = {'R_full', compress(R_full)};
        end
        if addBckgSamples
            bckgSamplesFields = {'bckgSamples', compress(bckgSamples)};
        end
        %%
        TP_fields = {'tf_Hz', tf_Hz};
        cellData.OSP = struct('R', {R}, 'ori', {oris}, 'sp', {sps}, 'ph', {phs}, TP_fields{:}, 'degPerPix', degreesPerPixel, 'stats', stats, ...
            OSP_full_fields{:}, OSP_standard_fields{:}, extraOSP_fields{:}, bckgSamplesFields{:});
    end
    
    
%         if ~exist('bckgRate', 'var')
%             if isfield(cellData, 'PSTH') && isfield(cellData.PSTH, 'bckgRate')
%                 bckgRate = cellData.PSTH.bckgRate;
%             else
%                 bckgRate = [];
%             end
%         end

%         tfs_Hz = (1000/frameLength_ms)./tps_frm;
%         tf_sec = 1/tf_Hz;%tps_frm*(frameLength_ms/1000);
%         TP_fields = {'tmp_frm', tps_frm, 'tp_sec', tps_sec};

        % old:  stats = calcStatsFromOSPfull(R, R_full, Gid, bckgRate);            


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%      3. CALCULATE  STA(s)              %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if calculateSTAs  && between(cellId, 0, 100)
        % STA options:
        calcSTAfromOSP = 1;
        useSavedData = 1;
        useEvenOddSTAs = 1;
        checksToDo = []; %'odd_OSP';
        
        % MID options
        includeMID = 0;
        includePredictedOSP = 0;
        addFittedGaborStats = 1;
        includeMIDfit = 1  +0;
        

        % optimal-window STA
        [L_bin, R_bin] = dealV( PSTHdata.timeWindow_bins );
        windowProfile = PSTHdata.windowProfile;
        timeWindow_ms = PSTHdata.timeWindow_ms;
        STA = getSTA(Gid, cellId, L_bin, R_bin, timeWindow_ms, windowProfile, 1, calcSTAfromOSP, useSavedData, checksToDo, R_full);
        
        % fixed-window STA
        if nFixedTimeWindows > 0
            STAs_extra = cell(1, nFixedTimeWindows);
            for sta_i = 1:nFixedTimeWindows                
                [L_bin_i, R_bin_i] = dealV( allTimeWindows_bins(sta_i,:) );                
                STAs_extra{sta_i} = getSTA(Gid, cellId, L_bin_i, R_bin_i, allTimeWindows_ms(sta_i,:), [], useEvenOddSTAs, calcSTAfromOSP, useSavedData, checksToDo, R_full_extra{sta_i});
            end            
            fixedWindowSTA_fields = {'fixedWindowSTAs', struct('windows_bins', allTimeWindows_bins, 'windows_ms', allTimeWindows_ms, 'STAs', {STAs_extra})};
        else
            fixedWindowSTA_fields = {};
        end
        
        
%         if (nFixedTimeWindows > 0) && strcmp(gratingType, 'flashed');
% 
%             binE = binCent2edge(PSTHdata.bins)';
%             allTimeWindows_bins = binarySearch(binE, fixedTimeWindows);
%             
%             PSTH_data_extra = PSTHdata;
%             PSTH_data_extra.windowProfile = [];            
%             
%             [fixedTimeWindows_C, R_full_extra] = deal( cell(1, nFixedTimeWindows) );
%             keepAllTrials_flag_extra = 0; % just want odd/even trials;
%             for extra_i = 1:nFixedTimeWindows
%                 fixedTimeWindows_C{extra_i} = fixedTimeWindows(extra_i, :);
%                 PSTH_data_extra.timeWindow_bins = allTimeWindows_bins(extra_i,:);
%                 [~, R_full_extra{extra_i}] = getOriSpfPhaseProfile_simple(Gid, cellId, PSTH_data_extra, keepAllTrials_flag_extra);
%             end
%             3;
%             extraOSP_fields = {'fixedWindowSTAs', struct('extraOSP_windows', {fixedTimeWindows_C}, 'extraOSPs', {R_full_extra})};
%         else
%             extraOSP_fields = {};
%         end
        
        %%
%         STA = getSTA(Gid, cellId, L_bin, R_bin, timeWindow_ms, windowProfile, 1, calcSTAfromOSP, useSavedData, checksToDo, R_full);
        MID_fields = {};   
        fixedWindowMID_fields = {};
        if includeMID
            trialMode = iff(useEvenOddSTAs, 'odd/even', 'all');
            
            MID_struct = getMIDstruct(Gid, cellId, 'best', trialMode);
            MID_fields = {'MID_s', MID_struct};
                    
            if nFixedTimeWindows > 0
                MIDs_extra = cell(1, nFixedTimeWindows);
                for mid_i = 1:nFixedTimeWindows                    
                    MIDs_extra{mid_i} = getMIDstruct(Gid, cellId, allTimeWindows_ms(mid_i, :), trialMode);                    
                end
                fixedWindowMID_fields = {'fixedWindowMIDs', struct('windows_bins', allTimeWindows_bins, 'windows_ms', allTimeWindows_ms, 'MIDs', {MIDs_extra})};
            end
            
        end
        
            
        %%
        %{
        MID_orig = S.MID;
        MID_90 = rot90(MID_orig, -1);
        MID_180 = rot90(MID_orig, -2);
        MID_270 = rot90(MID_orig, -3);

        MID_flip = fliplr(S.MID);
        MID_f90 = rot90(MID_flip, -1);
        MID_f180 = rot90(MID_flip, -2);
        MID_f270 = rot90(MID_flip, -3);
                
        R_pred_orig = predictOSPfromMID(Gid, MID_orig, R);   R_pred_orig = R_pred_orig/max(R_pred_orig(:));
        R_pred_90 = predictOSPfromMID(Gid, MID_90, R);       R_pred_90 = R_pred_90/max(R_pred_90(:));
        R_pred_180 = predictOSPfromMID(Gid, MID_180, R);     R_pred_180 = R_pred_180/max(R_pred_180(:));
        R_pred_270 = predictOSPfromMID(Gid, MID_270, R);     R_pred_270 = R_pred_270/max(R_pred_270(:));
        R_pred_flip = predictOSPfromMID(Gid, MID_flip, R);   R_pred_flip = R_pred_flip/max(R_pred_flip(:));
        R_pred_f90 = predictOSPfromMID(Gid, MID_f90, R);     R_pred_f90 = R_pred_f90/max(R_pred_f90(:));
        R_pred_f180 = predictOSPfromMID(Gid, MID_f180, R);   R_pred_f180 = R_pred_f180/max(R_pred_f180(:));
        R_pred_f270 = predictOSPfromMID(Gid, MID_f270, R);   R_pred_f270 = R_pred_f270/max(R_pred_f270(:));

        i90 = 36/2; 
        R_orig = R_pred_orig;  R_orig = R_orig/max(R_orig(:));
        R_90   = R_orig([i90:end, 1:i90-1],:,:);
        R_180  = flipdim(R_orig,3);
        R_270  = flipdim(R_orig([i90:end, 1:i90-1],:,:), 3);
        R_flip = flipdim(R_orig,1);
        R_f90  = R_flip([i90:end, 1:i90-1],:,:);
        R_f180 = flipdim(R_flip,3);
        R_f270 = flipdim(R_flip([i90:end, 1:i90-1],:,:), 3);        
        
        figure(1); clf; imagesc(MID_orig); axis equal tight;
        figure(2); clf; imagesc(MID_90); axis equal tight;
        figure(3); clf; imagesc(MID_180); axis equal tight;
        figure(4); clf; imagesc(MID_270); axis equal tight;
        figure(5); clf; imagesc(MID_flip); axis equal tight;
        figure(6); clf; imagesc(MID_f90); axis equal tight;
        figure(7); clf; imagesc(MID_f180); axis equal tight;
        figure(8); clf; imagesc(MID_f270); axis equal tight;
                
        figure(11); clf; imageOSP(R_pred_orig, 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(12); clf; imageOSP(R_pred_90,   'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(13); clf; imageOSP(R_pred_180, 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(14); clf; imageOSP(R_pred_270, 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(15); clf; imageOSP(R_pred_flip, 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(16); clf; imageOSP(R_pred_f90,   'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(17); clf; imageOSP(R_pred_f180, 'mean:sp', 'OSP', 'nolabels');caxis([0 1])
        figure(18); clf; imageOSP(R_pred_f270, 'mean:sp', 'OSP', 'nolabels');caxis([0 1])
        
        figure(21); clf; imageOSP(R_orig, 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(22); clf; imageOSP(R_90,   'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(23); clf; imageOSP(R_180, 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(24); clf; imageOSP(R_270, 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(25); clf; imageOSP(R_flip, 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(26); clf; imageOSP(R_f90,   'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(27); clf; imageOSP(R_f180, 'mean:sp', 'OSP', 'nolabels');caxis([0 1])
        figure(28); clf; imageOSP(R_f270, 'mean:sp', 'OSP', 'nolabels');caxis([0 1])
        
        figure(31); clf; imageOSP(abs(R_orig-R_pred_orig), 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(32); clf; imageOSP(abs(R_90-R_pred_90),   'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(33); clf; imageOSP(abs(R_180-R_pred_180), 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(34); clf; imageOSP(abs(R_270-R_pred_270), 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(35); clf; imageOSP(abs(R_flip-R_pred_flip), 'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(36); clf; imageOSP(abs(R_f90-R_pred_f90),   'mean:sp', 'OSP', 'nolabels'); caxis([0 1])
        figure(37); clf; imageOSP(abs(R_f180-R_pred_f180), 'mean:sp', 'OSP', 'nolabels');caxis([0 1])
        figure(38); clf; imageOSP(abs(R_f270-R_pred_f270), 'mean:sp', 'OSP', 'nolabels');caxis([0 1])

        %}
    
        
        3;
        
%         [R_pred, plaw_coef] = predictOSPfromMID(Gid, STA, R, 80, cellData.OSP);
        [xs, ys] = getStimulusXY(Gid);
        3;
%         [R_pred, plaw_coef] = predictOSPfromMID(Gid, S.MID, R, 90, cellData.OSP);
        if includePredictedOSP && MID_available
            
            if useSavedData
                R_pred_C = getOspDataForPsthWindow(Gid, cellId, [], [], LR_bins(1), LR_bins(2), windowProfile, {'osp_pred'});
                R_pred = R_pred_C{1};
                nlinFunc = R_pred_C{2};
                R_pred_cc = R_pred_C{3};
            else
%                 [R_pred, nlinFunc, R_pred_cc] = predictOSPfromMID(Gid, STA, R, 80, cellData.OSP);
                [R_pred, nlinFunc, R_pred_cc] = predictOSPfromMID(Gid, S.MID, R, 90, cellData.OSP);
%                 [R_pred, ~, R_pred_cc] = predictOSPfromMID(Gid, S.MID, R);
            end
%             [R_pred, nlin, cc] = predictOSPfromMID(Gid, S.MID, R);
%             [R_pred, ~, R_pred_cc] = predictOSPfromMID(Gid, S.MID, R, 90, cellData.OSP);
            
            cellData.OSP.R_pred = single(R_pred);
            cellData.OSP.R_pred_nlinFunc = nlinFunc;
            cellData.OSP.R_pred_cc = R_pred_cc;
        end        
                
        cellData.STAs = struct('STA', STA, fixedWindowSTA_fields{:}, ...
                               MID_fields{:}, fixedWindowMID_fields{:}, ...
            'timeWindow_ms', timeWindow_ms, 'degPerPix', degreesPerPixel );
        3;
    end
%         figure(44); 
%         subplot(1,2,1); imagesc(STA);
%         subplot(1,2,2); imagesc(STA2);

    
    if showResults
        plotCellData(cellData, 300);
        pause(1);
         3;
    end
        
    
%     if nargout == 0         % make sure the data is saved in the main file
%         varname = getName('celldata', Gid, cellId);
%         if exist(varname, 'var')
%             save('allCellResults.mat',varname, '-append');
%         end
%     end
%     disp('[Completed Cell!]');

end

function STA = getSTA(Gid, cellId, L_bin, R_bin, timeWindow_ms, windowProfile, odd_even_flag, calcFromOSP_flag, useSavedData_flag, checkSTAs_flag, R)

    calcOddEven = exist('odd_even_flag', 'var') && (isequal(odd_even_flag, 1) || strcmp(odd_even_flag, 'odd/even') ); % default: no
    useSavedData = ~exist('useSavedData_flag', 'var') || ~isequal(useSavedData_flag, 0); % default: yes
    calcFromOSP = ~exist('calcFromOSP_flag', 'var') || ~isequal(calcFromOSP_flag, 0); % default: yes
    
    checkEvenOdd = exist('checkSTAs_flag', 'var') && ~isempty(strfind(checkSTAs_flag, 'odd'));
    checkOSPmethod = exist('checkSTAs_flag', 'var') && ~isempty(strfind(checkSTAs_flag, 'OSP'));
    
    if ~exist('R', 'var')
        R = [];
    end
    
    if ~calcOddEven
        STA_field = 'STA';
    else
        STA_field = 'STA_oe';
    end
    
    if useSavedData
        STA = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin, R_bin, windowProfile, {STA_field});

    else % calculate from scratch
        if calcFromOSP
            if ~calcOddEven
                if size(R,4) > 1,
                    R = mean(R, 4); % average over trials
                end
                STA = getSTAfromOSP(Gid, R);
            else
                if size(R,4) > 2
                    nTrials = size(R, 4);
                    R_odd = mean(R(:,:,:,1:2:nTrials), 4) ;
                    R_even = mean(R(:,:,:,2:2:nTrials), 4) ;
                else
                    R_odd = R(:,:,:,1);
                    R_even = R(:,:,:,2);
                end                
                STA_odd = getSTAfromOSP(Gid, R_odd);
                STA_even = getSTAfromOSP(Gid, R_even);
                STA = cat(3, STA_odd, STA_even);
            end
        else
            %         relContrOfFrameToSpike = getParsedSpikes( 'frame', Gid, cellId, timeWindow_ms, windowProfile);
            relContrOfFrameToSpike = getSpkResponseToEachFrame(Gid, cellId, timeWindow_ms, windowProfile);            
            if ~calcOddEven            
                STA = getSTAforCell(Gid, cellId, timeWindow_ms, windowProfile, relContrOfFrameToSpike);                
            else
                [~, ~, ~, ~, stim_RepeatId] = getStimulusFrameSequence(Gid, 'OSP');  idx_odd = odd(stim_RepeatId);                
                relContrOfFrameToSpike_odd = relContrOfFrameToSpike; relContrOfFrameToSpike_odd(~idx_odd) = 0;
                relContrOfFrameToSpike_even = relContrOfFrameToSpike; relContrOfFrameToSpike_even(idx_odd) = 0;
                STA_odd = getSTAforCell(Gid, cellId, timeWindow_ms, windowProfile, relContrOfFrameToSpike_odd);
                STA_even = getSTAforCell(Gid, cellId, timeWindow_ms, windowProfile, relContrOfFrameToSpike_even);
                STA = cat(3, STA_odd, STA_even);                
            end
               
        end
    end

    dontRecurseCheck = [];
    if calcOddEven && checkEvenOdd  % compare odd/even with all-trial
        STA_all = mean(STA, 3); 
        STA_all_chk = getSTA(Gid, cellId, L_bin, R_bin,  timeWindow_ms, windowProfile, 0,             calcFromOSP,    0, dontRecurseCheck, R);
        max_diff1 = max(abs(normVec(STA_all) - normVec(STA_all_chk)));
        assert(max_diff1 < 1e-3);        
    end
    if calcFromOSP && checkOSPmethod  % compare OSP with individual trial
%         STA_all = mean(STA, 3); % 
        STA_indiv_chk = getSTA(Gid, cellId, L_bin, R_bin, timeWindow_ms, windowProfile, calcOddEven, 0,               0,        dontRecurseCheck, R);
        max_diff2 = max(abs(normVec(STA) - normVec(STA_indiv_chk)));
        assert(max_diff2 < 1e-4);    
    end
    if (calcOddEven && checkEvenOdd) && (calcFromOSP && checkOSPmethod)
        STA_indiv_chk = getSTA(Gid, cellId, L_bin, R_bin, timeWindow_ms, windowProfile, 0,           0,               0,        dontRecurseCheck, R);
        max_diff3 = max(abs(normVec(mean(STA,3)) - normVec(STA_indiv_chk)));
        assert(max_diff3 < 1e-4);            
    end
    
    
end

function MID_struct = getMIDstruct(Gid, cellId, timeWindow, trialMode)
    includeGaborFitParams = 1;
    includeMIDfit = 1;
        
    
    if nargin < 3 || isempty(timeWindow)
        timeWindow = 'best';
    end
    
    if nargin < 4 || isempty(timeWindow)
        trialMode = 'all';
    end
    
    if strcmp(trialMode, 'all')
        allTrialModes = {'all'};
    elseif strcmp(trialMode, 'odd/even')
        allTrialModes = {'odd', 'even'};
    end
        
    MID_filenames = cellfun(@(tm) mid_getPreferredMIDfile(Gid, cellId, timeWindow, tm), allTrialModes, 'un', 0);
            
    MID_available = all(cellfun(@(fn) exist(fn, 'file'), MID_filenames));                        

    if MID_available
        %% 
        N_mid = length(MID_filenames);               
         
        for mid_i = 1:N_mid
            
            S_mid = load(MID_filenames{mid_i});
            MID(:, :, mid_i) = S_mid.MID; %#ok<AGROW>
            [M,N,nJacks] = size(S_mid.v_MID{1}); assert(M == N);
            MID_jacks = reshape(S_mid.v_MID{1}, [M*N, nJacks]);
            
            [~, pairwiseJackCCs] = pearsonRm(MID_jacks);
            meanJackCC(mid_i) = mean(pairwiseJackCCs); %#ok<AGROW>
            if includeGaborFitParams
                [gaborParams(mid_i,:), rsqr(mid_i), MID_fit(:,:,mid_i)] = mid_getCellGaborParams(Gid, cellId, timeWindow, allTrialModes{mid_i}); %#ok<AGROW>
            else
                [gaborParams, rsqr, MID_fit] = deal([]);
            end
                          
        end
            
        MID_struct = struct('RF', MID, 'timeWindow', timeWindow, 'jackCC', meanJackCC, 'gaborParams', gaborParams, 'rsqr', rsqr);
        if includeMIDfit
            MID_struct.MID_fit = MID_fit;
        end     
        %%
    else
        MID_struct = [];
    end

end

function v = normVec(X)
    v = X(:);
    L = lims(v);
    v = (v-L(1))/(L(2)-L(1));
end


%         if strncmp(osp_method, 'spike', 3)
%             [R, R_full, oris, sps, phs, tps_frm] = getOriSpfPhaseProfile_simple_spike(Gid, [relContrOfFrameToSpike{:}]);            
%         elseif strncmp(osp_method, 'stimulus', 3)
%             [R, R_full, oris, sps, phs, tps_frm] = getOriSpfPhaseProfile_simple_stim(Gid, cellId, PSTHdata);
%         end
