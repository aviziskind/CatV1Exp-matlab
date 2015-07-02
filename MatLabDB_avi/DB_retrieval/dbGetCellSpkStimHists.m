function [bins, allHistVals, meanRate, bckgSamples] = dbGetCellSpkStimHists(Gid, cellId, opt)

    persistent  allSpkStimHists saveCount;    
    global psthStatsSettings

    if isempty(psthStatsSettings)
        setPsthGlobalSettings;
    end
    
    redo_all = false; 
    redo_current_cell = 0;
    saveCountSpacing = 30;    
    if nargin < 3
        opt = struct;
    end
    
    bkgrSkip_ms      = psthStatsSettings.bkgrSkip_ms; % 250;    

    % valid fields for 'opt' struct:
    %  psthMethod, psthWindow_ms, trialGrouping, splitWindowIfCph, keepFirstDriftingGratingCycle,

    % determine whether flashed or drifting gratings.
    if ischar(Gid) && strcmp(Gid, 'save')
        if isempty(allSpkStimHists) || (saveCount == 0)
            return;
        end        
        %%
        fprintf('[saving Cell Hists ...\n '); tic;
        fn = fieldnames(allSpkStimHists);        
        for i = 1:length(fn)
            psthType = fn{i};
            histsFileName = [CatV1Path 'MatLabDB_avi' filesep 'allSpkStimHists_' psthType '.mat'];
            tmp = orderfields(allSpkStimHists.(psthType));          %#ok<NASGU>
            save(histsFileName, '-struct', 'tmp', '-v6');                            
            fprintf('   %s (%d entries) ... ', psthType, length(fieldnames(allSpkStimHists.(psthType)))); tic;        
        end
        t_sec = toc; fprintf('done (%.2f s)]\n', t_sec);                
        saveCount = 0;
        return;
    end

        
    assert( isnumeric(Gid) )
    gratingType = flashedOrDrifting(Gid, 'str');  

    % PSTH method:
    psthMethod = 'stimulus';
    if isfield(opt, 'psthMethod')
        psthMethod = opt.psthMethod;
    end
    psthMethod = psthMethod(1:2);
    histFunction = switchh(psthMethod, {'st', 'sp'}, {@getIndividualPSTHs, @getIndividualPSSHs});

    % PSTH window
    if isfield(opt, 'psthWindow_ms')
        psthWindow_ms = opt.psthWindow_ms;
    else
        switch gratingType, 
            case 'flashed',  psthWindow_ms = [-300, 200];
            case 'drifting', psthWindow_ms = [0, 8+1/3];
        end    
    end
    
    
    % Whether to store all trials, or just save even/odd trials     
    groupingType = curGroupingType('');    
    switch groupingType
        case {'clustersPruned'},   % we never care about individual trials for this clustering
            trialGrouping_calc = 'odd/even';    
            keepAllBckgSamples = false;            
        case {'clusters', 'cells'},  % need for these clusterings                        
            trialGrouping_calc = 'individual';    
            keepAllBckgSamples = true;
    end    
    trialGrouping_calc = lower(trialGrouping_calc(1:3));        
    if ~any(strcmp(trialGrouping_calc, {'odd', 'ind'}))
        error('Invalid "trialGrouping" option');
    end    
    
    keepSpikeFirstFrameId = strcmp(gratingType, 'drifting') && strcmp(groupingType, 'clusters');
    
    trialGrouping_return = trialGrouping_calc;
    if isfield(opt, 'trialGrouping')
        trialGrouping_return = opt.trialGrouping;
    end
    
    
    % Whether to split the psth window if is counterphase flashed grating
    splitWindowIfCph = false; % default value
    if isfield(opt, splitWindowIfCph)
        splitWindowIfCph = opt.splitWindowIfCph;
    end
    cphSplitTime_ms  = psthStatsSettings.cphSplitTime_ms; % currently 70 ms;

    
    % Whether to keep first cycle of drifting gratings
    keepFirstDriftingGratingCycle = false;
    if isfield(opt, 'keepFirstDriftingGratingCycle')
        keepFirstDriftingGratingCycle = opt.keepFirstDriftingGratingCycle;
    end
            
    name_fun = @(Gid, cellId) sprintf('Gid_%04d_cell_%s_%s_%s', Gid, iff(cellId == -1, 'n1', num2str(cellId)), lower(psthMethod(1:2)), lower(trialGrouping_calc(1:3)) );
        
    % if not loaded, load/create it     
    if isempty(allSpkStimHists)
        allSpkStimHists = struct;
    end    
    fet_str = iff(curMatchDB, '_DB', ['_' curSortingFeatures('')]);
    
    psthType = [groupingType fet_str '_' gratingType(1) '_' psthMethod '_' trialGrouping_calc];
    
    histsFileName = [CatV1Path 'MatLabDB_avi' filesep 'allSpkStimHists_' psthType '.mat'];

    if ischar(Gid) && strcmp(Gid, 'save') && ~isempty(allSpkStimHists.(psthType))
        tmp = orderfields(allSpkStimHists.(psthType));         %#ok<NASGU>
        fprintf('\n[saving cell stim hists (%s) ... ', psthType); tic;
        save(histsFileName, '-struct', 'tmp', '-v6');                
        t_sec = toc; fprintf('done (%.2f s)]\n', t_sec);        
        saveCount = 0;        
        return;
    end    
    
    
    if ~isfield(allSpkStimHists, psthType) || redo_all        
        if ~exist(histsFileName, 'file') || redo_all
            allSpkStimHists.(psthType) = struct;
        else    
            method_str = iff(strcmp(psthMethod, 'st'), 'stimulus', 'spike');
            grp_str = iff(strcmp(trialGrouping_calc, 'odd'), 'odd/even', 'individual');            
            fprintf('Loading saved %s %s grating %s-trial psth''s using %s method  ... ', groupingType, upper(gratingType), upper(grp_str), upper(method_str)); tic; 
            allSpkStimHists.(psthType) = load(histsFileName);
            fprintf(' done. '); toc;
        end
    end

    if isempty(saveCount)
        saveCount = 0;
    end
    
    % if doesn't exist in loaded file, calculate it
    if ~exist('Gid', 'var') || isempty(Gid)  % do for all gids/cellids in flashed grating cells file
        
%         S1 = load('cellsGroups_movie_fg.mat');
%         st = {S1.movieGroups_fg.stimType};
%         grp_idx = cellfun(@(s) strncmp(s, 'Movie:Flashed_Gratings:36x10x8(4x1)C', 35), st);
%         allXGids = [S1.movieGroups_fg(grp_idx).Gid];            
% 
%         S = load('flashedGratingCells_all.mat');        
%         idx_select = find( arrayfun(@(g_id) any(allXGids == g_id), [S.allCells.Gid]) );        
%         Gid = [S.allCells(idx_select).Gid];
%         cellId = [S.allCells(idx_select).cellId];            
        fd = curGratingType;

        if fd == 1
            S = load('flashedGratingCells_all.mat');        
        else
            S = load('driftingGratingCells_all.mat');                    
        end
        Gid = [S.allCells.Gid];
        cellId = [S.allCells.cellId];                    

    end

    if length(Gid) > 1        
        progressBar('init-', length(Gid), 60);        
    end
    
    saveFile = false;    
    for i = 1:length(Gid)      
        
        if length(Gid) > 1
            progressBar;
        end
        psth_name = name_fun(Gid(i), cellId(i));
           
        fieldExists = isfield(allSpkStimHists.(psthType), psth_name);        
        if fieldExists && ~curMatchDB && strcmp(groupingType, 'cells') 
            curClustIds = getClusterCellAssignment(Gid, cellId);
            savedClustIds = allSpkStimHists.(psthType).(psth_name).clustIds;
            if ~isequal(savedClustIds, curClustIds)
                redo_current_cell = 1;
            end
        end
        
        if fieldExists  && anyNanInHists(allSpkStimHists.(psthType).(psth_name).histVals);
            redo_current_cell = 1;
        end
        
        if (~fieldExists || redo_current_cell) && (cellId ~= 100)...
                
            [bins, allHistVals_orig, meanRate, bckgSamples_orig] = histFunction(Gid(i), cellId(i), 'OSP', psthWindow_ms, bkgrSkip_ms, trialGrouping_calc);
                                    
            histVals_S = compress(allHistVals_orig); % compress vals (to conserve space);
%             assert(isequal(decompress(histVals_S), allHistVals_orig));
            if keepAllBckgSamples
                bckgSamples = compress(bckgSamples_orig); % compress vals (to conserve space);
            else
                bckgSamples = [mean(bckgSamples_orig), std(bckgSamples_orig)];
            end           
            
            s = struct('bins', bins, 'histVals', histVals_S, 'meanRate', meanRate, 'bckgSamples', bckgSamples);
            
            if ~curMatchDB && strcmp(groupingType, 'cells') 
                curClustIds = getClusterCellAssignment(Gid, cellId);
                s.clustIds = curClustIds;
            end
            allSpkStimHists.(psthType).(psth_name) = s;
            saveCount = saveCount + 1;   
            saveFile = true;
        end        
    end
    if length(Gid) > 1
        progressBar('done')
    end
    if saveFile && (saveCount > saveCountSpacing) 
        %%
        tmp = orderfields(allSpkStimHists.(psthType));         %#ok<NASGU>

        fprintf('\n[saving cell stim hists ... '); tic;
        save(histsFileName, '-struct', 'tmp', '-v6');                
        t_sec = toc; fprintf('done (%.2f s)\n', t_sec);        
        saveCount = 0;
    end                
    
        
    if length(Gid) == 1
        
%         curClustIds = getClusterCellAssignment(Gid, cellId);
%         allSpkStimHists.(psthType).(psth_name).clustIds = curClustIds;
        stimType = getGratingStimType(Gid);        
        sd = siteDataFor('Gid', Gid, 1);
        
        if (cellId ~= 100)
        
            psth_name = name_fun(Gid, cellId);
            s = allSpkStimHists.(psthType).(psth_name);        
            [bins, histVals_S, meanRate, bckgSamples_S] = deal(...
                s.bins, s.histVals, s.meanRate, s.bckgSamples);

            allHistVals = decompress( histVals_S );   
            
            if nargout >= 4
                if isstruct( bckgSamples_S )                    
                    bckgSamples = decompress( bckgSamples_S );
                elseif isnumeric( bckgSamples_S )
                    bckgSamples = bckgSamples_S;
                end                    
            end
            
        elseif (cellId == 100)
            if ~curMatchDB                
                S_sorting = load(getFileName(curGroupingType(''), Gid));
                allCellIds = S_sorting.uClustIds;
%                 allCellIds = allCellIds(allCellIds >= 0);
            else
                sd = siteDataFor('Gid', Gid, 1);
                allCellIds = sd.cellIds;
            end
%             allCellIds = allCellIds(allCellIds >= 0);         

            opts = struct('splitWindowIfCph', 0, 'keepFirstDriftingGratingCycle', 1);
            for ci = 1:length(allCellIds)
                if ci == 1
                    [bins, allHistVals_ci, meanRate_ci, bckgSamples] = dbGetCellSpkStimHists(Gid, allCellIds(ci), opts);
                    allHistVals = allHistVals_ci;
                    meanRate = meanRate_ci;                    
                else
                    [~, allHistVals_ci, meanRate_ci] = dbGetCellSpkStimHists(Gid, allCellIds(ci), opts);
                    allHistVals = allHistVals + allHistVals_ci;
                    meanRate = meanRate + meanRate_ci;                    
                end                                
            end                        
            
        end        
        
        if length(psthWindow_ms) == 1
            psthWindow_ms = [0 psthWindow_ms];
        end                       
        
        bin_idx = find( ibetween(bins, psthWindow_ms) );         
        if ~isequal(bin_idx(:)', [1:length(bins)])
            bins = bins(bin_idx);
            allHistVals = allHistVals(bin_idx,:,:);
        end        
        
        % remove first cycle of drifting grating response (with artifacts);
        if strcmp(gratingType, 'drifting') && strncmp(trialGrouping_calc, 'individual', 3) && ...
            ~keepFirstDriftingGratingCycle
            allHistVals = removeFirstCycleOfDriftingGratingResponse(Gid, allHistVals);            
        end
        
        % for counter-phase flashed gratings, estimating which bins are reproducible requires 
        % splitting the window like this, to remove the residual spiking effect of the 
        % previous/subsequent counterphase frame.
        isCphFlashed = stimType.isCphFlashed;
        if isCphFlashed && splitWindowIfCph
            bins_before =  ( bins <  cphSplitTime_ms );
            bins_after  =  ( bins >= cphSplitTime_ms );                
            allHistVals = cat(1, allHistVals(bins_before, :, [1 2], :), ...
                                 allHistVals(bins_after,  :, [3 4], :) );
        end
            
        nTrials = size(allHistVals, 3);
        if (~isCphFlashed || ~splitWindowIfCph) && strncmp(trialGrouping_return, 'odd/even', 3) && (nTrials > 2)                         
            idx_odd = 1:2:nTrials;
            idx_even = 2:2:nTrials;
            allHistVals = cat(3, mean(allHistVals(:,:,idx_odd),  3), ...
                                 mean(allHistVals(:,:,idx_even), 3) );
        end
            
        
    end
end
    


%     if ~exist('trialGrouping', 'var') || isempty(trialGrouping)
%         switch gratingType
%             case 'flashed',  trialGrouping = 'individual'; % 'odd/even' or 'individual'
%             case 'drifting', trialGrouping = 'individual'; % 'odd/even' or 'individual'
%         end            
%     end

function allPsthVals = removeFirstCycleOfDriftingGratingResponse(Gid, allPsthVals)
    % for drifting gratings, we want to discard first cycle which may contain artifacts (of sudden contrast change from mean luminance)
    % Additionally, the last frame (presenting the last phase) is occasionally truncated (was not shown
    % in the experiment). In these rare cases, to avoid having a nan in R_full (and instead of just putting a zero), we'll paste in the value of
    % the last phase from the first cycle (usually the artifacts will be at the beginning/middle of
    % the first cycle, so the end of the first cycle should be ok).

    gratingType = flashedOrDrifting(Gid, 'str');
    groupingType = curGroupingType('');
    assert(strcmp(gratingType, 'drifting'));
    [nOri, nSp, nPh, nCycles, nRep] = dbGetUniqueOriSpPh('Gid', Gid, 'length');
        
    if nCycles ~= round(nCycles)
        error('not implemented yet');
    end
    nCycles = ceil(nCycles); % for case in which is a fraction
        
    justOddEven = ~strcmp(groupingType, 'cells') && size(allPsthVals, 3) == 2;
    if justOddEven
        nCycles = 1;
        nRep = 2;
    end
%         nCycTotal = nCycles*nRep;

%         R_full = reshape(allPsthVals, [1, nOri*nSp*nPh, nCycTotal]);
    % compress from all trials to even/odd                           
    OSP_full = reshape(allPsthVals, [nOri, nSp, nPh, nCycles, nRep]);

    if (any(isnan(OSP_full(:))))            
        % find index of first nan in phase.                        
        error('Shouldn''t have this any more')
        OSP_full_ph_last = permute(OSP_full, [1, 2, 4, 5, 3]); % put in phase last
        OSP_full_vs_ph = reshape(OSP_full_ph_last, [nOri*nSp*nCycles*nRep, nPh]);
        idxFirstNan = find(any(isnan(OSP_full_vs_ph), 1), 1, 'first');
        idxsToReplace = idxFirstNan : nPh;

        % replace (missing) end of last cycle with end of first cycle (which will be discarded anyway).
        OSP_full(:, :, idxsToReplace, nCycles, :) = OSP_full(:, :, idxsToReplace, 1, :);                        
    end
    assert(~any(isnan(OSP_full(:))));

    % discard the first cycle of each repetition.
    okCycleIdxs = 2:nCycles;        
    OSP_full = OSP_full(:,:,:, okCycleIdxs, :);

    % reshape back into [1 x nStim x nTrials],  
    %    (the first dimension is the bin #, which is always of length 1 for drifting gratings)
    % then separate out even & odd trials
    nTrialsUsed = nRep*(nCycles-1);
    allPsthVals = reshape(OSP_full, [1, nOri*nSp*nPh, nTrialsUsed]);        
    
end


function tf = anyNanInHists(histVals)
    if isstruct(histVals)
        if isfield(histVals, 'uVals')
            tf = any(isnan(histVals.uVals(:)));
        elseif isfield(histVals, 'orig_vals')
            tf = any(isnan(histVals.orig_vals(:)));
        end
    else
        tf = any(isnan(histVals(:)));
    end
end


%             psth_name = name_fun(Gid, cellId);
%             s = allSpkStimHists.(psthType).(psth_name);        
%             [bins, histVals_S, meanRate2] = deal(...
%                 s.bins, s.histVals, s.meanRate);
% 
%             allHistVals2 = decompress( histVals_S );        
%             assert(isequal(allHistVals, allHistVals2));
%             assert(abs(meanRate-meanRate2)<1e-5);
