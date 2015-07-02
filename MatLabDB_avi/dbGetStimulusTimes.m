function R_full_t = dbGetStimulusTimes(Gid, opt, redo_flag)
    
    persistent allStimulusTimes saveCount
    
    if nargin < 2 || isempty(opt)
        opt = struct;
    end
    skipFirstCycleOfDriftingGratings = true && ~(isfield(opt, 'skipFirstCycleOfDriftingGratings') && isequal(opt.skipFirstCycleOfDriftingGratings, 0));
        
    stimulusTimes_file = [CatV1Path 'MatLabDB_avi' filesep 'allStimulusTimes.mat'];
    
    redo_all = 0;
    redo_current = 0 || (exist('redo_flag', 'var') && isequal(redo_flag, 1));
    saveCountSpacing = 50;        
    
    if isempty(allStimulusTimes)        
        if exist(stimulusTimes_file, 'file') && ~redo_all
            S_file = load(stimulusTimes_file);
            allStimulusTimes = S_file.allStimulusTimes;
        else
            allStimulusTimes = struct;
        end        
        saveCount = 0;        
    end

    if strcmp(Gid, 'save')
        save(stimulusTimes_file, 'allStimulusTimes', '-v6');                
        saveCount = 0;
        return;
    end        
    
    group_fld_name = sprintf('StimulusTimes_Gid_%d', Gid);            
    
%     group_fld_name_err = [group_fld_name '_err'];
    
    if (~isfield(allStimulusTimes, group_fld_name) || redo_current)
        R_full_t = calcStimulusTimes(Gid);

        allStimulusTimes.(group_fld_name) = R_full_t;
        saveCount = saveCount + 1;

        if saveCount > saveCountSpacing
            save(stimulusTimes_file, 'allStimulusTimes', '-v6');                
            saveCount = 0;
        end                
    end
        
    R_full_t = allStimulusTimes.(group_fld_name);
    
    if flashedOrDrifting(Gid) == 2 && skipFirstCycleOfDriftingGratings
%         stimType = getGratingStimType(Gid);
%         nOri, nSpf, nPh, nCycles, nRep]
        [nOri, nSpf, nPh, nCycles, nRep] = dbGetUniqueOriSpPh('Gid', Gid, 'length');

        % reshape so that cycles and repeats are in separate dimensions
        R_full_t = reshape(R_full_t, [nOri, nSpf, nPh, nCycles, nRep]);
        
        % discard the first cycle of each repetition.
        R_full_t = R_full_t(:,:,:, 2:nCycles, :);

        % reshape back 
        nTrialsUsed = nRep*(nCycles-1);
        R_full_t = reshape(R_full_t, [nOri, nSpf, nPh, nTrialsUsed]);        
    end


end


function [R_full_t, frameStartTimes] = calcStimulusTimes(Gid)
    
    [frameStimIds, ~, ~, ~, stim_RepeatId] = getStimulusFrameSequence(Gid, 'OSP', 1);
%     [frameStimIds, ~, ~, ~, stim_RepeatId_orig] = getStimulusFrameSequence(Gid, 'OSP');
        
    [uStimIds, stimIdLists] = uniqueList(frameStimIds);
    nStim = length(uStimIds); 
    sd = siteDataFor(Gid);
    nPres = length(sd.presIds);
    
    nStimTot = length(frameStimIds);
    nStimsPerPres = nStimTot/nPres;    
        
    [nOri, nSpf, nPh] = deal(length(sd.ori_deg), length(sd.spPeriod_pix), length(sd.spPh_deg));
    assert(nStim == nOri*nSpf*nPh);
    nRepTot = nStimTot/nStim;
    if (nRepTot ~= round(nRepTot))
        error('uneven: nRep = %.5f', nRepTot);
    end
    
    syncTimes0 = dbGetSyncs('Gid', Gid, 'sec');
    assert(length(syncTimes0) == (nStimsPerPres+1)*nPres);    
    syncTimes_sec = [0; syncTimes0];
        
    cumulativeFrames = [0, [1:nPres]*nStimsPerPres];    
%     frameIdForIntervalId  = zeros(1, length(syncTimes_ms)   );  

    frameStartTimes = zeros(1, length(frameStimIds));
    for iPres = 1:nPres
        offset = iPres;
        frame_idxs = cumulativeFrames(iPres)+1:cumulativeFrames(iPres+1);

        frameStartTimes(frame_idxs) = syncTimes_sec(frame_idxs+offset);        
    end    
    
    
    R_full_t = zeros(nOri*nSpf*nPh, nRepTot, 'single');
    for stim_i = 1:nStim
        idx_order = stim_RepeatId( stimIdLists{stim_i} );
        % this is important for the counter-phase flashed gratings
        % which are not stored in the original order --> we want
        % R_full_t to match the times of R_full, so we have to use the
        % correct idx_order (for other stimuli, this has no effect).

        R_full_t(stim_i,:) = frameStartTimes( stimIdLists{stim_i}(idx_order) );
    end
    R_full_t = reshape(R_full_t, [nOri, nSpf, nPh, nRepTot]);
    3;

end

%{
% allGids = setdiff(allGids, [513, 605, 607, 608 610 612, 615 616 619 620]);
allGids = getAllGids;
progressBar('init-', length(allGids), 60);
for i = 1:length(allGids)
    dbGetStimulusTimes(allGids(i));
    progressBar(i);
end
dbGetStimulusTimes('save');

%}