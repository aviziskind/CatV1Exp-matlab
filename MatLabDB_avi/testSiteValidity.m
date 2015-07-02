function [siteOk, msg, presOk] = testSiteValidity(Gid, exitFlag, redoFlag) 

    persistent allSiteTestData;
    maxGid = 5336;
%     allSiteTestData = [];
    
    redoNow = exist('redoFlag', 'var') && ~isempty(redoFlag) && redoFlag;
    
    siteTestFileName = [CatV1Path 'MatLabDB_avi' filesep 'siteTestData.mat'];
    if isempty(allSiteTestData)        
        if exist(siteTestFileName, 'file')
            S = load(siteTestFileName);
            allSiteTestData = S.allSiteTestData;
        else
            allSiteTestData = cell(maxGid, 1);
        end
    end
    siteTestData = allSiteTestData{Gid};
    calcNow = isempty(siteTestData) || redoNow;
            
    if calcNow
        [siteOk, msg, presOk] = testOneSite(Gid, 'details');        
        allSiteTestData(Gid) = {{siteOk, msg, presOk}};
        save(siteTestFileName, 'allSiteTestData');
    else
        siteData = allSiteTestData{Gid};
        [siteOk, msg, presOk] = deal(siteData{:});
    end

    
end



function [siteOk, msg, presOk] = testOneSite(Gid, exitFlag)
	acceptableErrorIds = {'db:noSid', 'db:noSyncs', 'db:badSyncs', 'db:flawedSyncs', 'db:invalidTbTe', ...
                          'db:syncsTbTeMismatch', 'db:syncNframesMismatch', ...
                          'stim:interrupted', 'stim:singleTrial', 'db:SinglePhaseMovie' }; % 
    
    global nFramesCmp;
    if nargin < 2
        exitFlag = 'report';
    end
    checkAllPres = (nargout > 2);
       
    switch exitFlag
        case 'error', exitWithError = true; fullDetails = false; % 
        case 'report', exitWithError = false; fullDetails = false;
        case 'details', exitWithError = false; fullDetails = true;
    end
    
    hnd = dbOpenExpDb;
    Did = dbLookup('Did',  'Gid', Gid);
    [stimulusType, stimulusSubType] = getStimulusTypeForDid(Did);
    stimTable = getDatabaseTableForDid(Did, stimulusType);
    stimPresIdField = [upper(stimulusType) '_PRES_ID'];
    presIds = getFieldsFromDatabaseTable(hnd, stimPresIdField, stimTable, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
    nPres = length(presIds);    
    addIdStr = false;
    presOk = false(1, nPres);
    
    idStr = sprintf(' {Did:%4d; Gid:%4d}%s', Did, Gid, iff(isempty(stimulusSubType), '', ['[' stimulusSubType ']']) );
    pct_s = '%%';
    
    if exitWithError
        [siteOk, msg] = doSiteTest(Gid, Did);
    else
        
        try
            [siteOk, msg] = doSiteTest(Gid, Did);
        catch Merr        
            if exitWithError      || ~any(strcmp(Merr.identifier, acceptableErrorIds))
                rethrow(Merr);
            else            
                [siteOk, msg] = deal(Merr.identifier, Merr.message);
            end
        end
        
    end
    idStr = iff(addIdStr, idStr, '');
    msg = [msg  idStr];
    if exitWithError && ~strcmp(siteOk, 'ok')
        error(siteOk, msg);
    end

    function fprintf2(varargin)
%         fprintf(varargin{:})
    end

    
    function [siteOk, msg] = doSiteTest(Gid, Did)

        
        errorForFgCompExps = true;
        [siteOk, msg] = deal([], []);

        brk = ' || ';
    
        % 1. Get spikes (will error if a problem with spikes) 
        allSpikes_ticks = dbGetSpikes(Gid, -1);
        
        % 2. Get Syncs (will error if a problem with syncs)
        syncTimes_ticks = dbGetSyncs('Gid', Gid);

        addIdStr = true;
        presOk = true(1, nPres);

        % 3. Get Database fields
        [nDisplayedFramesDB, nSustainedFrames, nPreBlankFrames, nPostBlankFrames] = getFieldsFromDatabaseTable(hnd, ...
            {'LNG_N_DISPLAYED_FRM', 'LNG_N_SUSTAINED_FRM', 'LNG_N_PRE_BLANK_FRM', 'LNG_N_POST_BLANK_FRM'}, stimTable, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
        if strcmpi(stimulusType, 'Grating')
            [nFadeInFrames, nFadeOutFrames] = getFieldsFromDatabaseTable(hnd, {'LNG_N_FADEIN_FRM', 'LNG_N_FADEOUT_FRM'}, stimTable, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
        else
            [nFadeInFrames, nFadeOutFrames] = deal(0);
        end

        
        % 3. Parse Syncs - check if presentations are missed / frames are unaccounted for 
        lastTick = allSpikes_ticks(end); % We assume that the recording ends around the last spike.  (We don't know exactly when the recording ends.)
        lastTick = max([lastTick, syncTimes_ticks(end)+dbConvertTimeMeasures(Did, 5, 'frame', 'tick')]); % In case no post-blank spikes, make sure there is a tick after the last frame, after the post blank frames
        syncTimes_ticks0 = [0; syncTimes_ticks; lastTick]; % add my own ticks at start & end of experiment.

        [nFramesPerPres_Runs, parseSyncsEndStatus, tmp, runs_presOk] = dbParseSyncs(Did, syncTimes_ticks0);
        switch parseSyncsEndStatus
            case -1, [siteOk, curmsg] = deal('db:badSyncs', 'Couldn''t parse syncs.');
            case 2,  [siteOk, curmsg] = deal('db:flawedSyncs', 'Parsed syncs, but missing some frames.');
            case 1,  curmsg = 'Syncs a little uneven, but complete.';
        end
        if parseSyncsEndStatus ~= 0
            msg = [msg, iff(~isempty(msg), brk), curmsg];
        end
        if parseSyncsEndStatus == -1
            presOk = false(1,nPres);
        else
            presOk = presOk & runs_presOk;
        end

        if fullDetails && (parseSyncsEndStatus~=0)
            fprintf2(['ParsedSyncStatus: ' curmsg '\n']);
            if length(nFramesPerPres_Runs) == 1
                fprintf2('Frames in pres : %d\n', nFramesPerPres_Runs)
            else
                m = median(nFramesPerPres_Runs);
                i = find(m ~= nFramesPerPres_Runs);
                if ~isempty(i)
                    fprintf2('Most presentations: %d frames.  %s\n', m, ...
                        sprintf('[Pres %d : %d frames]', i, nFramesPerPres_Runs(i)) );
                end
            end
        end
        if ~isempty(siteOk) && ~fullDetails && ~checkAllPres,   error(siteOk, msg);  end


        % 4. Check start/end ticks
        % 4a. Check that are valid (if syncs ok, then correctable)
        [beginTime_ticks, endTime_ticks] = dbGetTbTe(hnd, Did);
        presStartSyncIds = binarySearch(syncTimes_ticks0, beginTime_ticks);
        presEndSyncIds = binarySearch(syncTimes_ticks0, endTime_ticks);            
        if any(isnan([beginTime_ticks + endTime_ticks])) || isempty(beginTime_ticks) || isempty(endTime_ticks) %% maybe can get around this by using first/last syncs instead
            [siteOk, curmsg] = deal('db:invalidTbTe', 'Invalid begin and end ticks for experiment');
            msg = iff(isempty(msg), curmsg, [msg brk curmsg]);                    
        end        
        if ~isempty(siteOk) && ~fullDetails && ~checkAllPres,   error(siteOk, msg);  end

        % 4b. Check that Start/End Ticks align with syncs
        if any( syncTimes_ticks0(presStartSyncIds) ~= beginTime_ticks) || any( syncTimes_ticks0(presEndSyncIds) ~= endTime_ticks)
            [siteOk, curmsg] = deal('db:syncsTbTeMismatch', ['Begin/End Ticks do not align with syncs']);
            msg = iff(isempty(msg), curmsg, [msg brk curmsg]);        
            if fullDetails
                for i = find( syncTimes_ticks0(presStartSyncIds) ~= beginTime_ticks)';
                    fprintf2('[BeginTick.Vs.Sync] Begin tick #%d (%d) does not match closest sync (%d) [difference of %d]\n', ...
                        i, beginTime_ticks(i), syncTimes_ticks0(presStartSyncIds(i)), beginTime_ticks(i)-syncTimes_ticks0(presStartSyncIds(i)) );
                end                
                for i = find( syncTimes_ticks0(presEndSyncIds) ~= endTime_ticks)';
                    fprintf2('[EndTick.Vs.Sync] End tick #%d (%d) does not match closest sync (%d) [difference of %d]\n', ...
                        i, endTime_ticks(i), syncTimes_ticks0(presEndSyncIds(i)), endTime_ticks(i)-syncTimes_ticks0(presEndSyncIds(i)) );
                end                
            end
        end
        if ~isempty(siteOk) && ~fullDetails && ~checkAllPres,   error(siteOk, msg);  end

        
        % 5. Check that the number of frames from syncs matches the number of frames from DB fields 
        nFrames_fromDB    = nDisplayedFramesDB-nPreBlankFrames-nPostBlankFrames;  nTotFrames_fromDB = sum(nFrames_fromDB);    
        nFrames_fromSyncs = presEndSyncIds - presStartSyncIds;             nTotFrames_fromSyncs = sum(nFrames_fromSyncs);    
        nFramesDifference = nFrames_fromSyncs - nFrames_fromDB;
        [syncDB_discrepancy, indpres] = max(abs(nFramesDifference));

        syncDB_discrepancy_pct = syncDB_discrepancy / median(nFrames_fromDB) * 100;
                
        maxAllowedDiscrepancy = 0;         
        if (syncDB_discrepancy > 0) 
            if fullDetails            
                for i = find(nFramesDifference(:)' > 0)            
                    fprintf2('[sync.VS.db]  Pres# %d : NFrames_Syncs : %d. NFrames_DB  %d. Diff of %d frames (%.3f %%)]\n', ...
                        i, nFrames_fromSyncs(i), nFrames_fromDB(i), nFramesDifference(i), nFramesDifference(i)/nFrames_fromSyncs(i)*100);
                end                    
            end        

            cmp = [Gid, nTotFrames_fromSyncs, nTotFrames_fromDB, nTotFrames_fromSyncs-nTotFrames_fromDB syncDB_discrepancy];
            nFramesCmp = [nFramesCmp; cmp];
            [siteOk, curmsg] = deal('db:syncNframesMismatch', ...
                sprintf('[sync.VS.db] Nframes_from_syncs does not match Nframes_from_DB in %d (out of %d) pres''s. Max diff is at pres # %d: %d frames (%.3f %%)', ...
                    nnz(nFramesDifference), nPres, indpres, syncDB_discrepancy, syncDB_discrepancy_pct ));
            msg = iff(isempty(msg), [curmsg], [msg brk curmsg]);                
        end                
        presOk = presOk & (nFramesDifference(:)' <= maxAllowedDiscrepancy);
        
        if ~isempty(siteOk) && ~fullDetails && ~checkAllPres,   error(siteOk, msg);  end
        
        % 6. Check the consistency of the number of presentations recorded in the database (correctable)
        nSupposedToShowFrames = nSustainedFrames + nFadeInFrames + nFadeOutFrames ;
        nShowedFrames         = nDisplayedFramesDB - nPreBlankFrames - nPostBlankFrames;
        nFramesDifference = nShowedFrames - nSupposedToShowFrames;    
        [nFramesDiscrep, indpres] = max(abs(nFramesDifference));
        
        [gratingMaxNDiscrepPct, noiseMaxNDiscrepPct, movieMaxNDiscrepPct, mseqMaxNDiscrepPct] = ...
            deal(5,             30,                     5,                  0);  %#ok<NASGU> % can have a bit more leeway with the Grating stimuli b/c same stimulus is repeated across many frames.
                
        maxAllowedDiscrepPct = eval([lower(stimulusType) 'MaxNDiscrepPct']);
        
        if (nFramesDiscrep > 1)
            if fullDetails            
                for i = find(nFramesDifference(:)' > 0)
                    fprintf2('[sust.VS.disp] Pres# %d : planned: %d. shown: %d. Diff of %d frames (%.3f %%)]\n', ...
                        i, nSupposedToShowFrames(i), nShowedFrames(i), nFramesDifference(i), nFramesDifference(i)/nSupposedToShowFrames(i)*100);
                end                    
            end        
            pct = nFramesDiscrep/nSupposedToShowFrames(indpres)*100;
            msg = [msg, sprintf('[sust.VS.disp][%d (out of %d) pres''s had frame discreps. max diff was in pres# %d (planned: %d. shown: %d. Diff of %d frames (%.3f%s)', ...
                                    nnz(nFramesDifference), nPres, indpres, nSupposedToShowFrames(indpres), nShowedFrames(indpres), nFramesDifference(indpres), pct, pct_s)];
            if (pct > maxAllowedDiscrepPct) 
                siteOk = 'stim:interrupted';
            end
        end
        presOk = presOk & (nFramesDifference(:)' <= maxAllowedDiscrepPct);
        
        if ~isempty(siteOk) && ~fullDetails && ~checkAllPres,   error(siteOk, msg);  end
        
        % 7. Check that Syncs are evenly spaced within each presentation (more important for drifting gratings). 
    %     checkEvenlySpacedFrames(syncTimes_ms, presStartSyncIds, presEndSyncIds, th);
        maxStdDevDriftingGratings = .05;  % more important that drifting gratings be regular.
        maxStdDevOtherStimuli = .2;

        syncTimes_ms = dbConvertTimeMeasures(Did, syncTimes_ticks0, 'tick', 'ms');
        
        maxStdDevTh = iff( (strcmpi(stimulusType, 'Grating') && ~strcmp(stimulusSubType, 'Flashed Grating')), ...
            maxStdDevDriftingGratings, maxStdDevOtherStimuli);
        nUnevenFramesPctTh = 1;
        
    %     frameLength_ms = getFrameLength('Did', Did);
        stdevs = zeros(1,nPres);
        ndevs    = zeros(1,nPres);
        ndevs_pct = zeros(1,nPres);
        for pres_i = 1:nPres
            syncsForPresI = syncTimes_ms(presStartSyncIds(pres_i):presEndSyncIds(pres_i));
            Ls = diff(syncsForPresI);
            m = median(Ls); s = std(Ls);
            ndevs(pres_i) = nnz( abs(Ls-m)/m > maxStdDevTh );        
            ndevs_pct(pres_i) = ndevs(pres_i) / nShowedFrames(pres_i) * 100;
            stdevs(pres_i) = s/m;
            [maxStdDev, indMaxStdDev] = max(stdevs);
    %         [maxNDev, indMaxNDev] = max(ndevs);
        end        
        presOk = presOk & (ndevs_pct < nUnevenFramesPctTh);
                
        if fullDetails
            for i = find(stdevs > maxStdDevTh)            
                fprintf2('[UnevenFrames] In Pres# %d, there are %d frames (%.3f %%) not of correct length. StdDev of frame length = %.3f\n', ...
                    i, ndevs(i), ndevs_pct(i), stdevs(i) );
            end        
        end
        if (maxStdDev > maxStdDevTh)        
            curmsg = sprintf('[Frames.uneven]. %d (out of %d) pres have stddev_framelength > %.3f (max stdev is %.3f at pres # %d, with %d (%.3f %s) uneven frames)', ...
                nnz(stdevs > maxStdDevTh), nPres, maxStdDevTh, maxStdDev, indMaxStdDev, ndevs(indMaxStdDev), ndevs_pct(indMaxStdDev), pct_s);
            msg = iff(isempty(msg), [curmsg], [msg brk curmsg]);
            if fullDetails && (length(find(stdevs > maxStdDevTh)) > 1)
                fprintf2('[Frames.uneven] Max fluctuation was in Pres# %d (StdDev of %.3f, with %d uneven frames)\n', ...
                    indMaxStdDev, maxStdDev, ndevs(indMaxStdDev));
            end            
        end
        if ~isempty(siteOk) && ~fullDetails && ~checkAllPres,   error(siteOk, msg);  end
        
        % 8. Check that not too many frames were missed (not > 2)
        nMissedFrames = getFieldsFromDatabaseTable(hnd, 'LNG_N_MISSED_FRM', stimTable, {'DATAFILE_ID', Did});
        if any(nMissedFrames > 2)  % there are actually only 3 noise presentations that satisfy this criterion. For others, there are a bunch of single (or double) missed frames.
            pct = max(nMissedFrames)/median(nSupposedToShowFrames)*100;
            [siteOk, curmsg] = deal('db:MultipleMissedFrames', sprintf('[Missed.Frames]A total of %d / %d frames (ie. %.3f%s) were missed', nMissedFrames, nSupposedToShowFrames, pct, pct_s));
            msg = iff(isempty(msg), [curmsg], [msg brk curmsg]);            
        end
        presOk = presOk & (nMissedFrames(:)' <= 2);
        if ~isempty(siteOk) && ~fullDetails && ~checkAllPres,   error(siteOk, msg);  end

        % (9. Flashed grating movies: Check that are multiple phases.)
        if errorForFgCompExps
            fgCompGids = [ 2625 2665 2677 2723 2737 2747 2775 2791 2805 2887 2921 2957 2981 3047 3061 3077 3091];
            if any(Gid == fgCompGids)
                [siteOk, msg] = deal('db:SinglePhaseMovie', 'This movie does not have multiple phases for each ori/spfreq.');                                
                addIdStr = false;
                presOk = false(1,nPres);
            end            
        end
        if ~isempty(siteOk) && ~fullDetails && ~checkAllPres,   error(siteOk, msg);  end

        % (10. Flashed grating movies: Check that are multiple trials.)
        fgSingleTrialGids = [4191];
        if any(Gid == fgSingleTrialGids)
            [siteOk, msg] = deal('stim:singleTrial', 'Only 1 repetition of each grating stimulus was presented (need at least 2)');
            addIdStr = false;
            presOk = false(1,nPres);
        end            
        if ~isempty(siteOk) && ~fullDetails && ~checkAllPres,   error(siteOk, msg);  end

        
        if isempty(siteOk) && all(presOk == false)
            siteOk = 'no good presentations';
        end
        if ~isempty(siteOk) && ~fullDetails && ~checkAllPres,   error(siteOk, msg);  end
        
        if isempty(siteOk)
            if isempty(msg)
                siteOk = 'ok';
                addIdStr = false;
%                 if fullDetails
%                     disp('No problems found');
%                 end
            else
                siteOk = 'warning';
            end
        end
    
    end
end


% function checkEvenlySpacedFrames(syncTimes_ms, presStartSyncIds,
% presEndSyncIds, th);
