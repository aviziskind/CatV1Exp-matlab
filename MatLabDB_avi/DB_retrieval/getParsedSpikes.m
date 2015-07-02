function varargout = getParsedSpikes( calcMode_S, Gid, cellId, varargin )
    %    [spkTsRelToFrame, bckgSamples, meanStimFiringRate, spike_firstFrameId] = getParsedSpikes('timing', Gid, cellId, BkgrSkip_ms, extFrameLength_ms, framesToAnalyze)
    %    [relContrOfFrameToSpike, meanStimFiringRate] = getParsedSpikes('frame',  Gid, cellId, spkTimeWindow_ms, windowProfile )
    %    [framesTimesForEachSpike, framesIdsForEachSpike] = getParsedSpikes('spikes',  Gid, cellId, spkTimeWindow_ms, windowProfile )
    
    % Integrates the spike times and frame sync times for a particular
    % experiment, so that you can tell how many spikes there were in
    % response to each frame, and what the spike times were relative to the
    % beginning of the frame time (useful for generating PSTHs, OSPs & STAs).
    %
    % In general, you might call this function under one of two
    % circumstances, corresponding to the two 'calcModes' available
    % (a) PSTH generation (calcMode = 'timing') : to see what the spike
    %    profile ( P(spike)vs stimulus onset) is for a particular cell. 
    %    For this, a timeWindow input with zero delay is typically used,
    %    and if the frameLength is short (<50 ms), the profile over 
    %    many framelengths is typically considered. If you're interested in
    %    only certain stimuli, you can analyze only these stimuli by 
    %    providing the frames in which they occur in the optional parameter
    %    'framesToAnalyze'.  The spike profile can be computed by
    %    binning the spkTsRelToFrame output.
    %
    % (b) OSP / STA generation (calcMode = 'frame'): Finding how much each frame 
    %     contributed to each spike, usually using information about the
    %     PSTH / spike profile. A timeWindow is usually used with non-zero delay.
    %     In this case, the output provided is the '(mean) relative contribution
    %     of each frame to (any particular) spike', or 'relContrOfFrameToSpike'
    %   NOTE: This use of this code is deprecated. I found a much more efficient way to
    %    estimate the OSP profile. The PSTH data is saved (retrieved using dbGetCellSpkStimHists)
    %
    % (c) SPIKES - also for OSP generation.
    %
    % Notes on input parameters:
    % (3) spkTimeWindow_ms:  can be either 
    %        [delay_ms, windowEnd_ms]
    %    or 
    %        [delay_ms] in which case frameLength_ms is used as the windowSize)
    %    or   
    %        [] in which 0 delay, and frameLength_ms is used as the windowsize
    %                   
    % (5)  BkgrSkip_ms:  If a value for BkgrSkip_ms is provided,
    %   the mean (& std deviation) background rate (when no stimulus was presented) is
    %   calculated, (skipping the first 'BkgrSkip_ms' ms of each period of
    %   'downtime' between presentations).
    %
    % (6,7) extFrameLength_ms & framesToAnalyze are optional inputs for the 'timing' mode,
    %   for use with flashed grating movies where the length of each frame is
    %   too short to generate a PSTH. Then you can use these parameters to
    %   obtain the response to selected subsets of the frames to which the cell responded
    %   using an extended frameLength_ms of your choosing, and selecting
    %   the groups of frames which correspond to the same flash-grating stimulus
    %   using the framesToAnalyze parameter.

    doTest = strcmp(Gid, 'test');
    doSiteValidation = true;
    TIMING = 1;
    FRAME = 2;
    SPIKES = 3;

    % plotFramesVariables
    drawSpikeLabels = true;
    drawSpikes = true;    
    includeWindowEnd = false;
    drawGratingStimuli = true;
    
    %%%%%%%%%%%%%    INITIALIZE, PARSE INPUT PARAMETERS     %%%%%%%%%%%%% 
    
    if ~doTest
%         hnd = dbOpenExpDb;
        sd = siteDataFor(Gid);
        Did = sd.Did;        
        presIds = sd.presIds;
        nPres = length(presIds);
        frameLength_ms = sd.frameLength_ms;
    else        
%         rand('state', 0);
        nFrames = 30;
        frameLength_ms = 35;
        syncTimes_ms = [0:frameLength_ms:frameLength_ms*nFrames];
        nSyncs = length(syncTimes_ms);
        nSpikes = 10;
        spikeTimes_ms = sort( frameLength_ms*3 + rand(nSpikes, 1)*frameLength_ms*(nFrames-5) );
        
        spkTimeWindow_ms = [30 80];

        delay_ms = spkTimeWindow_ms(1);
        undelayedSpikeTimes_ms = spikeTimes_ms - delay_ms;        
        
        nPres = 1;
        presStartSyncIds = 1;
        presEndSyncIds = length(syncTimes_ms);
        doSyncValidation = false;
    end        
%     dbug = true;
       
    % Process Arguments
    calcMode_S = lower(calcMode_S); % (lower case)
    switch calcMode_S
        case 'timing',
            calcMode = TIMING;
            if (nargin >= 4), BkgrSkip_ms       = varargin{1}; end;
            if (nargin >= 5), extFrameLength_ms = varargin{2}; end;
            if (nargin >= 6), framesToAnalyze   = varargin{3}; end;
            
            extendedFrameMode = exist('extFrameLength_ms', 'var') && ~isempty(extFrameLength_ms);
            if extendedFrameMode                
                extWindowSize_ms = rectified(extFrameLength_ms(end)-frameLength_ms);  % = how much further beyond frameLength.
                if length(extFrameLength_ms) == 1
                    spkTimeWindow_ms  = [0 extWindowSize_ms];   % must look for spikes to include for each frame.
                elseif length(extFrameLength_ms) == 2
                    spkTimeWindow_ms = [extFrameLength_ms(1), extWindowSize_ms];
                end
            else
                spkTimeWindow_ms  = [0 0];  % frameLength_ms is far enough.
            end
%             windowProfile  = [1 1];
        case {'frame', 'spikes'}
            
            calcMode = iff(strcmp(calcMode_S, 'frame'), FRAME, SPIKES);
            if (nargin >= 4),
                spkTimeWindow_ms  = varargin{1}; 
            else
                spkTimeWindow_ms  = []; 
            end
            if (nargin >= 5), windowProfile     = varargin{2}; end
            
            if ~exist('windowProfile', 'var') || isempty(windowProfile);
                windowProfile  = [1 1];
            end
            extendedFrameMode = false;            
    
        otherwise
            error('Unknown mode (must be "Timing" or "Frame" or "Spikes")');
    end
    
    calcBkgr = exist('BkgrSkip_ms', 'var') && ~isempty(BkgrSkip_ms);    
    if isempty(spkTimeWindow_ms), spkTimeWindow_ms = [0 frameLength_ms]; end
    
    delay_ms = spkTimeWindow_ms(1);
    windowSize_ms = diff(spkTimeWindow_ms);
    
    useWindowProfile = (calcMode == FRAME) && (windowSize_ms > 0);
    makeFlatProfile = false;
    
    
    %%%%%%%%%%%%%    RETRIEVE SPIKES  &  SYNCS   %%%%%%%%%%%%% 

    if ~doTest
        
        spikeTimes_ms = getSpikes(Gid, cellId, 'ms');        
        
        zeroSpikeTimes = (spikeTimes_ms == 0);
        spikeTimes_ms(zeroSpikeTimes) = [];
%         if nnz(zeroSpikeTimes) > 10
%             warning('db:corruptedSpikeData', 'Spike data is slightly messed up (has lots of 0s)');
%         end
        undelayedSpikeTimes_ms = spikeTimes_ms - delay_ms;
        
%         [syncTimes_ticks0 = dbGetSyncs('Gid', Gid);    
%         [beginTime_ticks, endTime_ticks] = dbGetTbTe(hnd, Did);
        [syncTimes_ticks0, beginTime_ticks, endTime_ticks] = getSyncs(Gid, 'tick');
            
        dfInfo = sd.dataFileInfo;
        
        ticksPerMs = dfInfo.samplingRateHz/1000;
        lastTick = dfInfo.duration_sec * dfInfo.samplingRateHz;
        lastSpike_tk = spikeTimes_ms(end) * ticksPerMs;
        % sanity check:
        assert(lastTick > lastSpike_tk ); % make sure last tick is after last spike.
        assert(lastTick > syncTimes_ticks0(end)); % make sure last tick is after end of last frame.

        syncTimes_ticks = [0; syncTimes_ticks0; lastTick]; % add my own ticks at start & end of experiment.
        syncTimes_ms    = dbConvertTimeMeasures(Did, syncTimes_ticks, 'tick', 'ms');
        nSyncs = length(syncTimes_ticks);

        %%%%%%%%%%%%%    RETRIEVE  START/END TICKS OF PRESENTATIONS    %%%%%%%%%%%%% 
%         [beginTime_ticks, endTime_ticks] = dbGetTbTe(hnd, Did);
        
        presStartSyncIds = binarySearch(syncTimes_ticks, beginTime_ticks);
        presEndSyncIds = binarySearch(syncTimes_ticks, endTime_ticks);
    end
    
    nFramesInEachPres = presEndSyncIds - presStartSyncIds;
    cumulativeFrames = [0; cumsum(nFramesInEachPres)];
    nTotalFrames = sum(nFramesInEachPres);

    if ~exist('framesToAnalyze', 'var')
        framesToAnalyze = 1:nTotalFrames;
    end

    %%%%%%%%%%%%%    DO SOME VALIDATION OF SYNCS AND START/END TICKS    %%%%%%%%%%%%%     
    if doSiteValidation

        [siteOK, siteOKdetail, presOK] = testSiteValidity(Gid, 'report');

    end
        

    % Calculate mean firing rate across the entire experiment (for window profile processing, and for output) 
    presTime_sec = sum( diff([syncTimes_ms(presStartSyncIds), syncTimes_ms(presEndSyncIds)], 1, 2) / 1000 );
    nSpikesInPresentations = sum(elementsInRange(spikeTimes_ms, [syncTimes_ms(presStartSyncIds), syncTimes_ms(presEndSyncIds)], 'count'));
    meanStimFiringRate = nSpikesInPresentations / presTime_sec;
        
    
    % Window Profile processing (for FRAME mode)
    if useWindowProfile
        interpolate = false;
        nrep = 100;
        nSpikeProfileDivisions = iff(interpolate, 200, length(windowProfile)*nrep);
        timesRange = [0, windowSize_ms + frameLength_ms];
        if makeFlatProfile
            windowProfile = [1 1];
        else            
            windowProfile = rectified( fliplr(windowProfile(:)') - meanStimFiringRate ); 
        end
            
%         figure(5);
%         bar(windowProfile); drawHorizontalLine(meanStimFiringRate, 'Color', 'r');
        
        if all(windowProfile == 0) %ie. the profile was not actually that meaningful, as none of the values were above the mean firing rate, so can just use default (flat) profile.
            windowProfile = ones(size(windowProfile));
        end
        
        dT_ms = windowSize_ms / length(windowProfile);
        itpProfileTs = linspace(0, windowSize_ms, nSpikeProfileDivisions);
        profileTs_atCenters = linspace(dT_ms/2, windowSize_ms-dT_ms/2, length(windowProfile));
        if interpolate
            profileTs_atCenters = linspace(dT_ms/2, windowSize_ms-dT_ms/2, length(windowProfile));            
            itpProfileVals = rectified( interp1(profileTs_atCenters, windowProfile,  itpProfileTs, 'cubic') );
        else
            eps = 1e-5;
            itpProfileTs = bsxfun(@plus, profileTs_atCenters, linspace(-dT_ms/2+eps, dT_ms/2-eps, nrep)' );
            itpProfileTs = itpProfileTs(:)';
            
            itpProfileVals = repmat(windowProfile, nrep, 1);
            itpProfileVals = itpProfileVals(:)';
        end
        
        
        allTimesFromFrmStartToUDSpikeTimes = linspace( timesRange(1), timesRange(2), nSpikeProfileDivisions);
        frameProminenceInSpikeProfile = zeros(nSpikeProfileDivisions, 1);
        for di = 1:nSpikeProfileDivisions        
            spkWindowTimesForThisOffset = itpProfileTs - (windowSize_ms) + allTimesFromFrmStartToUDSpikeTimes(di);
            inds = elementsInRange(spkWindowTimesForThisOffset, [0 frameLength_ms], 'index', '[)');
            frameProminenceInSpikeProfile(di) = sum( itpProfileVals(inds)  );
        end
        frameProminenceInSpikeProfile = frameProminenceInSpikeProfile/(sum( itpProfileVals ));
    end
%{
     if dbug
        figure(701);
        subplot(3,1,1); bar(windowProfile, 1);
        subplot(3,1,2); plot(itpProfileTs, itpProfileVals, '-o');        
        subplot(3,1,3); plot(itpProfileTs, frameProminenceInSpikeProfile, 'g-.');            
     end
%}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% BEGIN PARSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Prepare variables

    intervalsContrToEachSpike = cell(1,nSpikesInPresentations);
    if calcMode == TIMING
        spkTsRelToFrame  = cell(1,nTotalFrames);         
    elseif calcMode == FRAME
        relContrOfFrameToSpike = arrayfun(@(n) zeros(1,n), nFramesInEachPres, 'UniformOutput', false);                    
    elseif calcMode == SPIKES
        [framesTimesForEachSpike, framesIdsForEachSpike] = deal( cell(1,nSpikesInPresentations) );                 
    end

    frameStarts_ms = syncTimes_ms(1:end-1);     % note: these 'frames' include the interspike intervals
    if ~extendedFrameMode   % note: I only actually use frameEnds_ms for plotting.
        frameEnds_ms   = syncTimes_ms(2:end);        
    else
        frameEnds_ms   = frameStarts_ms + extFrameLength_ms(end);  % (2nd element if 2-vector. 1st element if scalar.)
%         % special cases: end of blank frames:
%         frameEnds_ms(presStartSyncIds-1) = syncTimes_ms(presStartSyncIds);
    end
    
    % For each spike, find the indices of the first & last frame that
    % contributed to that spike 
    firstIntervalContainingSpikeWindow = binarySearch(frameStarts_ms, undelayedSpikeTimes_ms-windowSize_ms, 1, -1);
    lastIntervalContainingSpikeWindow  = binarySearch(frameStarts_ms, undelayedSpikeTimes_ms,               1, -1);


    % Use above variables to find out the spikes whose windows overlap with 
    % each interval (interval = frame or blanks between presentations.
    
    
    estPeakHz = 100;
    estNumSpkPerInterval = round(estPeakHz * diff (spkTimeWindow_ms)/1000)+1;    
    nSpkPerInterval = zeros(1,nSyncs-1);    
    spikesContributingToEachInterval = cell(1,nSyncs-1);
    spikesContributingToEachInterval(:) = {zeros(1, estNumSpkPerInterval)};    
    for spk_i = 1:length(spikeTimes_ms)
        intervalsForSpkI = firstIntervalContainingSpikeWindow(spk_i):lastIntervalContainingSpikeWindow(spk_i);
        if (windowSize_ms == 0) && length(intervalsForSpkI) > 1
            keyboard;
        end
        for int_j = intervalsForSpkI
            nSpkPerInterval(int_j) = nSpkPerInterval(int_j)+1;
            spikesContributingToEachInterval{int_j}(nSpkPerInterval(int_j)) = spk_i ;
        end
        intervalsContrToEachSpike{spk_i} = intervalsForSpkI;
    end    
    spikesContributingToEachInterval = cellfun(@(x,n) x(1:n), spikesContributingToEachInterval, num2cell(nSpkPerInterval), 'un', 0);
    
    
    
    syncIdForFrameStartId = zeros(1, sum(nFramesInEachPres) );  % find exactly which 'sync' corresponds to each frame.    
    frameIdForIntervalId  = zeros(1, length(syncTimes_ms)   );  % find exactly which 'frame' corresponds to each sync.    
%     intervalIsFrame       = false(1, length(syncTimes_ms)-1 );  % find which intervals are stimulus frames (vs inter-presentation intervals)
    for iPres = 1:nPres
        offset = iPres;
        frame_idxs = cumulativeFrames(iPres)+1:cumulativeFrames(iPres+1);
        syncIdForFrameStartId(frame_idxs) = frame_idxs+offset;        
        frameIdForIntervalId(frame_idxs+offset) = frame_idxs;        
    end
    
%     framesWithSpikes = ~cellfun(@isempty, spikesContributingToEachInterval);
    
    
%     framesToAnalyze_withSpikes = framesToAnalyze(framesWithSpikes);
%     framesToAnalyze = framesToAnalyze(framesToAnalyze_withSpikes);
%     plotFrames(1, 12, 2);
    
    %%% MAIN LOOP: process each frame
    whichPres = zeros(1, sum(nFramesInEachPres));
    
    if any(calcMode == [TIMING, FRAME])
        for iFrm = framesToAnalyze;             
                        
            iPres = find(iFrm <= cumulativeFrames, 1, 'first')-1;        
            iFrmInPres = iFrm-cumulativeFrames(iPres);
            
            whichPres(iFrm) = iPres;
            
%             frmSyncId = presStartSyncIds(iPres) + iFrmInPres -1;             
            frmSyncId = syncIdForFrameStartId(iFrm);
            frmStart = syncTimes_ms(frmSyncId);

            idx_spikesWithWindowInFrame = spikesContributingToEachInterval{frmSyncId};        

            if ~isempty(idx_spikesWithWindowInFrame)
                if calcMode == TIMING
                    spkTsRelToFrame{iFrm} = (spikeTimes_ms(idx_spikesWithWindowInFrame))' - frmStart;  % *not* undelayedSpikeTimes: want actual spike times.

                elseif calcMode == FRAME
                    if useWindowProfile
                        timesFromFrameStartToUDSpikes = undelayedSpikeTimes_ms(idx_spikesWithWindowInFrame) - frmStart;
                        indsHowMuchOfFrameContrToSpike = binarySearch(allTimesFromFrmStartToUDSpikeTimes, timesFromFrameStartToUDSpikes);
                        frameContributionToSpikes = frameProminenceInSpikeProfile(indsHowMuchOfFrameContrToSpike);
                        relContrOfFrameToSpike{iPres}(iFrmInPres) = sum(frameContributionToSpikes);
                    else
                        relContrOfFrameToSpike{iPres}(iFrmInPres) = length(idx_spikesWithWindowInFrame);
                    end                    
                end

            end

        end
        
    elseif (calcMode == SPIKES)
         
         for iSpk = 1:length(spikeTimes_ms);              
             spkTime = spikeTimes_ms(iSpk);
             intervalsContr = intervalsContrToEachSpike{iSpk};
             intervalsContr = [intervalsContr(end)+1, intervalsContr(end:-1:1)];
             
             frmTimes = syncTimes_ms(intervalsContr)'; 
             frameIds = frameIdForIntervalId(intervalsContr);
             
             framesTimesForEachSpike{iSpk} = spkTime - frmTimes;
             framesIdsForEachSpike{iSpk} = frameIds;             
         end
    end                
        
        
    if calcMode == FRAME
        3;
    end
    
    
    if calcBkgr

        % (1) Pre-blank frames        
    %    [relContrOfFrameToSpike, meanStimFiringRate] = getParsedSpikes('frame'
        
        % if drifting gratings, we know the temporal frequency of the cycles, so we can bin the
        % background spikes right now into the appropriate-sized bins.
%         binBkgrSpikes = (flashedOrDrifting(Gid) == 2);
        minChunkSize_ms = 4+(1/6);
    
        stimType = sd.stimType;
        msPerSec = 1000;
        % divide the interpresentation time into "chunks" so that can calculate mean & std deviation background spike rate.
        if strncmp(stimType, 'Grating:Orientation Batch', 25) || strncmp(stimType, 'Grating:Spatial Frequency Batch', 25)            
            % use the same binSize as for the stimulus.
            chunks_ms = sd.tempPeriod_sec * msPerSec;
        else % flashed gratings
            
            chunks_ms = minChunkSize_ms;
%             if extendedFrameMode
%                 chunks_ms = extFrameLength_ms(end);    
%             else
%                 [tmp1, tmp2, extFrameLength_ms] = getBinSizeParams(frameLength_ms);
%                 chunks_ms = extFrameLength_ms;
%             end            
            
        end
            
%         preBlankInterval        = [BkgrSkip_ms           syncTimes_ms(2)];        
%         nPreBlankSpikes        = elementsInRange(undelayedSpikeTimes_ms, preBlankInterval, 'count');                
%         bkgTime_ms             = diff(preBlankInterval);
%         nBkgSpks = nPreBlankSpikes;        
        allChunks_C = cell(1, nPres+1); % each inter-presentation interval, + beginning + end.
        
        preBlankInterval_chunks = [BkgrSkip_ms:chunks_ms:syncTimes_ms(2)];
        nPreBlankSpikes_chunks = elementsInRange(spikeTimes_ms, preBlankInterval_chunks, 'count');
        
        allChunks_C{1} = nPreBlankSpikes_chunks';        
        
        % (2) Inter-presentation intervals
        for iPres = 1:nPres-1            
            frmSyncId = presEndSyncIds(iPres);

%             interPresInterval        = [syncTimes_ms(frmSyncId)+BkgrSkip_ms,              syncTimes_ms(frmSyncId+1)];
%             if diff(interPresInterval) > 0
%                 nBkgSpikesInInterval = elementsInRange(spikeTimes_ms, interPresInterval, 'count');                
%                 nBkgSpks = nBkgSpks + nBkgSpikesInInterval;
%                 bkgTime_ms = bkgTime_ms + diff(interPresInterval);    
%             end            
            
            interPresInterval_chunks = [syncTimes_ms(frmSyncId)+BkgrSkip_ms : chunks_ms : syncTimes_ms(frmSyncId+1)];
            if length(interPresInterval_chunks) >= 2
                nBkgSpikesInInterval_chunks = elementsInRange(spikeTimes_ms, interPresInterval_chunks, 'count');                
                allChunks_C{iPres+1} = nBkgSpikesInInterval_chunks';
            end
            
        end        
        
        % (3) Post-blank frames
%         postBlankInterval = [syncTimes_ms(end-1) + BkgrSkip_ms,    syncTimes_ms(end)];
%         if diff(postBlankInterval) > 0
%             nPostBlankSpikes = elementsInRange(spikeTimes_ms, postBlankInterval, 'count');
%             nBkgSpks = nBkgSpks + nPostBlankSpikes;
%             bkgTime_ms = bkgTime_ms + diff(postBlankInterval);
%         end
        
        postBlankInterval_chunks = [syncTimes_ms(end-1) + BkgrSkip_ms : chunks_ms : syncTimes_ms(end)];        
        if length(postBlankInterval_chunks) >= 2
            nPostBlankSpikes_chunks = elementsInRange(spikeTimes_ms, postBlankInterval_chunks, 'count');
            allChunks_C{nPres+1} = nPostBlankSpikes_chunks';
        end
        
        allSpikesInEachChunk = [allChunks_C{:}];
        % calculate background rate
%         bckgRate = nBkgSpks/bkgTime_ms * 1000; % from ms to sec

        chunkPerSec = 1000 / chunks_ms;
        
        allSpikeRates_Hz_eachChunk = allSpikesInEachChunk * chunkPerSec;
        
        bckgRate_mean = mean(allSpikeRates_Hz_eachChunk);
        bckgRate_std  = std(allSpikeRates_Hz_eachChunk);
        
        bckgSamples =  allSpikeRates_Hz_eachChunk;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if calcMode == TIMING
        varargout{1} = spkTsRelToFrame;
        if calcBkgr
            varargout{2} = bckgSamples;
        end
        varargout{3} = meanStimFiringRate;
        nSpikesInPres_actual = nSpikesInPresentations; % sum(elementsInRange(spikeTimes_ms, [syncTimes_ms(presStartSyncIds), syncTimes_ms(presEndSyncIds)], 'count'));
        nSpikesInPres_calc = sum(cellfun(@length, spkTsRelToFrame));
        discrep = abs(nSpikesInPres_actual - nSpikesInPres_calc )/nSpikesInPres_actual * 100;
%         fprintf('Timing: Check: # spikes in presentations: %d. # spikes in output : %d.  (discrepancy of %.4f %%)\n', nSpikesInPres_actual, nSpikesInPres_calc , discrep);
        
        varargout{4} = whichPres;
                
        if nargout >= 5
            firstFrameContainingSpikeWindow = frameIdForIntervalId(firstIntervalContainingSpikeWindow);
            varargout{5} = firstFrameContainingSpikeWindow; % "spike_firstFrameId"
        end

    elseif calcMode == FRAME
        nSpikesInPres_calc = sum(cellfun(@sum, relContrOfFrameToSpike));

        varargout{1} = relContrOfFrameToSpike;
        varargout{2} = meanStimFiringRate;        
        
        meanStimFiringRate2 = nSpikesInPres_calc / (nTotalFrames * frameLength_ms /1000);        
        
        nSpkDiscrep = abs(nSpikesInPresentations-nSpikesInPres_calc)/nSpikesInPresentations * 100;
        rateDiscrep = abs(meanStimFiringRate-meanStimFiringRate2)/meanStimFiringRate * 100; % is usually only about .2602% 
        
%         fprintf('Frame: Check: # spikes in presentations: %d. # spikes in output : %d.  (discrepancy of %.4f %%)\n', nSpikesInPres_actual, round(nSpikesInPres_calc), nSpkDiscrep);
%         fprintf('Frame 2: Check: spkrate 1 : %.3f # : spkrate 2 : %.3f.  (discrepancy of %.4f %%)\n', meanStimFiringRate, meanStimFiringRate2, rateDiscrep);
%         disp(['nSpikes : '  num2str(nSpikes) '; calc: ' num2str(nSpkCalc) ' [ ' num2str(discrep, '%3.4f') '%]']);

     elseif calcMode == SPIKES
         varargout{1} = cellfun(@single, framesTimesForEachSpike, 'un', 0);
         varargout{2} = cellfun(@uint16, framesIdsForEachSpike, 'un', 0);
         
         
    end

    
    
    %---------------- Method for visualizing frames & spikes during debugging
    function h = drawSpikeWindow(spkT, m, d, col)
        l = m-d/2;
        u = m+d/2; 

        h(1) = line([0,0] + spkT,                         [l m], 'Color', col);
        if windowSize_ms > 0
            h(2) = line([spkT, spkT - spkTimeWindow_ms(1)],   [m m], 'Color', col, 'LineStyle', ':'); 
            h(3) = line([spkT - spkTimeWindow_ms],            [m m], 'Color', col); 
            h(4) = line([0,0] + [spkT - spkTimeWindow_ms(2)], [m u], 'Color', col); 
            h(5) = line([0,0] + [spkT - spkTimeWindow_ms(1)], [m u], 'Color', col); 
        end        
    end

    
    
    
    
    
    function plotFrames(startId, endId, locMode) % spkTimeWindow_ms, syncTimes_ms, spikeTimes_ms, windowSize_ms, frameLength_ms, extendedFrameMode, relContrOfFrameToSpike, itpProfileTs, itpProfileVals) 
        if nargin < 3
            locMode = 'sync';
        end
        mainAx = gca;
        switch locMode
            case 'sync',
                [syncStartId, syncEndId] = deal(startId, endId);
                [intStartId,  intEndId]  = deal(syncStartId, syncEndId-1);                
            case 'interval', 
                [intStartId,  intEndId] = deal(startId, endId);
                [syncStartId, syncEndId]  = deal(intStartId, intEndId+1);                
                        
            case 'frame',
                [frameStartId,  frameEndId] = deal(startId, endId);                
                presStartId = find(frameStartId <= cumulativeFrames, 1, 'first')-1;
                presEndId = find(frameEndId <= cumulativeFrames, 1, 'first')-1;                
                presFrameStartId = frameStartId - cumulativeFrames(presStartId);
                presFrameEndId   = frameEndId - cumulativeFrames(presEndId);                
                                                
            case 'pres/frame',                
                [presStartId, presFrameStartId, presEndId, presFrameEndId] = ...
                    deal(startId(1), startId(2), endId(1), endId(2));
                frameStartId = cumulativeFrames(presStartId)+presFrameStartId;                    
                frameEndId   = cumulativeFrames(presEndId)+presFrameEndId;                
        end
        if ~exist('intStartId', 'var') % ie. locMode == 'frame' or 'pres/frame':
            syncStartId = presStartSyncIds(presStartId) + presFrameStartId -1;
            syncEndId   = presStartSyncIds(presEndId)   + presFrameEndId ;
            [intStartId,  intEndId]  = deal(syncStartId, syncEndId-1);
        end
        
        if syncStartId < 1,       syncStartId = 1;       end;
        if syncEndId   > nSyncs,  syncEndId   = nSyncs;  end;
        if intStartId < 1,        intStartId = 1;        end;
        if intEndId   > nSyncs-1, intEndId   = nSyncs-1; end;
        
        syncIds = syncStartId:syncEndId;
        intIds = intStartId:intEndId;        

        [presentnIds, frameIds, presFrameIds] = deal(zeros(1, length(intIds)));
        for int_i = 1:length(intIds);
            if any(intIds(int_i) == presStartSyncIds-1)
                continue;
            end
            pres_i = find(intIds(int_i) >= presStartSyncIds-1, 1, 'first');

            if ~isempty(pres_i),
                presentnIds(int_i) = pres_i;
                frameIds(int_i) = intIds(int_i) - pres_i;
                presFrameIds(int_i) = frameIds(int_i)-cumulativeFrames(pres_i);
            end

        end        
        
        nfrm_cyc = 8;
        nspk_cyc = 5;
        nspkw_cyc = 6;
        
        Hspikelines = 7;
        HspikeWindows = 15;
        Hheader = 8;
        
        y_bottom = 0;
        y_spkLines_spkWindows = y_bottom + Hspikelines;
        y_spkWindows_header = y_spkLines_spkWindows + HspikeWindows;
        y_top = y_spkWindows_header + Hheader;
        ys = [y_bottom, y_spkLines_spkWindows, y_spkWindows_header, y_top];
        
        dy = (y_top-y_bottom)/20;
%         y0 = 0;   %bottom of plot
%         y1 = .7;  %y0-y1: spike lines
%         y2 = 1;   %y2-pres/presframe
%         y3 = 1.3; %y2-frame
%         y4 = 2.3;   %top of plot
        
        syncTimes = syncTimes_ms(syncIds);        
        syncTimes = syncTimes(:)';
        frameSpks = elementsInRange(spikeTimes_ms, [syncTimes(1), syncTimes(end)+includeWindowEnd*(frameLength_ms+windowSize_ms)], 'index');        
        range_cyc = @(a, b, n_cyc, N)  (mod([0:N-1],n_cyc)*((b-a)/(n_cyc-1)))+a;
        
        % 1. Set x/y limits ;
        xlim([syncTimes(1) syncTimes(end)+includeWindowEnd*frameLength_ms]);
        ylim([y_bottom y_top]); 
        set(mainAx, 'xtick', syncTimes, 'ytick', []);
        box on;
        htest = text(frameStarts_ms(syncIds(1)),y_bottom,'X', 'visible', 'off'); htest_loc = get(htest, 'extent');
        textHgt = htest_loc(4);        
        curYPos = y_top;
        
        % 1. Draw dotted lines demarking frame onset times
        line([1;1]*syncTimes(:)',   [y_spkLines_spkWindows; y_top] * ones(1,length(syncTimes)), 'Linestyle', ':', 'Color', 'k');
                 
        % 1b. Labels on each frame:
        if exist('relContrOfFrameToSpike', 'var')
            r = [relContrOfFrameToSpike{:}];
        end
        for i = 1:length(intIds)
            x = mean(syncTimes(i:i+1));
%             text( x, y2, num2str(intIds(i)), 'color', 'b', 'Ho', 'center');          
            text( x, curYPos-textHgt/2, num2str(intIds(i)), 'color', 'b', 'Ho', 'center');                      
            
            if (frameIds(i) > 0)                
%                 text( x, y3, {num2str(presentnIds(i)), num2str(presFrameIds(i)), ['[' num2str(frameIds(i)) ']']}, 'Ho', 'center');            

                if exist('relContrOfFrameToSpike', 'var')
%                     text( x, [y_top - textHgt*3.5], num2str(r(frameIds(i)), '%.2f') , 'Ho', 'center', 'color', 'r');
                end
            end
        end
        curYPos = curYPos-textHgt;
        if extendedFrameMode
            for frm_i = 1:length(frameIds)
%                 y = y4-(mod(frm_i,nfrm_cyc)+.5)/(nfrm_cyc/y1);
                line( [frameStarts_ms(frameIds(frm_i)) frameEnds_ms(frameIds(frm_i))], [y y], 'Color', color_s(frm_i));
                text( frameStarts_ms(frameIds(frm_i)), y + .05, num2str(frameIds(frm_i)), 'Color', color_s(frm_i) );
            end            
        end
        
        % 2. Spikes
        if drawSpikes
            nColors = 7;
            nSpikes = length(frameSpks);
            for ci = 1:nColors
                inds = ci:nColors:nSpikes;
                if ~isempty(inds)
                    spks = spikeTimes_ms(frameSpks(inds));
                    line( [1;1] * spks(:)', [y_bottom; y_spkLines_spkWindows]*ones(1,length(inds)), 'Color', color_s(ci));
                end
            end
        end
        
                
      % (Spike number labels)
        if drawSpikeLabels
            spikeLabelHeights = range_cyc(textHgt/2, y_spkLines_spkWindows-textHgt/2, nspk_cyc, length(frameSpks));
            for si = 1:length(frameSpks)
                text(spikeTimes_ms(frameSpks(si)), spikeLabelHeights(si), ['\leftarrow ' num2str(frameSpks(si))]);
            end
        end
        
        drawSpikeWindows = true;
        % 3. Spike-Window lines        
        plotSpikeProfiles = useWindowProfile && ~all(diff(itpProfileVals) == 0);
        
            rng = [y_spkLines_spkWindows+dy, y_spkWindows_header-dy];
%             if plotSpikeProfiles
                rng(2) = rng(2)- diff(rng)/nspkw_cyc;
%             end
            spikeWindowLineHeights = range_cyc(rng(1), rng(2), nspkw_cyc, length(frameSpks));
    %         heights = linspace( mod(  rng(1), rng(2), length(frameSpks)+2);

            d = diff(rng)/(length(frameSpks)+2);
        if drawSpikeWindows
            for si = 1:length(frameSpks)
                drawSpikeWindow(spikeTimes_ms(frameSpks(si)), spikeWindowLineHeights(si), d, color_s(si));        
            end
        end
%         disp( ['frame lengths (in seconds) = ' num2str(diff(syncTimes(2:3)/1000)) ] )
        
        % 4. Spike Profiles in windows
        if drawSpikeWindows && plotSpikeProfiles
            hold on        
            % spike profile stairs:
            skp = round(nSpikeProfileDivisions/30);
            
            if skp < 1, skp = 1; end
            scl = .8*(diff(rng)/(nspkw_cyc)) * (1/max(itpProfileVals));
            for si = 1:length(frameSpks)    
                   ts = fliplr(linspace(itpProfileTs(1), itpProfileTs(end), length(1:skp:length(itpProfileTs))+1));
                   vs = itpProfileVals(1:skp:end); vs = [vs, vs(end)]; %#ok<AGROW>
                   
                   x = spikeTimes_ms(frameSpks(si))  - delay_ms - ts;
                   y = spikeWindowLineHeights(si) + scl * vs;

    %             x = spikeTimes_ms(frameSpks(si)) + linspace( delay_ms 
    %             spkT + linspace(spkTimeWindow_ms(1):
    %             y = 
                stairs(x, y, color_s(si));
            end  
            hold off;
        end        
        
        if drawGratingStimuli
%             for i = 1:
            getFrame = @getMovieStimulusFrame;
            getFrame('load', Gid);
            for frm_id = frameIds(find(frameIds))
                x1 = frameStarts_ms(frm_id+1); x2 = frameEnds_ms(frm_id+1); dx = x2-x1;
                try
                    frm = getFrame(frm_id);                
                catch
                    frm = zeros(5);
                end
                ll_xy = [x1+dx/5,  y_top-textHgt*3];
                ur_xy = [x2-dx/5,  y_top-textHgt*.7];                
                stimBoxDs = corners2box(ll_xy, ur_xy);
                stimBoxPos = dsxy2figxy(mainAx, stimBoxDs);
                hStim = axes('Position', stimBoxPos);
                imagesc(frm);
                colormap('gray')
                axis equal tight;
                set(gca, 'xtick', [], 'ytick', []);
                
            end
            3;
            
            
        end
        
        
    end
    
    function viewExperiment
        clf;
        M = floor(sqrt(nPres));
        N = ceil(nPres/M);        
%         nFrm = max(nSustainedFrames);
        h = zeros(1,length(relContrOfFrameToSpike));
        for i = 1:nPres
            subplot(M,N,i);
            spks = relContrOfFrameToSpike{i}; nFrm = length(spks);            
            m = max(floor(sqrt(nFrm)), 1);
            n = max(ceil(nFrm/m),      1);
            
            spks(m*n) = 0;%-max(spks); %pad with negatives.            
            h(i) = imagesc(1:m,1:n, reshape(spks, [m n])');
            set(gca, 'xtick', [], 'ytick', [])  
            set(gca, 'position', get(gca, 'outerposition'))
            title(num2str(i));
        end
        matchAxes('C', h);
    end
    
    
    if ~doTest
%         edbClose(hnd);
    end

    if doTest
        figure(1); clf;
        plotFrames(1, nFrames+1)  
    
        if calcMode == TIMING
            for di = 1:nFrames
                disp([' frame ' num2str(di) ' : ' num2str(round(spkTsRelToFrame{di}))]);
            end
%             spkTsRelToFrame;
        elseif calcMode == FRAME
            for di = 1:nPres
                disp([' pres ' num2str(di) ' :   ' num2str(round(relContrOfFrameToSpike{1}))]);
            end
%             [1:nFrames; relContrOfFrameToSpike{1}]
%             relContrOfFrameToSpike;
%             nSpkCalc = sum(relContrOfFrameToSpike{1});
%             discrep = abs(nSpkCalc-nSpikes)/nSpikes*100;
            
            disp(['nSpikes : '  num2str(nSpikes) '; calc: ' num2str(nSpkCalc) ' [ ' num2str(discrep, '%3.4f') '%]']);
            
        end
    end
    
    
    
end








%         pres_i = find(frmSyncId >= [presStartSyncIds], 1, 'first');
%         if ~isempty(pres_i)            
%             assert(pres_i == iPres);
%         else
%             pres_i = 0;
%         end
%         movFrm_i = iFrm - pres_i+1;
%         assert(iFrm == movFrm_i);