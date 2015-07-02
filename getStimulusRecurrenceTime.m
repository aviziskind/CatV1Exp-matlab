function [stimPeriod_sec, stimPeriod_nFrames] = getStimulusRecurrenceTime(Gid, fracStimuli)

    if nargin < 2
        fracStimuli = 1; % all stimuli
    end

    % if isflashed -> each stimulus is in its own 'frame'/'presentation'
    % if drifting -> multiple cycles are presented together in a
    % 'presentation' - join these cycles to have one 'stimulus'
    isFlashed = (flashedOrDrifting(Gid) == 1);
    
    t = dbGetStimulusTimes(Gid);
    T = sort(t(:));
    median_dt = median(diff(T));
    mean_dt = mean(diff(T));
    if isFlashed
        dt_use = median_dt;
    else
        dt_use = mean_dt; % drifting grating inter-presentation gaps are important
    end

%%
    frameStimIds = getStimulusFrameSequence(Gid, 'OS');
    nUniqueStimuli = max(frameStimIds);
    
    
    if ~isFlashed
        idx_firstStimFrameEachPres =  find ( diff([frameStimIds;0]) ~= 0 );

        nFramesPerPres = median( diff(idx_firstStimFrameEachPres) );
        groupedFrameStimIds = frameStimIds(idx_firstStimFrameEachPres)';
        dt_use = dt_use * nFramesPerPres;        
    else
        groupedFrameStimIds = frameStimIds;
    end



    stimPeriod_nFrames = nUniqueStimuli * fracStimuli;
    
    stimPeriod_sec = stimPeriod_nFrames * dt_use;
     
    
    
    
    
    %{
    estimateTimeUntilSeeSameStimulusAgain = 0;
    
    if estimateTimeUntilSeeSameStimulusAgain
        meanRecurTime_frm_eachStim = zeros(1, maxStimId);
        for stim_i = 1:maxStimId
            idx_stim_i = find(groupedFrameStimIds == stim_i);
            if length(idx_stim_i) > 1 % multiple trials
                meanRecurTime_frm_eachStim(stim_i) = mean( diff(idx_stim_i));
            else % only 1 trial
                meanRecurTime_frm_eachStim(stim_i) = length(groupedFrameStimIds);
            end

        end

        meanRecurrencTime_frm = mean( meanRecurTime_frm_eachStim );
        meanRecurrencTime_sec = meanRecurrencTime_frm * dt_use;
        
        stimPeriod_sec = meanRecurrencTime_sec;
        stimPeriod_nFrames = meanRecurrencTime_frm;

    end
    
    
%%    
    estimateTimeToSeeAllStim = 0;
    if estimateTimeToSeeAllStim
        nSamplePoints = 10;

        meanFullSampleTime = zeros(1, nSamplePoints);
        startPoints = floor(linspace(1,  length(groupedFrameStimIds), nSamplePoints+5));
        startPoints = startPoints(1:end-1);
        for i = 1:nSamplePoints
            foundEachStim = false(1, maxStimId);
            j = startPoints(i);
            while j < length(groupedFrameStimIds) && ~all(foundEachStim)
                foundEachStim( groupedFrameStimIds(j) ) = true;
                j = j + 1;
            end
            meanFullSampleTime(i) = j - startPoints(i);

        end
    end
    
    estimateInverseOfStimChange = 0;
    if estimateInverseOfStimChange
        OS_types = {'OS', 'O', 'S'};
        for ti = 1:length(OS_types)
        %%
            [frameStimIds, uOri, uSp, uPh] = getStimulusFrameSequence(Gid, OS_types{ti});
            maxStimId = max(frameStimIds);

            % circDist is appropriate for orientations. abs(dist) is better for
            % spatial frequencies, but can use circDist for both (results in
            % slight overestimation of period for spf, but this is ok).
            allFrameStimIdChanges{ti} = circDist(frameStimIds(1:end-1), frameStimIds(2:end), maxStimId);
            allFrameStimIdChanges_ti = allFrameStimIdChanges{ti};
    %         allFrameStimIdChanges{ti} = abs( frameStimIds(1:end-1) - frameStimIds(2:end) );

            allFrameStimChanges{ti} = allFrameStimIdChanges{ti} / maxStimId;

            meanFrameStimChange(ti) = mean(allFrameStimChanges{ti});
            stimPeriod_nFrames(ti) = 1/meanFrameStimChange(ti);
            stimPeriod_sec(ti) =  stimPeriod_nFrames(ti) * dt_use;
        end    
    end
    
    %}
    
    
    
    
    
    %%
    
        % circDist is appropriate for orientations. abs(dist) is better for
        % spatial frequencies, but can use circDist for both (results in
%         % slight overestimation of period for spf, but this is ok).
%         allFrameStimIdChanges{ti} = circDist(frameStimIds(1:end-1), frameStimIds(2:end), maxStimId);
%         allFrameStimIdChanges_ti = allFrameStimIdChanges{ti};
% %         allFrameStimIdChanges{ti} = abs( frameStimIds(1:end-1) - frameStimIds(2:end) );
%         
%         allFrameStimChanges{ti} = allFrameStimIdChanges{ti} / maxStimId;
% 
%         meanFrameStimChange(ti) = mean(allFrameStimChanges{ti});
%         stimPeriod_nFrames(ti) = 1/meanFrameStimChange(ti);
%         stimPeriod_sec(ti) =  stimPeriod_nFrames(ti) * dt_use;
%     end    
    
%     if ~isFlashed
%         [frameStimIds, uOri, uSp, uPh] = getStimulusFrameSequence(Gid, 'OS');
%         idx_firstStimFrameEachPres =  find ( diff(frameStimIds) ~= 0 );
%         nFramesPerPres = median(diff(idx_firstStimFrameEachPres));
%         stimPeriod_nFrames = stimPeriod_nFrames / nFramesPerPres;
%     end
    
    
  %%  
    
    
    
end

%{

    OS_types = {'OS', 'O', 'S'};
    for ti = 1:length(OS_types)
    
        [frameStimIds, uOri, uSp, uPh] = getStimulusFrameSequence(Gid, OS_types{ti});
        maxStimId = max(frameStimIds);

        % circDist is appropriate for orientations. abs(dist) is better for
        % spatial frequencies, but can use circDist for both (results in
        % slight overestimation of period for spf, but this is ok).
        allFrameStimIdChanges{ti} = circDist(frameStimIds(1:end-1), frameStimIds(2:end), maxStimId);
%         allFrameStimIdChanges{ti} = abs( frameStimIds(1:end-1) - frameStimIds(2:end) );
        
        allFrameStimChanges{ti} = allFrameStimIdChanges{ti} / maxStimId;

        meanFrameStimChange(ti) = mean(allFrameStimChanges{ti});
        stimPeriod_nFrames(ti) = 1/meanFrameStimChange(ti);
        stimPeriod_sec(ti) =  stimPeriod_nFrames(ti) * dt_use;
    end    

%}


%{

  frameLength_sec = getFrameLength('Gid', Gid, 'sec');
  
    [allOri_deg, allSp_pix, allPhase_deg] = getOriSpPhaseForEachStimFrame(Gid);    
    
    frameStimIds0 = [0; frameStimIds; 0];
    idx_firstStimFrame = find ( diff(frameStimIds0) ~= 0 );
    nStimPres = length(idx_firstStimFrame);   
    stimIds = frameStimIds(idx_firstStimFrame(1:end-1));
    maxStimId = max(stimIds);
    
    
    allStimChanges = circDist(stimIds(1:end-1), stimIds(2:end), maxStimId) / maxStimId;
    meanStimChange_perStim = mean(allStimChanges);
    meanStimChange_stim = 1/meanStimChange_perStim;

    
    
    t = dbGetStimulusTimes(Gid);
    T = sort(t(:));
    stim_av_t = zeros(1, nStimPres);
    for i = 1:nStimPres-1
        stim_av_t(i) = mean( T( idx_firstStimFrame(i):idx_firstStimFrame(i+1)-1 ) );
    end
    
    
    t_perStim = median(diff(stim_av_t));
    
    
    
    %%
    
    uStims = unique(frameStimIds);
    

    %%
    nOri = length(uOri);
    nSp = length(uSp);
    nUStims = nOri*nSp;
    meanRecurrencTime_allStim = zeros(1, nUStims);

    %%
    [frameOriIds, frameSpfIds] = ind2sub([nOri, nSp], frameStimIds);
       
    %%
    if allowOriSpfRecurTimeSeparate
         for ori_i = 1:nOri
             for spf_i = 1:nSpf
                idx_stim_i = find( stimIds == stim_i);
            
            
            
            stimTs = stim_av_t(idx_stim_i); %#ok<FNDSB>
            meanRecurTime_i = mean(diff(stimTs));
             
            meanRecurrencTime_allStim(stim_i) = meanRecurTime_i;
            end
        end
        
    else
        
        for stim_i = 1:nUStims
            idx_stim_i = find(stimIds == stim_i);
            
            
            
            stimTs = stim_av_t(idx_stim_i); %#ok<FNDSB>
            meanRecurTime_i = mean(diff(stimTs));
            meanRecurrencTime_allStim(stim_i) = meanRecurTime_i;
        end
    end
   
    recur_t_allStim = meanRecurrencTime_allStim;
    recur_t_av = mean(meanRecurrencTime_allStim);
    

%}