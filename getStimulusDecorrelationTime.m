function [decor_time_sec, decor_time_frames] = getStimulusDecorrelationTime(Gid)

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

    [frameStimIds] = getStimulusFrameSequence(Gid, 'OS');
    nUniqueStimuli = max(frameStimIds);

    if ~isFlashed
        idx_firstStimFrameEachPres =  find ( diff([frameStimIds;0]) ~= 0 );

        nFramesPerPres = median( diff(idx_firstStimFrameEachPres) );
        groupedFrameStimIds = frameStimIds(idx_firstStimFrameEachPres)';
        dt_use = dt_use * nFramesPerPres;        
    else
        groupedFrameStimIds = frameStimIds;
    end

    
    % circDist is appropriate for orientations. abs(dist) is better for
    % spatial frequencies, but can use circDist for both (results in
    % slight overestimation of period for spf, but this is ok).
    allFrameStimIdChanges = circDist(groupedFrameStimIds(1:end-1), groupedFrameStimIds(2:end), nUniqueStimuli);
    
    allFrameStimChanges = allFrameStimIdChanges / nUniqueStimuli;

   
    meanFrameStimChange = mean(allFrameStimChanges);
    decor_time_frames = 1/meanFrameStimChange;
    decor_time_sec =  decor_time_frames * dt_use;

end

