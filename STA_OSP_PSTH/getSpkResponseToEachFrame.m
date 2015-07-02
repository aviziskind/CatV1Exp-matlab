function relContrOfFrameToSpike = getSpkResponseToEachFrame(Gid, cellId, timeWindow_ms, windowProfile, outputOrder )

    if nargin < 5
%         outputOrder = 'stimOrder';
        outputOrder = 'original';
    end

    windowProfile = windowProfile/sum(windowProfile);    
    
    histOpts = struct('psthWindow_ms', timeWindow_ms, 'trialGrouping', 'individual');
    [bins, allHistVals] = dbGetCellSpkStimHists(Gid, cellId, histOpts);
    [nBins, nStim, nTrials] = size(allHistVals);
    
    windowProfile = windowProfile(:);
    if ~isempty(windowProfile)
        allHistVals = bsxfun(@times, allHistVals, windowProfile(:));
    end
        
    response_stim_trials = reshape(sum(allHistVals, 1), [nStim, nTrials])';
    
%     imagesc(mean(reshape(sum(response_stim_trials, 1), [36, 10, 8]),3))    
    
    switch outputOrder
        case 'original'
            relContrOfFrameToSpike = zeros(1,nStim*nTrials);
            frameIdSequence = getStimulusFrameSequence(Gid);
            [uStim, stimIdx] = uniqueList(frameIdSequence(:));
            origIdx = [stimIdx{:}];
            relContrOfFrameToSpike(origIdx(:)) = response_stim_trials;                
        case 'stimOrder',
            relContrOfFrameToSpike = response_stim_trials;
        otherwise
            error('Unknown order');            
    end
    relContrOfFrameToSpike = double(relContrOfFrameToSpike);   
    
end