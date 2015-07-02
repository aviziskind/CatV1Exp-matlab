function PSTHvals = generatePSTHs(r, stimIds, nBinsPerFrame, nFramesPerExtFrame)

    nTotalPres = length(stimIds);
    nStimuli = length(unique(stimIds));

    nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;
    spikesPerBin = zeros(1, (nTotalPres+1)*nBinsPerFrame + nBinsPerExtFrame);

    for i = 1:nTotalPres
        inds = (i-1)*nBinsPerFrame +1:(i-1)*nBinsPerFrame  + nBinsPerExtFrame;            
        spikesPerBin(inds) = spikesPerBin(inds) + r(:,stimIds(i))';
    end        
    frame2BinInds = @(frm_i) (frm_i-1)*nBinsPerFrame  + 1 : frm_i*nBinsPerFrame ;

    PSTHvals = zeros(nFramesPerExtFrame*nBinsPerFrame , nStimuli);
    % get PSTH for each stim type
    for stim_i = 1:nStimuli
        frameIndsForStimI = find(stimIds == stim_i);

        for frm_i = 1:nFramesPerExtFrame
            relevantFrames = frameIndsForStimI+frm_i-1;
            for relv_frm = relevantFrames
                PSTHbins = frame2BinInds(frm_i);
                stimBins = frame2BinInds(relv_frm);
                PSTHvals( PSTHbins , stim_i) = PSTHvals( PSTHbins, stim_i) + spikesPerBin(stimBins)' / length(frameIndsForStimI);
            end
        end

    end
end
