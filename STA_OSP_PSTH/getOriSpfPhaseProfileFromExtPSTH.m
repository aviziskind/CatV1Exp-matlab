function OSP_PSTH = getOriSpfPhaseProfileFromExtPSTH(Gid, cellId)

    frameLength_ms = getFrameLength('Gid', Gid);    
    extFrameLengthApprox_ms = 90;
    nFramesPerExtFrame = floor(extFrameLengthApprox_ms/frameLength_ms);
    extFrameLength_ms = nFramesPerExtFrame * frameLength_ms;

    binSizeApprox_ms = 5;
    nBinsPerFrame = round(frameLength_ms / binSizeApprox_ms);
    nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;    

%     [relContrOfFrameToSpike, nSpikesEachFrame] = getParsedSpikes(Gid, cellId);
%     nSpikesEachFrame = [ nSpikesEachFrame{:} ];
%     meanSpikesPerFrame = sum(nSpikesEachFrame) / length(nSpikesEachFrame);
%     meanSpikesPerSecond = meanSpikesPerFrame / frameLength_sec;
    
    % get list of OSPs, and parse spikes
    [allOri_deg, allSp_pix, allPhase_deg] = getOriSpPhaseForEachStimFrame(Gid);
    allStimuliOSPs = [allOri_deg(:), allSp_pix(:), allPhase_deg(:)];
    uOriSpPhStims = unique(allStimuliOSPs, 'rows');
    nOriSpPhStims = size(uOriSpPhStims, 1);
    allStimuliOSPids = zeros( length(allOri_deg), 1);
    for stim_i = 1:nOriSpPhStims;
        allStimuliOSPids(  findRows(uOriSpPhStims(stim_i,:), allStimuliOSPs)  ) = stim_i;
    end

	spkTsRelToFrame_ms = getParsedSpikes('timing', Gid, cellId, [], extFrameLength_ms);

    % 1. Get PSTH vals
    PSTH_vals = zeros(nBinsPerExtFrame, nOriSpPhStims);
    for stim_i = 1:nOriSpPhStims
        frameIndsForStimI = find(allStimuliOSPids == stim_i);
        [PSTH_bins, PSTH_vals(:,stim_i)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI), extFrameLength_ms, length(frameIndsForStimI), [], nBinsPerExtFrame);
    end
    PSTH_means = mean(PSTH_vals, 1);
    stimIndsInOrder = ord(PSTH_means, 'descend');
    gridSubPlot(4,5, 500);
    for stim_i = 1:100;%length(stimIndsInOrder)
        gridSubPlot;
        plotThisPSTH(PSTH_bins, PSTH_vals(:,stimIndsInOrder(stim_i)), extFrameLength_ms);
        title(['(' outOf(stim_i, nOriSpPhStims) ') ' num2str(PSTH_means(stimIndsInOrder(stim_i)))]);
    end

    OSP_PSTH = reshape(PSTH_means, [nSp, nOri, nPh]);

%     figure(600);
%     subplot(1,2,1);
%     imagesc(OSP_PSTH);
% 
%     subplot(1,2,2);
%     imagesc(sum(OSP,3));

end