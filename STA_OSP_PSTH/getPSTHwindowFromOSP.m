function [timeWindow, windowProfile, PSTH_bins, PSTH_vals] = getPSTHwindowFromOSP(Gid, cellId,  OSP, oris, sps, phs,  extFrameLen_ms)
                
    % 1. use OSP to determine which frames to look at (in
    % extended fashion) -> get spkTsRelToFrame (good frames)
    bestOriSpPhs = findBestOriSpPhs(OSP, oris, sps, phs);

    stimulatoryFrameInds = findFrameIndsWithCertainOSP(Gid, bestOriSpPhs);
    spkTsRelToFrame_ms = getParsedSpikes('timing', Gid, cellId, [], extFrameLen_ms, stimulatoryFrameInds);

    % 2. use spkTsRelToFrame to get a better PSTH --> and a
    % better timeWindow.
    [PSTH_bins, PSTH_vals] = calcPSTH( spkTsRelToFrame_ms, extFrameLen_ms, length(stimulatoryFrameInds) );

    [timeWindow, windowProfile] = getBestTimeWindowFromPSTH(PSTH_bins(:), PSTH_vals(:));
end