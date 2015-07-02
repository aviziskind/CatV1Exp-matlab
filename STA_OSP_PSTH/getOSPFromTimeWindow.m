function [OSP, oris, sps, phs] = getOSPFromTimeWindow(Gid, cellId,   timeWindow, windowProfile)
                
    % 1. get estimate of relContrOfFrameToSpike using timeWindow
    estRelContrOfFrameToSpike = getParsedSpikes('frame', Gid, cellId, timeWindow, windowProfile);

    % 2. use estimated relContrOfFrameToSpike to get OSP
    [OSP, oris, sps, phs] = getOriSpfPhaseProfile(Gid, estRelContrOfFrameToSpike);

end