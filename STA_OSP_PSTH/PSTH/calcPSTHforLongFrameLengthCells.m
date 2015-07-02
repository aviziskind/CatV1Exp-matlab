function [PSTH_bins, PSTH_vals, spkTsRelToFrame_ms, bckgRate, meanFiringRate] = calcPSTHforLongFrameLengthCells(Gid, cellId, bkgrSkip_ms)
    if nargin < 3
        bkgrSkip_ms = 200;
    end
    frameLength_ms = getFrameLength('Gid', Gid, 'ms');

    % 1. Get Position of spike times relative to the frames
    [spkTsRelToFrame_ms, bckgRate, meanFiringRate] = getParsedSpikes('timing', Gid, cellId, bkgrSkip_ms );

    % 2. Use these spike times to calculate the PSTH
    [PSTH_bins, PSTH_vals] = calcPSTH( spkTsRelToFrame_ms, frameLength_ms);
    
end