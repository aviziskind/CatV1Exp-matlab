function spike_stimIds = getSpikesStimIds(Gid)
    
    fn = getFileName('spikeStimId', Gid);
    if exist(fn, 'file')
        S = load(fn);
        spike_stimIds = S.spike_stimIds;
    else            
        spike_stimIds = calculateSpikeFrameIds(Gid);
        save(fn, 'spike_stimIds');        
    end
    


end


function stimIdsForEachSpike = calculateSpikeFrameIds(Gid)

    psthWindow = [0 0];
    [framesTimesForEachSpike, framesIdsForEachSpike] = getParsedSpikes('spikes', Gid, 100, psthWindow); 
    spike_firstFrameId = cellfun(@(x) x(2), framesIdsForEachSpike);
%     spike_firstTime = cellfun(@(x) x(2), framesTimesForEachSpike);
    
    frameStimIds = getStimulusFrameSequence(Gid, 'OSPt');
%%
    idx = (spike_firstFrameId > 0);
    
%     stimIdsForEachSpike = frameStimIds(framesIdsForEachSpike);        
    stimIdsForEachSpike = zeros(size(spike_firstFrameId));
        
    stimIdsForEachSpike(idx) = frameStimIds(spike_firstFrameId(idx));               

end















