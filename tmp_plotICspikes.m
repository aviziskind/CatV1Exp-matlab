function tmp_plotICspikes(Gid)

    fn = getFileName('startEndIC', Gid);
    S_spk = load(fn);
    
    [channelMean, channelVar] = getIntraCellularChannelMeanAndVariance(Gid);
    channelStd = sqrt(channelVar);        
    
    figure(66);
    spkH = S_spk.spikeHeights;
    spkW = double(S_spk.spikeEnds-S_spk.spikeStarts);
    subplot(1,2,1); hist(spkH, 50);
    xlabel('# std deviations above background'); ylabel('Count'); 
    max_w = max(spkW); binC = 1:max_w; binE = binCent2edge(binC);
    subplot(1,2,2); hist(spkW, binC); xlim([0, max_w+1])
    xlabel('spike width'); 
    3;
    title(sprintf('Intracellular spikes for Gid = %d (# spikes = %d). Channel Std = %.2f', Gid, length(spkH), channelStd));

    


end