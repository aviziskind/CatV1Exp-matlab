Gid = 4470;
grp = siteDataFor('Gid', Gid, 1);
cellIds = grp.cellIds;
tWindow = [-.25, 1];
for i = 2;%2:length(cellIds)
    figure(i);
    [wvfm, t_ms] = getSpikeWaveforms(Gid, cellIds(i));
    nspk = size(wvfm,3);
    plotSpikeWaveforms(t_ms, wvfm);    
    title(sprintf('Gid = %d. cellId = %d (%d spikes)', Gid, cellIds(i), nspk))
end