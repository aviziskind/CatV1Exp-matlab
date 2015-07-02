function plotSyncsAndSpikes(Gid)
    Did = dbLookup('Did',  'Gid', Gid);
    N = 100;
    
    outputUnits = 'sec';
    figure(10);
    h1 = subplot(2,1,1);
    syncs = dbGetSyncs('Did', Did, outputUnits, 1);
    if isempty(syncs)
        disp('No Syncs');
        return;
    end
    hist(syncs, N);
%     hh = findobj(gca,'Type','patch');
%     set(hh,'FaceColor','w','EdgeColor','k')
    
    h2 = subplot(2,1,2);
    cellIds = dbLookupNumSpikes(Gid, [], outputUnits);
    if isempty(cellIds)
        return;
    end
    grpSpikes = dbGetSpikes(Gid);
    spks = cellfun(@(c_id) grpSpikes( grpSpikes(:,2)==c_id ,1), cellIds, 'un', 0);
%     if isempty(grpSpikes)
%         disp('No spikes')
%         return;
%     end
%     spks = cell(1,length(cellIds));
    ns = zeros(length(cellIds), N);
%     spks = cell(1,length(cellIds));
    l = 0; h = -inf;
    for i = 1:length(cellIds)
%         spks{i} = dbGetSpikes(grpSpikes, cellIds(i));                
%         spks{i} = dbConvertTimeMeasures(Did, spks{i}, 'tick', outputUnits);
        l = min(l, min(spks{i}));
        h = max(h, max(spks{i}));
    end
    for i = 1:length(cellIds)
        [ns(i,:), xout] = hist(spks{i},linspace(l,h,N));
    end    
    bar(xout, ns', 'stacked'); 
    xlim([0 h]);
    
    matchAxes('X', h1, h2);

end