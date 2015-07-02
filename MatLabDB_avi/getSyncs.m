function [syncs, startTicks, endTicks] = getSyncs(Gid, outputType)

    syncsFileName = getFileName('syncs', Gid);

    redo = 0;
    if exist(syncsFileName, 'file') && ~redo
        S = load(syncsFileName);
        syncs = S.syncs;
        startTicks = S.startTicks;
        endTicks = S.endTicks;

        if exist('outputType', 'var') && ~isempty(outputType) && ~strcmp(outputType, 'tick')
            sd = siteDataFor(Gid);
            Did = sd.Did;                
            syncs      = dbConvertTimeMeasures(Did, syncs, 'tick', outputType);
            startTicks = dbConvertTimeMeasures(Did, startTicks, 'tick', outputType);
            endTicks   = dbConvertTimeMeasures(Did, endTicks, 'tick', outputType);
        end
    else
        if nargin < 2 
            outputType = [];
        end            
        syncs = dbGetSyncs('Gid', Gid, outputType);
        hnd = dbOpenExpDb;
        
        sd = siteDataFor(Gid);
        Did = sd.Did;                
        [startTicks, endTicks] = dbGetTbTe(hnd, Did);
        
        S.syncs = syncs;
        S.startTicks = startTicks;
        S.endTicks = endTicks;
        
        assert(all(binarySearch(syncs, startTicks, 1, 0)));
        assert(all(binarySearch(syncs, endTicks, 1, 0)));
        assert(all(startTicks < endTicks));
        
        save(syncsFileName, '-struct', 'S', '-v6');
    end
        

end