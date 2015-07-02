function [bckgData, bckgMeanRate] = getBackgroundSpikes(Gid, cellId, bckgSkip_ms, nAv)

    persistent allBckgSpikes saveCount
        
    bckgSpikes_file = [CatV1Path 'MatLabDB_avi' filesep 'allBckgSpikes' curMatchDB('') '.mat'];
        
    redo_all = 0;
    redo_current = 0;
    saveCountSpacing = 50;
    chunkSize_ms = (4+(1/6));
    
    convertOutputToHz = 1;

    if ~exist('bckgSkip_ms', 'var') || isempty(bckgSkip_ms)
        bckgSkip_ms = 250;
    end
    
    if ~exist('nAv', 'var') || isempty(nAv)
        nAv = 0;
    end    
    
    if strcmp(Gid, 'save') 
        if (saveCount > 0)
            save(bckgSpikes_file, 'allBckgSpikes', '-v6');                
            saveCount = 0;
        end
        return;
    end        
    
    if isempty(allBckgSpikes)        
        if exist(bckgSpikes_file, 'file') && ~redo_all
            S_file = load(bckgSpikes_file);
            allBckgSpikes = S_file.allBckgSpikes;
        else
            allBckgSpikes = struct;
        end        
        saveCount = 0;
    end
    
    getAllGrpSpk = cellId == 100;
                
    group_fld_name = sprintf('BckgSpikes_Gid_%d_cellId_%d', Gid, cellId);
    group_fld_name = strrep(group_fld_name, '-1', 'n1');
    
    if (~isfield(allBckgSpikes, group_fld_name) || redo_current) && ~getAllGrpSpk
        allChunks_C = calcBackgroundSpikes(Gid, cellId, chunkSize_ms);
        allChunks_C_comp = compressCell(allChunks_C);
        
        allBckgSpikes.(group_fld_name) = allChunks_C_comp;
        saveCount = saveCount + 1;
        
        if saveCount > saveCountSpacing
            save(bckgSpikes_file, 'allBckgSpikes', '-v6');        
            saveCount = 0;
        end                
    end
    
    if getAllGrpSpk 

        if ~curMatchDB                
            S_sort = load(getFileName('cells', Gid));
            cellIds = S_sort.uClustIds;
        else
            sd = siteDataFor('Gid', Gid, 1);
            cellIds = sd.cellIds;
        end        
        bckgData_C = getBackgroundSpikes(Gid, cellIds(1), 0, []); % don't remove or average yet.
        for i = 2:length(cellIds)
%             b2 = getBackgroundSpikes(Gid, cellIds(i));
            bckgData_C = cellfun(@plus, bckgData_C, getBackgroundSpikes(Gid, cellIds(i), 0, []), 'un', 0);
        end             
    else
    
        bckgData_C = decompressCell( allBckgSpikes.(group_fld_name) ); 
        bckgData_C = cellfun(@double, bckgData_C, 'un', 0);
    end
    
    
    if bckgSkip_ms > 0
        binSize_ms = 4+1/6;
        bckgSkip_bins = round(bckgSkip_ms/binSize_ms);
        bckgData_C = cellfun(@(x) x(bckgSkip_bins+1:end), bckgData_C, 'un', 0);        
    end
%%    
    nBkgSpks = sum([bckgData_C{:}]);        
    bkgTime_chunks = sum(cellfun(@length, bckgData_C));
    bkgTime_sec = bkgTime_chunks * chunkSize_ms / 1000;
    bckgMeanRate = nBkgSpks / bkgTime_sec;
    
    if (nAv > 0)
        %%
        bckgData = [bckgData_C{:}];
        nSamples = length(bckgData);
        nSamples_new = floor(nSamples/nAv);
        new_idxs = arrayfun(@(i1) i1+[1:nAv], [0:nSamples_new-1]*nAv, 'un', 0);                
        bckgData = cellfun(@(idxs) mean(bckgData(idxs)), new_idxs);                
        
        chunkSize_sec = chunkSize_ms/1000; % don't multiply by nAv , since took mean on previous line.

        if convertOutputToHz
            % convert from # spikes to spike rate in spikes/sec.
            bckgData = bckgData / chunkSize_sec;                        
            
            meanBckgData = mean(bckgData);
            if (meanBckgData > 0) && (bckgMeanRate > 0)
                diff_ratio = abs( 1 - (meanBckgData / bckgMeanRate) );
                abs_diff = abs(meanBckgData - bckgMeanRate);                
                assert(diff_ratio < .05 || abs_diff < .5);
            end
        end
        
    else
        bckgData = bckgData_C;
    end

            
        
%         bkgTime_ms = length(allSpikesInEachChunk)*chunkSize_ms;
    
    
%     if ischar(cellIds) && strcmp(cellIds, 'all')
%         cellIds = grpData.cellIds;
%     elseif ischar(cellIds) && strcmp(cellIds, 'cells')
%         cellIds = grpData.cellIds;
%         cellIds = cellIds(cellIds > 0);
%     end
        
        
end    

function allChunks_C = calcBackgroundSpikes(Gid, cellId, chunkSize_ms)
    BkgrSkip_ms = 0; % get all spikes.
    
    spikeTimes_ms = getSpikes(Gid, cellId, 'ms');        
        
%     syncTimes_ms0 = dbGetSyncs('Gid', Gid, 'ms');    
    
    [syncTimes_ms0, presBeginTime_ms, presEndTime_ms] = getSyncs(Gid, 'ms');
    
    sd = siteDataFor(Gid);
    nPres = length(sd.presIds);
%     Did = sd.Did;
%     hnd = dbOpenExpDb;
    lastTick_ms = sd.dataFileInfo.duration_sec * 1000;
    syncTimes_ms = [0; syncTimes_ms0; lastTick_ms]; % add my own ticks at start & end of experiment.
%     nSyncs = length(syncTimes_ms);        

%     
%     [beginTime_ticks, endTime_ticks] = dbGetTbTe(hnd, Did);
%     presBeginTime_ms = dbConvertTimeMeasures(Did, beginTime_ticks, 'tick', 'ms');
%     presEndTime_ms = dbConvertTimeMeasures(Did, endTime_ticks, 'tick', 'ms');
    
    presStartSyncIds = binarySearch(syncTimes_ms, presBeginTime_ms, [], 2);
    presEndSyncIds = binarySearch(syncTimes_ms, presEndTime_ms, [], 2);

    
    
        % (1) Pre-blank frames                

        allChunks_C = cell(1, nPres+1); % each inter-presentation interval, + beginning + end.
        
        preBlankInterval_chunks = [BkgrSkip_ms:chunkSize_ms:syncTimes_ms(2)];
        nPreBlankSpikes_chunks = elementsInRange(spikeTimes_ms, preBlankInterval_chunks, 'count');
        
        allChunks_C{1} = nPreBlankSpikes_chunks';        
        
        % (2) Inter-presentation intervals
        for iPres = 1:nPres-1            
            frmSyncId = presEndSyncIds(iPres);

%             interPresInterval        = [syncTimes_ms(frmSyncId)+BkgrSkip_ms,              syncTimes_ms(frmSyncId+1)];
%             if diff(interPresInterval) > 0
%                 nBkgSpikesInInterval = elementsInRange(spikeTimes_ms, interPresInterval, 'count');                
%                 nBkgSpks = nBkgSpks + nBkgSpikesInInterval;
%                 bkgTime_ms = bkgTime_ms + diff(interPresInterval);    
%             end            
            
            interPresInterval_chunks = [syncTimes_ms(frmSyncId)+BkgrSkip_ms : chunkSize_ms : syncTimes_ms(frmSyncId+1)];
            if length(interPresInterval_chunks) >= 2
                nBkgSpikesInInterval_chunks = elementsInRange(spikeTimes_ms, interPresInterval_chunks, 'count');                
                allChunks_C{iPres+1} = nBkgSpikesInInterval_chunks';
            end
            
        end        
        
        % (3) Post-blank frames
%         postBlankInterval = [syncTimes_ms(end-1) + BkgrSkip_ms,    syncTimes_ms(end)];
%         if diff(postBlankInterval) > 0
%             nPostBlankSpikes = elementsInRange(spikeTimes_ms, postBlankInterval, 'count');
%             nBkgSpks = nBkgSpks + nPostBlankSpikes;
%             bkgTime_ms = bkgTime_ms + diff(postBlankInterval);
%         end
        
        postBlankInterval_chunks = [syncTimes_ms(end-1) + BkgrSkip_ms : chunkSize_ms : syncTimes_ms(end)];        
        if length(postBlankInterval_chunks) >= 2
            nPostBlankSpikes_chunks = elementsInRange(spikeTimes_ms, postBlankInterval_chunks, 'count');
            allChunks_C{nPres+1} = nPostBlankSpikes_chunks';
        end
        
        allSpikesInEachChunk = [allChunks_C{:}];
        
        assert(all( ibetween(allSpikesInEachChunk, 0, 10)));
        allChunks_C = cellfun(@uint8, allChunks_C, 'un', 0);
        3;

        %% calculate background rate
%         nBkgSpks = sum(allSpikesInEachChunk);
%         bkgTime_ms = length(allSpikesInEachChunk)*chunkSize_ms;
%         bckgRate = nBkgSpks/bkgTime_ms * 1000; % from ms to sec
% 
%         chunkPerSec = 1000 / chunkSize_ms;        
%         allSpikeRates_Hz_eachChunk = allSpikesInEachChunk * chunkPerSec;
%         
%         bckgRate_mean = mean(allSpikeRates_Hz_eachChunk);
%         bckgRate_std  = std(allSpikeRates_Hz_eachChunk);        
%         bckgSamples =  allSpikeRates_Hz_eachChunk;
end

   

function s = compressCell(C)
    cell_lengths = cellfun(@length, C);
    C_cat = [C{:}];
    s = compress(C_cat);
    s.uVals = s.uVals(:)';
    s.cell_lengths = cell_lengths;
end

function C = decompressCell(s)
    cell_lengths = s.cell_lengths;    
    cumLengths = [0, cumsum(cell_lengths)];    
    C_cat = decompress(s);
    C = arrayfun(@(i1, i2) C_cat(i1:i2), cumLengths(1:end-1)+1, cumLengths(2:end), 'un', 0);    
end




% Compute:

%{

[allGids, allCellIds] = getAllGids;
progressBar('init-', length(allGids), 60);
for i = 1:length(allGids)
    bckgData = getBackgroundSpikes(allGids(i), allCellIds(i));    
    progressBar(i);
end
getBackgroundSpikes('save');

allGids = getAllGids;
progressBar('init-', length(allGids), 60);
for i = 1:length(allGids)
    bckgData = getBackgroundSpikes(allGids(i), 100);    
    progressBar(i);
end
getBackgroundSpikes('save');


%}