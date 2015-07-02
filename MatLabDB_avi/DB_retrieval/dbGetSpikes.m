function spikes = dbGetSpikes(Gid, cellId, outputType, keepTetrodeDataFlag, keepAllSpikesFlag)
    % Can return either all spikes for a group of cells, or just 1 cell.
    %       groupSpikes = dbGetSpikes(Gid, [])  % all cells in the group
    %       groupSpikes = dbGetSpikes(Gid, cellId)  % only this cell
    %       groupSpikes = dbGetSpikes(..., outputType) %outputType = 'ms'/'sec'/'tick' 
    %       groupSpikes = dbGetSpikes(..., outputType, keepTetrodeDataFlag)
    
    %%% GET GROUP SPIKE DATA

    hnd = dbOpenExpDb;

    if isscalar(Gid)
%         Did = dbLookup('Did',  'Gid', Gid);

        groupSpikeFile = getFileName('spikes', Gid, 1);
        if exist(groupSpikeFile, 'file')
            load(groupSpikeFile);
            groupSpikes = double(groupSpikes); %#ok<NODEF>
        else
            spk_field = getFieldsFromDatabaseTable(hnd, 'MEM_SPIKETIMES_MTX', 'TBL_GROUPS', {'GROUP_ID', Gid});
            if ~iscell(spk_field) && isnan(spk_field);
                error('db:noSpikes', 'No spikes for this site');
            end
            spk_str = cat( 1,spk_field{:} );
            if isnan(spk_str)
                error('db:noSpikes', 'No spikes for this site');
            end
            groupSpikes = edbStr2Mtx(spk_str); %dim = [nspikes x 6] % column1: spike times.  column2: which cell spiked.  columns 3-6: voltage recordings on each of the 4 electrodes
        end
        
    else
        groupSpikes = Gid;
    end
        
    
    % 0. default : remove spikes that don't meet current detection
    % criteria, unless provide a flag to keep them.
    keepAllSpikes = exist('keepAllSpikesFlag', 'var') && ~isempty(keepAllSpikesFlag) && (keepAllSpikesFlag);     
    if ~keepAllSpikes
        deleteIdx = deletedSpikes(Gid);
        groupSpikes(deleteIdx,:) = [];
    end    
    
    %%% OUTPUT  GROUP/CELL DATA

    % some groups (Gid = 2953 have some spikeTimes equal to zero (t==0). (weird).
    % Not sure where they came from, but for now, we just remove these spikes.
    if any(Gid == [2951, 2953, 2955, 2965, 2975  4378])   % 2953 & 2955 are very weird - both are ~half for cell0, ~half for cell4, and 1 for cell 3 (???)
        spkT0 = (groupSpikes(:,1) == 0);
        groupSpikes(spkT0,:) = [];        
    end
        
    [nSpk, nCol] = size(groupSpikes);    
    
    % 1. Select which columns to keep
    if exist('keepTetrodeDataFlag', 'var') && ~isempty(keepTetrodeDataFlag)        
        col_idx = 1:nCol;
    else
        col_idx = 1:2;
    end        
    
    % 3. Select which spikes to keep (and whether to keep cellId column).    
    spkCellIds = groupSpikes(:,2);    
    if ~exist('cellId', 'var') || isempty(cellId) || (cellId < 0)      % retrieve spikes for *all* cells (& keep cellId column)
        spk_idx = 1:nSpk;        
    elseif (cellId == 100)  % retrieve spikes for *all* cells (& don't keep cellId column)
        spk_idx = 1:nSpk;
        col_idx(col_idx == 2) = [];
    else
        % retrieve spikes for 1 specific cell (so can remove cell_id column)
        spk_idx = spkCellIds ==cellId;  
        col_idx(col_idx == 2) = [];
    end    
    
    spikes = groupSpikes(spk_idx, col_idx);
    

    if isempty(spikes)
        error('db:noSpikes', 'No spikes for this cell');
    end

    if exist('outputType', 'var') && ~isempty(outputType) && ~strcmp(outputType, 'tick')
        if isscalar(Gid)
            Did = dbLookup('Did',  'Gid', Gid);
            spikes(:,1) = dbConvertTimeMeasures(Did, spikes(:,1), 'tick', outputType);
        else
            error('Must supply input Gid to convert spike times from ticks to %s', outputType);
        end
    end

end
