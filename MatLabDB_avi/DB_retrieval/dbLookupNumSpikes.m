function varargout = dbLookupNumSpikes(Gid, cellId)
    persistent spikeCounts;
    % a quick way to see how many spikes are in a given group/cell
    % The table 'spikeCounts' is loaded (or created and saved, if it doesn't exist)
    % so that the number of spikes in Group i, cell j, is found in
    % spikeCounts(i,j);
    
    if isempty(spikeCounts) 
        path = getName('MatlabDB_path');
        filename = [path 'dbNumSpikesInGroups.mat'];
        if exist(filename, 'file')
            S = load(filename);
            spikeCounts = S.spikeCounts;
        else
%             hnd = dbOpenExpDb;
%             allGroupIds = getFieldsFromDatabaseTable(hnd, 'GROUP_ID', 'TBL_GROUPS');
            allGroupIds = dbGetStimulusGids;

            maxGroupId = max(allGroupIds);
            nMaxCellsPerGroup = 12;
            spikeCounts = -ones(maxGroupId, nMaxCellsPerGroup);

            progressBar('init-', length(allGroupIds));
            for Gid_i = 1:length(allGroupIds);
                progressBar(Gid_i);

                Gid = allGroupIds(Gid_i);            
                
                
                try
                    groupSpikes = dbGetSpikes(Gid);
                catch
                    groupSpikes = [];
                end
                
                if ~isempty(groupSpikes)
                    [spikingCells, nSpikes] = uniqueCount(groupSpikes(:,2));
                    spikeCounts(Gid, [spikingCells+1] ) = nSpikes;
                end
                
            end
            
            save(filename, 'spikeCounts');
            
        end        
    end

    switch nargin
        case 0
            varargout = {spikeCounts};
        case 1
            groupSpikes = spikeCounts(Gid, :);
            cellIds = find(groupSpikes > 0);
            nSpikes = groupSpikes(cellIds);
            varargout = {cellIds-1, nSpikes};  %-1 b/c first one can be zero, whereas indexing is from 1
        case 2
            groupSpikes = spikeCounts(Gid, :);
            nSpikes = groupSpikes(cellId);
            if ~(nSpikes > 0)
                error(['Group ' num2str(Gid) ' does not have a cell # ' num2str(cellId) '.']);
            end            
            varargout = {nSpikes};
    end
            
    




end