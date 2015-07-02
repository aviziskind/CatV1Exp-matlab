function checkDataFileIds
%check that the datafile_ids listed in TBL_GROUPS match the ones in the
%separate tables
    hnd = dbOpenExpDb;
    stimTypes = {'movie', 'noise', 'grating', 'mseq'};

    allStimGidsC = cell(1,4);
    allStimDidsC = cell(1,4);
    for ti = 1:length(stimTypes)
        stimType = stimTypes{ti};
        S = load(['cellsGroups_' stimType '_all_DB']);
        groupData = S.([stimType 'Groups_all']);
        stimDids =  [groupData(:).Did];
        stimGids =  [groupData(:).Gid];
        
        allStimDidsC{ti} = stimDids;
        allStimGidsC{ti} = stimGids;
    end
    allStimGids = [allStimGidsC{:}];
    allStimDids = [allStimDidsC{:}];
    allStimGidDidPairs = unique([allStimGids(:), allStimDids(:)], 'rows');

    [gpGids, gpDids] = getFieldsFromDatabaseTable(hnd, {'GROUP_ID', 'DATAFILE_ID'}, 'TBL_GROUPS');

    allGpGidDidPairs = unique([gpGids(:), gpDids(:)], 'rows');

%     stimPairsNotInTable = setdiff(allStimGidDidPairs, allGpGidDidPairs, 'rows');   % none
    tablePairsNotInStims = setdiff(allGpGidDidPairs, allStimGidDidPairs, 'rows');

    %remove from TBL_GROUPS
    for i = 1:length(tablePairsNotInStims)
        Gid = tablePairsNotInStims(i,1);
        removeRecordFromDatabaseTable(hnd, 'TBL_GROUPS', {'GROUP_ID', Gid});
    end
    
    %remove from TBL_DATA_FILES
    for i = 1:length(tablePairsNotInStims)
        Did = tablePairsNotInStims(i,2);
        removeRecordFromDatabaseTable(hnd, 'TBL_DATA_FILES', {'DATAFILE_ID', Did});
    end

    
    end

%     extra GIDs / DIDs in the TBL_GROUPS / TBL_DATA_FILES tables that are not in the stimulus tables
%         GID          DID
%          500         485  crash
%          890        1069  'junk'
%          891        1070  'junk'
%          896        1080  'junk'
%          897        1081  'junk'
%          909        1104  'junk'
%          910        1105  'junk'
%          930        1141 ?
%         1280        1509  ?
%         1288        1524  ?
%         1446        1721 ?
%         1501        1776  crash
%         1509        1784  crash
%         1552        1829  crash
%         1561        1838  ?
%         1563        1839
%         1565        1840
%         1567        1841
%         1571        1842
%         1573        1843
%         1575        1844
%         1579        1846
%         1581        1848
%         1591        1857
%         1592        1858
%         1593        1859 ???
%         1870        2014  crash
%         1871        2014 crash
%         2364        2314 crash
%         2365        2314 crash
%         2641        2462 ?
%         2642        2462 ? 
%         2937        2610 ?
%         2938        2610 ?
%         3323        2803 ?
%         3324        2803
%         4286        3297
%         4287        3297  ?
%         5064        3686
%         5065        3686
%         5264        3786
%         5265        3786



















