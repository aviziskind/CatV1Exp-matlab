function dbRemoveEntriesWithNoSyncs
    
    hnd = dbOpenExpDb;
    tableName = 'TBL_SYNCHRO_GROUPS';
    [syncIds syncDids] = getFieldsFromDatabaseTable(hnd, {'SYNCHRO_GROUP_ID', 'DATAFILE_ID'}, tableName);

    % Remove Sids that have empty sync times
    SidsNoSyncs = false(size(syncIds));
    progressBar('init=', length(syncIds));
    for i = 1:length(syncIds)
        progressBar(i);
        Sid = syncIds(i);
        syncs_field = getFieldsFromDatabaseTable(hnd, 'MEM_SYNCHROTIMES_MTX', tableName, {'SYNCHRO_GROUP_ID', Sid});
        
        if ~iscell(syncs_field) && isnan(syncs_field)
            SidsNoSyncs(i) = true;
            disp(['Remove: Sid = ' num2str(Sid)]);
        end
    end
    
    SidsToRemove = find(SidsNoSyncs);
    disp('Removing following Sids:');
    disp(SidsToRemove(:)');
%     removeRecordFromDatabaseTable(hnd, tableName, criteria);


end
