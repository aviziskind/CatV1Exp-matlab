function dbRemoveGroupsWithNoSpikes

    hnd = dbOpenExpDb;
    
    allGroupIds = getFieldsFromDatabaseTable(hnd, 'GROUP_ID', 'TBL_GROUPS', [], 'GROUP_ID');
    nGroups = length(allGroupIds);
    grpsNoSpikes = false(1, nGroups);
    progressBar('init-', nGroups, 60);
    for i = 1:nGroups        
        progressBar;
        spk_field = getFieldsFromDatabaseTable(hnd, 'MEM_SPIKETIMES_MTX', 'TBL_GROUPS', {'GROUP_ID', allGroupIds(i)});        
        grpsNoSpikes(i) = ~iscell(spk_field) && isnan(spk_field);
    end
    
    fprintf('There are %d (out of %d) groups without any spikes.\n', nnz(grpsNoSpikes), nGroups);

    progressBar('init-', nnz(grpsNoSpikes), 60);
    for i = find(grpsNoSpikes);        
        progressBar;
        nCellSpks = dbLookupNumSpikes(allGroupIds(i));
        assert(isempty(nCellSpks));        
        removeRecordFromDatabaseTable(hnd, 'TBL_GROUPS', {'GROUP_ID', allGroupIds(i)})
    end
        
    
end