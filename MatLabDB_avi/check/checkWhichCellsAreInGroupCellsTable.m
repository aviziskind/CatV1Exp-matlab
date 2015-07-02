function checkWhichCellsAreInGroupCellsTable
    hnd = dbOpenExpDb;
    [GC_Gids, GC_cellIds, GC_nSpikes] = getFieldsFromDatabaseTable(hnd, {'GROUP_ID', 'LNG_CLUSTER_NO', 'LNG_N_SPIKES'}, 'TBL_GROUP_CELLS');
    
    G_Gids = getFieldsFromDatabaseTable(hnd, {'GROUP_ID'}, 'TBL_GROUPS');
    
    % total of 831 groups in GROUP_CELLS table
    % 602/831 are Grating  (total grating: 1379)
    % 226/831 are Noise    (total noise: 226)
    % 3  /831 are mseq     (total mseq: 7)
    
    
end