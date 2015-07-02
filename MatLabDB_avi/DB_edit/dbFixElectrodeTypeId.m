function dbFixElectrodeTypeId
    return;
    
    
    hnd = dbOpenExpDb;
%     fieldnames = {'TBL_DATA_FILES.DATAFILE_ID', 'TBL_GROUPS.GROUP_ID', 'TBL_ANIMALS.ANIMAL_ID', 'TBL_PENETRATIONS.PENETRATION_ID', 'TBL_LOCATIONS.LOCATION_ID', 'TBL_ELECTRODES.ELECTRODE_ID', 'TBL_ELECTRODE_TYPES.ELECTRODE_TYPE_ID'};
% 
%     TBL_ANIMALS, TBL_ELECTRODES, TBL_PENETRATIONS, TBL_LOCATIONS, TBL_LOCS_FILES_LINKS, TBL_LOC_DEPTHS, TBL_DATA_FILES, TBL_GROUPS
%     
    fieldnames = {'TBL_GROUPS.GROUP_ID', 'TBL_ELECTRODE_TYPES.ELECTRODE_TYPE_ID'};
    
    T1 = {'TBL_ELECTRODE_TYPES', 'ELECTRODE_TYPE_ID', 'TBL_ELECTRODES'};
    T2 = {T1, 'TBL_ELECTRODES', 'ELECTRODE_ID', 'TBL_AP_ML_ZERO'};
    T3 = {T2, 'TBL_AP_ML_ZERO', 'AP_ML_ZERO_ID', 'TBL_PENETRATIONS'};
    T4 = {T3, 'TBL_PENETRATIONS', 'PENETRATION_ID', 'TBL_LOCATIONS'};
    T5 = {T4, 'TBL_LOCATIONS', 'LOCATION_ID', 'TBL_LOC_DEPTHS'};    
    T6 = {T5, 'TBL_LOC_DEPTHS', 'LOC_DEPTH_ID', 'TBL_LOCS_FILES_LINKS'};    
    T7 = {T6, 'TBL_LOCS_FILES_LINKS', 'DATAFILE_ID', 'TBL_DATA_FILES' };
    T8 = {T7, 'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_GROUPS'};    
    joinedTables = T8;
% 
%     %             criterea = {'TBL_ELECTRODES.ELECTRODE_TYPE_ID', 2};
%     criterea = [];
    [Gids_tbl, electTypeIds_tbl] = getFieldsFromDatabaseTable(hnd, fieldnames, joinedTables);
% 

%     [uElectrodeId, idx_first] = unique(electId, 'first');
%     uGid = unique(Gid)
%     for gi = 1:length

    putativeSingleElectGids = Gids_tbl(electTypeIds_tbl==1);

    allGids = dbGetStimulusGids;
    allGids = putativeSingleElectGids(1:30);
    electrodeTypes_spk = zeros(size(allGids));
    electrodeTypes_db = zeros(size(allGids));
    rs = zeros(size(allGids));
    for i = 1:length(allGids)
        Gid = allGids(i);
        Gid_tbl_idx = find(Gids_tbl == Gid, 1);
        electrodeTypes_db(i) = electTypeIds_tbl(Gid_tbl_idx);
        groupSpikes = dbGetSpikes(Gid, [], [], 1);
        electrodeData = groupSpikes(:,3:end);
        nElectrodes = size(electrodeData,2);
        if nElectrodes == 1
            electrodeTypes_spk(i) =1;
        else
            withinElecVar = mean( std(electrodeData, [], 1) );
            betweenElecVar = mean( std(electrodeData, [], 2) );
            r = betweenElecVar / withinElecVar;
            rs(i) = r;
            if r < .05
                electrodeTypes_spk(i) =1;
            else
                electrodeTypes_spk(i) =2;
            end
        end
    end
        
    

  3;  
    
%     LocID_tbl =



    quickCorrect = true;  % set to true to just skip to the ones we know need to be corrected. false to recheck everything.

end