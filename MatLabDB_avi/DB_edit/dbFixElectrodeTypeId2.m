function dbFixElectrodeTypeId2
%     return;
    
    
    hnd = dbOpenExpDb;
%     fieldnames = {'TBL_DATA_FILES.DATAFILE_ID', 'TBL_GROUPS.GROUP_ID', 'TBL_ANIMALS.ANIMAL_ID', 'TBL_PENETRATIONS.PENETRATION_ID', 'TBL_LOCATIONS.LOCATION_ID', 'TBL_ELECTRODES.ELECTRODE_ID', 'TBL_ELECTRODE_TYPES.ELECTRODE_TYPE_ID'};
% 
%     TBL_ANIMALS, TBL_ELECTRODES, TBL_PENETRATIONS, TBL_LOCATIONS, TBL_LOCS_FILES_LINKS, TBL_LOC_DEPTHS, TBL_DATA_FILES, TBL_GROUPS
%     
%     fieldnames = {'TBL_GROUPS.GROUP_ID', 'TBL_ELECTRODE_TYPES.ELECTRODE_TYPE_ID'};
    old_value = 2;
    new_value = 1;
    [electrode_id] = getFieldsFromDatabaseTable(hnd, 'ELECTRODE_TYPE_ID', 'TBL_ELECTRODES', {'ELECTRODE_ID', 139});
    if electrode_id == old_value
    
        updateValueInDatabaseTable(hnd, new_value, 'ELECTRODE_TYPE_ID', 'TBL_ELECTRODES', {'ELECTRODE_ID', 139});
    end
    
  
end