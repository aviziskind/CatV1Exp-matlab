function dbDoSmallEdits
    hnd = dbOpenExpDb;

    % add 1 second to time-stamp on one experiment, so that matches standard DTM mask;
    updateValueInDatabaseTable(hnd, '9/25/2002 00:00:01', 'DTM_CREATED', 'TBL_GRATING_PRES', {'GRATING_PRES_ID', 166228}, 'DATE')
   
    % change name of stimulus type #3: 'Flash Grating Batch' ==> 'Flashed Grating Batch'. 
    updateValueInDatabaseTable(hnd, 'Flashed Grating Batch', 'TXT_STIMULUS_TYPE', 'TBL_STIMULUS_TYPES', {'STIMULUS_TYPE_ID', 3});    
    
end
