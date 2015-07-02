function [stimulusType subType] = getStimulusIdForDid(Did)
    hnd = dbOpenExpDb;

    stim_id = getFieldsFromDatabaseTable(hnd, 'STIMULUS_TYPE_ID', 'TBL_DATA_FILES', {'DATAFILE_ID', Did});
    stim_txt = getFieldsFromDatabaseTable(hnd, 'TXT_STIMULUS_TYPE', 'TBL_STIMULUS_TYPES', {'STIMULUS_TYPE_ID', stim_id});
    subType = stim_txt{1};
    
    switch subType
        case {'Single Grating', 'Flashed Grating Batch', 'Orientation Batch', 'Spatial Frequency Batch', 'Temporal Frequency Batch', 'Free Grating Batch'}
            stimulusType = 1;
        case 'Noise'
            stimulusType = 2;
        case 'M-sequence'
            stimulusType = 3;
        case {'Movie', 'Movie Batch'}
            stimulusType = 4;
    end
    
end