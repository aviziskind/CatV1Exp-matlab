function [stimulusType subType] = getStimulusTypeForDid(Did)
    hnd = dbOpenExpDb;

    stim_id = getFieldsFromDatabaseTable(hnd, 'STIMULUS_TYPE_ID', 'TBL_DATA_FILES', {'DATAFILE_ID', Did});
    stim_txt = getFieldsFromDatabaseTable(hnd, 'TXT_STIMULUS_TYPE', 'TBL_STIMULUS_TYPES', {'STIMULUS_TYPE_ID', stim_id});
    stimulusType = stim_txt{1};
    subType = '';
    
    switch stimulusType
        case {'Single Grating', 'Flashed Grating Batch', 'Orientation Batch', 'Spatial Frequency Batch', 'Temporal Frequency Batch', 'Free Grating Batch'}
            subType = stimulusType;
            stimulusType = 'Grating';
        case 'Noise'
            subType = '';
        case 'M-sequence'
            stimulusType = 'Mseq';
            subType = '';
        case {'Movie', 'Movie Batch'}
            % Note: There are only 28 'Movie' instances (vs. 1947 'Movie Batch'
            % instances.) all of which are Band_Noise movies. So I just rename all 'Movie
            % Batch' types to 'Movie' types, since 'Movie Batch' isn't a
            % very useful separate designation.
            stimulusType = 'Movie';
            if nargout > 1 % only bother doing this if requested
                subType = getMovieType('Did', Did);
            end                
    end
    
    % Types {SubTypes}
    %     Grating {Single / Flash / Orientation / Spatial Frequency / Temporal Frequency / Free Grating Batch}
    %     Noise {''}
    %     Mseq {''}
    %     Movie {'Oriented_Bars', 'Noise__Band', 'Noise__Sparse', 'Flashed_Gratings', 'Natural_Scenes'}
    
end