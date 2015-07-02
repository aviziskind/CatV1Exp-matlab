function degreesPerPixel = dbGetDegreesPerPixel(Did)
    hnd = dbOpenExpDb;
    
    stimType = getStimulusTypeForDid(Did);
    tableName = getDatabaseTableForDid(Did, stimType);

    degreesPerPixel = getFieldsFromDatabaseTable(hnd, 'DBL_DEGREES_PER_PIXEL', tableName, {'DATAFILE_ID', Did});
    if length(degreesPerPixel) > 1
        assert(max(abs(diff(degreesPerPixel))) < 1e-4);
    end
    if any(strcmp(stimType, {'Movie', 'Noise'}))
        pixelsPerBlock = getFieldsFromDatabaseTable(hnd, 'LNG_BLOCK_HEIGHT_PIX', tableName, {'DATAFILE_ID', Did}); %always have (block height = block width)
    else  % for grating (& mseq) stimuli, don't have this field in the table. I assume the value is 1.
        pixelsPerBlock = 1;
    end
    degreesPerPixel = degreesPerPixel(1)*pixelsPerBlock(1);

end
