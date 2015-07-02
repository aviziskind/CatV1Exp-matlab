function dbCorrectDegPerPixFields
    % the DegPerPix values for two sites are an order of magnitude larger than for the rest. 
    % (0.67 degPerPix instead of ~0.067 for most other sites.) .
    % I'm assuming this was just a typo and the correct number is just a factor of 10 lower.

    
    return; % update : it appears that these fields are actually correct... (!)
    
    Gids = [4946, 4950];
    Dids = [3627, 3629];
    
    hnd = dbOpenExpDb;
    

    for i = 1:length(Dids)
        Did = Dids(i);
        degPerPix_current = getFieldsFromDatabaseTable(hnd, {'DBL_DEGREES_PER_PIXEL'}, 'TBL_MOVIE_PRES', {'DATAFILE_ID', Did; 'LNG_PRESENT_NO', 1});
        
        if degPerPix_current < 0.1        
            degPerPix_amended = degPerPix_current*10;        
            updateValueInDatabaseTable(hnd, degPerPix_amended, 'DBL_DEGREES_PER_PIXEL', 'TBL_MOVIE_PRES', {'DATAFILE_ID', Did});
        end
        
    end        
        
    
end