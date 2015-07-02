function dbInsertMissingRecords
    
    %%% (1) what should have been moviePresId # 2751 is missing from tbl_movie_pres
    %%% it would be the same as #2749, except for a few details: (
    %%% presId, date/time, lng_present_no, and start/end ticks (which we can fix
    %%% afterwards). We put it at the end of the table, with presId 5057
    hnd = dbOpenExpDb;
    dateMaskOut = 'mm/dd/yyyy HH:MM:SS PM';
    dateMaskIn = 'mm/dd/yyyy HH:MM:SS';
    pres5057 = getFieldsFromDatabaseTable(hnd, 'MOVIE_PRES_ID', 'TBL_MOVIE_PRES', {'MOVIE_PRES_ID', 5057});
    if isempty(pres5057)
        pres5057 = getFieldsFromDatabaseTable(hnd, '*', 'TBL_MOVIE_PRES', {'MOVIE_PRES_ID', 2749});  % is more similar to pres#2 than pres#3
        pres5057{1} = 5057; % MOVIE_PRES_ID
        pres5057{5} = 4;    % LNG_PRESENT_NO

        dtm_num = datenum(pres5057{end}) + datenum(0,0,0,0,4,0); % 4 minutes after 2nd pres. (~2 minutes after 3rd pres)
        dtm_str = datestr(dtm_num, dateMaskIn);
        
        % remove date field, and add it later (cast into date format)
        pres5057 = pres5057(1:end-1);%{end}(1) = []; %skip date for now        
        pres5057 = pres5057(cellfun(@(x) ~isnan(x(1)), pres5057)); % remove ticks if are NaNs
        insertRecordIntoDatabaseTable(hnd, 'TBL_MOVIE_PRES', pres5057); 
        Did = 2639;
        dbMakeTbTe(hnd, Did);        
                        
        updateValueInDatabaseTable(hnd, dtm_str, 'DTM_CREATED', 'TBL_MOVIES', {'MOVIE_PRES_ID', 5057}, 'DATE');                
    end

    
    %%% (1) moviePresId # 34 missing from tbl_movie_pres
    %%% it should be the same as #35, except for a few details: (
    %%% presId, date/time, lng_present_no, and start/end ticks (which we can fix
    %%% afterwards). 
    %%% It is added to the table with presId #5058

    pres5058 = getFieldsFromDatabaseTable(hnd, 'MOVIE_PRES_ID', 'TBL_MOVIE_PRES', {'MOVIE_PRES_ID', 5058});
    if isempty(pres5058)
        pres5058 = getFieldsFromDatabaseTable(hnd, '*', 'TBL_MOVIE_PRES', {'MOVIE_PRES_ID', 35});        
        pres5058{1} = 5058; % MOVIE_PRES_ID
        pres5058{5} = 1;  % LNG_PRESENT_NO

        % remove date field, and add it later (cast into date format)
        dtm_num = datenum(pres5058{end}) - datenum(0,0,0,0,3,0); % 3 minutes before 2nd exp.
        dtm_str = datestr(dtm_num, dateMaskIn);

        pres5058 = pres5058(1:end-1); %skip date for now
        pres5058 = pres5058(cellfun(@(x) ~isnan(x(1)), pres5058)); % remove ticks if are NaNs
        insertRecordIntoDatabaseTable(hnd, 'TBL_MOVIE_PRES', pres5058);
        Did = 1622;
        dbMakeTbTe(hnd, Did); 
        
        updateValueInDatabaseTable(hnd, dtm_str, 'DTM_CREATED', 'TBL_MOVIE_PRES', {'MOVIE_PRES_ID', 5059});
        
        
    end
    
    
    updateValueInDatabaseTable(hnd, '9/25/2002 00:00:01', 'DTM_CREATED', 'TBL_GRATING_PRES', {'GRATING_PRES_ID', 166228}, 'DATE')
    
end

