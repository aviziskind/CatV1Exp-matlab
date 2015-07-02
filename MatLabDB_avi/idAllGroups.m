% function idAllGroups
    hnd = dbOpenExpDb;
%     allGids  = getFieldsFromDatabaseTable(hnd, 'GROUP_ID', 'TBL_GROUPS');
%     allGids2 = 1:5339;
    
    figure(1); clf; hold on;
    labels = {};
    for stim_id = 1:10
    	stim_txt = getFieldsFromDatabaseTable(hnd, 'TXT_STIMULUS_TYPE', 'TBL_STIMULUS_TYPES', {'STIMULUS_TYPE_ID', stim_id});
        disp(stim_txt{1});
        Dids = unique(getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_DATA_FILES', {'STIMULUS_TYPE_ID', stim_id}));

        switch stim_txt{1}
            case {'Single Grating', 'Flashed Grating Batch', 'Orientation Batch', 'Spatial Frequency Batch', 'Temporal Frequency Batch', 'Free Grating Batch'}
                col = 'b';
            case 'Noise'
                col = 'r';
            case 'M-sequence'
                col = 'm';
            case {'Movie', 'Movie Batch'}
                col = 'g';
        end        
        
        if ~isempty(Dids)
            G = getFieldsFromDatabaseTable(hnd, 'GROUP_ID', 'TBL_GROUPS', {'DATAFILE_ID', Dids});
            GidsThisStim{stim_id} = G; %#ok<AGROW>
            allGids = setdiff(allGids, G);
            plot(G, ones(size(G)), [col, marker(stim_id)]);
            labels = {labels{:}, stim_txt};
        end
                
    end

    for stim_id = 1:10
    	stim_txt = getFieldsFromDatabaseTable(hnd, 'TXT_STIMULUS_TYPE', 'TBL_STIMULUS_TYPES', {'STIMULUS_TYPE_ID', stim_id});        
        switch stim_txt{1}
            case {'Single Grating'} %, 'Flashed Grating Batch', 'Orientation Batch', 'Spatial Frequency Batch', 'Temporal Frequency Batch', 'Free Grating Batch'}
                dataFileTable = 'TBL_GRATING_PRES';
                col = 'b';
            case 'Noise'
                dataFileTable = 'TBL_NOISE_PRES';
                col = 'r';
            case 'M-sequence'
                dataFileTable = 'TBL_MSEQ_PRES';
                col = 'm';
            case {'Movie'}%, 'Movie Batch'}
                dataFileTable = 'TBL_MOVIE_PRES';
                col = 'g';
            otherwise
                continue;
        end
        tic;
        datafilesFromDFtable = unique(getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', dataFileTable));
        toc;
        G2 = getFieldsFromDatabaseTable(hnd, 'GROUP_ID', 'TBL_GROUPS', {'DATAFILE_ID', datafilesFromDFtable});
        plot(G2,2*ones(size(G2)), [col, marker(stim_id)]);
        allGids2 = setdiff(allGids2, G2);
    end


    plot(allGids, zeros(size(allGids)), 'k.');
    plot(allGids2, .2+zeros(size(allGids2)), 'b.');
    ylim([-1 5]);
    legend([labels{:}]);
    
    
    
    load movieCellRecords
    inds = findInStructArray(movieCellRecords, 'movieCells', [], @isnotempty);
    [GidsC{1:length(inds)}] = movieCellRecords(inds).Gid;
    GidsWorked = cell2mat(GidsC);
    plot(GidsWorked, 2.2 + zeros(size(GidsWorked)), 'g.');

    
    
    
    
    
    
    
