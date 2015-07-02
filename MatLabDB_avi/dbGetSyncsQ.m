function syncs = dbGetSyncsQ(Did)
    persistent allSyncs;
    % a quick way to see how many syncs are in a given group/cell
    % The table 'allSyncs' is loaded (or created and saved, if it doesn't exist)
    % so that the number of syncs in Datafile i in allSyncs(i);
    redo = true;
    if isempty(allSyncs) || redo
        path = getName('MatlabDB_path');
        filename = [path 'dbAllSyncs.mat'];
        if exist(filename, 'file')
            S = load(filename);
            allSyncs = S.allSyncs;
        else
            hnd = dbOpenExpDb;
            allDatafileIds = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_DATA_FILES');

            maxDatafileId = max(allDatafileIds);            
            allSyncs = cell(maxDatafileId, 1);

            progressBar('init-', length(allDatafileIds));
            for Did_i = 1:length(allDatafileIds);
                progressBar(Did_i);
                Did = allDatafileIds(Did_i);            
                syncs = dbGetSyncs('Did', Did, 'tick', 1);
                allSyncs{Did} = syncs;
            end            
            save(filename, 'allSyncs');            
        end        
    end

    syncs = allSyncs{Did};


end