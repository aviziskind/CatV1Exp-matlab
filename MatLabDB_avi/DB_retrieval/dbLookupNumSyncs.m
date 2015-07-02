function varargout = dbLookupNumSyncs(Did)
    persistent syncCounts;
    % a quick way to see how many syncs are in a given group/cell
    % The table 'syncCounts' is loaded (or created and saved, if it doesn't exist)
    % so that the number of syncs in Datafile i in syncCounts(i);
    
    if isempty(syncCounts) 
        path = getName('MatlabDB_path');
        filename = [path 'dbNumSyncsInDatafile.mat'];
        if exist(filename, 'file')
            S = load(filename);
            syncCounts = S.syncCounts;
        else
            hnd = dbOpenExpDb;
            allDatafileIds = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_DATA_FILES');

            maxDatafileId = max(allDatafileIds);            
            syncCounts = -ones(maxDatafileId, 1);

            progressBar('init-', length(allDatafileIds));
            for Did_i = 1:length(allDatafileIds);
                progressBar(Did_i);
                Did = allDatafileIds(Did_i);            
                nSyncs = length(dbGetSyncs('Did', Did, 'tick', 1));                
                syncCounts(Did) = nSyncs;
            end            
            save(filename, 'syncCounts');            
        end        
    end

    switch nargin
        case 0
            varargout = {syncCounts};
        case 1
            nSyncs = syncCounts(Did);
            if (nSyncs == -1)
                error(['Datafile id # ' num2str(Did) ' doesn''t exist']);
            end
            varargout = {nSyncs};
    end
            

end