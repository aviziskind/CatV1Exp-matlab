function toValues = dbLookup(toName,   fromName, fromValues, oneToOneFlag)
%     eg. Gid = lookupId('Gid',  'Did', 45); === "look up Gid where Did == 45;

    if strcmp(toName, fromName)
        toValues = fromValues;
        return;
    end
    
    if ~isOldXPS
        assert(length(fromValues) == 1)
        assert( any( strcmp(fromName, {'Did', 'Gid'}) ) )
        sd = siteDataFor(fromName, fromValues);
        toValues = sd.(toName);
        return
    end
    
    hnd = dbOpenExpDb;
    
 	DatafileTable = 'TBL_DATA_FILES';
    SynchroTable = 'TBL_SYNCHRO_GROUPS';
	GroupTable = 'TBL_GROUPS';
    MovieTable = 'TBL_MOVIE_PRES';
    
    if exist('oneToOneFlag', 'var') && ~isempty(oneToOneFlag) && (length(fromValues) > 1)
        toValues = zeros(size(fromValues));
        for i = 1:numel(fromValues)
            tv = dbLookup(toName,  fromName, fromValues(i));
            toValues(i) = tv(1);
        end
        return;
    end
    
    switch titleCase(fromName)
        case {'Did', 'Datafile_id'}
            Dids = fromValues;
        
        case {'Gid', 'Group_id'}
            Dids = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', GroupTable, {'GROUP_ID', fromValues});

        case {'Sid', 'Sync_id', 'Synchro_id'}
            Dids = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', SynchroTable, {'SYNCHRO_GROUP_ID', fromValues});

        case {'Mid', 'Movie_id'}
            Dids = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', MovieTable, {'MOVIE_ID', fromValues});            

        case {'Stimid', 'Stimulus_id'}
            Dids = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', DatafileTable, {'STIMULUS_TYPE_ID', fromValues});            
            
        otherwise
            error('Unknown id group type');            
    end

    if isempty(Dids)
        error(['No Did for ' fromName ' == ' num2str(fromValues) ]);
    end
%     if length(Did) > 1
%         warning('lookup:multipleDids', 'More than one datafile for specified ID tag');
%     end
    
    switch titleCase(toName)
        case {'Did', 'Datafile_id'}
            toValues = Dids;
        
        case {'Gid', 'Group_id'}
            toValues = getFieldsFromDatabaseTable(hnd, 'GROUP_ID', GroupTable, {'DATAFILE_ID', Dids});

        case {'Sid', 'Sync_id', 'Synchro_id'}
            toValues = getFieldsFromDatabaseTable(hnd, 'SYNCHRO_GROUP_ID', SynchroTable, {'DATAFILE_ID', Dids});

        case {'Mid', 'Movie_id'}
            toValues = getFieldsFromDatabaseTable(hnd, 'MOVIE_ID', MovieTable, {'DATAFILE_ID', Dids});
            
        case {'Stimid', 'Stimulus_id'}
            toValues = getFieldsFromDatabaseTable(hnd, 'STIMULUS_TYPE_ID', DatafileTable, {'DATAFILE_ID', Dids});            

        otherwise
            error('Unknown id group type');            
    end

    if isempty(toValues);
        warning('lookup:Empty', ['No ' toName ' for ' fromName '=' num2str(fromValues)]);
    end

end