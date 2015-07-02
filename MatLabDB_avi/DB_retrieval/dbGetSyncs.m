function syncs = dbGetSyncs(idType, id, outputType, returnEmpty_flag, dontAmend_flag)

    returnEmptyIfError = exist('returnEmpty_flag', 'var') && ~isempty(returnEmpty_flag);
    dontAmendSyncs = exist('dontAmend_flag', 'var') && ~isempty(dontAmend_flag);

    if isXPS15
        assert(strcmp(idType, 'Gid'))
        syncs = getSyncs(id, outputType);
        return;
    end
    
    
    Sid = dbLookup('Sid', idType, id);
    if isempty(Sid)
        if returnEmptyIfError
            syncs = [];
            return;
        else
            error('db:noSid', 'No Synchro-group Id');
        end
    end

    hnd = dbOpenExpDb;
    syncs_field = getFieldsFromDatabaseTable(hnd, 'MEM_SYNCHROTIMES_MTX', 'TBL_SYNCHRO_GROUPS', {'SYNCHRO_GROUP_ID', Sid});
    
    if isempty(syncs_field) || (~iscell(syncs_field) && isnan(syncs_field))
        if returnEmptyIfError
            syncs = [];
            return;
        else
            error('db:noSyncs', 'No sync pulses');
        end
    end
    
    syncs_str = [ syncs_field{:} ];
    if any( ~ischar(syncs_str) )
        if returnEmptyIfError
            syncs = [];
            return;            
        else
            error('Invalid syncs (~char)');
        end
    end

    syncTimes_ticks = edbStr2Mtx(syncs_str);
    if any( isnan(syncTimes_ticks) )
        if returnEmptyIfError
            syncs = [];
            return;            
        else
            error('Invalid syncs (NaN)');
        end
    end
        
    syncTimes_ticks = syncTimes_ticks(:, 1);  % second column is not used

    if ~dontAmendSyncs
        syncTimes_ticks = amendSyncs(Sid, syncTimes_ticks);        
    end
    
    if exist('outputType', 'var') && ~isempty(outputType)
        if strfind(outputType, 'frame'), error('Recursive call'); end
        Did = dbLookup('Did',  idType, id);
        syncs = dbConvertTimeMeasures(Did, syncTimes_ticks, 'tick', outputType);
    else
        syncs = syncTimes_ticks;
    end
        
end

function syncTimes_ticks = amendSyncs(Sid, syncTimes_ticks)
    S_amend = dbGetSyncAmendments('Sid', Sid);

    if ~isempty(fieldnames(S_amend))
        if isfield(S_amend, 'syncsToRemove')
            syncTimes_ticks = setdiff(syncTimes_ticks, S_amend.syncsToRemove(:) );           
            syncTimes_ticks = syncTimes_ticks(:);
        end
        if isfield(S_amend, 'syncsToAdd')
             syncTimes_ticks = unique([syncTimes_ticks(:); S_amend.syncsToAdd(:)]);
        end
        
        if isfield(S_amend, 'cutInHalf')
            L = length(syncTimes_ticks);
            syncTimes_ticks = syncTimes_ticks(1:L/2);
        end                
    end
        
end