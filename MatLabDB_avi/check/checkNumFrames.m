function checkNumFrames
% compare nSustained, nDisplayed, nSyncs between start & end

    function fprintf2(s)
        fprintf(fid1, s);
        fprintf(fid2, s);
    end
    function xout = iff(x)
        if ~isempty(x),
            xout = x;
        else
            xout = 0;
        end
    end


    hnd = dbOpenExpDb;

    allGids = getFieldsFromDatabaseTable(hnd, 'GROUP_ID', 'TBL_GROUPS', [], 'GROUP_ID');
    allDids = dbLookup('Did',  'Gid', allGids);
    
    fieldsToGet = {'LNG_N_PRE_BLANK_FRM', 'LNG_N_POST_BLANK_FRM', 'LNG_N_SUSTAINED_FRM', 'LNG_N_DISPLAYED_FRAMES', 'LNG_N_MISSED_FRM', 'LNG_START_TICK', 'LNG_END_TICK'};
	
    
    for Gid_i = 1:length(allGids)
        Gid = allGids(Gid_i);
        Did = allDids(Gid_i);
        stimType = getStimulusTypeForDid(Did);
        tableName = getDatabaseTableForDid(Did);
        presIdField = [upper(stimType) '_PRES_ID'];

        presIds = getFieldsFromDatabaseTable(hnd, presIdField, tableName, {'DATAFILE_ID', Did});
        
        [nPreBlank, nPostBlank, nSustainedFrm, nDisplayedFrm, nMissedFrm, startTicks, endTick] = ...
            getFieldsFromDatabaseTable(hnd, stimTableName, fieldsToGet, {presIdField, presIds});
        
        Sid  = dbLookup('Sid',  'Did', Did);
        if ~isempty(Sid), 
            syncs = edbGetSyncs(hnd, Sid);
            if isempty(syncs)
                continue;
            end
        else
            continue;
        end;

        fprintf2(['Gid = ' pad(num2str(Gid), 4),  ': ']);
        for p_i = 1:length(presIds)
            st = startTicks(p_i);
            [nSyncStartId nSyncEndId] = elements( binarySearch(syncs, [startTicks(p_i) endTicks(p_i)]);
            if ((syncs(nSyncStartId) == startTicks(p_i)) && syncs(nSyncEndId) == endTicks(p_i))
                nFrames_syncs = nSyncEndId - nSyncStartId +1;
                
                matchDisp = (nFrames_syncs == nDisplayedFrm(p_i));
                matchSust = (nFrames_syncs == nSustainedFrm(p_i));
                fprintf2(' presId = ' num2str(presIds(p_i) ' N_Displayed=' num2str(nDisplayedFrm(p_i)) '. N_Sustained=' num2str(nSustainedFrm(p_i)) '. N_frames=' num2str(nFrames_syncs)
                
                
            end
            
            
            
        end
        
        if ~isempty(frameLen_ms_Sync) && ~isempty(frameLen_ms_Defn)
            if (frameLen_ms_Sync == frameLen_ms_Defn)
                fprintf2(['Good: Both are ' num2str(frameLen_ms_Defn) '.\n' ]);
            else
                discrep = abs(frameLen_ms_Defn - frameLen_ms_Sync) / (frameLen_ms_Sync+eps) * 100;
                fprintf2([' -- ERR: ( ' pad( num2str(discrep, '%3.2f'),5) ' %% )     From syncs: ' pad( num2str(frameLen_ms_Sync, '%4.2f'), 7) '.  From Defn: ' pad(num2str(frameLen_ms_Defn, '%10.3f'),7) '\n']);
            end

            allDiscreps(Gid_i) = discrep;
        else
            fprintf2([' -- ERR: ( ' pad( num2str(0, '%3.2f'),5) ' %% )     From syncs: ' pad( num2str(iff(frameLen_ms_Sync), '%4.2f'), 7) '.  From Defn: ' pad(num2str(iff(frameLen_ms_Defn), '%10.3f'),7) '\n']);
            allDiscreps(Gid_i) = NaN;

        end
        
        
    end
    
    
    
    
end