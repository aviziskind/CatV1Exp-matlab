%%% ALL ARE CORRECT
function correctStimTypeId

    hnd = dbOpenExpDb;
    stimTypes = {'movie', 'noise', 'grating', 'mseq'};
    stimTypeIds = {[8, 10], 4, [1,2,3,6,7,9], 5};

    useMatFiles = true;
    
    %gather all correct GroupIds, DatafileIds;
    allStimGids = cell(1,4);
    allStimDids = cell(1,4);
    for ti = 1:length(stimTypes)
        stimType = stimTypes{ti};
        S = load(['cellsGroups_' stimType '_all_DB']);
        groupData = S.([stimType 'Groups_all']);
        stimDids = unique( [groupData(:).Did] );
        stimGids = unique( [groupData(:).Gid] );
        
        allStimDids{ti} = stimDids;
        allStimGids{ti} = stimGids;
    end
%     allGoodGids = [allStimGids{:}];
%     allGoodDids = [allStimDids{:}];

    % go through each datafile table, and correct incorrect stimPresId entries 
    % for movies - just put as 'movie' (8), not 'movie batch'.
    % for gratings ....
    
    for ti = 1:length(stimTypes)
        stimType = stimTypes{ti};
        disp(['Correcting stimPresId for ' stimTypes{ti} ' DatafileIds'])
        
        % are there any datafiles of this stimulus type in TBL_DATA_FILES
        % with another stimulus_type_id?
        
        
%         [dids, stimTypeIds] = getFieldsFromDatabaseTable(hnd, {'DATAFILE_ID', 'STIMULUS_TYPE_ID'}, 'TBL_DATA_FILES', {'STIMULUS_TYPE_ID', stimTypeIds{ti}});
        [dids, stimT_Ids] = getFieldsFromDatabaseTable(hnd, {'DATAFILE_ID', 'STIMULUS_TYPE_ID'}, 'TBL_DATA_FILES', {'DATAFILE_ID', allStimDids{ti}});
        %%%% DATAFILE_IDS 
        varBreakdown(stimT_Ids);
        
        DidIndsToCorrect = find( ~any (crossOp(stimT_Ids, @eq, stimTypeIds{ti} ) , 2) );
        [dids(DidIndsToCorrect) stimT_Ids(DidIndsToCorrect)]
        
            
        % now go to Datafile Table, and remove extra entries 
        
%         end
        3;
    end
    

end