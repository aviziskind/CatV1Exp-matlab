function dbRemoveOrphanedGidsDids

    hnd = dbOpenExpDb;
    stimTypes = {'movie', 'noise', 'grating', 'mseq'};

    useMatFiles = true;

    removeOrphansFromDataFilesTable = true;
    removeOrphansFromGroupsTable = true;

    %gather unorphaned GroupIds, DatafileIds; These are the groups/data
    %files that are in the stim_pres tables.
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
    allGoodGids = [allStimGids{:}];
    allGoodDids = [allStimDids{:}];
    

    doEachStimulusSeparately = false;
    
    %%%% DATAFILE_IDS
    if removeOrphansFromDataFilesTable

        disp([' DATAFILE_IDs ']);                    
        dfDids = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_DATA_FILES');
        orphanDids = setdiff(dfDids, allGoodDids);

        disp(['Found ' num2str(length(allGoodDids)) ' Datafile Ids in the presentation tables.']);
        disp(['Found ' num2str(length(dfDids)) ' Datafile Ids in the TBL_DATA_FILES  table.']);
        if ~isempty(orphanDids)
            disp(['There are the thus ' num2str(length(orphanDids)) ' entries in the TBL_DATA_FILES that will be removed: ']);
            yn = input('Proceed? [Y] ', 's');
            if isempty(yn) || strncmpi(yn, 'Y', 1);
                progressBar('init-', length(orphanDids), 30);
                for i = 1:length(orphanDids)
                    progressBar;
                    removeRecordFromDatabaseTable(hnd, 'TBL_DATA_FILES', {'DATAFILE_ID', orphanDids(i)})
                end
                progressBar('done');
            end
        end
    end                


    %%%% GROUP_IDS
    
    % troubleSomeGids = [134 470 622 801 802]; % these don't fully remove - leave "#deleted" entries in DB
    
    if removeOrphansFromGroupsTable
        disp([' GROUP_IDs ']);    
        gpGids = getFieldsFromDatabaseTable(hnd, 'GROUP_ID', 'TBL_GROUPS');        
        orphanGids = setdiff(gpGids, allGoodGids);
        
        disp(['Found ' num2str(length(allGoodGids)) ' Group Ids in the presentation tables.']);
        disp(['Found ' num2str(length(gpGids)) ' Group Ids in TBL_GROUPS.']);
        if ~isempty(orphanGids)
            disp(['There are thus ' num2str(length(orphanGids)) ' in the TBL_GROUPS that will be removed: ']);
%             disp(orphanGids(:)');
            yn = input('Proceed? [Y] ', 's');            
            if isempty(yn) || strncmpi(yn, 'Y', 1);
                progressBar('init-', length(orphanGids), 30);
                for i = 1:length(orphanGids)
                    progressBar;
                    removeRecordFromDatabaseTable(hnd, 'TBL_GROUPS', {'GROUP_ID', orphanGids(i)})                    
                end
                progressBar('done');
            end
            fprintf('\nDone\n\n');
        end
    end

    if removeOrphansFromPresTable % are only 2 movie grp ids of these kinds
        for ti = 1;%1:length(stimTypes)
            stimType = stimTypes{ti};
            disp([' *** ' upper(stimTypes{ti}) ' ***']);
            stimTableName = ['TBL_' upper(stimType) '_PRES'];

            tic;
            Dids = unique ( getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', stimTableName) );
            toc;
            
            existingGids = dbLookup('Gid',  'Did', Dids);
            validDids = unique( dbLookup('Did',  'Gid', existingGids) );
            orphanedDids = setdiff(Dids, validDids);
%             extraDids = setdiff(validDids, Dids);
            
            disp(['Found ' num2str(length(Dids)) ' Datafile Ids for ' stimType ' stimuli in ' stimTableName ' table.']);
            if isempty(orphanedDids)
                disp('all are valid');
            else
                disp(['However, only ' num2str(length(validDids)) ' of these have Group Id entries.']);
                disp(['These are the ' num2str(length(orphanedDids)) ' entries that are orphaned, and will be removed: ']);
                disp(orphanedDids(:)');

                for i = 1:length(orphanedDids)
                    removeRecordFromDatabaseTable(hnd, stimTableName, {'DATAFILE_ID', orphanedDids(i)})
                end
            end
            3;
        end
    end    
    
    end


    
%         if removeOrphansFromDataFilesTable || removeOrphansFromGroupsTable
%         for ti = 2:length(stimTypes)
%             stimType = stimTypes{ti};
%             disp([' *** ' upper(stimTypes{ti}) ' ***']);
%             stimTableName = ['TBL_' upper(stimType) '_PRES'];
% 
%             if useMatFiles
%                 stimDids = allStimDids{ti};
%                 stimGids = allStimGids{ti};
%             else
%                 stimDids = unique(getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', stimTableName)); %#ok<*UNRCH>
%                 stimGids = unique( dbLookup('Gid',  'Did', stimDids) );
%             end
% 
%             %%%% DATAFILE_IDS
%             if removeOrphansFromDataFilesTable
%                 
%                 disp([' DATAFILE_IDs ']);
%                 [dfDids] = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_DATA_FILES', {'STIMULUS_TYPE_ID', stimTypeIds{ti}});
%                 orphanDids = setdiff(dfDids, stimDids);
% 
%                 disp(['Found ' num2str(length(dfDids)) ' Datafile Ids for ' stimType ' stimuli in TBL_DATA_FILES table.']);
%                 disp(['Found ' num2str(length(allStimDids{ti})) ' Datafile Ids for ' stimType ' stimuli in ' stimTableName ' table.']);
%                 disp(['These are the extra ' num2str(length(orphanDids)) ' entries in the TBL_DATA_FILES that will be removed: ']);
%                 if isempty(intersect(orphanDids, allGoodDids))
%                     disp('none are accounted for for other stimulus types')
%                 else
%                     keyboard;
%                 end
%                 disp(orphanDids(:)');
% 
%                 for i = 1:length(orphanDids)
%                     removeRecordFromDatabaseTable(hnd, 'TBL_DATA_FILES', {'DATAFILE_ID', orphanDids(i)})
%                 end
% 
%                 fprintf('\nDone\n\n');
%             end
% 
%             %%%% GROUP_IDS
%             if removeOrphansFromGroupsTable
%                 gpGids = dbLookup('Gid',  'Did', dfDids);
%                 orphanGids = setdiff(gpGids, stimGids);
% 
%                 disp([' GROUP_IDs ']);
%                 disp(['Found ' num2str(length(gpGids)) ' Group Ids for ' stimType ' stimuli in TBL_GROUPS table.']);
%                 disp(['Found ' num2str(length(allStimGids{ti})) ' Group Ids for ' stimType ' stimuli in ' stimTableName ' table.']);
%                 if ~isempty(orphanGids)
%                     disp(['These are the extra ' num2str(length(orphanGids)) ' in the TBL_GROUPS that will be removed: ']);
%                     disp(orphanGids(:)');
%                     for i = 1:length(orphanGids)
%                         removeRecordFromDatabaseTable(hnd, 'TBL_GROUPS', {'GROUP_ID', orphanGids(i)})
%                     end
%                 end
%                 fprintf('\nDone\n\n');
%             end
% 
%         end
%     end
