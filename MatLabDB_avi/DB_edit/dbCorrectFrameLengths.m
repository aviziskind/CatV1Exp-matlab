function dbCorrectFrameLengths

    quickCorrect = true;  % set to true to just skip to the ones we know need to be corrected. false to recheck everything.

    function frameLen_ms = getFramesLengthsFromSyncs(Did)
        syncTimes_ms = dbGetSyncs('Did', Did, 'ms', 1);
        if ~isempty(syncTimes_ms)
            frameLen_ms = median(diff(syncTimes_ms)); % remove inter-presentation intervals (outliers)            
        end
    end

    promptForEachEdit = true;

    hnd = dbOpenExpDb;
    xifne = @(x) iff(isempty(x), x);
    pad = @(s, N) [repmat(' ', 1, N-length(s)), s];

    stimTypes = {'movie', 'noise', 'grating', 'mseq'};
    if quickCorrect
        typesToDo = 2;
    else
        typesToDo = [1 2 3 4];
    end

    for ti = typesToDo
        stimType = stimTypes{ti};
        fprintf(['DOING ' upper(stimType) ' GROUPS\n-----------------------\n']);

        if quickCorrect
            allStimGids = [136, 786];  % only these 2 need to be corrected.
        else
            S = load(['cellsGroups_' stimType '_all_DB']);
            groupData = S.([stimType 'Groups_all']);
            eval([stimType 'Groups = groupData;']);
            allStimGids = [groupData.Gid];
        end

        tableName = getDatabaseTableForDid(dbLookup('Did',  'Gid', allStimGids(1) ));

        allFrameLen_ms_Sync = zeros(size(allStimGids));
        allFrameLen_ms_DB = zeros(size(allStimGids));
        allDiscreps = zeros(size(allStimGids));

        for Gid_i = 1:length(allStimGids);
            Gid = allStimGids(Gid_i);
            Did = dbLookup('Did',  'Gid', Gid);

            fprintf(['(' outOf(Gid_i, length(allStimGids)) ') Gid = ' pad(num2str(Gid), 4),  '. Did = ' num2str(Did) ' : ']);

            % Movie/Noise/Grating IDs,  begin/end ticks of presentations
            frameLen_ms_Sync = getFramesLengthsFromSyncs(Did);
            frameLen_ms_DB   = dbConvertTimeMeasures(Did, 1, 'frame', 'ms');

            allFrameLen_ms_Sync(Gid_i) = frameLen_ms_Sync;
            allFrameLen_ms_DB(Gid_i)   = frameLen_ms_DB;

            if ~isempty(frameLen_ms_Sync) && ~isempty(frameLen_ms_DB)
                discrepPercent = abs(frameLen_ms_DB - frameLen_ms_Sync) / (frameLen_ms_Sync+eps) * 100;
                if (discrepPercent < 2) 
                    fprintf(['Good: Both are ~ ' num2str(frameLen_ms_DB) '.\n' ]);
                else
                    
                    fprintf([' -- discrepancy: ( ' pad( num2str(discrepPercent, '%3.2f'),5) ' %% )     From syncs: ' pad( num2str(frameLen_ms_Sync, '%4.2f'), 7) '.  From DB: ' pad(num2str(frameLen_ms_DB, '%10.3f'),7) '\n']);
                    fpuSync = round(frameLen_ms_Sync/(1000/120));
                    fpuDB   = round(frameLen_ms_DB/(1000/120));
                    if (fpuSync ~= fpuDB)
                        fprintf('    This is likely due to an incorrect fpu value in the DB: %d instead of %d.', fpuDB, fpuSync);
                        if promptForEachEdit
                            yn = input(' Make correction? [Y]', 's');
                            if isempty(yn)
                                yn = 'Y';
                            end
                        else
                            yn = 'Y';
                        end
                        if strcmpi(yn, 'Y')
                            updateValueInDatabaseTable(hnd, fpuSync, 'LNG_FRAMES_PER_UPDATE', tableName, {'DATAFILE_ID', Did});
                            fprintf('Corrected\n');
                        end
                        
                    end
                    
                end

                allDiscreps(Gid_i) = discrepPercent;
            else
                fprintf([' -- ERR: ( ' pad( num2str(0, '%3.2f'),5) ' %% )     From syncs: ' pad( num2str(xifne(frameLen_ms_Sync), '%4.2f'), 7) '.  From DB: ' pad(num2str(xifne(frameLen_ms_DB), '%10.3f'),7) '\n']);
                allDiscreps(Gid_i) = NaN;

            end
    
        end

    end
    
end