function generateGroupData_fromDBgroups(stimTypesToDo, groupingTypesToDo)
    
    % takes the groupData files generated from generateCellGroupData_DB
    % and changes the 'cellId' and 'nSpikes' fields appropriately
    % to match the results from the automatic sorting.
    mergeDiscardedAndMU = 1;
    curMatchDB(0);
    
    if nargin < 1
        stimTypesToDo = {'movie_fg', 'grating_dOr', 'grating_dSf'};
%         stimTypesToDo = {'movie_fg'};
%         stimTypesToDo = {'grating_dSf'};

    end            
    
    if nargin < 2
        groupingTypesToDo = {'cells'};
%         groupingTypesToDo = {'clusters', 'clustersPruned', 'cells'};
    end
    
    for grpType_i = 1:length(groupingTypesToDo)
        groupingType = groupingTypesToDo{grpType_i};
        
        
        for stim_i = 1:length(stimTypesToDo);
            gratType = stimTypesToDo{stim_i};

            [destPath, destFile] = getFileName('Groups', gratType);
            fprintf('Doing %s %s groups (%s)\n', gratType, groupingType, destFile);
            
            S = load(['cellsGroups_' gratType '_DB.mat']);
            fld_nm = fieldnames(S); fld_nm = fld_nm{1};
            allGroups = S.(fld_nm);
            
            
            allGids = [allGroups.Gid];
            progressBar('init-', length(allGids));
            for i = 1:length(allGids)
                progressBar;
                [cellIds, nSpikes] = getCellIds(allGids(i), groupingType, mergeDiscardedAndMU);
                
                allGroups(i).cellIds = cellIds;
                allGroups(i).nSpikes = nSpikes;
                
                % add decor fields:
                [~, stimRecur_frm] = getStimulusRecurrenceTime(allGids(i));
                [~, stimDecor_frm] = getStimulusDecorrelationTime(allGids(i));
                stimOrdered = isequalToPrecision(stimRecur_frm, stimDecor_frm, 1e-3);
                allGroups(i).stimOrdered = stimOrdered;
                
                if strcmp(groupingType, 'cells')
                    cellClusters = getClusterCellAssignment(allGids(i));
                    allGroups(i).cellClusters = cellClusters;
                end
            end
            progressBar('done');
            
            S.(fld_nm) = allGroups;
            S_old = load([destPath destFile]);
            save([destPath destFile], '-struct', 'S');
                        
        end
        
    end
        
                                

    
end

function [cellIds, nSpikes] = getCellIds(Gid, grpType, mergeDiscardedAndMU)        
    S_sorting = load(getFileName(grpType, Gid), 'uClustIds','clustCounts');
    uClustIds = S_sorting.uClustIds;
    clustNSpikes = S_sorting.clustCounts;
    if mergeDiscardedAndMU && any(uClustIds == -1)
        idxDiscarded = find(uClustIds == -1, 1);
        clustNSpikes(find(uClustIds == 0, 1)) = sum(clustNSpikes(uClustIds <= 0));
        uClustIds(idxDiscarded) = [];            
        clustNSpikes(idxDiscarded) = [];
        assert(sum(clustNSpikes) == sum(S_sorting.clustCounts));
    end                    
    cellIds = uClustIds;
    assert(isequal(cellIds, cellIds(1):cellIds(end)))
    nSpikes = clustNSpikes;
end


%         if all( dbDoesFieldExist(hnd, frameNumberFields, tableName ) );
%             
%             nSustainedFrames = stimTable.LNG_N_SUSTAINED_FRM(DidIdxs);
%             nDisplayedFrames = stimTable.LNG_N_DISPLAYED_FRM(DidIdxs);
%             preBlankFrames   = stimTable.LNG_N_PRE_BLANK_FRM(DidIdxs);
%             postBlankFrames  = stimTable.LNG_N_POST_BLANK_FRM(DidIdxs);
%                         
%             presOK = [(nMissedFrames == 0) & (nDisplayedFrames == nSustainedFrames + preBlankFrames + postBlankFrames) ]';
%         else
%             presOK = (nMissedFrames == 0)';
%         end



%         % compile flashed grating movies, with extra fg-specific data:
%         fprintf('Compiling flashed grating groups... '); tic;
%         fg_inds = findInStructArray(movieGroups_all, 'stimType', [], @(s) strncmp(s, 'Movie:Flashed_Grating', 19));
%         for i = 1:length(fg_inds)
%             s = [fnames, struct2cell(movieGroups_all(fg_inds(i)))]';
%             indMovieFileNames = strcmp(fnames, 'movieFiles');
%             indGid = strcmp(fnames, 'Gid');
%             Gid = s{2, indGid};
%             movieFileNames = s(2, indMovieFileNames);
%             s(:,indMovieFileNames) = [];
% %             s = s';
%             [allOri_deg, allSp_pix, allPhase_deg, uori, usp, uph] = getOriSpPhaseForEachStimFrame(Gid, 'Movie', 'Planned');
%             movieGroups_fg_all(i) = struct(s{:}, 'movieFiles', movieFileNames, 'ori_deg', uori, 'spPeriod_pix', usp, 'spPh_deg', uph, 'tempPeriod_sec', inf, 'tempPh_deg', 0);            
%         end
%         fprintf('done \n'); toc;
%         indsWithSiteOK = findInStructArray(movieGroups_fg_all, 'siteOK', [], @(s) strcmp(s, 'ok') );                
%         movieGroups_fg = movieGroups_fg_all(indsWithSiteOK);


%         vname_all = ['grating' GroupingType 'Groups_all'];
%         vname     = ['grating' GroupingType 'Groups'];
%         vname_fOrSph = ['grating' GroupingType 'Groups_fOrSph'];
%         vname_dSng     = ['grating' GroupingType 'Groups_dSng'];
%         vname_dOr     = ['grating' GroupingType 'Groups_dOr'];
%         vname_Sf     = ['grating' GroupingType 'Groups_dSf'];
%         


%         movieGids = [ 1194 1195 1201 1207 1215 1432 1433 1441 1442 1444 1474 1491 1493 1495 1499 1515 1521 1531 1534 1535 ...
%                         1536 1662 1709 1790 1861 1983 2132 2138 2139 2158 2159 2162 2165 2166 2170 2178 2238 2250 2460 2601 ...
%                         2625 2631 2655 2657 2665 2677 2715 2723 2737 2747 2775 2791 2805 2855 2887 2921 2957 2981 3035 3047 ...
%                         3051 3059 3061 3073 3077 3091 3093 3125 3151 3215 3373 3397 3415 3467 3637 3709 3719 3729 3731 3775 ...
%                         3841 3997 4023 4051 4117 4127 4312 4344 4362 4768 5100 5146 ];



% allSortFet = {{'Neg'}, {'Neg', 'Pos', 'Ptp'}, {'Neg', 'GLFs', 'GLFc'}, {'Neg', 'GLFs'} };
% nSets = length(allSortFet)
% for i = 1:nSets
%     curFets = allSortFet{i};    
%     curSortingFeatures(curFets{:});    
%     curFet_str = cellstr2csslist(curFets, '_');
%     
%     S = load(['clustersGroups_' curFet_str '_movie_fg']);
%     NClusts(:,i) = arrayfun(@(s) length(s.cellIds), S.movieGroups_fg);        
%     NClusts_C{i} = arrayfun(@(s) length(s.cellIds), S.movieGroups_fg);        
% end