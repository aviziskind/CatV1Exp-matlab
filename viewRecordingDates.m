function viewRecordingDates(idType, idVal)
 3;
    idType = 'PenId';
    idVal = 173;
%%
%     idType = 'LocId';
%     idVal = 507;

    dateMask = 'mm/dd/yyyy HH:MM:SS PM';
    
    loc_data = dbGetLocationData;
    loc_data = loc_data(~cellfun(@isempty, {loc_data.Gid}));
    
    n = arrayfun(@(s) length( s.PenetrId), loc_data);
    idx2 = find(n==2);
    Gids_double = [loc_data(idx2).Gid];
    %%
    allAnimalId = [loc_data.AnimalId];
    allGid = [loc_data.Gid];
    allDid = [loc_data.Did];
    allPenId1   = arrayfun(@(s) s.PenetrId(1), loc_data);
    allPenId2   = arrayfun(@(s) s.PenetrId(end), loc_data);
    allLocId1   = arrayfun(@(s) s.LocId(1), loc_data);
    allLocId2   = arrayfun(@(s) s.LocId(end), loc_data);

    allPenDates1 = arrayfun(@(s) s.Pen_date(1), loc_data);
    allPenDates2 = arrayfun(@(s) s.Pen_date(end), loc_data);
    allLocDates1 = arrayfun(@(s) s.Loc_date(1), loc_data);
    allLocDates2 = arrayfun(@(s) s.Loc_date(end), loc_data);
    allGrpDates = arrayfun(@(s) s.Grp_date(1), loc_data);

    %%
    date2num = @(date_str) datenum(date_str, dateMask);
    allPenDateNums1 = cellfun(date2num, allPenDates1);
    allPenDateNums2 = cellfun(date2num, allPenDates2);
    allLocDateNums1 = cellfun(date2num, allLocDates1);
    allLocDateNums2 = cellfun(date2num, allLocDates2);
    allGrpDateNums  = cellfun(date2num, allGrpDates);
    
    
    
    3;
    %%
    switch idType
        case 'LocId', idx_use = find( (allLocId1 == idVal) | (allLocId2 == idVal) );
        case 'PenId', idx_use = find( (allPenId1 == idVal) | (allPenId2 == idVal) );    
    end
        
    %%
    iter = 1;
    V = {};
    V_srt = {};
    for i = 1:length(idx_use)  
        j = idx_use(i);
        if length(loc_data(j).PenetrId) == 1
            V(iter,:) = {i, ' ', allPenId1(j), allPenDates1{j}, allLocId1(j), allLocDates1{j}, allGid(j), allGrpDates{j} };
            V_srt(iter,:) = {j, 0, allPenId1(j), allPenDateNums1{j}, allLocId1(j), allLocDateNums1{j}, allGid(j), allGrpDateNums{j} };
            iter = iter+1;
            
        else
            V(iter,:) = {i, 'a', allPenId1(j), allPenDates1{j}, allLocId1(j), allLocDates1{j}, allGid(j), allGrpDates{j} };
            V_srt(iter,:) = {j, 1, allPenId1(j), allPenDateNums1(j), allLocId1(j), allLocDateNums1(j), allGid(j), allGrpDateNums(j) };
            iter = iter+1;
            V(iter,:) = {i, 'b', allPenId2(j), allPenDates2{j}, allLocId2(j), allLocDates2{j}, allGid(j), allGrpDates{j} };
            V_srt(iter,:) = {j, 2, allPenId2(j), allPenDateNums2(j), allLocId2(j), allLocDateNums2(j), allGid(j), allGrpDateNums(j) };
            iter = iter+1;
        end        
    end
       
    %%
    [~, idx] = sortrows(V_srt, [7]);
    V = V(idx,:);
    for i = 1:size(V,1)
        fprintf('(%d %s) Pen %4d. (%s) |  Loc %4d. (%s)  | Gid %4d (%s)  |\n', ...
        V{i,:});
3;        
    end
    
    3;
    
    
end
    

%     CatId: 36
%     PenId: 173
%     LocId: 507    

    
    
    
%     
% 
%     get
% 
% 
%             fieldnames = {'TBL_DATA_FILES.DATAFILE_ID', 'TBL_GROUPS.GROUP_ID', 'TBL_ANIMALS.TXT_LAB_NAME', ...
%                 'TBL_ANIMALS.ANIMAL_ID', 'TBL_PENETRATIONS.PENETRATION_ID', 'TBL_LOCATIONS.LOCATION_ID', ...
%                 'TBL_LOCATIONS.TXT_LOCATION_NAME', ...
%                 'TBL_ELECTRODE_TYPES.ELECTRODE_TYPE_ID', 'TBL_PENETRATIONS.BRAIN_STRUCT_TYPE_ID', ...
%                 'TBL_GROUPS.DTM_CREATED', 'TBL_DATA_FILES.DTM_CREATED', 'TBL_PENETRATIONS.DTM_CREATED', 'TBL_LOCATIONS.DTM_CREATED', 'TBL_AP_ML_ZERO.DTM_CREATED'};
% 
%             T1 = {'TBL_AP_ML_ZERO', 'ELECTRODE_ID', 'TBL_ELECTRODES'};
%             T2 = {T1, 'TBL_AP_ML_ZERO', 'AP_ML_ZERO_ID', 'TBL_PENETRATIONS'};               
%             T3 = {'TBL_ANIMALS', 'ANIMAL_ID', T2, 'TBL_ELECTRODES'};
%             T4 = {T3, 'TBL_PENETRATIONS', 'PENETRATION_ID', 'TBL_LOCATIONS'};
%             T5 = {'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_LOCS_FILES_LINKS'};
%             T6 = {'TBL_LOC_DEPTHS', 'LOC_DEPTH_ID',  T5, 'TBL_LOCS_FILES_LINKS'};
%             T7 = {T4, 'TBL_LOCATIONS', 'LOCATION_ID', T6, 'TBL_LOC_DEPTHS'};
%             T8 = {T7, 'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_GROUPS'};
%             T9 = {'TBL_ELECTRODE_TYPES', 'ELECTRODE_TYPE_ID', T8, 'TBL_ELECTRODES'};
%             joinedTables = T9;
%             
%                         
% %             criterea = {'TBL_ELECTRODES.ELECTRODE_TYPE_ID', 2};
%             criterea = {};
%             
%             [Did_tbl, Gid_tbl, CatName_tbl, AnimalID_tbl, PenetrID_tbl, LocID_tbl, LocName_tbl, elect_tbl, brainStruct_tbl, ...
%                 date_grp, date_df, date_pen, date_loc, date_apml] = ...
%                 getFieldsFromDatabaseTable(hnd, fieldnames, joinedTables, criterea);    
% 
% 




%{
    for i = 1:length(idx_use)  
        j = idx_use(i);
        if length(loc_data(j).PenetrId) == 1
            fprintf('(%d) Pen %4d. (%s) |  Loc %4d. (%s)  | Gid %4d (%s)  |\n', ...
                j, allPenId1(j), allPenDates1{j}, allLocId1(j), allLocDates1{j}, allGid(j), allGrpDates{j} );
            
        else            
            fprintf('(%da)Pen %4d. (%s) |  Loc %4d. (%s)  | Gid %4d (%s)  |\n', ...
                j, allPenId1(j), allPenDates1{j}, allLocId1(j), allLocDates1{j}, allGid(j), allGrpDates{j});
            fprintf('(%db)Pen %4d. (%s) |  Loc %4d. (%s)  | Gid %4d (%s)  |\n', ...
                j, allPenId2(j), allPenDates2{j}, allLocId2(j), allLocDates2{j}, allGid(j), allGrpDates{j});
        end
        
        
    end
%}









