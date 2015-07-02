function discLoc_pix = getOpticDiskLocation(idType, idVal, leftRight)
%     S_amend = dbGetOpticDiscLocations(idType, idVal)        
%     S_amend = dbGetOpticDiscLocations(Gid)

    persistent allOpticDiscLocations allLocIds
    
    opticDiscLocations_file = [CatV1Path 'MatLabDB_avi' filesep 'allOpticDiscLocations.mat'];
    
    redo = 1;

    if nargin == 1 ; %
        idVal = idType; 
        idType = 'Gid';
    end    

    
    if isempty(allOpticDiscLocations)        
        if exist(opticDiscLocations_file, 'file') && ~redo
            S_file = load(opticDiscLocations_file);
            allOpticDiscLocations = S_file.allOpticDiscLocations;
            allLocIds = S_file.allLocIds;
        else
            [allOpticDiscLocations, allLocIds] = calcOpticDiskLocations;
            save(opticDiscLocations_file, 'allOpticDiscLocations', 'allLocIds', '-v6');                
        end                
    end
    
    if strcmp(idType, 'Gid')
        sd = siteDataFor('Gid', Gid);
        locId = sd.locationData.LocId;
    elseif strcmpi(idType, 'LocId')
        locId = idVal;
    end

    idx = find(allLocIds ==locId);
    if isempty(idx)
        error('No LocID = %d found', locId);
    end
    
    discLoc_pix = allOpticDiscLocations(idx).LeftDisc;
    

end




function [opticDiscLocations, allLocIds] = calcOpticDiskLocations
%     allLocationData = repmat(struct, maxGid,1);
    
    hnd = dbOpenExpDb;
    
%     fieldNames_loc = {'LOCATION_ID', 'PENETRATION_ID', 'TXT_LOCATION_NAME', ...
%         'DBL_LEFT_OPT_DISK_X', 'DBL_LEFT_OPT_DISK_Y', 'DBL_RIGHT_OPT_DISK_X', 'DBL_RIGHT_OPT_DISK_Y', 'DTM_CREATED'};
%     [LocID, PenId, LocName, LeftDiskX, LeftDiskY, RightDiskX, RightDiskY, dtm] = ...
%         getFieldsFromDatabaseTable(hnd, fieldNames_loc, 'TBL_LOCATIONS');            
% %     allLR = [LeftDiskX,LeftDiskY,RightDiskX,RightDiskY];
%     
%     fieldNames_pen = {'BLN_HEMISPHERE_R', 'BLN_HEMISPHERE_L', 'PENETRATION_ID', 'TXT_PENETRATION_NAME'};
%     [tf_PenLeft, tfPenRight, PenIds_key, PenName] = ...
%         getFieldsFromDatabaseTable(hnd, fieldNames_pen, 'TBL_PENETRATIONS');

%%
    tbl_merge = {'TBL_LOCATIONS', 'PENETRATION_ID', 'TBL_PENETRATIONS'};
    fieldnames_loc_pen = {'LOCATION_ID', 'TXT_LOCATION_NAME', ...
        'DBL_LEFT_OPT_DISK_X', 'DBL_LEFT_OPT_DISK_Y', 'DBL_RIGHT_OPT_DISK_X', 'DBL_RIGHT_OPT_DISK_Y', 'TBL_LOCATIONS.DTM_CREATED', ...
        'TBL_PENETRATIONS.PENETRATION_ID', 'BLN_HEMISPHERE_R', 'BLN_HEMISPHERE_L'};
    
    [LocID, LocName, LeftDiskX, LeftDiskY, RightDiskX, RightDiskY, loc_dtm, ...
        PenIds, tf_PenLeft, tfPenRight] = ...
        getFieldsFromDatabaseTable(hnd, fieldnames_loc_pen, tbl_merge);
%     allLR = [LeftDiskX,LeftDiskY,RightDiskX,RightDiskY];
        
    
    pen_hemiLR = repmat('L', length(tfPenRight), 1);
    pen_hemiLR(tfPenRight ~= 0) = 'R';
    %%
    
    nLocs = length(LocID);
    [uPenIds, uPenList] = uniqueList(PenId);
    3;
    assert( length(LocID) == length(unique(LocID)) );
    locationData(nLocs) = struct('LocID', [], 'PenID', [], 'LocName', [], 'LeftDisc', [], 'RightDisc', [], 'Hemisphere', []);
    
%     jj = find(LeftDiskX > 1);
%     figure(44); clf;
%     plot(LeftDiskX(jj),LeftDiskY(jj), 'b.', RightDiskX(jj),RightDiskY(jj), 'r.')
        
    fix_str = @(str) strrep(strrep(strrep(str, '200', '0'), ' AM', 'AM'), ' PM', 'PM');
    dtm_fix = cellfun(fix_str, dtm, 'un', 0);
    
    for p_i = 1:length(uPenIds)
        pen_idx = uPenList{p_i};       
        
        left_x = LeftDiskX( pen_idx );
        left_y = LeftDiskY( pen_idx );
        right_x = RightDiskX( pen_idx );
        right_y = RightDiskY( pen_idx );
        
        missing_tf = isnan( left_x ) | ( left_x == 0) |  ( left_x == 1);        
        idx_first_ok = find(~missing_tf, 1);
        
%         left_x_before = left_x;
        
        hemiLR = pen_hemiLR( pen_idx );
        
        idx_src = zeros(size(pen_idx));
        
        PenId_i = uPenIds(p_i);

        if all(missing_tf)
            [left_x, left_y, right_x, right_y] = deal(nan(size(pen_idx)));
%             continue;
        else

%             if all(~missing_tf)
%                 continue;
%             end

            if missing_tf(1)
    %             left_x(1) = left_x( find( ~missing_tf, 1) );                        
                idx_src(1) = find( ~missing_tf, 1);                     
                missing_tf(1) = 0;
            end

            idx_still_missing = find( missing_tf' );
            for i = idx_still_missing(1:end)
    %             left_x(i) = left_x(i-1);
                if i < idx_first_ok      
                    idx_src(i) = idx_first_ok;
                else
                    idx_src(i) = find( ~missing_tf(1:i), 1, 'last');
                end
            end

            missing_tf_after = (isnan( left_x ) | ( left_x == 0)) & (idx_src == 0);
            assert(~any(missing_tf_after));

            idx_dest = find(idx_src);
            left_x(idx_dest)  = left_x(idx_src(idx_dest));
            left_y(idx_dest)  = left_y(idx_src(idx_dest));
            right_x(idx_dest) = right_x(idx_src(idx_dest));
            right_y(idx_dest) = right_y(idx_src(idx_dest));        

            C = [LocName(pen_idx), num2cell([LocID(pen_idx), PenId(pen_idx), LeftDiskX( pen_idx ), left_x]), dtm_fix(pen_idx)]
            
        end
        
        for i = 1:length(pen_idx)
            pen_idx_i = pen_idx(i);
            loc_i = struct('LocID', LocID(pen_idx_i), 'PenID', PenId(pen_idx_i), 'LocName', LocName(pen_idx_i), ...
                'LeftDisc', [left_x(i), left_y(i)], 'RightDisc', [right_x(i), right_y(i)], 'Hemisphere', hemiLR(pen_idx_i) );
            
            locationData(pen_idx_i) = loc_i;
        end
        
%         [uV, vC] = uniqueCount(left_x);
%         if length(uV) > 1
%             fprintf('*******************\n');
%         end
        3;
    end
               3;         
            
            
end


%{
if missing_tf(1)
    idx_src(1) = find( ~missing_tf, 1);                     
    missing_tf(1) = 0;
end

idx_still_missing = find( missing_tf' );
for i = idx_still_missing(1:end)
    left_x(i) = left_x(i-1);
end

missing_tf_after = isnan( left_x ) | ( left_x == 0);
assert(~any(missing_tf_after));


        left_x_tmp = left_x_before;        
        idx_dest = find(idx_src);
        left_x_tmp(idx_dest) = left_x_tmp(idx_src(idx_dest));
        
        assert(isequal(left_x_tmp, left_x))


%}                

%{
    for i = 1:length(LocName)
        LN = LocName{i};
        PN = PenName{i};
        idx = strfind(LN, '.');
        LN = LN(1:idx(end)-1);
        assert(strcmp(LN, PN))
        fprintf('%s = %s\n', LN, PN)
    end
    %}