function dbCorrectLeftRightHemisphere
    hnd = dbOpenExpDb;

    [pen_right, pen_left, pen_name, pen_id, ap, ml, ap_ml0_id] = getFieldsFromDatabaseTable(hnd, ...
        {'BLN_HEMISPHERE_R', 'BLN_HEMISPHERE_L', 'TXT_PENETRATION_NAME', 'PENETRATION_ID', 'DBL_AP', 'DBL_ML', 'AP_ML_ZERO_ID'}, 'TBL_PENETRATIONS');
    idx_not_ok = find( ~xor(pen_right , pen_left) );
    
    [ap_ml0_id_key, ml0, ap0] = getFieldsFromDatabaseTable(hnd, ...
        {'AP_ML_ZERO_ID', 'DBL_ML_ZERO', 'DBL_AP_ZERO'}, 'TBL_AP_ML_ZERO');
    
    
    %%
    iscat = cellfun(@(s) s(1) == 'K', pen_name);
    pen_right = pen_right(iscat);
    pen_left = pen_left(iscat);
    pen_name = pen_name(iscat);
    pen_id = pen_id(iscat);
    ap = ap(iscat);
    ml = ml(iscat);
    ap_ml0_id = ap_ml0_id(iscat);

    cat_id = cellfun(@(s) str2double( strtok(s(2:end), '.') ), pen_name);
    
    [ucatid, cat_idx] = uniqueList(cat_id);
    
    %%
    
    
    for i = 1:length(ucatid)
        pen_right_i = pen_right(cat_idx{i});
        pen_left_i = pen_left(cat_idx{i});
%         if ~ all( xor(pen_right_i, pen_left_i))
%             continue;
%         end
       
        
        if 1 %% any(pen_right_i) && any(pen_left_i)
%             figure(1); clf; hold on;            

            figure(i); clf; hold on;
            
            mls = ml( cat_idx{i} );
            aps = ap( cat_idx{i} );
            
%             idx_R = (pen_right_i ~= 0) ;
%             idx_L = (pen_left_i ~= 0);

            idx_R = (pen_right_i) & (~pen_left_i);
            idx_L = (pen_left_i) & (~pen_right_i);
            
            
            ids = binarySearch(ap_ml0_id_key, ap_ml0_id( cat_idx{i} ));
            
            ml0s = ml0(ids); 
            ap0s = ap0(ids);

            idx_ok = between(ml0s, 1, 1000) & between(ap0s, 1, 1000);

            if ~isUnique(idx_ok)
                3;
            end
            
            mls = mls(idx_ok);
            aps = aps(idx_ok);
            ml0s = ml0s(idx_ok);
            ap0s = ap0s(idx_ok);
            idx_L = idx_L(idx_ok);
            idx_R = idx_R(idx_ok);
            
            
            mls = mls-ml0s;
            aps = aps-ap0s;
            
%             plot(ml0s, ap0s, 'k*');
            3;
            plot(mls(idx_L), aps(idx_L), 'bo', mls(idx_R), aps(idx_R), 'ro'); hold on;
            
            
            3;
            
            drawHorizontalLine(0, 'color', 'g', 'linewidth', 3); 
            drawVerticalLine(0, 'color', 'g', 'linewidth', 3);
            drawHorizontalLine(0, 'color', 'g', 'linewidth', 3); 
            drawVerticalLine(0, 'color', 'g', 'linewidth', 3);
            xlabel('ML'); ylabel('AP');
            
            3;
        end
        
        
    end

%             drawHorizontalLine(0); 
%             drawVerticalLine(0);
    
            3;
    
    3;
    %%
    if ~isempty(idx_not_ok)
       
        for i = 1:length(idx_not_ok)
            switch pen_id(i)
                case 32, % (neither selected) first 9/10 are in L hemisphere. assume this one is left as well
                    L = -1;
                    R =  0;
                case 163 % (L and R selected) all 11 are on the L hemisphere. but for #5, both are selected)
                    L = -1;
                    R =  0;                    
                case 178  % (L and R selected) all 3 are on the L hemisphere. but for 35, both are selected)
                    L = -1;
                    R =  0;                                        
                case 293 % (L and R selected) all 8 are on the L hemisphere. but for #6, both are selected)
                    L = -1;
                    R =  0;                                        
            end
        
%             updateValueInDatabaseTable(hnd, L, 'BLN_HEMISPHERE_L', 'TBL_PENETRATIONS', {'PENETRATION_ID', pen_id(i)} );
%             updateValueInDatabaseTable(hnd, R, 'BLN_HEMISPHERE_R', 'TBL_PENETRATIONS', {'PENETRATION_ID', pen_id(i)} );            
        end
    
    end
end


