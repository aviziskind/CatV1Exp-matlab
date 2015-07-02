function [Gid, cellId] = getGidCellIdFromVarname(varname)

    if strncmp(varname, 'celldata', 8)
        Gid_nm_start_ind = strfind(varname, 'Group')+6;
        varname2 = varname(Gid_nm_start_ind:end);
        Gid_nm_end_ind   = strfind(varname2, '_')-1;
        Gid = str2double(varname2(1:Gid_nm_end_ind));

        cellId_nm_start_ind = strfind(varname, 'Cell')+5;
        varname2 = varname(cellId_nm_start_ind:end);
        cellId = str2double(varname2);
        
        assert( strcmp(getName('celldata', Gid, cellId), varname) );    
    
    elseif  strncmp(varname, 'siteTest', 8)
        Gid_nm_start_ind = strfind(varname, 'Group')+6;
        varname2 = varname(Gid_nm_start_ind:end);        
        Gid = str2double(varname2);
        cellId = 0;
    end    
    
end