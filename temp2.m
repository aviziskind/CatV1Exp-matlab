for i = 1:length(gratingGroups)
    grp = gratingGroups(i);
    Gid = grp.Gid;
    g_type = flashedOrDrifting(Gid);
    if g_type == 1
        nph = length(grp.spPh_deg);
    elseif g_type == 2
        nph = round( grp.tempPeriod_sec / (grp.frameLength_ms/1000) );
    end
    dims = [length(grp.ori_deg), length(grp.spPeriod_pix), nph];
    cellIds = grp.cellIds;
    for j = 1:length(cellIds)
        nm = getName('celldata', Gid, cellIds(j));
        if exist(nm, 'var')
            vr = eval(nm);
            sz = [size(vr.OSP.R, 1), size(vr.OSP.R, 2), size(vr.OSP.R, 3)];
%             if isequal(sz, dims);
%                 fprintf('[%d, %d]: ok\n', Gid, cellIds(j));
%             end
            if ~isequal(sz, dims);
                if length(dims) ~= 3
                    3;
                end
                fprintf('[%d, %d]: should be: [%d x %d x %d], but is [%d x %d x %d]\n', Gid, cellIds(j), dims, sz);
            end
        end
    end
    
    
end