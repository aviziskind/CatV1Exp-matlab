s = who;
r = {'R', 'ori', 'sp', 'ph', 'tp_frm', 'tp_sec', 'stats', 'degPerPix', 'R_full', 'R_full_factor'};
for i = 1:length(s)
    nm = s{i};
    v = eval(nm);
    if isstruct(v) && isfield(v, 'OSP') && ~isfield(v.OSP, 'tp_sec')
        Gid = v.Gid;
        OSP = v.OSP;
        frameLength_ms = getFrameLength('Gid', Gid, 'ms');
        if strcmp(flashedOrDrifting(Gid, 's'), 'drifting')
            tps_frm = OSP.tp;
            tps_sec = tps_frm*(frameLength_ms/1000);
        else
            tps_frm = Inf;
            tps_sec = Inf;
        end
        if isfield(OSP, 'tp')
            OSP = rmfield(OSP, 'tp');
        end
        OSP.tp_frm = tps_frm;
        OSP.tp_sec = tps_sec;
        OSP = orderfields(OSP, r);
        v.OSP = OSP;
        eval([nm ' = v;']);
    end
    
end