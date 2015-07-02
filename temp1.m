s = who('celldata*');
for i = 1:length(s)
    v = eval(s{i});
    
    if isstruct(v) && isfield(v, 'OSP') && isfield(v.OSP, 'R_full') && strcmp(class(v.OSP.R_full), 'double')
        R_full = v.OSP.R_full;
        
        f = 255/max(R_full(:));
        R_full_int = uint8(R_full*f);
                        
        v.OSP.R_full = R_full_int;
        v.OSP.R_full_factor = f;
        eval([s{i} ' = v;']);
    end
end

% save('indivCells_grating_full1_int', 'celldata*')