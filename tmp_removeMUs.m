
fn_d = fieldnames(allSpkStimHists.cells_GLFcuw8_d_st_ind);
fn_f = fieldnames(allSpkStimHists.cells_GLFcuw8_f_st_ind);
%%
t = 'f';
if t == 'f'
    fn_use = fn_f; 
elseif t == 'd'
    fn_use = fn_d; 
end

for i = 1:length(fn_use)
    cellId = nan;
    f_i = fn_use{i};
    A = sscanf(f_i, 'Gid_%d_cell_%d_st');
    if length(A) == 1
        A = sscanf(f_i, 'Gid_%d_cell_n%d_st');
        cellId = -A(2);
    elseif length(A) == 2
        cellId = A(2);
    end
    
    if cellId == 0 || cellId == 100 || cellId == -1
        if t == 'f'
            allSpkStimHists.cells_GLFcuw8_f_st_ind = rmfield(allSpkStimHists.cells_GLFcuw8_f_st_ind, f_i);
        elseif t == 'd'
            allSpkStimHists.cells_GLFcuw8_d_st_ind = rmfield(allSpkStimHists.cells_GLFcuw8_d_st_ind, f_i);
        end
        fprintf('Removed %s\n', f_i);
    end
end
