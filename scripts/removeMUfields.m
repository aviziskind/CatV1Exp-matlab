function S = removeMUfields(S)
    %%
%     redoBefore = 735445.462858;
    
    if nargin == 1
        S = removeMU_flds(S);
    elseif nargin == 0
        filenames = getAllFileNames;
        for i = 1:length(filenames)
            filename_i_full = filenames{i};
            [~, filename_i, ext_i] = fileparts(filename_i_full); 
%                         
            tic;
            fprintf('Processing %s%s ...', filename_i, ext_i);            
            file_s = load(filename_i_full);
            if length(fieldnames(file_s)) == 1
                fn1 = fieldnames(file_s);
                [file_s.(fn1{1}), n_removed, n_remaining] = removeMU_flds(file_s.(fn1{1}));
            else
                [file_s, n_removed, n_remaining] = removeMU_flds(file_s); %#ok<ASGLU>
            end
            save(filename_i_full, '-struct', 'file_s', '-v6');
            fprintf('done (removed %d fields; %d remaining)', n_removed, n_remaining);            
            toc;
        end
        
    end
    


end

function [S, n_removed, n_remaining] = removeMU_flds(S)
    fn = fieldnames(S);
    n_orig = length(fn);
    is_mu = @(str) ~isempty(strfind(lower(str), 'cell_0')) || ~isempty(strfind(lower(str), 'cell_n1')) || ~isempty(strfind(lower(str), 'cell_100')) || ...
                   ~isempty(strfind(lower(str), 'cellid_0')) || ~isempty(strfind(lower(str), 'cellid_n1')) || ~isempty(strfind(lower(str), 'cellid_100'));
    idx_fn_to_remove = cellfun(is_mu, fn);
    fn_to_remove = fn(idx_fn_to_remove);
    
    n_removed = length(fn_to_remove);
    
    S = rmfield(S, fn_to_remove);

    n_remaining = length(fieldnames(S));
    assert(n_orig - n_removed == n_remaining);
end

function s = getAllFileNames()
   

s = {[CatV1Path 'MatlabDB_avi' filesep 'allBckgSpikes.mat'];
     [CatV1Path 'MatlabDB_avi' filesep 'allPsthWindowData_cells_GLFcuw8_d_st_mean.mat'];
    [CatV1Path 'MatlabDB_avi' filesep 'allPsthWindowData_cells_GLFcuw8_f_st_mean.mat'];
    [CatV1Path 'MatlabDB_avi' filesep 'allWindowOspData_cells_GLFcuw8_d_st_mean.mat'];
    [CatV1Path 'MatlabDB_avi' filesep 'allWindowOspData_cells_GLFcuw8_f_st_mean.mat'];
    [CatV1Path 'MatlabDB_avi' filesep 'allCellPSTHs_GLFcuw8.mat'];    
    [CatV1Path 'indivCells_GLFcuw8_movie_fg.mat'];
    [CatV1Path 'indivCells_GLFcuw8_grating_dOr.mat'];
    [CatV1Path 'indivCells_GLFcuw8_grating_dSf.mat'];
    };
    ...'allSpkStimHists_cells_GLFcuw8_d_st_ind.mat';
    ...'allSpkStimHists_cells_GLFcuw8_f_st_ind.mat';    

end