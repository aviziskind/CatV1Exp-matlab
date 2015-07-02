
% edits_pruning_auto --> edits_prunedClusts_states
% 
% clusterPruningsAuto --> prunedClustersStats

function tmp_renameFiles
    basepath = 'C:\ExperimentDB\Spikes\clusters\';
%     basepath = 'C:\ExperimentDB\Spikes\clusters\GLFcuw4\';

    
    chk_folderAndSub(basepath, 1);
    
    return;
    s = dir(basepath);
    s = s([s.isdir]);
    
    is_ok = cellfun(@(nm) ~any(strncmp(nm, {'.', '~'}, 1)), {s.name});
    s = s(is_ok);
    folders = {s.name};

    
    for i = 1:length(folders)
        f_i = [basepath folders{i} '\'];
        
%         orig_folder_name = [f_i 'edits_pruning_auto'];
%         new_folder_name = [f_i 'edits_prunedClusts_stats'];
        if exist(orig_folder_name, 'dir') || exist(new_folder_name, 'dir')
            if exist(orig_folder_name, 'dir')
                movefile(orig_folder_name, new_folder_name);
            end
            
            s_i = dir([new_folder_name]);
            s_i = s_i(~[s_i.isdir]);
            allFileNames = {s_i.name};
            idx_to_change = cellfun(@(nm) ~isempty(strfind(nm, 'clusterPruningsAuto')), allFileNames);
            %%
            for file_i = find(idx_to_change)
                new_file_name_i = strrep(allFileNames{file_i}, 'clusterPruningsAuto', 'prunedClustersStats');
                movefile([new_folder_name '\' allFileNames{file_i}], [new_folder_name '\' new_file_name_i]);           
            end
            fprintf('*')
        end
        
        
    end
    
end


function chk_folderAndSub(folderName, flag)
    
    s_i = dir([folderName '*allFalsePosCrossPruned*']);
    for i = 1:length(s_i);
        orig_name = s_i(i).name;
        new_name = strrep(orig_name, 'allFalsePosCrossPruned_', 'allGroupsCrossPrunedStats_');
        movefile([folderName  orig_name], [folderName  new_name]);           
        fprintf('%s => %s\n', orig_name, new_name);
    end
    
    s_i = dir([folderName '*allFalsePosPruned*']);
    for i = 1:length(s_i);
        orig_name = s_i(i).name;
        new_name = strrep(orig_name, 'allFalsePosPruned_', 'allGroupsPrunedStats_');
        movefile([folderName  orig_name], [folderName new_name]);           
        fprintf('%s => %s\n', orig_name, new_name);
        3;
    end    
    
    s = dir(folderName);
    subdirs = s([s.isdir] & cellfun(@(str) ~strncmp(str, '.', 1), {s.name}) );
    
    for j = 1:length(subdirs)        
        nm = subdirs(j).name;
        if flag == 1
            fprintf('--%s-- \n', nm)
        end
        chk_folderAndSub([folderName nm '\'], 0);
    end
    


end



% edits_pruning_auto --> edits_prunedClusts_states
% 
% clusterPruningsAuto --> prunedClustersStats