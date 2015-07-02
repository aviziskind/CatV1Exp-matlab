function tmp_testClust1
%     operateOnCurDir;	

%     fixKNN;
%     fixGLF;
    dir_name = 'C:\ExperimentDB\Spikes\clusters\';

%         case 'clusterData1',      fileTypeSubdir = ['clusters\' clustFet '\data1\'];              fileTypeName = ['clusters_data1_' clustFet];
%         case 'clusterData1_IC',   fileTypeSubdir = ['clusters\' clustFet '\data1_IC\'];           fileTypeName = ['clusters_data1_IC_' clustFet];

    subdirNames = subdirs(dir_name);
    idx_use = ~strncmp(subdirNames, '~', 1)  & ~strcmp(subdirNames, 'GLFcuw8');
    subdirNames = subdirNames(idx_use);
    
    progressBar('init-', length(subdirNames))
    for i = 1:length(subdirNames)
        data1_folder = [dir_name subdirNames{i} '\data1\'];
        s = dir([data1_folder 'Group_6*.mat']);
        for j = 1:10:length(s);
            S = load([data1_folder s(j).name]);
            assert(isempty(S.data.isiData.isis_mb))
        end
        
        progressBar(i);
    end
    
end


function subdirNames = subdirs(arg1)

    if ischar(arg1)
        dir_name = arg1;
        s = dir(dir_name);
    elseif isstruct(arg1)
        s = arg1;
    end
        
    isDir = [s.isdir] & ~arrayfun(@(si) any(strcmp(si.name(1), {'.', '~'})), s)';
    s = s(isDir);
    subdirNames = {s.name};
end
        
