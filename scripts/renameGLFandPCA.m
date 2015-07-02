function renameGLFandPCA

    doOnAllSubdirs = 1;
    curDir = cd;
    if doOnAllSubdirs
        s_all = dir;
        s_all = s_all(3:end);
        s_all = s_all([s_all.isdir]);

        if ~isempty(s_all)
            for dir_i = 1:length(s_all)
                cd([curDir '\' s_all(dir_i).name]);
                renameGLFandPCA;
            end            
        end
    end
    cd(curDir);
    
    fprintf('Operating on dir %s... \n', curDir);
    s = dir;
    
    fets = {'GLFc', 'GLFs', 'PCAc', 'PCAs'};
    nFiles = length(s);
%     progressBar('init', nFiles, 30);
    for i = 1:nFiles        
        nm = s(i).name;
        
        for j = 1:length(fets)
            nm = strrep(nm, [fets{j} '_'], [fets{j} '2_']);
            nm = strrep(nm, [fets{j} '.'], [fets{j} '2.']);
        end
        if ~strcmp(s(i).name, nm)
            fprintf('%s --> %s\n', s(i).name, nm);
            3;
            movefile(s(i).name, nm);
        end
        
%         progressBar;
    end

end