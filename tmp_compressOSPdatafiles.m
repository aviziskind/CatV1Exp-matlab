function tmp_compressOSPdatafiles

%     allWindowOspData_cells_DB_d_st_mean
    cd([CatV1Path 'MatLabDB_avi' filesep]);
    s_dir = dir('allWindowOspData*.mat');
    filenames = {s_dir.name};
    
    for fi = 1:length(filenames)
        filename = filenames{fi};

        S1 = load(filename);
%         S1 = load('allWindowOspData_cells_GLFcuw8_d_st_mean');
        q1 = whos('S1');
        size_before = q1.bytes;
    
        fn = fieldnames(S1);

        progressBar('init-', length(fn));
        for i = 1:length(fn)
            v = S1.(fn{i});

            flds = fieldnames(v);
            osp_fields_idx = cellfun(@(s) ~isempty(strfind(s, 'osp')), flds);

            flds = flds(osp_fields_idx);
            for j = 1:length(flds)
                fld_j = flds{j};
                idx_notEmpty = find(~cellfun(@isempty, v.(fld_j)));
                for k = idx_notEmpty
                    v.(fld_j){k} = compress(v.(fld_j){k});
                end
            end        

            S1.(fn{i}) = v;
            progressBar;
        end
    
        q2 = whos('S1');
        size_after = q2.bytes;
        fprintf('\n%s : %.1f MB -> %.1f MB (decrease of %.2f%%)\n', ...
            filename, size_before/1024^2, size_after/1024^2, (size_before-size_after)/size_before * 100);

        3;

        save(filename, '-struct', 'S1');

    end
end


