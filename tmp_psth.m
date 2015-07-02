% psthWindowDataFile = [CatV1Path 'MatLabDB_avi' filesep 'psthWindowData_st_compressed'];
% psthWindowDataFile2 = [CatV1Path 'MatLabDB_avi' filesep 'psthWindowData_st'];
% 
S1 = load(psthWindowDataFile);
cell_names = fieldnames(S1);
% name_fun = @(s) sprintf('Gid_%04d_cell_%d_st', Gid, cellId);            

for i = 1:length(cell_names)
    v = S1.(cell_names{i});
    stat_names = fieldnames(v);
    
    ispval = find( cellfun(@(s) ~isempty(strfind(s, '_p')), stat_names ) );
    for j = ispval(:)'
        st = v.(stat_names{j});
        st(1:end-2) = pval2NegLogPval(  st(1:end-2)  );
        
        v.(stat_names{j}) = st;
        if any(v.(stat_names{j})(:) == 0)
            3;
        end
    end    
    v_decomp = structfun(@decompressOutNans, v, 'un', 0);
    
%     v_sparse = structfun(@(x) sparse(double(x)), v_decomp, 'un', 0);
    v_sparse = structfun(@compressOutNans, v_decomp, 'un', 0);
    
    S1.(cell_names{i}) = v_sparse;    
end

