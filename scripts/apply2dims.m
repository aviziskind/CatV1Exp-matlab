function X = apply2dims(func, X, dims, dim_arg_id, concatFirst_flag)
    if ~exist('dim_arg_id', 'var') || isempty(dim_arg_id)
        dim_arg_id = 1;
    end
    concatDimsFirst = exist('concatFirst_flag', 'var') && isequal(concatFirst_flag, 1);
    
    %%
    if concatDimsFirst
                
        all_dims = 1: max(length(size(X)), max(dims));        
        all_dim_sizes = size(X);
        all_dim_sizes(end+1:length(all_dims)) = 1;
                            
        dims_av = dims;
        dims_leave = setdiff(all_dims, dims_av);
        
        % put all dims to average over at the end:
        X = permute(X, [dims_leave, dims_av]);
        
        dims_leave_sizes = all_dim_sizes(dims_leave);
        dims_av_size     = prod(all_dim_sizes(dims_av));
        if isempty(dims_leave_sizes)
            dims_leave_sizes = 1;
        end
        X = reshape(X, [dims_leave_sizes, dims_av_size]);
        last_dim_idx = length(dims_leave_sizes)+1;
        
        if dim_arg_id == 1
            X = func(X, last_dim_idx);
        elseif dim_arg_id == 2
            X = func(X, [], last_dim_idx);
        end                
        
    else
        for di = 1:length(dims)
            if size(X, dims(di)) > 1
                if dim_arg_id == 1
                    X = func(X, dims(di));
                elseif dim_arg_id == 2
                    X = func(X, [], dims(di));
                end
            end
        end
    end
    
end
