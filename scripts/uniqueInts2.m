function u = uniqueInts2(X)

    if iscell(X)
        X = catAll(X);
    end
    
    i_max = max(X);
    tf = false(i_max, 1);
    tf(X) = true;
    
    u = find(tf(:));
end
    
function C = catAll(C)
    % concatenates all elements of the cell array into a long vector.
    if iscell(C)
        idx_sub_cell = cellfun(@iscell, C);
        C(idx_sub_cell) = cellfun(@catAll, C(idx_sub_cell), 'un', 0);        
    end
    C = cellfun(@(x) x(:), C, 'un', 0); % make all into column vectors;
    C = cat(1, C{:});
end
