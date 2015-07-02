function u = uniqueInts(X)

    i_max = maxAll(X);
        
    tf = false(i_max, 1);
    
    tf = setIdxTrue(tf, X);
        
    u = find(tf(:));
end
    
function mx = maxAll(C)
    % concatenates all elements of the cell array into a long vector.
    
    if iscell(C)
        idx_sub_cell = cellfun(@iscell, C);
        C(idx_sub_cell) = cellfun(@maxAll, C(idx_sub_cell), 'un', 0);        
        idx_notEmpty = ~cellfun(@isempty, C);       
        mx = max(cellfun(@max, C(idx_notEmpty)));
    else
        mx = max(C);
    end
        
end

function tf = setIdxTrue(tf, C)

    
    if iscell(C)
        for i = 1:length(C)
            if iscell(C{i})
                tf = setIdxTrue(tf, C{i});
            else
                tf(C{i}) = true;
            end
        end                        
    else
        tf(C) = true;
    end


end