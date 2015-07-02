function Y = insertNansAt(X, idx_nan)    
%     idx_inOrig = arrayfun(@(i) find(idx_or find(idx_nan
    Y = [];
    idx_nan = [1, idx_nan];
    for i = 1:length(idx_nan)-1
        Y = [Y; X(idx_nan(i):idx_nan(i+1)-1, :);  nan(1, size(X,2)) ];         %#ok<AGROW>
    end
    Y = [Y; X(idx_nan(end):end, :)];
end
