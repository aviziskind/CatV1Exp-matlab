function i = indmax_tmp(x, dim)
    if isvector(x)
        [val, i] = max(x);
        if nnz(x == val) > 1
            i = nan;
        end
    elseif isamatrix(x)
        if nargin < 2
            dim = 1;
        end
        if dim == 2
            x = x';
        end

        [vals, i] = max(x, [], 1);
        for j = 1:size(x, 2)
            if nnz( (x(:,j) == vals(j)) ) > 1
                i(j) = nan;
            end
        end
    end

end