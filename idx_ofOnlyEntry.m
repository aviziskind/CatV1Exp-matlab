function i = idx_ofOnlyEntry(x, dim)
    if isvector(x)
        i = nan;
        if nnz(x) == 1
            i = find(x);
        end
    elseif isamatrix(x)
        if nargin < 2
            dim = 1;
        end
        if dim == 2
            x = x';
        end

        i = nan(1, size(x,2));
        for j = 1:size(x, 2)
            if nnz( x(:,j) ) == 1 
                i(j) = find( x(:,j) );
            end
        end
    end

end
