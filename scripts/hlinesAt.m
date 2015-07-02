function hlinesAt(x1, x2, y, varargin)
    nx = length(x1);
    assert(length(x2) == nx);
    line([x1(:)'; x2(:)'], [y;y]*ones(1, nx), 'color', 'k', varargin{:});

end