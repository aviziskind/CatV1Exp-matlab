function graySquareImage(varargin)
    % ways to call this function:
    %  graySquareImage(M)
    %  graySquareImage(x, y, M)
    %  graySquareImage(..., titleText);

    if isamatrix(varargin{1})
        M = varargin{1};
        [Nx, Ny] = size(M);
        xs = 1:Nx;
        ys = 1:Ny;
    elseif isamatrix(varargin{3})
        [xs, ys, M] = elements(varargin(1:3));
    end
    if ischar(varargin{end})
        titleText = varargin{end};
    end
    
    mx = maxElements(abs(M), 1);
%     eps = 1e-3;
    h = imagesc(xs, ys, M);
%     h = imagesc(xs, ys, M, [-mx -eps, mx +eps]);
    axis equal tight xy;
    colormap('gray');
%     set(gca, 'Position', [0 0 .95 1]);


    if exist('titleText', 'var')
        title(titleText);
    end
end