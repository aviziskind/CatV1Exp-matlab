function [crossCorr, idxX, idxY] = getSpikeCrossCorrelogram_Matlab(X, Y, maxDist)
    % computes the (unbinned) cross-correlogram between spike trains X and Y.

    autoCross = ischar(Y) && strcmp(Y, 'auto');        
    if autoCross
        Y = X;
    end
    if (isempty(X)) || (isempty(Y))
        crossCorr = [];
        return;
    end
    outputIdxs = nargout > 1;    
    
    if ~issorted(X)
        error('X must be sorted');
%         X = sort(X, 'ascend');
    end
    if ~autoCross && ~issorted(Y)
        error('Y must be sorted');
%         Y = sort(Y, 'ascend');
    end
    
    relevant_x_idxs =  ibetween(X, Y(1)-maxDist, Y(end)+maxDist) ;    
    relevant_y_idxs =  ibetween(Y, X(1)-maxDist, X(end)+maxDist) ;    
    
    x_idx_offset = find(relevant_x_idxs, 1)-1;
    y_idx_offset = find(relevant_y_idxs, 1)-1;
    
    X = X(relevant_x_idxs);
    Y = Y(relevant_y_idxs);
    
    % algorithm is optimized for when x is longer than y.
    % if nx < ny, it makes sense to just switch x and y so that algorithm 
    % will be faster. (but make sure to flip the sign, b/c have switched X and Y, and C(X,Y) = -C(Y,X) )
    sgn = 1;
    if length(X) < length(Y);
        [X, Y] = deal(Y, X);
        sgn = -1;
    end
    X = X(:);
    Y = Y(:);
    nx = length(X);
    ny = length(Y);
    
    
    Cy = cell(1,ny);
    idxX_C = cell(1,ny);
    idxY_C = cell(1,ny);
    
    idx_x_1 = 1;
    
    
    for yi = 1:ny
        % 1. find the first X in range of current Y
        while (idx_x_1 <= nx) && (X(idx_x_1) < Y(yi)-maxDist) 
            idx_x_1 = idx_x_1+1;
        end
        if (idx_x_1 > nx)
            break;
        end            
        if (X(idx_x_1) > Y(yi)+maxDist) 
            continue;  % no Xs in range of current Y
        end
                
        % 2. find the last X in range of current Y
        idx_x_2 = idx_x_1;
        while (idx_x_2 < nx) && (X(idx_x_2+1) < Y(yi)+maxDist)
            idx_x_2 = idx_x_2+1;
        end            
        idx_x = idx_x_1:idx_x_2;
        
        if autoCross 
            if (length(idx_x) <= 1)
                continue;
            else
                idx_x(idx_x == yi) = [];
            end
        end

        % 3. find the differences, and put into C.
        Cy{yi} = Y(yi)-X(idx_x);
        if outputIdxs
            idxX_C{yi} = idx_x;
            idxY_C{yi} = ones(1, length(idx_x))*yi;
        end
%         if ~isempty(idx_x)
%             diffs = ;
            
%             idx_rm = find( abs(diffs) > maxDist );
%             if ~isempty(idx_rm)
%                 diffs(idx_rm) = [];
%             end            
%         end
        
    end
    
    doSorting = 0;
    
    crossCorr = sgn*vertcat(Cy{:});    
    if doSorting
        [crossCorr, newOrder] = sort(crossCorr);    
    end
    assert(all(abs(crossCorr) <= maxDist))

    if outputIdxs
        idxX = [idxX_C{:}]';
        idxY = [idxY_C{:}]';
        if doSorting
            idxX = idxX(newOrder);
            idxY = idxY(newOrder);
        end        
        
        if (sgn == -1)
            [idxX, idxY] = deal(idxY, idxX);
        end
        
        idxX = idxX + x_idx_offset;
        idxY = idxY + y_idx_offset;        
    end    
    
end



%     tic;
%     Cy = cell(1,ny);
%     for yi = 1:ny
%         idx_x_1 = binarySearch(X, Y(yi)-maxDist, -1, 1);
%         idx_x_2 = binarySearch(X, Y(yi)+maxDist, 1, -1);
%         idx_x = idx_x_1:idx_x_2;
%         if ~isempty(idx_x)
%             diffs = Y(yi)-X(idx_x);
%             diffs = diffs(abs(diffs) < maxDist);
%             Cy{yi} = diffs;
%         end
%         
%     end
%     Cy = sort(sgn*[Cy{:}]);
%     crossCorr = Cy;
%     t1 = toc;


%     tic;
%     for B = 1:100
%         Cx = cell(1,nx);
%         for xi = 1:nx
%             idx_y_1 = binarySearch(Y, X(xi)-maxDist, -1, 1);
%             idx_y_2 = binarySearch(Y, X(xi)+maxDist, 1, -1);
%             idx_y = idx_y_1:idx_y_2;
%             if ~isempty(idx_y)
%                 diffs = Y(idx_y)-X(xi);
%     %             if (idx_y_1 == idx_y_2) && (idx_y_1 == 1) || (idx_y_1 == ny)
%                     diffs = diffs(abs(diffs) < maxDist);
%     %             end
%                 Cx{xi} = diffs;
%             end
%         end
%         Cx1 = sort(sgn*[Cx{:}]);    
%     end
%     t2 = toc;
    
%     fprintf('%.2f\n', t1/t2);    

%     assert( isequal (Cx1, Cy1));


%     [X,Y, nx, ny, sgn] = deal(Y,X, ny, nx, -sgn);
%     tic;
%     Cy2 = cell(1,ny);
%     idx_x_1 = 1;
%     idx_x_2 = 1;        
%     for yi = 1:ny
%         while (idx_x_1 < nx) && (X(idx_x_1) < Y(yi)-maxDist)
%             idx_x_1 = idx_x_1+1;
%         end
%         idx_x_2 = idx_x_1;
%         while (idx_x_2 < nx) && (X(idx_x_2+1) < Y(yi)+maxDist)
%             idx_x_2 = idx_x_2+1;
%         end            
%         idx_x = idx_x_1:idx_x_2;
%         if ~isempty(idx_x)
%             diffs = Y(yi)-X(idx_x);
%             diffs = diffs(abs(diffs) < maxDist);
%             Cy2{yi} = diffs;
%         end
%         
%     end
%     Cy2 = sort(sgn*[Cy2{:}]);
%     t2 = toc;
