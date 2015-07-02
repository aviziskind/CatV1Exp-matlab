function [minVal, xMin] = findInterpMin(x, y, nItp, minRange)

    itpFirst = false;  % true==> first interpolate in the range, then look for minimum
                       % false==> first look for minimum, then interpolate around it and look for minimum again. 
    show = false;

    if ~exist('nItp', 'var') || isempty(nItp)
        nItp = 5;
    end    
    if ~exist('minRange', 'var')  || isempty(minRange) 
        minRange = [x(1), x(end)];
    end
%     get_xMin = nargout > 1;
    
    nX = length(x);
    idxLookForMin = find( x >= minRange(1) &  x <= minRange(2) );    
    if ~any(nX == size(y))
        error('length of x must be the same as the number of rows (or columns) of y');
    end
    if (nX ~= size(y,1))
        y = y';
    end
    nY = size(y,2);
    
    if length(idxLookForMin) < nX
        x_inRange = x(idxLookForMin);
        y_inRange = y(idxLookForMin, :);
    else
        x_inRange = x;
        y_inRange = y;
    end

    if show
        figure(423); clf; hold on;
        cols = ['bmkc'];
        for yi = 1:nY
            plot(x, y(:,yi), [cols(yi) 'o']); 
        end
    end

    
    if itpFirst
        x_inRange_itp = linspace(x_inRange(1), x_inRange(end), (length(x_inRange)-1)*nItp +1);
        y_itp = interp1(x_inRange, y_inRange, x_inRange_itp, 'spline');
        [minVal, idx_min] = min(y_itp);
        xMin = x_inRange_itp(idx_min);        
            
        if show            
            plot(x_inRange_itp, y_itp, '.')
            plot(xMin,  minVal, 'r*');            
        end
        
    else
        % 1. find minimum of Y in range.
        y_inRange = y(idxLookForMin,:);    
        [minVal_tmp, idx_min_tmp] = min(y_inRange);
        
        % 2. interpolate n points around the minimum. 
        % Look for minimum 1 point around the minimum
        nItpAroundMin = 10;
        nLookAroundMin = 1;

        x_itpAround = [-nItpAroundMin:nItpAroundMin]';        
        x_lookAround = [-nLookAroundMin*nItp:nLookAroundMin*nItp]';        
        
        idxNearMin_forItp = bsxfun(@plus, idxLookForMin(idx_min_tmp), x_itpAround);
        idxInBounds = idxNearMin_forItp >= 1 & idxNearMin_forItp <= nX;
        
        if any(~idxInBounds(:))  % have to do all interpolations separately.
        
            minVal = zeros(1,nY);
            xMin   = zeros(1,nY);
            
            for yi = 1:nY
                idxNearMin_i = idxNearMin_forItp(idxInBounds(:,yi), yi); % make sure don't look out of range.
                x_nearMin = x(idxNearMin_i);
                y_nearMin = y(idxNearMin_i,yi);

                x_nearMin_itp = linspace(x_nearMin(1), x_nearMin(end), (length(x_nearMin)-1)*nItp +1);
                y_itp = interp1(x_nearMin, y_nearMin, x_nearMin_itp, 'spline');

                % 3. find minimum of interpolated values within range.
                [tmp, idx_ofMin] = min(abs(x_nearMin_itp - x_inRange(idx_min_tmp(yi))));
                idx_lookForMin = idx_ofMin+x_lookAround;                
                idx_lookForMin = idx_lookForMin( idx_lookForMin >= 1 & idx_lookForMin <= length(y_itp));
                
                [minVal(yi), idx_min] = min(y_itp(idx_lookForMin) );
                idx_min = idx_min + idx_lookForMin(1) - 1;
                xMin(yi) = x_nearMin_itp(idx_min);                    
                
                if show                    
                    plot(x_nearMin_itp, y_itp, [cols(yi) '.']);
                    plot(x_inRange(idx_min_tmp), minVal_tmp, 'ro', 'markersize', 10);
                    plot(xMin(yi), minVal(yi), 'rs')
                end         
                
            end
            
        else
            % can do all interpolations together.            
            x_nearMin = x(idxNearMin_forItp(:,1));
            x_offsets = x(idxNearMin_forItp(1,:))-x_nearMin(1);
            
            idxNearMin_y = bsxfun(@plus, idxNearMin_forItp, [0:nY-1]*nX);
            y_nearMin = y(idxNearMin_y);
            
%             y2 = [y(idxNearMin_forItp(:,1),1), y(idxNearMin_forItp(:,2),2)];

            x_itp = linspace(x_nearMin(1), x_nearMin(end), (length(x_nearMin)-1)*nItp +1);            
            y_itp = interp1(x_nearMin, y_nearMin, x_itp, 'spline');

            idx_mins = ceil(length(x_itp)/2);
            idx_lookForMin = idx_mins+x_lookAround;
            
            % 3. find minimum of interpolated values
            [minVal, idx_min] = min(y_itp(idx_lookForMin,:));
            idx_min = idx_min + idx_lookForMin(1)-1;
            xMin = x_itp(idx_min)+x_offsets;        
            
            if show
                for yi = 1:nY
                    plot(x_itp+x_offsets(yi), y_itp(:,yi), [cols(yi) '.']);
                    plot(x_inRange(idx_min_tmp(yi)), minVal_tmp(yi), 'ro', 'markersize', 10)
                    plot(xMin(yi), minVal(yi), 'rs');
                end
                
            end         
            
        end
            
            
        
    end
        
    
        
    
    
end