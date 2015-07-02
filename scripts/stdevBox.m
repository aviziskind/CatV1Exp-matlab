function hout = stdevBox(h, x, y, data, w1, w2, col, hv)

    if isempty(h) || all(h == 0)
        h = zeros(5,1);
    end
    
    if iscell(data)
        mu = data{1};
        stdev = data{2};
        sem = data{3};
    elseif isnumeric(data)        
        mu = mean(data);
        stdev = std(data);
        sem = stderr(data);
    end
    if sem == 0
        sem = eps;
    end

    isDegreeData = mu > 10;  %if so, convert from 0:180 ==> -1:1;
        
    if isDegreeData
        shiftToRange = @(x) (1/90)*x -1;
        rescaleToInterval = @(x) x/90;
        mu = shiftToRange(mu);
        stdev = rescaleToInterval(stdev);
        sem = rescaleToInterval(sem);
    end
    stdev = 1;
    
    xc  = x+mu;
    xl2 = x-stdev;
    xr2 = x+stdev;
        
    lineX1 = [x; x];            % center line
    lineY1 = [y-w2; y+w2];

    lineX2 = [xl2    xl2  xr2;  % outside lines
              xr2    xl2  xr2];
    lineY2 = [y,  y-w2,  y-w2;
              y,  y+w2,  y+w2];
    rect_pos = [xc-sem, y-w1, 2*sem, 2*w1];
    if any(isnan(rect_pos))
        rect_pos = [-100 -100 .1 .1 ];
    end
          
    if (h(1) == 0)
        h(1) = rectangle('position', rect_pos, 'facecolor', col);
    else
        set(h(1), 'position', rect_pos, 'facecolor', col);
    end
    
    h(2)   = lineOrUpdateLine(h(2),   lineX1, lineY1, 'color', 'k', 'linestyle', ':');
    h(3:5) = lineOrUpdateLine(h(3:5), lineX2, lineY2, 'color', 'k');
%     h(1) = line([x; x], [y+w1, y-w1], lineparams{:});
%     h(2) = line([xl2; xr2], [y, y], lineparams{:});
%     h(3) = line([xl2; xl2], [y-w2, y+w2], lineparams{:});
%     h(4) = line([xr2; xr2], [y-w2, y+w2], lineparams{:});
   
    if nargout > 0
        hout = h;
    end
   
    
end