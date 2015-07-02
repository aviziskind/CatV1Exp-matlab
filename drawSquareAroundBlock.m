function [x1, x2, y1, y2, h] = drawSquareAroundBlock(ax, xs, ys, x_i, y_i, varargin)
%         x_i = x_i + 1;
        xlims = get(ax, 'xlim');
        ylims = get(ax, 'ylim');        
        
        xbnd = linspace(xlims(1), xlims(2), length(xs)+1); % x boundaries
        ybnd = linspace(ylims(1), ylims(2), length(ys)+1); % y boundaries
                
        dx = diff(xbnd(1:2));
        dy = diff(ybnd(1:2));
        r = .05;
        
        x1 = xbnd(x_i)  - dx*r;
        x2 = xbnd(x_i+1)+ dx*r;    
        y1 = ybnd(y_i)  - dy*r;
        y2 = ybnd(y_i+1)+ dy*r;

%         x1 = mean(xbnd(x_i-1:x_i))  - dx*r;
%         x2 = mean(xbnd(x_i  :x_i+1))+ dx*r;    
%         y1 = mean(ybnd(y_i-1:y_i) ) - dy*r ;
%         y2 = mean(ybnd(y_i  :y_i+1)) + dy*r;
        lb = [x1 y1]; ur = [x2 y2];
        
        h = drawSquare(lb, ur, varargin{:});

end