function h_label = axlabel(whichLabel, h_label, h_ax, varargin)
    % Adds an xlabel or ylabel to a graph, without adjusting the
    % size of the axes.
    
    firstTime = isempty(h_label) || (h_label == 0) || ~ishandle(h_label);

    if firstTime
        switch lower(whichLabel)
            case 'x', labelfun = @xlabel;
            case 'y', labelfun = @ylabel;
            case 'z', labelfun = @zlabel;
            otherwise, error('First parameter must be "x", "y" or "z"');
        end
        
        if nargin < 3 || isempty(h_ax)
            h_ax = gca;
        end        
        curPos = get(h_ax, 'position');
        h_label = labelfun(h_ax, varargin{:});
        if all(curPos(3:4) > 0)
            set(h_ax, 'position', curPos)
        end

    else
        if ischar(varargin{1})
            set(h_label, 'String', varargin{:})
        else
            set(h_label, varargin{:});
        end
    end

end