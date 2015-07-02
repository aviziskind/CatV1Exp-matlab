function varargout = curWindowOffset(windowOffset_in)
    
    persistent windowOffset

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curWindowOffset.mat'];

    setType = (nargin >= 1) && ~isempty(windowOffset_in);
    
    if setType   % set current Comparison type 
                        
        if ( windowOffset_in < -2 || windowOffset_in > 2 )
            error('must be in range [-2, 2]');
        end
        windowOffset = windowOffset_in;
                    
        save(filename, 'windowOffset');
    end
    
    % retrieve current Comparison type
    if isempty(windowOffset)
        if exist(filename, 'file')
            load(filename);
        else
            error('WindowOffset file does not exist');
        end
    end

    windowOffset_str_out = iff(windowOffset == 0, '', strrep(sprintf('_W%d', windowOffset), '-', 'n'));

    if nargin == 1 && isempty(windowOffset_in)
        varargout = {windowOffset_str_out, windowOffset};
    else
        varargout = {windowOffset};
    end
            
%     elseif nargout == 2
%         varargout = {windowOffset, windowOffsetMode};
        
    
end


%{



function [windowOffset_out, windowOffset_str_out] = curWindowOffset(windowOffset_in)
    
    persistent windowOffset

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curWindowOffset.mat'];

    setType = (nargin == 1) && ~isempty(windowOffset_in);
    
    if setType   % set current Comparison type 
    
        if ~isnan(windowOffset_in) && ( windowOffset_in < 0 || windowOffset_in > 1 )
            error('must be in range [0 1]');
        end
        windowOffset = windowOffset_in;
        
        save(filename, 'windowOffset');
        
    else   % retrieve current Comparison type
        if isempty(windowOffset)
            if exist(filename, 'file')
                load(filename);
            else
                error('MaxPoe file does not exist');
            end
        end
        windowOffset_out = windowOffset;
        if ~isnan(windowOffset)
            windowOffset_str_out = sprintf('_MaxP_%.3f', windowOffset);
        else
            windowOffset_str_out = '';
        end        
        
        if (nargin==1) && isempty(windowOffset_in)
            windowOffset_out = windowOffset_str_out;
        end        

        
    end

end
%}