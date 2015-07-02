function varargout = curMaxP_oe(maxP_oe_in, maxP_oeMode_in)
    
    persistent maxP_oe maxP_oeMode

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curMaxP_oe.mat'];

    setType = (nargin >= 1) && ~isempty(maxP_oe_in);
    
    if setType   % set current Comparison type 

        if isnan(maxP_oe_in)
            maxP_oe = nan;
            maxP_oeMode = '';
        else
                        
            if ( maxP_oe_in < 0 || maxP_oe_in > 1 )
                error('must be in range [0 1]');
            end
            if nargin == 1 || isempty(maxP_oeMode_in) || ~any(strcmp(maxP_oeMode_in, {'oe', 'hoe', 'fs'}))
                error('must have second input (maxP_oeMode): "oe" or "hoe"')
            end
            maxP_oe = maxP_oe_in;
            maxP_oeMode = maxP_oeMode_in;
            
        end
        
        save(filename, 'maxP_oe', 'maxP_oeMode');
    end
    
    % retrieve current Comparison type
    if isempty(maxP_oe)
        if exist(filename, 'file')
            load(filename);
        else
            error('MaxPoe file does not exist');
        end
    end

    if nargout <= 1            
        if ~isnan(maxP_oe)
            maxP_oe_str_out = sprintf('_MaxP_%s_%.2f', maxP_oeMode, maxP_oe);
        else
            maxP_oe_str_out = '';
        end

        varargout = {maxP_oe_str_out};
    elseif nargout == 2
        varargout = {maxP_oe, maxP_oeMode};
    end            
        
    
end


%{



function [maxP_oe_out, maxP_oe_str_out] = curMaxP_oe(maxP_oe_in)
    
    persistent maxP_oe

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curMaxP_oe.mat'];

    setType = (nargin == 1) && ~isempty(maxP_oe_in);
    
    if setType   % set current Comparison type 
    
        if ~isnan(maxP_oe_in) && ( maxP_oe_in < 0 || maxP_oe_in > 1 )
            error('must be in range [0 1]');
        end
        maxP_oe = maxP_oe_in;
        
        save(filename, 'maxP_oe');
        
    else   % retrieve current Comparison type
        if isempty(maxP_oe)
            if exist(filename, 'file')
                load(filename);
            else
                error('MaxPoe file does not exist');
            end
        end
        maxP_oe_out = maxP_oe;
        if ~isnan(maxP_oe)
            maxP_oe_str_out = sprintf('_MaxP_%.3f', maxP_oe);
        else
            maxP_oe_str_out = '';
        end        
        
        if (nargin==1) && isempty(maxP_oe_in)
            maxP_oe_out = maxP_oe_str_out;
        end        

        
    end

end
%}