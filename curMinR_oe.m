function varargout = curMinR_oe(minR_oe_in, minR_oeMode_in)
    
    persistent minR_oe minR_oeMode

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curMinR_oe.mat'];

    setType = (nargin >= 1) && ~isempty(minR_oe_in);
    
    if setType   % set current Comparison type 

        if isnan(minR_oe_in)
            minR_oe = nan;
            minR_oeMode = '';
        else
                        
            if ( minR_oe_in < -1 || minR_oe_in > 1 )
                error('must be in range [-1 1]');
            end
            if nargin == 1 || isempty(minR_oeMode_in) || ~any(strcmp(minR_oeMode_in, {'oe', 'hoe', 'fs'}))
                error('must have second input (minR_oeMode): "oe" or "hoe" or "fs"')
            end
            minR_oe = minR_oe_in;
            minR_oeMode = minR_oeMode_in;
            
        end
        
        save(filename, 'minR_oe', 'minR_oeMode');
    end
    
    % retrieve current Comparison type
    if isempty(minR_oe)
        if exist(filename, 'file')
            load(filename);
        else
            error('MinRoe file does not exist');
        end
    end

    if nargout <= 1            
        if ~isnan(minR_oe)
            minR_oe_str_out = sprintf('_MinR_%s_%.2f', minR_oeMode, minR_oe);
        else
            minR_oe_str_out = '';
        end

        varargout = {minR_oe_str_out};
    elseif nargout == 2
        varargout = {minR_oe, minR_oeMode};
    end            
        
    
end
