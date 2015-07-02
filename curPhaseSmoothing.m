function varargout = curPhaseSmoothing(smoothVal_in, smoothMethod_in)
    
    persistent smoothVal smoothMethod

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curSmoothing.mat'];

    setType = (nargin >= 1) && ~isempty(smoothVal_in);
    
    if setType   % set current Comparison type 

        if isnan(smoothVal_in)
            smoothVal = nan;
            smoothMethod = '';
        else
                        
            if ( smoothVal_in < 0 || smoothVal_in > 100 )
                error('must be in range [0 100]');
            end
            if nargin == 1 || isempty(smoothMethod_in) || ~any(strcmp(smoothMethod_in, {'Gauss', 'Fermi'}))
                error('must have second input (smoothMethod): "Gauss" or "Fermi"')
            end
            smoothVal = smoothVal_in;
            smoothMethod = smoothMethod_in;
            
        end
        
        save(filename, 'smoothVal', 'smoothMethod');
    end
    
    % retrieve current Comparison type
    if isempty(smoothVal)
        if exist(filename, 'file')
            load(filename);
        else
            error('Smoothing file does not exist');
        end
    end

    if nargout <= 1            
        if ~isnan(smoothVal)
            if strcmp(smoothMethod, 'Gauss')
                if smoothVal > 0
                    smoothVal_str_out = sprintf('_G%.1f', smoothVal);
                else
                    smoothVal_str_out = '';
                end
            elseif strcmp(smoothMethod, 'Fermi')
                smoothVal_str_out = sprintf('_F%d', smoothVal);
            end
        else
            smoothVal_str_out = '';
        end

        varargout = {smoothVal_str_out};
    elseif nargout == 2
        varargout = {smoothVal, smoothMethod};
    end            
        
    
end
