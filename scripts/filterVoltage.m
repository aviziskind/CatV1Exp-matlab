function [x, filt_state] = filterVoltage(b, a, x, filt_state, dim, nonlin, bartlettN)
    if ~isempty(dim)
        dim_arg = {dim};
    else
        dim_arg = {};
    end
    
    doRMS = 0;
    doNEO = 0;
        
    if exist('nonlin', 'var') 
        if isstruct(nonlin)
            doRMS = nonlin.useVoltageRMS;
            doNEO = nonlin.useVoltageNEO;
            bartlettN = nonlin.bartlettN;
        elseif ischar(nonlin)
            switch upper(nonlin)
                case 'RMS', doRMS = 1;
                case 'NEO', doNEO = 1;
            end        
        end
    end

    if ~exist('bartlettN', 'var')
        bartlettN = 0;
    end
    
    if doRMS
        x = x.^2;        
    end    
        
    if ~isequal([a b], [1 1])
        [x, filt_state] = filter(b, a, x, filt_state, dim_arg{:});
    end

    if doNEO
        x = nonlinEnergyOp(x, dim, bartlettN);
    end        
    
    if doRMS
        x = sqrt(x);
    end

    


    
    

end