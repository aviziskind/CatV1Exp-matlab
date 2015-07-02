function varargout = getGlobals(varargin)

    persistent globalVars
    
    if isempty(globalVars)        
    
        globalVars.isi_allRanges_ms = [1, 2, 4, 8, 16, 32, 64, 128, 256];
        globalVars.isi_allNbins     = [10, 20, 40, 80, 160];        
        globalVars.refrPeriod_ms_range = [.05 : .05 : 8.0];
        globalVars.maxRefrPeriod_ms = 8;
        globalVars.minRefrPeriod_ms = 0.85;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:length(varargin)
        if ~isfield(globalVars, varargin{i})
            error('Unknown global variable "%s"', varargin{i})
        end
    end
        
    varargout = cellfun(@(s) globalVars.(s), varargin, 'un', 0);


end