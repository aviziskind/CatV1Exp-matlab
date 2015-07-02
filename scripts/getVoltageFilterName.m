function filt_name = getVoltageFilterName(filter_opt, nameMode, Gid)
    filt_name = '';
        
    if ~exist('nameMode', 'var') || isempty(nameMode)
        nameMode = 'dfile';
    end
    
    
    switch nameMode
        case 'dfile',   doFilterStrs = 1; doDetectStrs = 0; doMeasureStrs = 0;
        case 'detect',  doFilterStrs = 1; doDetectStrs = 1; doMeasureStrs = 0;
        case 'measure', doFilterStrs = 1; doDetectStrs = 1; doMeasureStrs = 1;
    end
    
    add_dfFilterType = nargin == 3 && ~isempty(Gid);
    
    if isfield(filter_opt, 'filterVoltageForDetection') && filter_opt.filterVoltageForDetection == 0
        return;
    end
    
    if doFilterStrs % add_dfFilterType
        
        if add_dfFilterType
            if Gid < 6000
                filtName_default = 'butter';
                freq_default = 300;
            else
                filtName_default = 'median';
                freq_default = 600;
            end
            if ~strcmp(filter_opt.filterName, filtName_default) || (filter_opt.highPass_freq_Hz ~= freq_default)
                df_filt_str = sprintf('%s_%d', filter_opt.filterName, filter_opt.highPass_freq_Hz);
                filt_name = addToFiltName(filt_name, df_filt_str);
            end

        end
    
       
    
        nzWindow = isfield(filter_opt, 'window_ms') && (filter_opt.window_ms) > 0;
        if nzWindow
            filt_name = [filt_name strrep(sprintf('wind_%.1f_ms', filter_opt.window_ms), '.', '_')];        
        end

    %     if nzWindow && isfield(filter_opt, 'window_shape') && ~isempty(filter_opt.window_shape)
    %         filt_name = [filt_name, '_' filter_opt.window_shape];        
    %     end
        if isfield(filter_opt, 'useVoltageRMS') && filter_opt.useVoltageRMS;
            if nzWindow 
                filt_name = [filt_name '_'];
            end            
            filt_name = [filt_name, 'rms'];
        end

        useVoltageNEO = isfield(filter_opt, 'useVoltageNEO') && filter_opt.useVoltageNEO;
        if useVoltageNEO        
            if nzWindow 
                filt_name = [filt_name '_'];
            end
            filt_name = [filt_name, 'neo'];
        end    
        if useVoltageNEO && isfield(filter_opt, 'bartlettN') && filter_opt.bartlettN > 0;
            filt_name = [filt_name, sprintf('_%d', filter_opt.bartlettN)];
        end        

    end
    
    
    
    if doDetectStrs 
        threshold_default = 8;

        if isfield(filter_opt, 'threshold') && (filter_opt.threshold ~= threshold_default);
            th_str = strrep(sprintf('th_%.1f', filter_opt.threshold), '.', '_');        
            filt_name = addToFiltName(filt_name, th_str);
        end
    end
        
    
    
    if doMeasureStrs
        itpMethod_default = 'sinc';
        interpN_default = 10;
        
        if filter_opt.interpN == 1
            filt_name = addToFiltName(filt_name, 'noInterp');
            
        elseif ~strcmp(filter_opt.interpMethod, itpMethod_default) || (filter_opt.interpN ~= interpN_default)
            itp_filt_str = sprintf('%s_%d', filter_opt.interpMethod, filter_opt.interpN);
            filt_name = addToFiltName(filt_name, itp_filt_str);
        end        
    end
    
    
    
end

function filt_name = addToFiltName(filt_name, extra_filt_str)
    if isempty(filt_name)
        filt_name = extra_filt_str;
    else
        filt_name = [filt_name '__' extra_filt_str];
    end        
    
end
