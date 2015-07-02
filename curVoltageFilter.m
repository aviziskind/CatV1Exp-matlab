function voltageFilterStruct_out = curVoltageFilter(varargin)
    persistent voltageFilterStruct  voltageFilterStruct_default

%     useVoltageRMS  useVoltageNEO  window_ms bartlettN threshold
    useDefaultForCatV1 = 1;

    setType = (nargin >= 1) && ~isempty(varargin{1});
    
    filename = [CatV1Path 'MatLabDB_avi' filesep 'curVoltageFilter.mat'];
    
    if isempty(voltageFilterStruct_default)
        if useDefaultForCatV1
            EC_highPass_freq_Hz_default = 300;
            EC_filterOrder_default = 1;
            EC_filterName_default = 'butter';
        else
            EC_highPass_freq_Hz_default = 600;
            EC_filterOrder_default = 1;
            EC_filterName_default = 'median';
        end
        
        window_ms_default = 0;
        doRMS_default = false;
        doNEO_default = false;
        bartlettN_default = 10;
        spikeDetectionThreshold_default = 8;
        
        interpMethod_default = 'sinc';
        interpN_default = 10;
        
  %     filter_opt_voltage = struct('highPass_freq_Hz', 600, 'filterOrder', 1, 'filterName', 'median');        

        allFieldNames = {'highPass_freq_Hz', 'filterOrder', 'filterName', ...
            'window_ms', 'useVoltageRMS', 'useVoltageNEO', 'bartlettN', 'threshold', 'interpMethod', 'interpN'};
        allFieldDefaults = {EC_highPass_freq_Hz_default, EC_filterOrder_default, EC_filterName_default, ...
            window_ms_default, doRMS_default, doNEO_default, bartlettN_default, spikeDetectionThreshold_default, ...
            interpMethod_default, interpN_default};    
       
        S = [allFieldNames; allFieldDefaults];
        voltageFilterStruct_default = struct(S{:});        
    end
        
    
    if isempty(voltageFilterStruct)            
        if exist(filename, 'file')
            S = load(filename);
            voltageFilterStruct = S.voltageFilterStruct;
        else
            voltageFilterStruct = voltageFilterStruct_default;
        end
    end
    
    
    if setType   % set current pair types                            
        if isstruct(varargin{1})
            S_input = varargin{1};
            fn = fieldnames(S_input);
            for i = 1:length(fn)
                field_nm_i = fn{i};
                if ~isfield(voltageFilterStruct, field_nm_i)
                    error('Invalid field name : %s', field_nm_i); 
                end
                voltageFilterStruct.(field_nm_i) = S_input.(field_nm_i);                
            end
        
        elseif ischar(varargin{1})
               
            if strcmp(varargin{1}, 'default')
                voltageFilterStruct = voltageFilterStruct_default;
            else
                            
                if (mod(length(varargin),2) ~= 0), error('Number of arguments must be a multiple of 2'); end            
                for i = 1:nargin/2
                    field_nm =  varargin{i*2-1};
                    field_val = varargin{i*2};
                    if ~isfield(voltageFilterStruct, field_nm)
                        error('Invalid field name : %s', field_nm); 
                    end
                    voltageFilterStruct.(field_nm) = field_val;
                end
                
            end
            
        else
            error('Invalid input syntax');
        end
            
        save(filename, 'voltageFilterStruct');        
        
    end
        
    voltageFilterStruct_out = voltageFilterStruct;        

end
