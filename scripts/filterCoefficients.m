function [b, a] = filterCoefficients(n, Wn, ftype, funcName)

    try
        switch funcName            
            case 'butter',
                [b,a] = butter(n, Wn, ftype);
            case 'fir1',            
                b = fir1(n, Wn, ftype);
                a = 1;
        end
                            
    catch
        filterCoefficients_saved(n, Wn, ftype, funcName);
    end

end



function [b_out, a_out] = filterCoefficients_saved(n, Wn, ftype, funcName)
    % currently don't have the Signal Processing toolbox installed. this is
    % the result of calling 
    %    [b,a] = fir1(n, Wn, ftype)
    % for the specific parameters that I use for filtering the datafile.
    %  (n = 1, Wn = 300/10000, ftype = 'high')

    isLaptop = strcmp(getenv('computername'), 'AVI-PC');    
    redo = 0 && ~isLaptop;
            
    all_filterOrders = [1, 10, 20, 30, 40, 50];
    singleFreqs = [[400, 434, 476, 526, 588, 666, 769, 909, 1000], 300:50:1000];
    singleFreqs_C = num2cell(singleFreqs);
    all_filterFreqs = [singleFreqs_C, [100 250], [150, 250]];
    all_filterFreqs_types = cellfun(@(f) iff(length(f)==1, 'high', 'bandpass'), all_filterFreqs, 'un', 0);
    all_samplingFreqs = [20000, 10000];    
    all_fnames = {'butter', 'fir1'};
    
    filterCoef_file = [CatV1Path 'scripts\filterCoefs.mat'];
    
    if ~exist(filterCoef_file, 'file') || redo
        if 00 && isLaptop
            error('File not available');
        else  % create file.
            
            allCoefs = struct;
            for ord_i = 1:length(all_filterOrders)
                for sfreqs_i = 1:length(all_samplingFreqs);
                    for f_i = 1:length(all_filterFreqs);                    
                        
                        for name_i = 1:length(all_fnames);
                            n_i = all_filterOrders(ord_i);
                            Wn_i = all_filterFreqs{f_i} / (all_samplingFreqs(sfreqs_i)/2) ;
                            ftype_i = all_filterFreqs_types{f_i};
                            fname_i = all_fnames{name_i};
                            fn = getFilterFieldName(n_i, Wn_i, ftype_i, fname_i);
                            %                                 try
                            switch fname_i
                                case 'butter',
                                    [b,a] = butter(n_i, Wn_i, ftype_i);
                                    3;
                                case 'fir1',
                                    n_i_tmp = iff(odd(n_i), n_i+1, n_i);                                    
                                    [b] = fir1(n_i_tmp, Wn_i, ftype_i);
                                    a = 1;
                            end
                            %                                 catch
                            %                                     [b,a] = deal(0);
                            %                                 end
                            allCoefs.(fn).b = b;
                            allCoefs.(fn).a = a;
                            
                        end
                    end
                end
            end
                           
            save(filterCoef_file, '-struct', 'allCoefs');            
            fprintf('Generated filter coefficients file\n');
        end
        
    end
    
    if nargin < 1
        return;
    end
    allCoefs = load(filterCoef_file);
    
    fn = getFilterFieldName(n, Wn, ftype, funcName);
    b_out = allCoefs.(fn).b;
    a_out = allCoefs.(fn).a;

    % [b,a] = butter(1, 300/10000, 'high');

end


function s = getFilterFieldName(n, Wn, ftype, fname)
    if length(Wn) == 1 % 
        Wn_str = sprintf('%.5g', Wn);
    elseif length(Wn) == 2  % bandpass
        Wn_str = sprintf('%.5g_%.5g', Wn(1), Wn(2) );
    end    
    Wn_str = strrep(Wn_str, '.', '_');
    
    s = sprintf('%s__O%d__Wn%s__%s', fname, n, Wn_str, ftype);

end

