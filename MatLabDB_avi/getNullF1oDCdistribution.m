function [dist_m, dist_s] = getNullF1oDCdistribution(nbins, nsamples, Fn, f1Type, reflect)

    persistent FNoDCdata
    
    % f1Type --> Fn --> nbins --> nsamples
    
    redo = false;
    nBoot = 2000;

    minNSamples = 8;
    allNSamples_default = [minNSamples:10, 12:2:30, 35:5:200, 210:10:500, 520:20:1500, 1600:100:5000, 6000:500:30000];        
    
    if nargin < 3
        Fn = 1;
    end
    if nargin < 4
        f1Type = 'cos';
    end
    if nargin < 5
        reflect = false;
    end
        
    dataFile = [CatV1Path 'MatLabDB_avi' filesep 'FNoDCdata.mat'];    
    
    % load if empty
    if isempty(FNoDCdata) 
        if ~exist(dataFile, 'file') || redo        
            FNoDCdata = struct;
        else
            S = load(dataFile);
            FNoDCdata = S.FNoDCdata;
        end           
    end
    
    
    fld_name = sprintf('%s_%d%s', f1Type, Fn, iff(reflect, '_reflect', ''));
    
    % sub1: f1Type
    if isfield(FNoDCdata, fld_name) && any( [FNoDCdata.(fld_name).nbin] == nbins)
        id = find( [FNoDCdata.(fld_name).nbin] == nbins, 1 );
        dataStruct = FNoDCdata.(fld_name)(id);        
    else        
        [samp_means, samp_stds] = calcNullF1oDCdistribution(nbins, allNSamples_default, Fn, f1Type, reflect, nBoot);
        dataStruct = struct('nbin', nbins, 'Nsamples', allNSamples_default, 'means', samp_means, 'stds', samp_stds);
        if ~isfield(FNoDCdata, fld_name)
            FNoDCdata.(fld_name) = dataStruct;
        else
            FNoDCdata.(fld_name)(end+1) = dataStruct;
        end
        save(dataFile, 'FNoDCdata');
    end
            
    samp_ns = dataStruct.Nsamples;
    samp_means = dataStruct.means;
    samp_stds = dataStruct.stds;
    
    assert(nsamples < max(samp_ns))
    
    if ~ibetween(nsamples, min(samp_ns), max(samp_ns))
        dist_m = nan;
        dist_s = nan;
    else    
        dist_m = interp1(samp_ns, samp_means, nsamples);
        dist_s = interp1(samp_ns, samp_stds, nsamples);    
    end
    

end


