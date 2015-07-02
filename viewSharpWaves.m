function viewSharpWaves

    redo = 1;
    allSharpWavesFile = [RatHcPath 'allSharpWaves.mat'];
    if ~exist(allSharpWavesFile, 'file') || redo;
        Gids = ratDatafilesAvailable;
        for i = length(Gids):-1:1
            sd = siteDataFor(Gids(i));
            dfInfo = sd.dataFileInfo;
            
            s = load(getFileName('sharpWaves', Gids(i)));             
            s.swLengths_tk = double(s.sharpWaveEnds-s.sharpWaveStarts+1);
            s.swLengths_ms = s.swLengths_tk / dfInfo.samplingRateHz * 1000;
            s.Gid = Gids(i);
            
            allSW.allSharpWaves(i) = s;        
        end                
        save(allSharpWavesFile, '-struct', 'allSW');
    end

    S = [allSW.allSharpWaves];
    allWidths = cat(1, S.swLengths_ms);
    allHeights = cat(1, S.sharpWaveHeights);
    3;

    
    nInEach = arrayfun(@(s) length(s.sharpWaveStarts), S);
    nTot = sum(nInEach);
    
    filter_opt = struct('highPass_freq_Hz', 600, 'filterOrder', 1, 'filterName', 'median');
    filtered_file_name_ext = getName('dfFiltExt', filter_opt);

    
    for gi = 1:nTot
        sd = siteDataFor(Gids(gi));
        dfInfo = sd.dataFileInfo;
        file_info = getOpenFileHandle(dfInfo); %, struct('fileExt', filtered_file_name_ext));
        file_info_M = getOpenFileHandle(dfInfo, struct('fileExt', 'M'));
        sampRateHz = dfInfo.samplingRateHz;
        freq_sharp = [100, 250];
        
        
        for wi = 1:nInEach(gi)
                        %%

            rng = round(lims(double([S(gi).sharpWaveStarts(wi), S(gi).sharpWaveEnds(wi)]), .05));            
            smp = readSamples(file_info, dfInfo.channelIds, rng(1), diff(rng));
            smp = bsxfun(@minus, smp, smp(:,1));
            smp_mean = mean(smp, 1);
            smp_filt = fft_bandpass(smp_mean, freq_sharp/(sampRateHz/2), 2);
            
%             smp = abs(mean(smp, 1));
            
            med = readSamples(file_info_M, 1, rng(1), diff(rng));
            
            figure(1); clf;
            plot(smp'); hold on;
            plot(smp_filt, 'r-');
            plot(med, 'k-');
            title(sprintf('Total Length: %d ms', diff(rng) / sampRateHz * 1000))
            3;
        
        end
        fclose(file_info.fid);
        fclose(file_info_M.fid);
        
    end
    
    




end