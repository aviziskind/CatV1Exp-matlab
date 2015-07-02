function test_spikerAmps


    Gid = 4470;

    global spikeWaveforms t_ms
    global channelMeans CovMtx
    if isempty(spikeWaveforms)
        tic;
        [spikeWaveforms, t_ms] = retrieveSpikeWaveForms(Gid, 3, [-2.5, 2.5]);
        toc;
    end
    [nT, nChannels, nSpks] = size(spikeWaveforms);    
    spikeWaveforms = double(spikeWaveforms);
    
    if isempty(channelMeans)
        [channelMeans, CovMtx] = estimateChannelMeansAndCovariance(Gid);
    end    
    
    
%     troughWind_ms = [-.05, .3];
    troughWind_ms = [-.5, .5];
    idx_lookForMin = find( t_ms > troughWind_ms(1) & t_ms < troughWind_ms(2) );    
    
    spk = dbGetSpikes(Gid, 3, [], 1);
    spk = spk(:,2:end);
    
    
%     spikeWaveforms1 = double( spikeWaveforms(:,:,1) ); 
    
    amps = zeros(nSpks, nChannels);
    
    chOffsets = [41, -21, 6, -32];
    
    nItp = 10;
    
    progressBar('init-', nSpks*4)
    for spk_i = 1:nSpks
        wvfms_spk_i = spikeWaveforms(:,:, spk_i);
        wvfms_spk_i_ms = bsxfun(@minus, wvfms_spk_i, chOffsets);
        wvfms_mag_i = sum(wvfms_spk_i_ms, 2);
        [tmp1, tmp2, wvfms_mag_i_itp] = fourierInterp(wvfms_mag_i, [], nItp, [], t_ms);

        for ch_i = 1:4            
            progressBar;
            spk_wvfm = spikeWaveforms(:,ch_i, spk_i);
            [mn_val, idxs_min] = min(spk_wvfm(idx_lookForMin));
            idxs_min = idxs_min + idx_lookForMin(1)-1;                
            [tmp, t_ms_itp, spkWaveform_itp] = fourierInterp(spk_wvfm, [], nItp, [], t_ms);
%             idx0_itp = indmin(abs(t_ms_itp));
            idx0 = indmin(abs(t_ms));

            methodId = 2;
            
            % generate 'idx_min_itp_best'
            switch methodId
                case 1 % go back until find minimum. if negative, go forward until find minimum
                       % or: go forward/backward until find both minimums.
                    diffs = diff(spk_wvfm);
                    idx_min_best1 = find(diffs' < 0 & 1:length(diffs) < idx0, 1, 'last');
                    idx_min_best2 = find(diffs' > 0 & 1:length(diffs) > idx0, 1, 'first');
                    
                    if diffs(idx0) > 0
                    elseif diffs(idx0) < 0
                    end
                    
                    idx_lookForMinItp = idx_min_best*nItp + [-nItp:nItp];
                    idx_min_itp_best = indmin( spkWaveform_itp(idx_lookForMinItp)) + idx_lookForMinItp(1)-1;

                    
                case 2,  % look for minima in interpolated wvfm
                    
                    idx_lookForMin_itp = find( t_ms_itp > troughWind_ms(1) & t_ms_itp < troughWind_ms(2) );    
        %             [mn_val_itp, idxs_min_itp] = min(spkWaveform_itp(idx_lookForMin_itp));
        %             idxs_min_itp = idxs_min_itp + idx_lookForMin_itp(1)-1;

                    idx_minima = findLocalMinima(spkWaveform_itp(idx_lookForMin_itp), 7);
                    idx_minima = idx_minima + idx_lookForMin_itp(1)-1;
                    idx_minima = idx_minima(spkWaveform_itp(idx_minima) < 0);  % exclude ones that are above 0
        %             idx_closestTo0 = indmin(abs(t_ms_itp(idx_minima)));

                    vec = [-spkWaveform_itp(idx_minima)'].^2 ./ exp( abs( t_ms_itp(idx_minima)) );
                    idx_min_itp_best = indmax ( vec  );

                    idx_minima = idx_minima(spkWaveform_itp(idx_minima) < 0);  % exclude ones that are above 0
                    
    %             idx_closestTo0 = indmin(abs(t_ms_itp(idx_minima)));
    %             vec = [-spkWaveform_itp(idx_minima)'].^2 ./ exp( abs( t_ms_itp(idx_minima)) );
    %             idx_min_itp_best = indmax ( vec  );
                    
                case 3  % compare the minima before and after 0. take one with lower value.
            
        %             [mn_val_itp, idxs_min_itp] = min(spkWaveform_itp(idx_lookForMin_itp));
        %             idxs_min_itp = idxs_min_itp + idx_lookForMin_itp(1)-1;

                    idx_minima = findLocalMinima(spk_wvfm(idx_lookForMin), 1);
                    idx_minima = idx_minima + idx_lookForMin(1)-1;
                    % find the two that are closest to zero (before/after)
                    idx_min1 = find( (t_ms(idx_minima) <= 0), 1, 'last');
                    idx_min2 = find( (t_ms(idx_minima) > 0), 1, 'first');
                    mins = idx_minima([idx_min1, idx_min2]);
                    if length(mins) == 1
                        idx_min = mins;
                    elseif length(mins) == 2
                        idx_whichMin = indmin(spkWaveform_itp(mins));
                        idx_min = mins(idx_whichMin);
                    elseif isempty(mins)
                        3;
                    end

            end                            
            
%             idxs_min_itp = idx_minima(idx_min_itp_best);
%             mn_val_itp = spkWaveform_itp(idxs_min_itp);
            mn_val_itp = spkWaveform_itp(idx_min_itp_best);
            
            recordedAmp = spk(spk_i,ch_i);
            diffRec = (mn_val_itp - recordedAmp) - chOffsets(ch_i);
            
            x = spk(spk_i, ch_i);
            y = mn_val_itp;
            if isempty(y)
                y = 0;
                diffRec = 0;
                mn_val_itp = 0;
            end
            if y > 0
               3; 
            end
            
            if (norm([x;y] - [-162; 17.42]) < 1) %|| abs(diffRec) > 10

                
                showWorking = 1;
                if showWorking
                    figure(1); clf;
                    plot(t_ms, spk_wvfm, 'bo'); hold on
                    plot(t_ms_itp, spkWaveform_itp, 'b.');
                    plot(t_ms, wvfms_mag_i, 'go'); 
                    plot(t_ms_itp, wvfms_mag_i_itp, 'g.');
                    plot(t_ms_itp(idx_min_itp_best), mn_val_itp, 'rs', 'markersize', 12)

                    drawHorizontalLine(x+ chOffsets(ch_i), 'linestyle', ':', 'color', 'r')
                    drawVerticalLine(0, 'linestyle', ':', 'color', 'b')
%                     plot(t_ms(idxs_min), mn_val, 'go')
%                     plot(t_ms_itp(idxs_min_itp), mn_val_itp)
%                     plot(t_ms_itp(idx_minima), spkWaveform_itp(idx_minima), 'go', 'markerfacecolor', 'g')

%                     plot(t_ms(idxs_min), mn_val, 'go')

%                     drawHorizontalLine(median(spkWaveform_itp), 'linestyle', ':')
%                     drawVerticalLine(troughWind_ms, 'linestyle', ':')

%                     closestIdx = find( abs(spkWaveform_itp-recordedAmp) < 10);
%                     plot(t_ms_itp(closestIdx), spkWaveform_itp(closestIdx), 'r*')
                    xlim([-1 1])
                    3;
                    
                end
            end
        
            amps(spk_i,ch_i) = mn_val_itp;
            
            
        end
    end
    3;
    progressBar('done')

    assert( all ( size(spk) == size(amps) ));
    figure(2); clf;
    plot(spk, amps, '.');
    a = min(axis); b = max(axis);
    hold on;
    plot([a b], [a b], 'k:')
    axis equal square
    xlabel('recorded'); ylabel('just computed')
3;




end