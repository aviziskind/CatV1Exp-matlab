function test_midStimSpkFiles(Gid, cellId, all_nReps)
    
    N = length(all_nReps);
    sd = siteDataFor(Gid);
%     stimType = getGratingStimType(sd.stimType);
    npix = [sd.stimulusInfo.nrows];
    nStim = length(sd.ori_deg)*length(sd.spPeriod_pix)*length(sd.spPh_deg);

    s = calculatePSTH_STAs_OSP_ForOneCell(Gid, cellId);
    figure(100);
    subplot(1, N+1, 1);
    imagesc(s.STAs.STA);
    for i = 1:N        
        frameMode = sprintf('%drep', all_nReps(i));
        spikeFileName = mid_getSpikeFileName(Gid, cellId, [], [], frameMode);
        movieFileName = mid_getMovieStimFileName(Gid, frameMode);
        nFrames = nStim*all_nReps(i);
        STA_i = getSTA(spikeFileName, movieFileName, npix, nFrames);

        subplot(1,N+1, i+1);
        imagesc(STA_i);
        3;
    end
    3;

end

function STA = getSTA(spikeFileName, movieFileName, npix, nFrames)
    fid_spk = fopen(spikeFileName, 'r');
    fid_mov = fopen(movieFileName, 'r');
        
    nspks = double(fread(fid_spk, nFrames));
    nspks = nspks/max(nspks);
    
    STA = zeros(npix);
    for frm_i = 1:nFrames
        frm = double(fread(fid_mov, [npix, npix]));        
        frm = (frm/255) - .5;
        STA = STA + frm*nspks(frm_i);
    end
    fclose(fid_spk);
    fclose(fid_mov);

end