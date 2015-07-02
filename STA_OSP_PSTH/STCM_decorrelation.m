function rf = STCM_decorrelation(Gid, cellId)

    % fieldName/Value: 'Gid', or 'Did', with corresponding Gid/Did value.
    % frameIds : index of frames 
    timeWindow_ms = [30 60];
    filename = ['STCMtmp_' num2str(Gid) '_' num2str(cellId)];
    
    relContrOfFrameToSpike = getParsedSpikes('frame', Gid, cellId, timeWindow_ms );
    Pspike_s = [relContrOfFrameToSpike{:}];
    Pspike_s = Pspike_s / 1.5*max(Pspike_s);
    
    getFrame = getFrameRetrieverFunction(Gid);
    getFrame('load', Gid, 'scaling', 'aroundZero');
    [nrows ncols] = elements(getFrame('size'));
    nPixels = nrows*ncols;
    
    if ~exist(filename, 'file');
        mean_SST_sp = zeros(nPixels, nPixels);
        mean_SST     = zeros(nPixels, nPixels);

        meanFrame = zeros(nrows, ncols);
        meanFrame_sp = zeros(nrows, ncols);

        nFramesEachPres  = cellfun(@length, relContrOfFrameToSpike);
        nTotalFrames     = sum(nFramesEachPres);

        progressBar('init', nTotalFrames, 30);
        for fi = 1:nTotalFrames
            progressBar(fi);
            S = getFrame(fi);
            meanFrame = meanFrame + S;
    %         currFrame = coarseGrain(currFrame, 2);

            SST = S(:) * S(:)';
            mean_SST = mean_SST + SST;
            if Pspike_s(fi) > 0
                mean_SST_sp = mean_SST_sp + Pspike_s(fi) * SST;
                meanFrame_sp = meanFrame_sp + S;
            end
        end
        mean_SST = mean_SST / nTotalFrames;
        mean_SST_sp = mean_SST_sp / nnz(Pspike_s);

        meanFrame = meanFrame / nTotalFrames;
        meanFrame_sp = meanFrame_sp/nnz(Pspike_s);

        save(filename, 'mean*');

    else
        load(filename);
    end

    A = mean_SST_sp - mean_SST;
    [v, e] = eig(A); e = diag(e);
    [topEigs, inds] = sort(e, 'descend'); 
%     keyboard;
    eigsToPlot = 3;
    figure(1); graySquareImage(meanFrame, 'mean Frame'); colorbar;
    figure(2); graySquareImage(meanFrame_sp,'mean Frame|spike'); colorbar;    
    mx = max(abs([meanFrame(:); meanFrame_sp(:)]));
    for ei = 1:eigsToPlot
        figure(ei+2); 
%         subplot(5,6, ei);
%             subplot(1,3,ei);
        vi = reshape(v(:,inds(ei)), nrows, ncols);
        mx = max([mx; abs(vi(:))]);
%         graySquareImage(vi); colorbar;
        graySquareImage(vi, [' eig = ' num2str(e(inds(ei)))]); colorbar;
    end
    for ei = 1:eigsToPlot;
        figure(ei+2); 
        caxis([-mx mx]);
    end

    
    
    % progressBar('init', nFramesToRetrieve, 50);
    getFrame('close');    
    


end


%         if showSTAs
%             subplot(1, length(delays_ms), si);
%             imagesc(STA);
%             axis equal tight xy; colormap('gray');
%             title([' STA: delay = ' num2str(delays_ms(si)) ' ms']);
%             % title([' STA: window = ' num2str(delays_ms(si)-windowSizes_ms(si)) ' - ' num2str(delays_ms(si)) ' ms']);
%             % toc;
%             fprintf('\n');
%         end
% 

