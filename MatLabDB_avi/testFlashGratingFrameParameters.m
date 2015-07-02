function testFlashGratingFrameParameters

    Gid = 4462;
    [allOri_deg, allSp_pix, allPhase_deg] = getOriSpPhaseForEachStimFrame(Gid);
    
    getFrame = getFrameRetrieverFunction(Gid);    
    getFrame('load', Gid);
    frameDims = getFrame('size');
    startFrame = 1;
    nframes = 8;
    specificFrameIds = find( (allOri_deg==0)&(allSp_pix ==42) );%startFrame + [1:2:nframes*2];
    frameIds = specificFrameIds(startFrame:startFrame+nframes-1);
    figure(242);
    for fi = 1:length(frameIds)
        frameId = frameIds(fi);
        ori = allOri_deg(frameId);
        sp  = allSp_pix(frameId);
        sph = allPhase_deg(frameId);
        
        subplot(2,nframes, fi);
        graySquareImage( getFrame(frameId) );
        axis xy;
%         title({['Ori = ' num2str(ori) '. Sp = '  num2str(sp)]});
        title({['Ori = ' num2str(ori) ], [ 'Sp = '  num2str(sp) ], ['Ph = ' num2str(sph)]});

        subplot(2,nframes, nframes + fi);                
        graySquareImage ( generateGratingFrame(frameDims, ori,  sp, sph) );
        axis xy;
        
    end
    
    


end