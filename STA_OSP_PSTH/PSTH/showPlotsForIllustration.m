
frmLength_ms = 16;
nBinsPerFrame = 3;
nFramesPerExtFrame = 5;
nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;

t = ([0.5:1:nBinsPerExtFrame] * (frmLength_ms/nBinsPerFrame));

nStimuli = 5;

r =   abs(randn(nBinsPerExtFrame, nStimuli)) * 1 ; %noise

    m = nBinsPerExtFrame/2;
    s = nBinsPerExtFrame/12;

    %add some peaks
    for i = 1:nStimuli
        thisPeak = histc( s*randn(1,5000)+m, 1:nBinsPerExtFrame )';
        thisPeak = thisPeak / max(thisPeak) * 15 * ((nStimuli-i)/nStimuli)^3;
        r(:,i) = r(:,i) + thisPeak;
    end
%     r = r(:,randperm(nStimuli));   % scramble the order

    figure(1);
    bar( t, r(:,1), .9 );
    xlims = xlim;
    set(gca, 'xtick', 0:frmLength_ms:xlims(2))
    
    Gid = 4470;
    gf = getFrameRetrieverFunction(Gid);
    gf('load', 'Gid', Gid);
    frm1 = gf(3);
    
    figure(2);
    imagesc(frm1);
    
    
    
    
    

