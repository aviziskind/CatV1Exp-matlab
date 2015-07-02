function testDecorrelatePSTH_basic

    showPSTHs = true; 
    nBinsPerFrame  = 3;
    nFramesPerExtFrame = 4;
    nStimuli = 6; 
    nPresEachStim = 20;
    plotStyle = '3D'; % '2D' or '3D'
        
    function s = normSum(r)
        r_sum = sum(r,1);        
        s = r_sum - min(r_sum);  
        s = s/max(s);
    end
    
%     function err = doDecorrelationTest
        
        nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;
        RandStream.setDefaultStream(RandStream('mt19937ar','seed',0));

        % Make up a frame sequence.
        stimIds = repmat(  1:nStimuli , 1, nPresEachStim);
        stimIds = stimIds(randperm(length(stimIds)));
        
        %****  Make up some r values; ****        
        r = randn(nBinsPerExtFrame, nStimuli) * 1 ; % take some base noise level 
                
        m = nBinsPerExtFrame/3;                     %add some peaks
        s = nBinsPerExtFrame/12;
        for i = 1:nStimuli
            thisPeak = histc( s*randn(1,5000)+m, 1:nBinsPerExtFrame )';
            thisPeak = thisPeak / max(thisPeak) * 15 * ((nStimuli-i)/nStimuli)^3;
            r(:,i) = r(:,i) + thisPeak;
        end
        
        r = r(:,randperm(nStimuli));   % scramble the order
        
%         xlims = [0 nBinsPerFrame *(nFramesPerExtFrame+1)];
%         ylims = [0 max(r(:))*1.1];
        
%         N = 12; 
%         showMoreThanN = false;        
%         nStimuliToShow = iff(showMoreThanN, nStimuli, min(nStimuli, N));
                
        % generate PSTH values
        R = generatePSTHs(r, stimIds, nBinsPerFrame, nFramesPerExtFrame);
                
        % Decorrelate PSTHs
        r2 = decorrelatePSTHs(R, stimIds, nBinsPerFrame, nFramesPerExtFrame);

        stimIndsInOrder = ord( sum(r,1), 'descend' );
        marg_r =  normSum(r);
        marg_R =  normSum(R);
        marg_r2 = normSum(r2);
        
        if showPSTHs
            plotPSTHseries({r, r2, r-r2, []}, plotStyle, {stimIndsInOrder});
            figure;
            plot(1:nStimuli, marg_r(stimIndsInOrder), 'bo-', ...
                 1:nStimuli, marg_R(stimIndsInOrder), 'rs-', ...
                 1:nStimuli, marg_r2(stimIndsInOrder), 'g.-');
            legend('r', 'R', 'r2');
        end
            
        err = mean(  abs( r(:)' - r2(:)') );
                
%     end

    
%     err = doDecorrelationTest;    
    disp(err);


end

