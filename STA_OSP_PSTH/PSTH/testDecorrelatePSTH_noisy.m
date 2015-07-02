function testDecorrelatePSTH_noisy

    showPSTHs = true; 
    nBinsPerFrame  = 3;
    nFramesPerExtFrame = 6;
    nStimuli = 100; 
    nPresEachStim = 30;
    plotStyle = '3D';
        
    stimIds = repmat(  1:nStimuli , 1, nPresEachStim);
    stimIds = stimIds(randperm(length(stimIds)));

    
    % make up some r values        
    nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;
    nStd = 1;
    frmWindow = [2 3];
    binWindow = [ [(frmWindow(1)-1) * nBinsPerFrame]+1, frmWindow(2) * nBinsPerFrame];
    r = makeUpPSTHs(nBinsPerExtFrame, nStimuli, binWindow, nStd);
        
    % generate PSTH values
    R = generatePSTHs(r, stimIds, nBinsPerFrame, nFramesPerExtFrame);
%     plotPSTHseries({r, R}, plotStyle);

    % RandStream.setDefaultStream(RandStream('mt19937ar','seed',0));
    
    function s = normSum(r)
        r_sum = sum(r,1);        
        s = r_sum - min(r_sum);  
        s = s/max(s);
    end


    function err = doDecorrelationTest(r, R, noiseStd, amtToSkip, eigThreshold)
        
        % ADD NOISE
        R = R + randn(size(R)) * noiseStd;
        
        % 3. Decorrelate PSTHs
        r2 = decorrelatePSTHs(R, stimIds, nBinsPerFrame, nFramesPerExtFrame, [], amtToSkip, eigThreshold);
        
        r3 = decorrelatePSTHs(R, stimIds, nBinsPerFrame, nFramesPerExtFrame, frmWindow, amtToSkip, eigThreshold);

        marg_r =  normSum(r);
        marg_R =  normSum(R);
        marg_r2 = normSum(r2);
        marg_r3 = normSum(r3);

        err_R  = sum( (marg_r - marg_R).^2 );
        err_r2 = sum( (marg_r - marg_r2) .^2);
        err_r3 = sum( (marg_r - marg_r3) .^2);

        err = [err_R, err_r2, err_r3];
        
        if showPSTHs                
            plotPSTHseries({r, R, r2}, plotStyle, {1:nStimuli, 441}, {'"hidden"', 'observed', 'calculated (all)'});            
            suptitle(['Noise: ' num2str(noiseStd) '. Skip: ' num2str(amtToSkip) '. Eig Threshold ' num2str(eigThreshold) ]);
            
            figure(442);
            plot(1:nStimuli, marg_r,  'bo-', ...
                 1:nStimuli, marg_R,  'r*-', ...
                 1:nStimuli, marg_r2, 'gs-', ...
                 1:nStimuli, marg_r3, 'mo-');
            title(['noise = ' num2str(noiseStd)]);
            legend('r', ['R : ' num2str(err_R)], ['r2 : ' num2str(err_r2)], ['r3 : ' num2str(err_r3)]);
        end
        
    end

    doDecorrelationTest(r, R, .5, 0, 0);
    return;
    % see what effect adding a little noise to each bin does.
    
    nTrials = 5;
    noiseStds  = logspace(-1.5, 0.5, 12);
    nsToIgnore = [0, 1, 2, 4];
%     eigThresholds = [0; .01];   
    vars = {noiseStds, nsToIgnore}; % {noiseStds, nsToIgnore, eigThresholds}
    labels = {'noise', {'R', 'r2', 'r3'}, 'N skipped'}; % {'noise', 'ns', 'th'}
    
    coreFunction = @(noiseS, nIg) doDecorrelationTest(r, R, noiseS, nIg, 0);    
    
    errs = iterateOverValues( coreFunction, vars, nTrials );
    plotAllResults(errs, vars, labels, {'semilogx', 'sub'});

    
    
end    
        

%     nFramesPerExtFrame = 4;
%     nStimuli = 5;
%     nPresEachStim = 10;
%     nTotalPres = nStimuli*nPresEachStim;
%     r = round( rand(nFramesPerExtFrame, nStimuli)*10 );
% %     r = [ [1;3;5], [2;1;2], [1;8;4] ];
% 
%     xlims = [0 nFramesPerExtFrame+1];
%     ylims = [0 9];
%     
%     stimIds = repmat(  [1:nStimuli], 1, nPresEachStim);
%     stimIds = stimIds(randperm(length(stimIds)));

    %     stimIds = [1 3 2 3 1 2 1 3 2 1 2 3 1 2 1 2 2 3 1];
