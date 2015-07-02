% testMIDdecorrelation
    global Xc dbug;
    hnd = dbOpenExpDb;

    stimulusOptions = {'noiseMovie', 'gratingMovie', 'naturalMovie', 'noise'};
    responseOptions = {'generatedResponse', 'actualResponse'};

    stimulusOption = stimulusOptions{2};
    responseOption = responseOptions{2};
    
    fig_id = 101;
	figure(fig_id); clf(fig_id);
    hrz_plots = 5;
    
    switch stimulusOption

        case 'naturalMovie'
            Gid = 3293;  % movieGroups( 540   542   552   556   562   570   580   586   594   598 );
            cellId = 0; %movieFileName = 'WALK1_IEEE_128X128X16384.RAW';
        
        case 'noiseMovie'
            Gid = 4620; % movieGroups(   806   810  820   824   828   836   848   852   856   874);
            cellId = 0; 
            
        case 'gratingMovie'
            Gid = 4796;
            cellId = 3;
                        
    end

    getFrame = getFrameRetrieverFunction(Gid);
    getFrame('load', Gid);
	[nrows ncols] = elements(getFrame('size'));
    getFrame('close');
    xs = 1:ncols;
    ys = 1:nrows;
	[xs_grid, ys_grid] = meshgrid(xs, ys);
    STA_timeWindow_ms = [30 60];
    
    switch responseOption
        case 'generatedResponse'

            % 1. Create a model Neuron with gabor function filter
            disp_debug('Creating new gabor filter...');

%             xlims = [-5 5]; Nx = 128;
%             ylims = [-5 5]; Ny = 128;
%             xs = linspace(xlims(1), xlims(2), Nx);
%             ys = linspace(ylims(1), ylims(2), Ny);

            Xc = [0;0];
            A = 1;
            x_mean = 1;
            y_mean = .7;
            x_std =  1;
            y_std =  1.4;
            k = 2;%rand*pi;
            phi = .2; 
            theta = pi/4;
            const = [];%0 %2*rand - 1;

            actualGaborParams = [A; x_mean; x_std; y_mean; y_std; k; phi; theta; const];
            gaborFunction = @(X) gabor(actualGaborParams, X);
            Nrep = 100;
    
            figure(fig_id);
            subplot(2,hrz_plots,1);
            zs = gaborFunction({xs_grid, ys_grid});
            graySquareImage(xs, ys, zs, 'Actual Filter');

            % 2. Generate spike train using some sample movie 
%             movieFileName = 'WALK1_IEEE_128X128X16384.RAW';
            disp_debug('Generating Spike Train ...');
            [spikeTrain, nrows, ncols] = generateSpikeTrainFromMovie(gaborFunction, simulusFileName, Nrep);

            % 3b. Get a noisy verstion of the filter for easy solving purposes
%             noisyGaborParams = actualGaborParams + (actualGaborParams/10) .* randn(size(actualGaborParams));    
%             noisyGaborFunc = @(X) gabor(noisyGaborParams, (X));

        case 'actualResponse'

%             assert(strcmp(stimulusOption, 'gratingStimulus'));
            disp_debug('Getting Parsed Spikes...');
            relContrOfFrameToSpike = getParsedSpikes('frame', Gid, cellId, STA_timeWindow_ms );
            relContrOfFrameToSpike = [ relContrOfFrameToSpike{:} ];
            Pspike_s = relContrOfFrameToSpike/(1.5*max(relContrOfFrameToSpike)); % convert to a pseudo-probability.
    end

%     disp_debug('Finding Important Pixels ...');
    mostImportantPixels = getMostImportantPixels('Gid', Gid, relContrOfFrameToSpike ) ;


    % 3a. Calculate the STA.
    disp_debug('Calculating STA...');
    N_lag = 1;
    
    STA = getSTAforCell( Gid, cellId, {STA_delay_ms, 'ms'}, {STA_windowSize_ms, 'ms'}, relContrOfFrameToSpike);

    
	figure(fig_id);
    subplot(2,hrz_plots,1);
    graySquareImage(xs, ys, STA, 'STA');
    
%     subplot(1, 6, 2);
%     graySquareImage(xs, ys, noisyGabor, 'Noisy Gabor');
    
    
    % Find MID using Gradient Search method.
    disp_debug('Calculating MID...');
    cellInfo = struct('Gid', Gid, 'cellId', cellId, 'Pspike_s', Pspike_s, 'importantPixels', {  });

    switch responseOption
        case 'generatedResponse'
            gb = findMostInformativeGabor(cellInfo, noisyGaborParams);

        case 'actualResponse'
            gaborParamsEstimate = estimateGaborParameters(xs, ys, STA);
            subplot(2,hrz_plots,hrz_plots+1);
            graySquareImage(xs, ys, gabor(gaborParamsEstimate, {xs_grid, ys_grid}) );

            gb = findMostInformativeGabor(cellInfo, gaborParamsEstimate);
            
    end

%     subplot(1,hrz_plots,3);
%     imagesc(MIfilter);
%     colormap('gray');
%     axis equal tight xy;
%     title('Most informative Gabor');

    