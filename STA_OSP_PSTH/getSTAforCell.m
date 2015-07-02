function STAs = getSTAforCell(Gid, cellId, timeWindowsV, windowProfilesC, relContrOfFrameToSpike, fig_id)

    if ~exist('timeWindowsV', 'var') || isempty(timeWindowsV)
        timeWindowsV = {[0 1], 'frame'};
    end
    [timeWindows, units] = elements(timeWindowsV);
    if ~exist('windowProfilesC', 'var') || isempty(windowProfilesC)
        windowProfilesC = cell(size(timeWindows,1),1);
    end    
    if ~iscell(windowProfilesC) % ok if one input
        windowProfilesC = {windowProfilesC};
    end
    givenRelContrOfFrameToSpike = (exist('relContrOfFrameToSpike', 'var') && ~isempty(relContrOfFrameToSpike));
    showSTAs = exist('fig_id', 'var');

    Did = dbLookup('Did',  'Gid', Gid);
    switch units
        case 'frame',  timeWindows_ms = dbConvertTimeMeasures(Did, timeWindows, 'frame', 'ms');
        case 'ms',     timeWindows_ms = timeWindows;
    end
        

    getFrame = getFrameRetrieverFunction(Gid);
    getFrame('load', Gid, 'scaling', 'aroundZero');
    [nrows ncols] = elements(getFrame('size'));

    STAs = zeros(nrows, ncols, size(timeWindows,1));
    if showSTAs
        figure(fig_id); clf;
    end
    for wi = 1:size(timeWindows,1)
        % tic;
        % fprintf(['Calculating STA # ' outOf(wi, length(delays_ms)) ])
        if ~givenRelContrOfFrameToSpike
            relContrOfFrameToSpike = getParsedSpikes('frame', Gid, cellId, timeWindows_ms(wi,:), windowProfilesC{wi} );
        end
        if iscell(relContrOfFrameToSpike)
            relContrOfFrameToSpike = [relContrOfFrameToSpike{:}];
        end
        relContrOfFrameToSpike = relContrOfFrameToSpike(:)';

        STA = zeros(nrows, ncols);
        significantFrames = find(relContrOfFrameToSpike);
        % progressBar('init', length(significantFrames), 30);

        for fi = significantFrames
            % progressBar;
            STA = STA + relContrOfFrameToSpike(fi) * getFrame(fi);
        end
%         STA = STA / sum(relContrOfFrameToSpike);


        if showSTAs
            subplot(1, length(delays_ms), wi);
            imagesc(STA);
            axis equal tight xy; colormap('gray');
            title([' STA: delay = ' num2str(delays_ms(wi)) ' ms']);
            % title([' STA: window = ' num2str(delays_ms(wi)-windowSizes_ms(wi)) ' - ' num2str(delays_ms(wi)) ' ms']);
            % toc;
            fprintf('\n');
        end

        STAs(:,:,wi) = STA;

    end

    getFrame('close');
end
% 
% 
%     
% end
% 
% 
% % good examples: 
% % (ncols = nrows), no missed frames, no preblank frames
% % Did: 1018. Gid = 865. cell 0
% 
% 

%{
    if ~exist('timeWindowsV', 'var') || isempty(timeWindowsV)
        timeWindowsV = {[0 1], 'frame'};
    end
    [timeWindows, units] = elements(timeWindowsV);
    if ~exist('windowProfilesC', 'var') || isempty(windowProfilesC)
        windowProfilesC = cell(size(timeWindows,1),1);
    end    
    if ~iscell(windowProfilesC) % ok if one input
        windowProfilesC = {windowProfilesC};
    end
    givenRelContrOfFrameToSpike = (exist('relContrOfFrameToSpike', 'var') && ~isempty(relContrOfFrameToSpike));
    showSTAs = exist('fig_id', 'var');

    Did = dbLookup('Did',  'Gid', Gid);
    switch units
        case 'frame',  timeWindows_ms = dbConvertTimeMeasures(Did, timeWindows, 'frame', 'ms');
        case 'ms',     timeWindows_ms = timeWindows;
    end
%}


%     if ~exist('delaysV', 'var') || isempty(delaysV)
%         delaysV = {0, 'ms'};            % [60 90 120];
%     end
%     if ~exist('windowSizeV', 'var') || isempty(windowSizeV)
%         windowSizeV = {1, 'frame'};     % [30 30 30];
%     end
