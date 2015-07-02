function [spikesForEachFrame nrows ncols] = generateSpikeTrain(spatialFilter, Did, Nrep)
    % spatial filter is either:
    % - a matrix representing which is dotted to each frame;
    % - a handle to a scalar 2D function

    %     xrange = [-5 5];   % dilemma: what xs & ys to use:
    %     yrange = [-5 5];

    %     maxLength_sec = 60;
    noise_sigma = 500;
    theta = 5000;
    %     latency_ms = 0;
    
    getFrame = getFrameRetrieverFunction(Gid);
    getFrame('load', Gid, 'scaling', 'aroundZero');
    nTotalFrames = getFrame('nTotalFrames');

    [nrows, ncols] = getFrame('size');
    %     fps = dbConvertTimeMeasures(Did, 1, 'sec', 'frame');
    %     npixels = nrows * ncols;

    spikesForEachFrame = zeros(1,nTotalFrames);

    if isnumeric(spatialFilter)
        zfilter = spatialFilter;
    elseif isa(spatialFilter, 'function_handle')
        xs = linspace(xrange(1), xrange(2), ncols);
        ys = linspace(yrange(1), yrange(2), nrows);
        [xs_grid,ys_grid] = meshgrid( xs , ys);
        zfilter = spatialFilter({xs_grid, ys_grid});
    end

    %     currents = zeros(1,nframes_max);

    %         fseek(fid, 0, 'bof');
    for frm_i = 1:nTotalFrames
        frm = getFrame(frm_i);
        %  movieFrame = fread(fid, [nrows ncols]);
        %  [movieFrame, nBytesRead] = fread(fid, [nrows ncols]);
        %  if( nBytesRead ~= npixels), error(['Error in reading file at frame: ' num2str(frm_i)]); end

        f = zfilter .* (frm);
        current = sum( f (:) );

        zeta = randn(1,Nrep);% might as well save computational time and do the repetitions in parallel:
        spikesForEachFrame(frm_i) =  sum ( stepFunction( current - theta + (noise_sigma * zeta) ) );
        %  currents(frm_i) = sum( f(:) ) - theta + (noise_sigma * r);
    end

%     spikesForEachFrame = spikesForEachFrame / Nrep;
%         figure(2); plot(currents);
%         figure(3); plot(spikes, '.-');

    getFrame('close');
end


%         for ri = 1:Nrep    % might as well save computational time and do the repetitions in parallel:
%             r = randn;
% %           currents(frm_i) = sum( f(:) ) - theta + (noise_sigma * r);
%             spikesForEachFrame(frm_i) = spikes(frm_i) + stepFunction( current - theta + (noise_sigma * r) );
%         end
