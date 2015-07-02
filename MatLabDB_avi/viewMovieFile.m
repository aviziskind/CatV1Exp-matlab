function viewMovieFile(filename)

    if nargin < 1
        filename = 'fgcomp_2x64x64x8192_1.raw';
    end

    hnd = dbOpenExpDb;
    [movieId, nCols, nRows, nFramesTotal] = getFieldsFromDatabaseTable(hnd, {'MOVIE_ID', 'LNG_N_COLUMNS', 'LNG_N_ROWS', 'LNG_N_FRAMES'}, 'TBL_MOVIES', ...
        {'TXT_MOVIE_FILE_NAME', upper(filename)});
    
    [arg1, arg2, arg3] = getFlashGratingOSPs_fgmovie(hnd, movieId);
    Gid = 3091;
    getFrame = getFrameRetrieverFunction(Gid);
    getFrame('load', Gid);
    uArg1 = unique(arg1);
    uArg2 = unique(arg2);
    uArg3 = unique(arg3);
    
    frms = find( (arg2 == uArg2(4)) & (arg3 == uArg3(1)) );
%     arg2valInds = ord ( arg2(frms) );
%     frms = frms(arg2valInds);

    inds = ord ( arg1(frms) );
    frms = frms(inds);

    
    numFramesX = 4;
    numFramesY = 2;
    nFramesToShow = 100;
    side_len = unique([nCols nRows]);
    fid = fopen(filename);
    figure(1); clf;
    xy_i = 1;
    gridSubPlot(numFramesY, numFramesX, [1 10]);
    for i = 1:length(frms)
        frmId = frms(i);
        gridSubPlot;

        thisFrame = getFrame(frmId);% fread(fid, [side_len side_len]);        
        imagesc(thisFrame); colormap('gray');
        title(sprintf('(%d/%d)Frame #%d [%d, %d, %d]', i, length(frms), frmId, arg1(frmId), arg2(frmId), arg3(frmId) ));
        set(gca, 'xtick', [], 'ytick', []);
        axis square
        if ~mod(i, numFramesX *numFramesY)
            input('');
        end
    end
    fclose(fid);
    
end