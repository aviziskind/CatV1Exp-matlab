function dbGetStimulusCenter

    hnd = dbOpenExpDb;
    % noise presentations
    
        
    stimFieldNames = {'DATAFILE_ID', 'LNG_ORIGIN_X', 'LNG_ORIGIN_Y', 'LNG_N_COLUMNS', 'LNG_N_ROWS', 'LNG_BLOCK_WIDTH_PIX', 'LNG_BLOCK_HEIGHT_PIX', 'DBL_DEGREES_PER_PIXEL'};    
    [Did, orig_x, orig_y, n_cols, n_rows, block_w, block_h, degPerPix] = ...
        getFieldsFromDatabaseTable(hnd, stimFieldNames, 'TBL_MOVIE_PRES');
    
    frame_w = n_cols .* block_w;
    frame_h = n_rows .* block_h;
    
    stim_UL = [orig_x, orig_y];
    stim_LR = [orig_x + frame_w, orig_y + frame_h];    
    stim_Center = [orig_x + (frame_w)/2, orig_y + (frame_h)/2];
    
    figure(1); clf;
    for i = 1:length(Did)
        rectangle('position', [orig_x(i), orig_y(i), frame_w(i), frame_h(i)] );
    end
    
    figure(2); clf;
    plot(stim_Center(:,1), stim_Center(:,2), '.');

    3;
    
    % movie presentations




end