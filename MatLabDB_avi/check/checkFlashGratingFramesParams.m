function checkFlashGratingFramesParams

    % check that occurences phase, sp_pix, ori_deg are evenly distributed
    % in flashed grating files

    S = load('cellsGroups_movie_fg');
    movieGroups = S.movieGroups_fg;
    movieFileNames = {movieGroups.movieFiles};
    movieFileNamesRAW = unique( [movieFileNames{:}]' );
    hnd = dbOpenExpDb;
    for mi = 1:length(movieFileNamesRAW)
        movieIds(mi) = getFieldsFromDatabaseTable(hnd, 'MOVIE_ID', 'TBL_MOVIES', {'TXT_MOVIE_FILE_NAME', movieFileNamesRAW{mi}});  %#ok<AGROW>
    end

    for mi = 1:length(movieIds)
        [phase_deg, sp_pix, ori_deg] = edbGetFlashGratIndex(hnd, movieIds(mi)); %#ok<NASGU>
        [vals, count] = uniqueCount(phase_deg);
        figure(mi);
        bar(1:length(vals), count);
        title(sprintf(movieFileNamesRAW{mi}), 'Interpreter', 'none');
        set(gca, 'xtick', 1:length(vals), 'xticklabel', cellfun(@(x) num2str(x), num2cell(vals), 'UniformOutput', false) );
        
    end

end


%     movieFileNames = {movieGroups.movieFiles};
% 
%     movieFileNamesRAW = unique( [movieFileNames{:}]' );
%     movieFileNamesTXT = strrep(movieFileNamesRAW, '.RAW', '.TXT');


%     allGids = [movieGroups.Gid];
%     for gi = 1:length(allGids)
%         movieGroups(gi).
% 
%     movieFileNamesRAW = unique( [movieFileNames{:}]' );
%     movieFileNamesTXT = strrep(movieFileNamesRAW, '.RAW', '.TXT');
%     
%     [allOri_deg, allSp_pix, allPhase_deg] = getOriSpPhaseForEachStimFrame(Gid)
