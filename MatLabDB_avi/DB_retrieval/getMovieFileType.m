function varargout = getMovieFileType(movieFilename)
%     if iscell(movieFilename)
%         GidCell = m2cell( repmat(Gid, size(movieFilename)));
%         [varargout{:}] = cellfun(@getMovieFileType, movieFilename, GidCell, 'UniformOutput', false);
%         return;
%     end

    movieTypes = {'Noise', 'Flashed_Gratings', 'Natural_Scenes'};
    
    Noise_names = {'bn', 'mcli', 'frm', 'medmovie_', 'mm', 'mv_', 'smtw_',  'mspn_', 'spn_', 'bar_'}; %#ok<NASGU>
        Noise_SubTypes = {'Band_Noise',    {'bn', 'frm', 'mcli', 'medmovie_', 'mm', 'mv', 'smtw'};
                          'Sparse_Noise',  {'mspn_', 'spn_'}; 
                          'Oriented_Bars', {'bar_'} };

    Flashed_Gratings_names = {'cph_', 'fg', 'mring', 'ring_'}; %#ok<NASGU>
        
    Natural_Scenes_names = {'photo_', 'vid_', 'walk'}; %#ok<NASGU>

    
    function tf = isType(filename, movieType)
        typeNames = eval([ movieType '_names']);
        minNameLength = min(cellfun(@length, typeNames));
        tf = any( strncmpi(filename, typeNames, minNameLength) ) ;
    end

    % Determing Movie type 
	movieType = 'Unknown';
    for ti = 1:length(movieTypes)
        if isType(movieFilename, movieTypes{ti})
            movieType = movieTypes{ti};
            break;
        end
    end
    varargout{1} = movieType;

    if nargout > 1
        if nargout == 3
            varargout{3} = 0;
        end
        % Determing Movie Subtype for Noise & Flashed Grating movies
        subType = '';
        if strcmp(movieType, 'Noise')        
            for st_i = 1:size(Noise_SubTypes,1)
                if any( strncmpi(movieFilename, Noise_SubTypes{st_i,2}, 2) )             
                    subType = Noise_SubTypes{st_i,1};
                    break;
                end
            end

        elseif strcmp(movieType, 'Flashed_Gratings')    
            hnd = dbOpenExpDb;
            movieId = getFieldsFromDatabaseTable(hnd, 'MOVIE_ID', 'TBL_MOVIES', {'TXT_MOVIE_FILE_NAME', movieFilename});
            [uOri, uSp, uPh, nRep] = dbGetUniqueOriSpPh('Mid', movieId);
            nOri = length(uOri); nSp = length(uSp); nPh = length(uPh);
%             try
%                 [phase_deg, sp_pix, ori_deg] = getFlashGratingOSPs_fgmovie(hnd, movieId);
%             catch
%                 [phase_deg, sp_pix, ori_deg] = deal(0);
%             end
%             nOri = length(unique(ori_deg));
%             nSp = length(unique(sp_pix));
%             nPh = length(unique(phase_deg));
%             nRep = length(phase_deg) / (nOri * nSp * nPh);
            subType = sprintf('%dx%dx%d', nOri, nSp, nPh);
            varargout{3} = nRep;
        end
        varargout{2} = subType;
    end
%     warning('Movie:unknownType', 'Unknown Movie File entered');
end


% 	subType = '';
%     if exist([movieType '_SubTypes'], 'var')
%         subTypes = eval([movieType '_SubTypes']);
%         for sti = 1:size(subTypes,1)
%             if any( strncmpi(movieFilename, subTypes{sti,2}, 2) )             
%                 subType = subTypes{sti,1};
%                 break;
%             end
%         end
%     end


% function movieType = getMovieFileType(movieFilename)
%     if iscell(movieFilename)
%         movieType = cellfun(@getMovieFileType, movieFilename, 'UniformOutput', false);
%         return;
%     end
% 
%     movieTypes = {'Oriented_Bars', 'Noise__Band', 'Noise__Sparse', 'Flashed_Gratings', 'Natural_Scenes'};
%     
%     Oriented_Bars_names = {'bar_'}; 
%     Flashed_Grating_names = {'cph_', 'fg', 'mring', 'ring_'}; 
%     Natural_Scenes_names = {'photo_', 'vid_', 'walk'}; 
%     Noise__Band_names = {'bn', 'frm', 'mcli', 'medmovie_', 'mm', 'mv', 'smtw'}; %#ok<*NASGU>
%     Noise__Sparse_names = {'mspn_', 'spn_'}; 
%     
%     function tf = isType(filename, movieType)
%         typeNames = eval([ movieType '_names']);
%         minNameLength = min(cellfun(@length, typeNames));
%         tf = any( strncmpi(filename, typeNames, minNameLength) ) ;
%     end
%     
% 	movieType = 'Unknown';
%     for ti = 1:length(movieTypes)
%         if isType(movieFilename, movieTypes{ti})
%             movieType = movieTypes{ti};
%             break;
%         end
%     end
% 
% %     warning('Movie:unknownType', 'Unknown Movie File entered');
% end
% 
