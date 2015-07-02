function [ori_deg, sp_pix, phase_deg] = getFlashGratingOSPs_fgmovie(movieId)

    if length(movieId) > 1
        [uMovieIds, movieIdx] = uniqueList(movieId);        
        [ori_deg_C_all, sp_pix_C_all, phase_deg_C_all] = deal(cell(1, length(movieId)));
        
        [ori_deg_C, sp_pix_C, phase_deg_C] = arrayfun(@getFlashGratingOSPs_fgmovie, uMovieIds, 'un', 0);                
        for mi = 1:length(uMovieIds)
            ori_deg_C_all(movieIdx{mi})   = ori_deg_C(mi);
            sp_pix_C_all(movieIdx{mi})    = sp_pix_C(mi);
            phase_deg_C_all(movieIdx{mi}) = phase_deg_C(mi);            
        end
        
        ori_deg = cat(1, ori_deg_C_all{:});
        sp_pix = cat(1, sp_pix_C_all{:});
        phase_deg = cat(1, phase_deg_C_all{:});
        return;
    end
    
    
    filename = getName('fgMovie_movieId', movieId);
    
    if ~exist(filename, 'file')
        [S.ori_deg, S.sp_pix, S.phase_deg] = retrieveFlashedGratingMovie_OSPs(movieId);
        save(filename, '-struct', 'S');        
    end    
    S = load(filename);
    ori_deg = double(S.ori_deg);
    sp_pix = double(S.sp_pix);
    phase_deg = double(S.phase_deg);        

end

function [ori_deg, sp_pix, phase_deg] = retrieveFlashedGratingMovie_OSPs(movieId)

    hnd = dbOpenExpDb;
    fullpath = dbGetMoviePathAndFilename(hnd, movieId);
    
    [path, filename, ext] = fileparts(fullpath); 
%     movieFileName = [filename, ext];
    
    idxFilename = fullfile(path, [filename, '.txt']);

    param_mtx = load(idxFilename);     % spatphase_deg, spatperiod_pix, orientation_deg
    if size(param_mtx,2) ~= 3
        error('Text file does not have 3 columns');
    end
    
    nFrames_DB = getFieldsFromDatabaseTable(hnd, 'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', movieId});
    if (size(param_mtx, 1) ~= nFrames_DB)
        error(['Text file does not have the correct number of rows (' num2str(nFrames_DB) ' instead of ' num2str(size(param_mtx, 1)) ')' ]);
    end
        
    phase_deg = param_mtx(:,1);
    sp_pix    = param_mtx(:,2);
    ori_deg   = param_mtx(:,3);

    if any(movieId == [112,  135:136,  137:138])
        phase_deg = round( rad2deg( phase_deg ) );
    end

    if any(movieId == [213:228])
        ori_deg = rem(360 + 90 - ori_deg,180);
    end
    
    if (movieId == 202)          
    %{
        This movie type was generated using a fourier transform method
        (See mm_fgcomplx2_1.m). The x and y parameters determine the x and y 
        components of the spatial frequency vector k, and A determines the 
        phase (-1 or 1 --> which are 180 degrees apart).
        For our purposes, we want movies with at least 4 (or 8) phases.        
        So we usually return an error if this kind of movie is accessed.
        Gids: [ 2625 2665 2677 2723 2737 2747 2775 2791 2805 2887 2921 2957 2981 3047 3061 3077 3091]
    %}        
        convertToOriSpPh = true;
        
        x = param_mtx(:,1); % [1..64]        x-component of spatial frequency vector
        y = param_mtx(:,2); % [1..64]        y-component of spatial frequency vector   
        A = param_mtx(:,3); % [either -1, 1]; controls phase /counterphase.

        if convertToOriSpPh
            ori_deg = rad2deg(atan2(x-1,y-1));        
            sp_pix = 64 ./ normV([x-1, y-1], 2);
            phase_deg = (A+1)*90;   % [0,2]*90 --> either 0 or 180    
        else        
            ori_deg = x;
            sp_pix  = y;
            phase_deg = A;
        end
                
    end
        
    ori_deg = single(ori_deg);
    sp_pix = single(sp_pix);
    phase_deg = single(phase_deg);    

end
% **************************************
% *** [s_ph_deg,  sp_pix,  ori_deg] ****
% **************************************
% cph_cir_10_2x64x64x5760_*  [230-245]
% fg_cir_10_2x64x64x2880_*   [213-228]	
% FGCIR10_2X64X64X2880_*	 [200-201]	
% MRINGLOGSP_2X64X64X28800_* [144-145]
% 
% 
% **************************************
% *** [s_ph_rad,  sp_pix,  ori_deg] ****
% **************************************
% mringlogsp_2x64x64x14400_* [137-138]
% mring_2x64x64x14400_*	     [135-136]	
% ring_2x128x128x14400@1     [112]
% 
% 
% ***********************************
% *********** [x, y, A] *************
% ***********************************
% fgcomp_2x64x64x8192_1   [202]



