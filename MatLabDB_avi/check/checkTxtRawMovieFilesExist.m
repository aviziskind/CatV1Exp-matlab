function checkTxtRawMovieFilesExist

    function s = nameOfFile(fullFileName)
        ind = strfind(fullFileName, '.');
        s = fullFileName(1:ind-1);
    end

    function doCheckFunction(ext)


        % get names of files in directory
        filenames_S = dir('C:\ExperimentDB\Movies\Database');
        filenames = {filenames_S.name};
        filenames = cellfun(@lower, filenames, 'UniformOutput', false)';

        inds = cellfun(@(s) ~isempty( regexp(s, ['.' ext], 'once')), filenames);
        filesInDir = filenames(inds);
        filesInDirNames = cellfun(@nameOfFile, filesInDir, 'Uniformoutput', false);

        % get names of files in group data file.
        S = load('cellsGroups_movie');
        movieGroups = S.movieGroups;
        flashGratingInds = strcmp('Flashed_Gratings', {movieGroups.moviesType});
        flashGratingGroups = movieGroups(flashGratingInds);
        filesNeeded = {flashGratingGroups.movieFiles};
        filesNeeded = unique( [filesNeeded{:} ]');
        filesNeededNames = cellfun(@nameOfFile, filesNeeded, 'Uniformoutput', false);
        filesNeededNames = cellfun(@lower, filesNeededNames, 'UniformOutput', false)';

        disp(['******  ' upper(ext) ' files not used : ******  ']);
        setdiff(filesInDirNames, filesNeededNames)' %#ok<NOPRT>

        disp(['******  ' upper(ext) ' files that need but don''t have: ******  ']);
        setdiff(filesNeededNames, filesInDirNames)' %#ok<NOPRT>
    end

    doCheckFunction('txt');
    doCheckFunction('raw');
    
end


% RESULTS:
%
% ******  TXT files not used : ******  
% ans = 
%     'bn_64x64x4096_1'
%     'cph_cir_10_2x64x64x5760_10'
%     'cph_cir_10_2x64x64x5760_11'
%     'cph_cir_10_2x64x64x5760_12'
%     'cph_cir_10_2x64x64x5760_13'
%     'cph_cir_10_2x64x64x5760_14'
%     'cph_cir_10_2x64x64x5760_15'
%     'cph_cir_10_2x64x64x5760_16'
%     'cph_cir_10_2x64x64x5760_5'
%     'cph_cir_10_2x64x64x5760_6'
%     'cph_cir_10_2x64x64x5760_7'
%     'cph_cir_10_2x64x64x5760_8'
%     'cph_cir_10_2x64x64x5760_9'
%     'fgcir10_2x64x64x2880_1_donotdel'
%     'fgcir10_2x64x64x2880_2_donotdel'
%     'mclibw_2x64x64x2048_01-05-1_1'
%     'mclibw_2x64x64x2048_01-05-1_2'
%     'mring_2x64x64x14400_1'
%     'mring_2x64x64x14400_2'
%     'mringlogsp_2x64x64x28800_1_donotdel'
%     'mringlogsp_2x64x64x28800_2_donotdel'
%     'spn_16x16x10240_1'
%     'spn_32x32x16384_1'
%     'spn_32x32x16384_2'
% ******  TXT files that need but don't have: ******  
% ans = 
%    Empty cell array: 1-by-0


% ******  RAW files not used : ******  
% ans = 
%     'bn32_128x128x32768_4'
%     'bn32_128x128x32768_5'
%     'bn32_128x128x32768_6'
%     'bn32_128x128x32768_loop6'
%     'bn64_128x128x32768_4'
%     'cph_cir_10_2x64x64x5760_10'
%     'cph_cir_10_2x64x64x5760_11'
%     'cph_cir_10_2x64x64x5760_12'
%     'cph_cir_10_2x64x64x5760_13'
%     'cph_cir_10_2x64x64x5760_14'
%     'cph_cir_10_2x64x64x5760_15'
%     'cph_cir_10_2x64x64x5760_16'
%     'cph_cir_10_2x64x64x5760_5'
%     'cph_cir_10_2x64x64x5760_6'
%     'cph_cir_10_2x64x64x5760_7'
%     'cph_cir_10_2x64x64x5760_8'
%     'cph_cir_10_2x64x64x5760_9'
%     'mcli_2x64x64x2048_01-05-1_1'
%     'mcli_2x64x64x2048_01-05-1_2'
%     'mcli_2x64x64x2048_01-05-1_3'
%     'mcli_2x64x64x2048_01-05-1_4'
%     'mcli_2x64x64x2048_01-1-2_1'
%     'mcli_2x64x64x2048_01-1-2_2'
%     'mcli_2x64x64x2048_01-1-2_3'
%     'mcli_2x64x64x2048_01-1-2_4'
%     'mclibw_2x64x64x2048_01-05-1_1'
%     'mclibw_2x64x64x2048_01-05-1_2'
%     'medmovie_m2'
%     'mm_2x64x64x2048_01-05-1@1'
%     'mm_2x64x64x2048_01-05-1@2'
%     'mm_2x64x64x2048_01-05-1@3'
%     'mm_2x64x64x2048_01-05-1@4'
%     'mm_4x32x32x8192_01-05-1@1'
%     'mm_4x32x32x8192_01-05-2@1'
%     'mm_4x32x32x8192_01-1-2@1'
%     'mmc_2x64x64x2048_01-05-1_1'
%     'mmc_2x64x64x2048_01-05-1_2'
%     'mmc_2x64x64x2048_01-05-1_3'
%     'mmc_2x64x64x2048_01-05-1_4'
%     'mmc_2x64x64x2048_01-1-2_1'
%     'mmc_2x64x64x2048_01-1-2_2'
%     'mmc_2x64x64x2048_01-1-2_3'
%     'mmc_2x64x64x2048_01-1-2_4'
%     'mmcw_2x64x64x2048_01-05-1_1'
%     'mmcw_2x64x64x2048_01-05-1_2'
%     'mmcw_2x64x64x2048_01-05-1_3'
%     'mmcw_2x64x64x2048_01-05-1_4'
%     'mring_2x64x64x14400_1'
%     'mring_2x64x64x14400_2'
%     'mspn_16x16x10240_1'
%     'mspn_32x32x16384_1'
%     'mspn_32x32x16384_2'
%     'mv_1_3std'
%     'mv_2_3std'
%     'mv_3_3std'
%     'mv_4_3std'
%     'photo_128x128x16384_1'
%     'photo_128x128x16384_2'
%     'photo_128x128x16384_3'
%     'photo_128x128x24576_1'
%     'smtw_2x64x64x2048_01-05-1_1'
%     'smtw_2x64x64x2048_01-05-1_2'
%     'smtw_2x64x64x2048_01-05-1_3'
%     'smtw_2x64x64x2048_01-05-1_4'
%     'spn_16x16x10240_1'
%     'spn_32x32x16384_1'
%     'spn_32x32x16384_2'
%     'vid_a01_average'
%     'vid_a02_average'
%     'vid_a03_average'
%     'vid_a04_average'
%     'walk1_ieee_128x128x16384'
%     'walk5_ieee_128x128x16384'
%     'walk6_bn32_6_halfs_128x128x32768'
%     'walk6_bn32_6test_128x128x32768_4096'
%     'walk6_ieee_128x128x16384'
%     'walk6_loop_128x128x330x55'
%     'walk6a_bn32_6test_128x128x32768_4096'
% ******  RAW files that need but don't have: ******  
% ans = 
%    Empty cell array: 1-by-0

