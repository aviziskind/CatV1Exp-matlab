function [MID, gparams] = mid_findMostInfoDim(Gid, cellId, opt)

%     baseDir = 'C:\Users\Avi\Documents\MATLAB\columbia\miller\gabor\MID\';
%     baseDir = 'C:\Users\Avi\Documents\MATLAB\columbia\miller\MID\';
    redo = 1;
    downSampFactor_default = 1;
    frameMode_default = '1rep';
    timeWindow_default = 'best';
    trialMode_default = 'all';
    calcGaborParams_default = 0 || (nargout > 1);
    responseType_default = curResponseType(''); % 'raw';
    
    if (nargin < 3) || isempty(opt)
        opt = struct;
    end
    if ~isfield(opt, 'downSampFactor'),  opt.downSampFactor  = downSampFactor_default;  end
    if ~isfield(opt, 'frameMode'),       opt.frameMode       = frameMode_default;       end
    if ~isfield(opt, 'timeWindow'),      opt.timeWindow      = timeWindow_default;      end
    if ~isfield(opt, 'trialMode'),       opt.trialMode       = trialMode_default;       end        
    if ~isfield(opt, 'calcGaborParams'), opt.calcGaborParams = calcGaborParams_default; end        
    if ~isfield(opt, 'responseType'),    opt.responseType    = responseType_default;    end        
    

    
%     cellName = sprintf('Group_%d__cell_%d', Gid, cellId);
%     cellName_detailed = sprintf('%s_d%d_%s', cellName, opt.downSmpFactor, opt.frameMode);
%     baseDir = ['C:\ExperimentDB\MID\'];
%     
%     mid_fileName = sprintf('%sCells\\%s\\MID_%s.mat', ...
%         baseDir, cellName, cellName_detailed);

    [mid_fileName, cellName] = getName('MID_file', Gid, cellId, opt.downSampFactor, opt.frameMode, opt.timeWindow, opt.trialMode, opt.responseType );
    
    if ~exist(mid_fileName, 'file') || redo
        [MID, v_MID, dimInfo, mid_resultFile, t_calc] = calculateMostInfoDim(Gid, cellId, opt, cellName); %#ok<NASGU>
        save(mid_fileName, 'MID', 'v_MID', 'dimInfo', 'mid_resultFile', 'opt', 't_calc');
        
    else
        S = load(mid_fileName);
        MID = S.MID;    
    end
    
    if opt.calcGaborParams
        gparams = mid_getCellGaborParams(Gid, cellId, opt.timeWindow, opt.trialMode);
    end
    
end

function [MID, v_MID, dimInfo, MIDfilename, t_calc] = calculateMostInfoDim(Gid, cellId, opt, cellName)
                        
    redoStimFile = 0;
    redoSpikeFile = 0;    
    chkFiles = 0;
    
    baseDir = [experimentDBPath 'MID' filesep];
    funcName_base = 'do_MaxInfoDim_Nvec';
    frameMode = opt.frameMode;
    trialMode = opt.trialMode;
    timeWindow = opt.timeWindow;
    responseType = opt.responseType;

    sd = siteDataFor(Gid);
    stimInfo = sd.stimulusInfo;
    nrows = stimInfo.nrows;
    ncols = stimInfo.ncols;

    x0 = 1; y0 = 1;
    dx = nrows; dy = ncols;
    cx = opt.downSampFactor;
    cy = cx;
    
    nLags = 1;
    nBins = 15;

%     onNYUserver = ~isempty(strfind(getenv('hostname'), 'nyu.edu'));
    
%     if onNYUserver
%         dataDir = ['~/MID/Cells/' cellName filesep ];
%     else
        dataDir = [baseDir 'Cells' filesep cellName filesep ];
%     end
        
    if ~exist(dataDir, 'dir')
        mkdir(dataDir);
    end

    stimFileName = mid_getMovieStimFileName(Gid, frameMode);                         
    if ~exist(stimFileName, 'file') || redoStimFile
        mid_generateStimulusFile(Gid, frameMode);
    end

    spikeFileName = mid_getSpikeFileName(Gid, cellId, timeWindow, trialMode, frameMode, responseType);
    if ~exist(spikeFileName, 'file') || redoSpikeFile
%         if ischar(timeWindow) &&  strcmp(timeWindow, 'best')
%             if exist('celldata', 'var')
%                 if isfield(celldata, 'PSTH')
%                     PSTH = celldata.PSTH;
%                 elseif isfield(celldata, 'timeWindow')
%                     PSTH = celldata;
%                 end
%                 timeWindow_ms = PSTH.timeWindow_ms;
%                 windowProfile_ms = PSTH.windowProfile;
%             else
%                 error('need celldata or pre-existing spikefile');
%                 timeWindow_ms = [];
%                 windowProfile_ms = [];
%             end
%         elseif isnumeric(timeWindow)
%             timeWindow_ms = timeWindow;
%             windowProfile_ms = [];
%         end
%         mid_generateSpikeFile(Gid, cellId, timeWindow_ms, windowProfile_ms, frameMode, timeWindow_mode);
        mid_generateSpikeFile(Gid, cellId, timeWindow, trialMode, frameMode, responseType);
        
    end

    if chkFiles
       % check stimulus file
       %%
       frameDims = [64, 64];
       verify_stim_file(stimFileName, frameDims(1), frameDims(2));

       % check spike file
       verify_spike_file(spikeFileName);
       3;
    end

    fileDetails = dir(stimFileName);
    fileSizeBytes = fileDetails.bytes;
    nFrames = fileSizeBytes / (nrows*ncols);                
    arch = switchh(computer, {'PCWIN', 'PCWIN64', 'GLNXA64'}, {'win32', 'win64', 'linux64'});
    funcName = [baseDir 'bin' filesep arch filesep funcName_base];
    
    if strcmp(getenv('hostname'), 'XPS') && strcmp(arch, 'linux64')
        funcName = ['~/Local/MID/' funcName_base];
    elseif onNYUserver
        funcName = ['~/MID/bin/' funcName_base];
    end
    
        
    findDimCmd = sprintf('%s %s %s %s %s %d %d %d %d %d %d %d %d %d %d %d', funcName, cellName, ...
        stimFileName, spikeFileName, dataDir, nrows, ncols, nFrames, x0, y0, dx, dy, cx, cy, nLags, nBins);
    
    dimInfo = [round(x0+dx-1)/cx, round(y0+dy-1)/cy, nLags, cx];
    name_suffix_str = sprintf('%ux%ux%u_%u',dimInfo);    
    
    redo = false;
    nJacks = 4;
    
%     nCPUs = 2;
%     started = [];
%     completedFileNames = cellfun(@(i) [dataDir 'done' num2str(i)], num2cell(1:nJacks), 'Un', false );        
    
    if redo
        for jack_i = 1:nJacks
            if exist(completedFileNames{jack_i}, 'file')
                system(['del ' completedFileNames{jack_i}]);
            end
        end
    end
    
    
    tStart_all = tic;    
    t_calc = zeros(1,nJacks);
    
    fprintf('\n**Starting MID search at %s\n', datestr(clock))
    for jack_i = 1:nJacks        
        jack_i_name = sprintf('%sTest_v_%s_%s_jack%d.dat', dataDir, cellName, name_suffix_str, jack_i);
        jack_i_name_clock = sprintf('%sClock_data_%s_%s_jack%d.mat', dataDir, cellName, name_suffix_str, jack_i);
        % the 'clock data' file stores how much time it took to calculate
        % the MID. More importantly, its existence indicates that the MID
        % calculation completed successfully.

        
        completedAlready = exist(jack_i_name, 'file') && exist(jack_i_name_clock, 'file');
        if completedAlready && ~redo            
           
            S_i = load(jack_i_name_clock);
            t_calc(jack_i) = S_i.timeElapsed_sec;
            fprintf('Already completed jack # %d previously -- (took %s)\n', jack_i, sec2hms(t_calc(jack_i)));
            
            
        else

%         findDimCmdFull = cmdstart( findDimCmd, 'wtitle', funcName, 'priority', 'belownormal');
%         dos(findDimCmdFull);
            date_start = now;
            tStart = tic;      
            fprintf('Starting jack %d at %s\n', jack_i, datestr(clock));

            sys_cmd = sprintf('%s %d', findDimCmd, jack_i);
            [status] = system(sys_cmd);
            
            if status ~= 0
                if status == 130
                    error('Ctrl-C: Operation terminated by user during during MID search')
                else
                    error('Unknown error: %d', status);
                end
            end

            date_end = now;
            t_calc(jack_i) = toc(tStart);                
            fprintf('(Took %s to do jack %d)\n', sec2hms(t_calc(jack_i)), jack_i);        

            S_i = struct('t_start', date_start, 't_end', date_end, 'timeElapsed_sec', t_calc(jack_i)); %#ok<NASGU>
            save(jack_i_name_clock, '-struct', 'S_i');
        end
        
    end        
%     t_tot = toc(tStart_all);
    fprintf('\n(Completed all %d jacks [total time : %s])\n', nJacks, sec2hms(sum(t_calc)));
    
    
%     figure(56);
%     STAfilename = [dataDir 'STA_' cellName];
%     v_STA = read_vec_pxpxt_file (STAfilename,dimInfo, nBins);
%     figure(57);
    
    
    MIDfilename  = [dataDir 'Test_v_' cellName];    
    v_MID = read_vec_pxpxt_file(MIDfilename,dimInfo);
        
    MID = mean(v_MID{1}, 3);
    
    
%         MID = mean(v_MID{1}, 3);
%         figure(10); imagesc(MID);
%         3;
%         figure(102);
%         for i = 1:4
%             subplot(2,3,i+2); imagesc(v_MID{1}(:,:,i)); axis square;
%         end
%         subplot(2,3,1); imagesc(mean(v_MID{1}, 3)); axis square;
%     end

%     dimInfo = [dx/cx, dy/cx, nLags, cx];
%     Test_v_Group4470_Cell0_16x16x1_4_jack1.dat
%     Test_v_Group4470_Cell0_64x64x1_4_jack1.dat


end


%{
    code for doing jacks in parallel using multiple CPUs
        jack_i = 1;
        nCompleted = @() nnz( cellfun(@(s) exist(s, 'file'),  completedFileNames) );
        while (nCompleted() < nJacks) && (jack_i <= nJacks)
            nCurrentlyDoing = length(started) - nCompleted();
            if ( nCurrentlyDoing < nCPUs)
                started = [started, jack_i]; 
                disp(['Starting # ' num2str(jack_i)]);
                findDimCmdFull = cmdstart( [findDimCmd num2str(jack_i)], 'wtitle', funcName, 'priority', 'belownormal');
                eval(['! ' findDimCmdFull ' & createfile ' completedFileNames{jack_i} ' & exit &']);
                jack_i = jack_i + 1;
            end
            pause(10);
        end
%}