function [name, name2] = getName(purpose, varargin)
%   Example syntax:
%   To retrieve variable names:
%         name = getName('celldata', Gid, cellId)
%
%   To retrieve path names:
%         name = getName('movieBasePath');
%         name = getName('MatlabDB_Path');
%
%   To retrieve stimulus file names:
%         name = getName('movieStimFile', Gid)
%         name = getName('noiseStimFile', idType, idval)  % (ie. 'Gid' or 'Did')
%             or getName('noiseStimFile', randSeed, nGrad)  % (ie. 'Gid' or 'Did')
%         name = getName('mseqStimFile', idType, idval)
%             or getName('mseqStimFile', tapRegister, nMseqBits)
%         name = getName('stimFile', idType, idval)  % (ie. 'Gid' or 'Did')
%
%   To retrieve response/data file names:
%         name = getName('spikefile', Gid, cellId)
%
    
    expDB_path = experimentDBPath;
    
    switch purpose
        % VARIABLES
        case 'celldata'
            [Gid, cellId] = elements(varargin);
            name = sprintf('celldata__Group_%d__Cell_%d', Gid, cellId);

        case 'siteTest'
            [Gid] = elements(varargin);
            name = ['siteTest__Group_' num2str(Gid) ];
            
        % PATHS
        case 'movieBasePath'   % this is the root directory where all movie stimuli are stored
            name = [expDB_path 'Stimuli' filesep 'Movies' filesep];

        case 'fgMoviePath',
            name = [expDB_path 'Stimuli' filesep 'Movies' filesep 'Flashed_Gratings' filesep];
            
        case 'fgMovie_movieId'
            movieId = varargin{1};            
            movie_path = [expDB_path 'Stimuli' filesep 'Movies' filesep 'Flashed_Gratings' filesep];
            movieOSP_name_mat = sprintf('FgMovie_%d.mat', movieId);
            name = [movie_path movieOSP_name_mat];            
            
        case 'MatlabDB_path'   % this where all DB .m files are stored            
            name = [CatV1Path 'MatLabDB_avi' filesep];

            
        case 'dfFiltExt',            
            if (length(varargin) == 1) && ischar(varargin{1})
                dfname = varargin{1};                
                idx_f = strfind(dfname, 'F');    
                name = dfname(idx_f(1):end);
            else

                if length(varargin) == 3
                    [highPass_freq_Hz, filterOrder, filterName] = deal(varargin{:});
                elseif length(varargin) == 1
                    opt = varargin{1};
                    doRawFile = isempty(opt);
                    if ~doRawFile
                        [highPass_freq_Hz, filterOrder, filterName] = deal(...
                            opt.highPass_freq_Hz, opt.filterOrder, opt.filterName);
                    end
                end              
                if ~doRawFile
                    name = sprintf('F_%d_%d_%s', highPass_freq_Hz, filterOrder, filterName);                        
                else
                    name = '';
                end
            end
            
        % FILE NAMES: STIMULI
        case 'compiledMoviePath'
%             [Gid, movieFileName] = elements(varargin);
            
            name = [expDB_path 'MID' filesep 'CompiledMovies' filesep];
            
        case 'compiledSpikesPath'
            name = [expDB_path 'MID' filesep 'CompiledSpikes' filesep];
            
            
        case 'noiseStimFile'            
            if ischar(varargin{1})
                [idtype, idval] = elements(varargin);
                Did = dbLookup('Did',  idtype, idval);
                [randSeed, nGradations] = getFieldsFromDatabaseTable(dbOpenExpDb, {'LNG_RAND_SEED', 'LNG_N_GRADATIONS'}, 'TBL_NOISE_PRES', {'DATAFILE_ID', Did}, [], 1, 1);
            else
                [randSeed, nGradations] = elements(varargin);
            end
            path = [expDB_path 'Stimuli' filsep 'Noise' filesep];
            filename = ['noise' num2str(randSeed) '_' num2str(nGradations)];            
            ext = '.raw';
           
            name = [path filename ext]; 

        case 'mseqStimFile'            
            if ischar(varargin{1})
                [idtype, idval] = elements(varargin);
                Did = dbLookup('Did',  idtype, idval);
                [tapRegister, nMseqBits] = getFieldsFromDatabaseTable(dbOpenExpDb, {'LNG_TAP_REGISTER', 'LNG_MSEQ_BITS'}, 'TBL_NOISE_PRES', {'DATAFILE_ID', Did}, [], 1, 1);
            else
                [tapRegister, nMseqBits] = elements(varargin);
            end            
            path = [expDB_path 'Stimuli' filesep 'Mseq' filesep];
            filename = ['mseq' num2str(tapRegister) '_' num2str(nMseqBits)];            
            ext = '.raw';
            name = [path filename ext]; 

        case 'MID_file'
            [Gid, cellId, downSmpFactor, frameMode, timeWindow, trialMode, responseType] = deal(varargin{:});
            if strcmp(frameMode, 'all')
                gt = getGratingStimType(Gid);
                frameMode = sprintf('%drep', gt.nTrials);
            end             
            
            if strcmp(timeWindow, 'best') || isempty(timeWindow)
                timeWindow_str = '';
            elseif isnumeric(timeWindow) && length(timeWindow) == 2        
                timeWindow_str = sprintf('__%d-%d', timeWindow(1), timeWindow(2));
            end
            
            switch trialMode
                case 'all', trialMode_str = '';                
%                 case {'odd/even', 'oe'}, trialMode_str = '__oe';
                 case {'odd', 'even'}, trialMode_str = ['__' trialMode];
                otherwise, error('Unknown trial mode')
            end
                        
            switch responseType
                case 'raw', responseType_str = '';
                case 'gainCorrected', responseType_str = '_GC';
            end
            
            cellName = sprintf('Group_%d__cell_%d_%s%s%s', Gid, cellId, frameMode, timeWindow_str, trialMode_str);
            cellName_detailed = sprintf('%s_d%d', cellName, downSmpFactor);
            
            cellName          = [cellName          responseType_str];
            cellName_detailed = [cellName_detailed responseType_str];
            
            baseDir = [expDB_path 'MID' filesep];
    
%             name = sprintf('%sCells\\%s\\MID_%s.mat', ...
%                 baseDir, cellName, cellName_detailed);
            name = [baseDir 'Cells' filesep 'MID_' cellName_detailed '.mat'];

            name2 = cellName;            
        % FILE NAMES: DATA
%         case 'spikeFile'
%             [Gid, cellId, frameMode] = elements(varargin);            
%             path = [expDB_path '\MID\CompiledSpikes\'];
%             filename = sprintf('Group_%d__cell_%d__%s', Gid, cellId, frameMode);
%             ext = '.spk';            
%             name = [path filename ext];
                    
    end
    
end


