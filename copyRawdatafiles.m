function copyRawdatafiles(Gids)

    if nargin < 1  % copy datafiles for *all* groups
        
        gratingType = curGratingType;

        if gratingType == 1
            load('cellsGroups_movie_fg_DB.mat');
            allGroups = movieGroups_fg;

        elseif gratingType == 2
            load('cellsGroups_grating_dSf_DB');
            allGroups = gratingGroups_dSf;

            load('cellsGroups_grating_dOr_DB');
            allGroups = [allGroups; gratingGroups_dSf];            
            
        end
        dfs = [allGroups.dataFileInfo];    
        Gids = [allGroups.Gid];    
        filenames = {dfs.dataFileName};    
        
        nFiles = length(filenames);
    else
        nFiles = length(Gids);
        
        for i = 1:nFiles
            sd = siteDataFor(Gids(i));
            dfs(i) = sd.dataFileInfo; %#ok<AGROW>
        end        
        filenames = {dfs.dataFileName};    
    end
    
    doBar = nFiles > 10;
    if doBar
        progressBar('init-', nFiles);
    end
    for i = 1:nFiles
        if doBar
            progressBar;
        end
        dfname = filenames{i};
        catName = strtok(dfname, '_');        
        dest_dir = ['C:\ExperimentDB\RawData\' catName '\'];
        if ~doBar
            fprintf('File %d : %s (for Gid = %d) [%d MB]... ', i, dfname, Gids(i), round( dfs(i).filesize/(1024^2))  )            ;
        end

        if ~exist(dest_dir, 'file')
            mkdir(dest_dir)
        end
        if ~exist([dest_dir dfname], 'file')
            source_dir = rawDataDir(dfname);
            copyfile([source_dir dfname], [dest_dir dfname]);    
            fprintf('copied \n')
        else
            fprintf('(already present) \n')
        end
        
        
    end


end

%{
    dir1 = 'C:\ExperimentDB\RawData\FlashedGratings\';
    dir2_base = 'C:\ExperimentDB\RawData\'
    s = dir(dir1)
    s = s(3:end);
    
    
    for i = 1:length(s)        
        fprintf('*');
        nm = s(i).name;
        catName = strtok(nm, '_');
        dir2 = [dir2_base catName '\'];
        if ~exist(dir2, 'dir')
            mkdir(dir2);
        end
        n1 = [dir1 nm];
        n2 = [dir2 nm];
        assert(exist(n1, 'file')>0);
        assert(~exist(n2, 'file'));
        try
            movefile(n1, n2);
        catch
        end
    end

%}