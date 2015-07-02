function siteData = siteDataFor(idType, idVal, includeCellData)
    persistent matchDB groupingType sortingFeatures
    persistent movieSites gratingSites movieGroups gratingGroups 
    persistent Gid_stimTypes Did_stimTypes  Gid_stimIdx Did_stimIdx
    persistent hcSites hcGroups
    
%     Gids = dbLookup('Gid',  idtype, idval); Gid = Gids(1);
%     Did = dbLookup('Did',  idtype, idval);    
%     stimulusType = getStimulusTypeForDid(Did);    
    useOnlyRelevantGratingStimuli = 1;

    if isnumeric(idType) && ((nargin < 2) || isempty(idVal))
        idVal = idType;
        idType = 'Gid';        
    end

    includeCellData = exist('includeCellData', 'var') && ~isempty(includeCellData);    
    
    reloadGroupData = includeCellData && ...
                    (~isequal(matchDB, curMatchDB) || ~isequal(groupingType, curGroupingType) || ~isequal(sortingFeatures, curSortingFeatures('')) );                
                
    if reloadGroupData
        movieGroups = [];
        gratingGroups = [];        
        
        matchDB = curMatchDB;
        groupingType = curGroupingType;
        sortingFeatures = curSortingFeatures('');
    end

    loadMovieGroups   = (includeCellData && isempty(movieGroups)) || (~includeCellData && isempty(movieSites));
    loadGratingGroups = (includeCellData && isempty(gratingGroups)) || (~includeCellData && isempty(gratingSites));
    loadHippocampusGroups = false; %isempty(hcSites);

    
    if loadMovieGroups || loadGratingGroups
        if includeCellData
            matchDB_str = curMatchDB('');
            if matchDB
                groupingType_str = 'cells';
                fetStr = '';
            else
                groupingType_str = curGroupingType('');
                fetStr = [curSortingFeatures('') '_'];
            end
        else
            matchDB_str = '_DB'; % just get the original database site info. don't need the cell/cluster info.
            groupingType_str = 'cells';
            fetStr = '';
        end
    end
    
    
    if loadMovieGroups
        
        if useOnlyRelevantGratingStimuli
            flashed_file_name = [CatV1Path groupingType_str 'Groups_' fetStr 'movie_fg' matchDB_str '.mat'];
        else
            flashed_file_name = [CatV1Path groupingType_str 'Groups_' fetStr 'movie_all' matchDB_str '.mat'];
        end
        S = load(flashed_file_name); fn = fieldnames(S);
        movieGroups_tmp = S.(fn{1});
        movieGids = [movieGroups_tmp.Gid];
        movieDids = [movieGroups_tmp.Did];
        
        Gid_stimTypes(movieGids) = 'm';        
        Did_stimTypes(movieDids) = 'm';                        
        Gid_stimIdx(movieGids) = 1:length(movieGids);
        Did_stimIdx(movieDids) = 1:length(movieDids);
        
        if isempty(movieSites)
            movieSites = arrayfun(@(s) rmfield(s, {'cellIds', 'nSpikes'}), movieGroups_tmp);            
        end
        if includeCellData
            movieGroups = movieGroups_tmp; % store the movie groups, with the cell data.
        end        
        
    end
    
    if loadGratingGroups
        
        if useOnlyRelevantGratingStimuli
            drifting_file_names = {[CatV1Path groupingType_str 'Groups_' fetStr 'grating_dSf' matchDB_str '.mat'], ...
                                   [CatV1Path groupingType_str 'Groups_' fetStr 'grating_dOr' matchDB_str '.mat']};
        else
            drifting_file_names = {[CatV1Path groupingType_str 'Groups_' fetStr 'grating_all' matchDB_str '.mat']};
        end
        
        gratingGroups_tmp = [];
        for di = 1:length(drifting_file_names)
            S = load(drifting_file_names{di}); fn = fieldnames(S);  
            if isempty(gratingGroups_tmp)
                gratingGroups_tmp = S.(fn{1});
            else
                gratingGroups_tmp = [gratingGroups_tmp; S.(fn{1})]; %#ok<AGROW>
            end
        end                            
            
        gratingGids = [gratingGroups_tmp.Gid];
        gratingDids = [gratingGroups_tmp.Did];

        Gid_stimTypes(gratingGids) = 'g';
        Did_stimTypes(gratingDids) = 'g';     
        Gid_stimIdx(gratingGids) = 1:length(gratingGids);
        Did_stimIdx(gratingDids) = 1:length(gratingDids);        
        
        if isempty(gratingSites)
            gratingSites = arrayfun(@(s) rmfield(s, {'cellIds', 'nSpikes'}), gratingGroups_tmp);            
        end
        if includeCellData
            gratingGroups = gratingGroups_tmp; % store the grating groups, with the cell data.
        end        
    end

    if loadHippocampusGroups 
%         hcSites hcGroups
        hc_file_name = [CatV1Path 'RatHC' filesep 'allHcDatafiles.mat'];
        S = load(hc_file_name);
        hcGroups = S.allDataFiles;
        hcSites = S.allDataFiles;
        hcGids = [hcGroups.Gid];
        hcDids = [hcGroups.Did];

        Gid_stimTypes(hcGids) = 'h';  
        Did_stimTypes(hcDids) = 'h';   
        Gid_stimIdx(hcGids) = 1:length(hcGids); 
        Did_stimIdx(hcDids) = 1:length(hcDids);    
    end        
    
    if isa(Gid_stimTypes, 'double')
       Gid_stimTypes = uint8(Gid_stimTypes);
       Did_stimTypes = uint8(Did_stimTypes);
        
       Gid_stimIdx = uint16(Gid_stimIdx);
       Did_stimIdx = uint16(Did_stimIdx);
    end    
    
    
    switch idType
        case 'Gid', 
            idx = Gid_stimIdx(idVal);
            if any(idx == 0)
                error('No group with Gid == %d', idVal(idx == 0));
            end
            stimType = char( Gid_stimTypes(idVal) );
            
        case 'Did',
            idx = Did_stimIdx(idVal);
            if any(idx == 0)
                error('No group with Did == %d', idVal(idx == 0));
            end
            stimType = char( Did_stimTypes(idVal) );
            
    end
    
    for ii = 1:length(idVal)
        
        if includeCellData
            switch stimType(ii)
                case 'g', 
                    siteData(ii) = gratingGroups(idx(ii)); %#ok<AGROW>
                case 'm', 
                    siteData(ii) = movieGroups(idx(ii)); %#ok<AGROW>
                case 'h', 
                    siteData(ii) = hcGroups(idx(ii)); %#ok<AGROW>
            end

        else
            switch stimType(ii)
                case 'g', 
                    siteData(ii) = gratingSites(idx(ii)); %#ok<AGROW>
                case 'm', 
                    siteData(ii) = movieSites(idx(ii)); %#ok<AGROW>
                case 'h', 
                    siteData(ii) = hcSites(idx(ii)); %#ok<AGROW>
            end        

        end
    
        if ii == 1
            siteData(2:length(idVal)) = repmat(blankStruct(siteData), 1, length(idVal)-1);
        end    
    end
        
end