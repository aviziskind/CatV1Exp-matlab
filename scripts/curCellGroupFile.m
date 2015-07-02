function fname = curCellGroupFile(gratType, fet_str, grpType, matchDB)

    if (nargin < 1) || isempty(gratType)
        gratType = curGratingType('');
    end
    if (nargin < 2) || isempty(fet_str)
        [~,~,fet_str] = curSortingFeatures('');
    end
    if (nargin < 3) || isempty(grpType)
        grpType = curGroupingType('');
        grpType = grpType(1:end-1);
    end
    if (nargin < 4) || isempty(matchDB)
        matchDB = curMatchDB('');
    end
    
    if matchDB        
        suffix = '_DB';
    else
        suffix = ['_' fet_str];
    end
    
    switch gratType(1)
        case 'f', grating_str = 'movie_fg';
        case 'd', grating_str = 'grating_dSf';
    end
        
    fname = [grpType 'Groups_' grating_str suffix];
    
end
