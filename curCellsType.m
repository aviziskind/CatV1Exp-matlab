function s = curCellsType
    if curMatchDB
        s = '_DB';
    else
        s = ['_' curSortingFeatures('')];
    end
end