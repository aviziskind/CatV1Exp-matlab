function catIds = getCatId(animalIds)
    
    persistent allCatIds;
    if isempty(allCatIds)
    
        allCatIds([1  2  4  5  9  10 11 12 14 16 18 19 26 32 33 34 35 36 37 41 42 43 44 47 48 50 51 52 53 54 ]) = ...
                  [12 13 14 15 16 17 18 19 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 ];
    end
    catIds = allCatIds(animalIds);

end