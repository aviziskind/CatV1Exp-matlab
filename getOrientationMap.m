function [Z_map, x_map, y_map] = getOrientationMap(nModes, n_wav, lambda, nPointsPerWavelength, rand_seed)

    persistent allMaps nDone;
    margin = 0;

    redo = 0;
    ori_map_file = [CatV1Path 'OriMaps.mat'];
    nBeforeSave = 5;
    
    
    get_field_name = @(nMod, n_wav, lambda, nPoints, rand_seed) ...
        fix_fld( sprintf('Map__Mod%d__Wav%.1f__lam%.1f__nPPw%d_seed%d', nMod, n_wav, lambda, nPoints, rand_seed) );

    if isempty(allMaps)
        if ~exist(ori_map_file, 'file') || redo
            allMaps = struct;        
        else
            allMaps = load(ori_map_file);
        end
        nDone = 0;
    end
    
    fld_name = get_field_name(nModes, n_wav, lambda, nPointsPerWavelength, rand_seed);
    if isfield(allMaps, fld_name)
        Z_map = allMaps.(fld_name).map;
        x_map = allMaps.(fld_name).x; 
        y_map = x_map;
        return;
    else
        rng(rand_seed);
        [Z_map, x_map, y_map] = generateOrientationMap(nModes, n_wav, lambda, margin, nPointsPerWavelength);
        assert(isequal(x_map, y_map));
        map_struct = struct('map', Z_map, 'x', x_map);
        allMaps.(fld_name) = map_struct;
        nDone = nDone + 1;
        if (nDone >= nBeforeSave)
            save(ori_map_file, '-struct', 'allMaps', '-v6')
            nDone = 0;
        end
        
    end

end

function s = fix_fld(s)
    s = strrep(s, '.', '_');
end
