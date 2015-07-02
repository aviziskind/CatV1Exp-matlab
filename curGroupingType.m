function [groupingType_id_out, groupingType_out] = curGroupingType(s)
    
    persistent groupingType

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curGroupingTypeFile.mat'];
    allGroupingTypes = {'clusters', 'clustersPruned', 'cells'};
    
    setType = (nargin == 1) && ~isempty(s) && (isnumeric(s) || strncmp(s, 'c', 1));
    
    if setType   % set current grouping type 

        if isnumeric(s)
            s = switchh(s, [1, 2, 3 4], {'clusters', 'clustersPruned', 'cells', 'cells_onlyIC', 'cells_onlyPrunedIC'});
        end
        
        switch lower(s)
            case {'clusters', 'cl', 'c'}
                groupingType = 'clusters';
            case {'clusterspruned', 'cp'}
                groupingType = 'clustersPruned';
            case {'cells'}
                groupingType = 'cells';
            case {'cells_onlyic', 'ic'}
                groupingType = 'cells_onlyIC';
            case {'cells_onlyprunedic', 'pic'}
                groupingType = 'cells_onlyPrunedIC';
            otherwise
                error('invalid grouping type: enter "clusters"("c"), "clustersPruned" ("cp")  "cells" ("ce"), cells_onlyIC ("ci"), or cells_onlyPrunedIC ("pic")');
        end
        save(filename, 'groupingType');
        
    else
        if isempty(groupingType)
            if exist(filename, 'file')
                load(filename);
            else
                error('Grouping Type file does not exist');
            end
        end
        groupingType_out = groupingType;
        groupingType_id_out = find(strcmp(groupingType, allGroupingTypes)); 
        
        if (nargin==1) && isempty(s)
            groupingType_id_out = groupingType;
%         elseif (nargin==1) && ~isempty(s)
%             error('check');
%             groupingType = strrep(groupingType, 'cells', 'cell');
%             groupingType = strrep(groupingType, 'clusters', 'cluster');
%             groupingType_id = groupingType;
        end
    end

end
