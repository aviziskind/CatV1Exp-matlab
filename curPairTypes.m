function [pairType_ids, pairTypes] = curPairTypes(varargin)
%     % set pair type syntax:
%     curPairTypes('Wcc', 'Bcc');
%     
%     % retrieve pair type syntax
%     [pts_i, pts] = curPairTypes;

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curPairTypesFile.mat'];
    allPairTypes = {'Wcc', 'Wrcc', 'Bcc', 'Wcm', 'Wrcm', 'Bcm', 'Wscc', 'Bmm'}; 

    setType = (nargin >= 1) && ~isempty(varargin{1});
    
    if setType   % set current pair types
        n = nargin;
        if iscellstr(varargin{1})
            varargin = varargin{1};
            n = length(varargin);
        end
        allIdxs = zeros(1,n);
        for i = 1:n
            pt_i = varargin{i};
            idx = find(strcmpi(pt_i, allPairTypes), 1);
            if isempty(idx)
                error([pt_i ' is not a valid pairType']);
            end                        
            allIdxs(i) = idx;
        end
            
        allIdxs = sort(allIdxs);
        pairTypes = allPairTypes(allIdxs);
        
        save(filename, 'pairTypes');
        
    else 
        if exist(filename, 'file')
            load(filename);
        else
            error('Pair Type file does not exist');
        end
        
        pairType_ids = cellfun(@(pt) find(strcmp(pt, allPairTypes), 1), pairTypes); %#ok<NODEF>
        if (nargin>=1) && isempty(varargin{1})
            pairType_ids = pairTypes;
        end        

        
    end

end
