function [featureNames_out, featureNs_out, feature_str_out, featureNamesAndN_out] = curPruningFeatures(varargin)

    persistent allFeatureNames allFeatureDefaultNs allFeatureChannelFactor
    persistent featureNames featureNs feature_str featureNamesAndN
%     % set feature type syntax:
%     curfeatures('Wcc', 'Bcc');
%     
%     % retrieve feature syntax
%     [pts_i, pts] = curfeatures;

    nPCs  = 2; % default value;
    nGLFs = 2;
    nCh = 4; % nChannels;

    channelActionDflt = 'c'; % == concatenated
    normChannelDflt = 'u'; % == false
    whitenFeaturesDflt = 'w'; % == true

    coefFeatures = {'PCA', 'GLF'};
    coefDfltNs = [nPCs, nGLFs];
    
    channelActionOptions = {'s', 'c'}; %separate, concatenated    
    channelActionFactors = [nCh,  1];
    
    normOptions = {'u', 'n'}; %unnormalized, vs. normalized 
    whitenOptions = {'r', 'w'}; %cross-channel whitened, vs. raw.
    
    
    filename = [CatV1Path 'MatLabDB_avi' filesep 'curPruningFeaturesFile.mat'];
    if isempty(allFeatureNames)
        allFeatureNames         = {'Neg', 'Pos', 'Ptp', 'Chd', 'Egy'}; %, 'PCAs', 'PCAc', 'GLFs', 'GLFc'};
        allFeatureDefaultNs     = [ 1,     1,     1,     1,     1   ]; %,  nPCs,   nPCs,  nGLFs,   nGLFs ];
        allFeatureChannelFactor = [ nCh,   nCh,   nCh,   nCh,   nCh ]; %,   nCh,    1,     nCh,     1,    ];

        nCoefs = length(coefFeatures);
        nChan = length(channelActionOptions);    
        nNorm = length(normOptions);
        nWhiten = length(whitenOptions);

        allCoefNames = cell(nCoefs, nChan, nWhiten, nNorm);
        allCoefNs = zeros(nCoefs, nChan, nWhiten, nNorm);
        allCoefChFactors = zeros(nCoefs, nChan, nWhiten, nNorm);
        for i = 1:nCoefs
            allCoefNs(i,:,:,:) = coefDfltNs(i);        
            for j = 1:nChan
                allCoefChFactors(i,j,:,:) = channelActionFactors(j);
                for k = 1:nNorm
                    for l = 1:nWhiten                    
                        allCoefNames{i,j,k,l} = [coefFeatures{i} channelActionOptions{j}, normOptions{k} whitenOptions{l} ];
                    end
                end
            end
        end    
        allFeatureNames         = [allFeatureNames, allCoefNames(:)'];
        allFeatureDefaultNs     = [allFeatureDefaultNs, allCoefNs(:)'];
        allFeatureChannelFactor = [allFeatureChannelFactor, allCoefChFactors(:)'];                
    end
        
    setType = (nargin >= 1) && ~isempty(varargin{1});
    
    if setType   % set current feature
        newFeatures = varargin;
        if length(newFeatures) == 1 && iscellstr(newFeatures{1})
            newFeatures = newFeatures{1};
        end
        nFeatures = length(newFeatures);
        
        allIdxs = zeros(1,nFeatures);
        featureNs  = zeros(2,nFeatures);
        isCoefType = false(1,nFeatures);
        for fet_i = 1:nFeatures
            
            fet_i_name = newFeatures{fet_i};            
            [fet_str, nCoefs] = extractDigits(fet_i_name);
            newFeatures{fet_i} = fet_str;
            
            validPrefix = ~isempty( find(strncmpiC(fet_str, allFeatureNames, 3), 1));
            isCoefType(fet_i) = any(strncmpiC(fet_str, coefFeatures, 3));
            if ~validPrefix
                error([fet_i_name ' is not a valid pruning feature']);
            end                        
            
            fet_idx = find(strncmpiC(fet_str, allFeatureNames, 6), 1); % find exact match
            if isempty(fet_idx) 
                if isCoefType(fet_i)  % didn't have exact match - use defaults.
                    if (length(fet_str) < 4)
                        fet_str(4) = channelActionDflt;
                    elseif ~any(strcmp(fet_str(4), channelActionOptions))
                        error('featureName(4) must be "s" or "c"');
                    end

                    if (length(fet_str) < 5)
                        fet_str(5) = normChannelDflt;
                    elseif ~any(strcmp(fet_str(5), normOptions))
                        error('featureName(5) must be "u" or "n" (not %s)', fet_str(5));
                    end
                    
                    if (length(fet_str) < 6)
                        fet_str(6) = whitenFeaturesDflt;
                    elseif ~any(strcmp(fet_str(6), whitenOptions))
                        error('featureName(5) must be "u" or "n" (not %s)', fet_str(6));
                    end                    
                    
                    fet_idx = find(strncmpiC(fet_str, allFeatureNames, 6), 1); % find exact match
                    assert(~isempty(fet_idx));                    
                else                
                    error([fet_i_name ' is not a valid pruning feature']);
                end
            end                        
                    
            fet_str = allFeatureNames{fet_idx}; % make the same upper/lower case as original 
            
            if isCoefType(fet_i) 
                if isempty(nCoefs)
                    nCoefs = switchh(fet_str(1:3), {'PCA', 'GLF'}, [nPCs, nGLFs]);                    
                end
            else 
                nCoefs = allFeatureDefaultNs(fet_idx); 
            end 
            featureNs(:,fet_i) = [nCoefs; allFeatureChannelFactor(fet_idx)];                    
            
            allIdxs(fet_i) = fet_idx;
        end
            
        if length(unique(allIdxs)) < length(allIdxs)
            [uIdx, uCnt] = uniqueCount(allIdxs);
            error('Duplicate feature "%s" entered', allFeatureNames{uIdx(find(uCnt>1,1))});
        end
        
        sort_idx = ord(allIdxs);
        featureNames = allFeatureNames(allIdxs(sort_idx)); 
        featureNs = featureNs(:,sort_idx);
        isCoefType = isCoefType(sort_idx);
                
        featureNamesAndN = arrayfun(@(s,n,tf) iff(tf, [s{1}, num2str(n)], s{1}), featureNames, featureNs(1,:), isCoefType, 'un', 0);
        
        feature_str = cellstr2csslist(featureNamesAndN, '_');
                
        save(filename, 'featureNames', 'featureNs', 'feature_str', 'featureNamesAndN');
        
    else
        if isempty(featureNames)            
            if exist(filename, 'file')
                S = load(filename);
                featureNames = S.featureNames;
                featureNs = S.featureNs;
                feature_str = S.feature_str;
                featureNamesAndN = S.featureNamesAndN;                
            else
                error('Feature file does not exist');
            end
        end
    end
        
    featureNames_out = featureNames;
    featureNs_out = featureNs;
    feature_str_out = feature_str;
    featureNamesAndN_out = featureNamesAndN;

%         feature_ids = cellfun(@(fet) find(strcmp(fet, allFeatureNames), 1), featureNames); %#ok<NODEF>
%         fet_str = cellstr2csslist(featureNames, '_'); %#ok<NODEF>
    if (nargin == 1) && isempty(varargin{1}) && (nargout < 2)
        featureNames_out = feature_str; 
    end
        
    

end


function [fet_str, n] = extractDigits(fet_str)
    n = [];
    idx_digits = isstrprop(fet_str, 'digit');
    if nnz(idx_digits) > 0
        n = str2double(fet_str(idx_digits));
        fet_str = fet_str(~idx_digits);
    end
end


function tf = strncmpiC(s,C,n)
    tf = cellfun(@(c) strncmpi(s,c,n), C);
end