function X = decompress(S)
    
    if isfield(S, 'uValsLists') || isfield(S, 'uValsListSkips')
        uVals = S.uVals;
        
        if isfield(S, 'uValsLists') 
            uValLists = S.uValsLists;
            nTot = sum( cellfun(@length, uValLists) );
        elseif isfield(S, 'uValsListSkips')            
            uValLists = { cumsum(double(S.uValsListSkips)) };            
            nTot = length(uValLists{1})+1;
        end        

        if ~isfield(S, 'defaultVal') || (S.defaultVal == 0)
            if isfield(S, 'inputIsSparse') && S.inputIsSparse  
                X = spalloc(S.dims(1), S.dims(2), nTot);
                if islogical(uVals(1))
                    X = logical(X);
                end                            
            else
                if islogical(uVals)
                    X = false(S.dims);
                else
                    X = zeros(S.dims, class(uVals));
                end
            end
            
        else            
            X = S.defaultVal(ones(S.dims));
        end        
        
        for i = 1:length(S.uVals)
            X(uValLists{i}) = uVals(i);
        end            
                
    elseif isfield(S, 'vals_idx')        
        uVals = S.uVals;
        X = uVals(S.vals_idx);            
    elseif isfield(S, 'orig_vals')        
        X = S.orig_vals;
    else
        X = S;
    end
    
    if isfield(S, 'cell_sizeX')
        %%
        sizeX = size(X);        
        n_cellElements = sizeX(S.cell_dimCompress);
        sizeX(S.cell_dimCompress) = [];
%         sizeX_C = num2cell(sizeX);
        X_C = cell(S.cell_sizeX);
        idx_c = repmat({':'}, 1, length(sizeX));
        for j = 1:n_cellElements
            X_C{j} = X(idx_c{:}, j);
        end
        3;
%         X_C2 = reshape( mat2cell(X, sizeX_C{:}, ones(1, n_cellElements) ), S.cell_sizeX);
%         assert(isequal(X_C, X_C2));
        X = X_C;

                
    end
    
    
end


%     if ~isstruct(S)
%         X = S;
%         return;
%     end
