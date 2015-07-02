function Z = intersectq(varargin)
    
    if (nargin == 2)
        X = varargin{1};
        Y = varargin{2};
        nx = length(X);
        ny = length(Y);
        if (nx > ny)
            indsYinX = binarySearch(X,Y, 1, 0);
            Z = X( nonzeros(indsYinX) );
        else
            indsXinY = binarySearch(Y,X, 1, 0);
            Z = Y( nonzeros(indsXinY) );
        end
        
    else
       Z = intersectq( varargin{1},  ...
                intersectq( varargin{2:end} ) ) ;        
    end
            
end