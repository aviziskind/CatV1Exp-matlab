function C = conv2_circ(A, B, varargin)

    % *** can make more efficient ***

    if nargin == 3
        [circ_dims] = elements(varargin);
    elseif nargin == 4
        [shape, circ_dims] = elements(varargin);
    end
    
    wrap1 = any(circ_dims == 1);
    wrap2 = any(circ_dims == 2);
        
    [mA, nA] = size(A);
    if wrap1        
        A = [A; A; A];       
        m_inds = mA+1:2*mA;
    else
        m_inds = 1:mA;
    end
    
    if wrap2
        A = [A, A, A];
        n_inds = nA+1:2*nA;
    else
        n_inds = 1:nA;
    end    
    
    if exist('shape', 'var')
        C = conv2(A, B, shape);
    else
        C = conv2(A, B);        
    end
    C = C(m_inds, n_inds);
        
    
end