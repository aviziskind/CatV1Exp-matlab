function X = msequenceGenerator(varargin)
    persistent uTapRegister nMseqBits uMseqRegister;

    if (nargin > 0) && ischar(varargin{1}) && strcmp(varargin{1}, 'init')

        [uTapRegister, nMseqBits] = elements(varargin(2:3));
        uTapRegister = uint16(uTapRegister);

        % init mseq generating register so that lower mseq_bits  are equal to 1 and higher bits are 0
        uMseqRegister=bitshift(1, nMseqBits)-1;
        return;
            
    elseif (nargin == 0) || isnumeric(varargin{1})
        if nargin == 0                        % func;
            dims = [1 1];
        elseif nargin == 1
            if numel(varargin{1}) > 1         % func([m,n]); func([m,n,p,...]);
                dims = varargin{1};
            elseif numel(varargin{1}) == 1    % func(m); --> func(m,m);
                dims = [varargin{1}, varargin{1}];
            end
        else
            dims = cell2mat(varargin);        % func(m,n); func(m,n,p,...);
        end
    end

    nBitsToDo = prod(dims);
    X = zeros(dims, 'uint8');
    progressBar('init-', nBitsToDo);
    for i = 1:nBitsToDo
        progressBar(i);
        % get the current bit
        iMseqValue=bitget(uMseqRegister, 1);  

        % use the tap register to mask bits
        uMaskedMseqReg= bitand(uMseqRegister, uTapRegister);

        % count number of bits equal to 1 in the uMaskedMseqReg
        bit_cnt = nnz( dec2bin(uMaskedMseqReg) == '1');

        % shift right 1 bit [the highest bit is 0 after this]
        uMseqRegister= bitshift(uMseqRegister, -1);

        % if number of set bits is odd, put 1 in the highest bit
        if mod( bit_cnt, 2 )
            uMseqRegister = bitset( uMseqRegister, nMseqBits, 1);
        end

        X(i) = uint8( iMseqValue == 1 );
    end
end

