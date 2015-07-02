% function s = stackMatrices(dim, varargin)
%     nmat = length(varargin);
%     % check matrices are all same size
%     sz = size(varargin{1});
%     for i = 2:nmat
%         if ~all( size(varargin) == sz)
%             error('All matrices must be the same size');
%         end
%     end
%     
%     nstacks = zeros(1,nmat);
%     for mi = 1:nmat
%         nstacks(mi) = size(varargin{mi}, 3);
%     end
%     cumnstacks = cumsum(nstacks);
% 
%     [nrows, ncols] = elements(size(varargin{1}), [1,2] );
%     s = zeros(nrows, ncols, sum(nstacks) );
%     for mi = 1:nmat
%         indices = [cumnstacks(mi)-nstacks(mi)+1:cumnstacks(mi)];
%         s(:,:,indices) = varargin{mi};
%     end
% end
use CAT instead !