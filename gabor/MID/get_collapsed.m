function vec = get_collapsed (mtx,Nbins,Nvec)

% Given a mtx (size of Nbins^Nvec), return the collapsed vec
% which is indexed along each dimension.  Note that this function
% should correspond with ind_convert_tbl in the C codes for writing
% px and pxt variables. 
% From compute_prob.c
%  int ind_convert_tbl[Nvec];
%  for (j=1;j<=Nvec;j++) ind_convert_tbl[j]=(int)pow(Nbins,j-1);
% e.g., for Nvec=2 and Nbins=5, 
%   [1 3] = (1-1)*5^(1-1) + (3-1)*5^(2-1) + 1 = 11.
%
% mtx is a joint prob distributions (px or pxt), 
% and returned vec (size of Nbins*Nvec) is a collapsed prob
% distribution (px_collapse or pxt_collapse).
%
% To test: 
%   mtx = zeros(1,5^2);
%   mtx(11) = 1;
%   vec = get_collapsed(mtx,5,2);
%   It should give vec(1,1)=1 and vec(3,2)=1, otherwise 0.

if Nbins^Nvec ~= length(mtx)
  error('size mismatch of Nbins^Nvec');
end

ind_convert_tbl = Nbins.^([1:Nvec]-1);

vec = zeros(Nbins,Nvec);

for i = 1:Nbins^Nvec
  num_tmp = i;
  idx = zeros(1,Nvec);

  for n = Nvec:-1:1
    remainder = rem(num_tmp,Nbins^(n-1));
    if remainder > 0 
      idx(n) = floor(num_tmp/(Nbins^(n-1))) + 1;
    else
      idx(n) = (num_tmp/(Nbins^(n-1)));
    end
    num_tmp = num_tmp - (idx(n)-1)*(Nbins^(n-1));

    vec(idx(n),n) = vec(idx(n),n) + mtx(i);
  end
end

