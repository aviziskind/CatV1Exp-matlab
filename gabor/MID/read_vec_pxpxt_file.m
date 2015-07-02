function [v_all,out_vec] = read_vec_pxpxt_file (filename,stim_dim,varargin)

% read_vec_pxpxt_file (filename,stim_dim,varargin)
% 
% Reads and displays the results of do_MaxInfoDim_Nvec_verX codes.
%
% filename = path/cellname.
% stim_dim = [dx,dy,nlags,cx].
% (filename will be appended with the appropriate stim_dim info.)
% varargin{1} = Nbins (default is 15).
% varargin{2} = Nvec (default is 1).
% varargin{3} = Nparts (default is 4).
% 
% Note these values should match the parameters used in finding the
% vectors.
makePlots = (nargin == 0);

if length(stim_dim)~=4
  error('Wrong input format.');
end
Nn = prod(stim_dim(1:3)); % Size of each vector.
nlags = stim_dim(3);
name_str = sprintf('%ux%ux%u_%u',stim_dim);

Nbins = 15; % Number of bins used for probability distribution.
Nparts = 4; % Number of jack-knife splits.
Nvec = 1; % Number of dimensions.
if ~isempty(varargin)
  Nbins = varargin{1};
  if length(varargin) > 1
    Nvec = varargin{2};
    if length(varargin) > 2
      Nparts = varargin{3};
    end
  end
end

v_all = [];

% Initialize.
x = zeros(1,Nbins*Nvec);
px = zeros(1,Nbins^Nvec);
pxt = zeros(1,Nbins^Nvec);
% ior = zeros(1,Nbins^Nvec);
% rbar = 0;
% info = 0;
x_mtx = zeros(Nparts,Nbins*Nvec); 
px_mtx = zeros(Nparts,Nbins,Nvec); % Collapsed probabilities.
pxt_mtx = zeros(Nparts,Nbins,Nvec);
% norm_spike = 1;

% Read vectors.
v_all = cell(1,Nvec);
for j = 1:Nparts
  fullname = sprintf('%s_%s_jack%u.dat',filename,name_str,j);
  if ~exist(fullname,'file')
    disp(['No such file: ' fullname]);
    v = zeros(1,Nn);
    %return;
  else
    fid = fopen(fullname,'rb');
    v = fread(fid,'double');
    fclose(fid);
  end

  for n = 1:Nvec
    v_tmp = v([1:Nn]+(n-1)*Nn);
    v_all{n}(:,:,j) = reshape(v_tmp,stim_dim(1),stim_dim(2)*nlags);
  end
end

% Read P(x) and P(x|spikes).
for j = 1:Nparts
  fullname = sprintf('%s_%s_pxpxt_jack%u.dat',filename,name_str,j);
  fullname = strrep(fullname,'_v_','_');
  if ~exist(fullname,'file')
    disp(['No such file: ' fullname]);
    %return;
  else
    [x,px,pxt,rbar,info] = read_pxpxt_Nvec (fullname,Nbins,Nvec);
    ior = pxt./(px+eps);
  end
  
  % Note these _mtx matrices will have Nparts rows.
  x_mtx(j,:) = x;
  px_mtx(j,:,:) = get_collapsed(px,Nbins,Nvec);
  pxt_mtx(j,:,:) = get_collapsed(pxt,Nbins,Nvec);
end

% Some MID vectors may be polarity-reversed wrt others. flip_sign keeps
% track of which ones should be flipped. This vector is based on the
% dot product with the first vector.
flip_sign = ones(Nparts,Nvec);
for n = 1:Nvec
  for j = 2:Nparts
    dot_val = sum(sum(v_all{n}(:,:,j).*v_all{n}(:,:,1)));
    if dot_val < 0
      flip_sign(j,n) = -1;
      v_all{n}(:,:,j) = -1*v_all{n}(:,:,j);
      
      px_mtx(j,:,n)  = px_mtx(j,end:-1:1,n);
      pxt_mtx(j,:,n) = pxt_mtx(j,end:-1:1,n);
    end
  end
end

ior_mtx = pxt_mtx./(px_mtx+eps);
out_vec = compute_average_ptx(x_mtx,px_mtx,pxt_mtx,ior_mtx,Nparts,Nbins,Nvec);

if makePlots
    plot_vec_pxpxt_Nvec (v_all,out_vec,x_mtx,px_mtx,ior_mtx,Nvec,Nparts);
end

%%-------------------------------------------------------
function out_vec = compute_average_ptx (x_mtx,px_mtx,pxt_mtx,ior_mtx,Nparts,Nbins,Nvec)

% Find average probability distributions along Nparts jack-knives.
% This function was not tested for Nvec>1.

if Nparts ~= size(x_mtx,1)
  error('dimension mismatch for x_mtx (Nparts).');
end
if Nbins*Nvec ~= size(x_mtx,2)
  error('dimension mismatch for x_mtx (Nbins,Nvec).');
end

Nbins_short = 14;
new_Nbins = Nbins_short-1;

npoints = zeros(Nparts,new_Nbins,Nvec);
px_rescaled = zeros(Nparts,new_Nbins,Nvec);
pxt_rescaled = zeros(Nparts,new_Nbins,Nvec);
ior_rescaled = zeros(Nparts,new_Nbins,Nvec);

for n = 1:Nvec
  xmin(n) = min(min(x_mtx(:,[1:Nbins]+(n-1)*Nbins)));
  xmax(n) = max(max(x_mtx(:,[1:Nbins]+(n-1)*Nbins)));
  edges{n} = linspace(xmin(n),xmax(n),Nbins_short);

  for i = 1:new_Nbins
    for j = 1:Nparts
      curr_x = x_mtx(j,[1:Nbins]+(n-1)*Nbins);
      ind = find( curr_x>=edges{n}(i) & curr_x<edges{n}(i+1) );
      
      if ~isempty(ind)
	if (min(ind)<1) || (max(ind)>Nbins)
	  error('Index outside of Nbins limit.');
	end
	npoints(j,i,n) = 1; % Number of jack-knife estimates.
	px_rescaled(j,i,n) = px_rescaled(j,i,n)+sum(px_mtx(j,ind,n));
	pxt_rescaled(j,i,n) = pxt_rescaled(j,i,n)+sum(pxt_mtx(j,ind,n));
      end
    end
    x_mean(i,n)=0.5*(edges{n}(i)+edges{n}(i+1));
  end
end

% Take the average across jack-knife dimension (1; indexed by j previously).
ior_rescaled = pxt_rescaled./(px_rescaled+eps);
nsamples = sum(npoints,1);
px_mean = sum(px_rescaled,1)./(nsamples+eps);
ior_mean = sum(ior_rescaled,1)./(nsamples+eps);
px_std  = sqrt(var(px_rescaled,[],1)./(nsamples+eps)).*(nsamples-1); 
ior_std = sqrt(var(ior_rescaled,[],1)./(nsamples+eps)).*(nsamples-1); 
% Note the (N-1)/N factor, yielding standard errors from jackknife estimates.

out_vec = zeros(new_Nbins,Nvec,5); 
out_vec(:,:,1) = x_mean;
out_vec(:,:,2) = px_mean;
out_vec(:,:,3) = ior_mean;
out_vec(:,:,4) = px_std;
out_vec(:,:,5) = ior_std;

%%--------------------------------------
function plot_vec_pxpxt_Nvec (v_all,out_vec,x_mtx,px_mtx,ior_mtx,Nvec,Nparts)

% This version will plot Nparts together, and the average
% separately.

if length(v_all) ~= Nvec
  error('Size of v_all is not equal to Nvec.');
end
num_row = 2;
num_col = Nvec;

color_vector1 = [0,0,1];
color_vector2 = [1,0,1];

px_scale = 10;

set(gcf,'Name','Average Vector and Nonlinaerity');
for n = 1:Nvec
  subplot(num_row,num_col,n+num_col);
  % Plot x_mean against ior_mean and ior_std.
  errorbar(out_vec(:,n,1),out_vec(:,n,3),out_vec(:,n,5), ...
      'o-','Color',color_vector1,'Linewidth',2)
  hold on;
  % Plot x_mean against px_mean.
  px_scale = max(out_vec(:,n,3))/max(out_vec(:,n,2))*0.75;
  plot(out_vec(:,n,1),out_vec(:,n,2)*px_scale,'.-','Color',color_vector2);
  hold off;
  
  xlabel('Proj. Value (per \sigma)');
  ylabel('Prob. of spike');
  
  % Set axis properties.
  xmax = max(out_vec(:,n,1))*1.1;
  xmin = min(out_vec(:,n,1));
  if (xmin<0)
    xmin=xmin*1.1;
  else
    xmin=xmin*0.9;
  end
  
  if (xmin==xmax)
    xmin = 0;
    xmax = 1;
  end
  %xlim([xmin xmax]);
  xlim([-1 1]*max([xmin xmax]));
  ylim([0 max([1, max(out_vec(:,n,3)*1.2)])]);
  axis square;
  box off;
  
  % Plot average vectors.
  subplot(num_row,num_col,n);
  clim = [-1 1]*max(max(abs(mean(v_all{n},3))));
  imagesc(mean(v_all{n},3),clim);
  axis image;
  set(gca,'XTick',[-2 -1]);
  set(gca,'YTick',[-2 -1]);
  colormap hot;
end

%set(gcf,'Position',[0 0 500 400]+100);

% Plot all jack-knife vectors with the same clim.
if 0
  num_row = Nvec;
  num_col = Nparts+1;
  figure(1001);
  set(gcf,'Name','All Jackknives');
  for n = 1:Nvec
    
    clim = [-1 1]*max(abs(v_all{n}(:)));
    if max(clim) < eps
      clim = [-1 1];
    end
    v_display = [];
    for j = 1:Nparts
      v_display = [v_display; v_all{n}(:,:,j)];    
    end
    v_display = [v_display; mean(v_all{n},3)];
    subplot(num_row,num_col,1+(n-1)*num_col);
    imagesc(v_display,clim);
    %colorbar;
    axis image;
    set(gca,'XTick',[-2 -1]);
    set(gca,'YTick',[-2 -1]);
    title(['Receptive Field (average at the bottom)']);
    
    for j = 1:Nparts
      subplot(num_row,num_col,1+(n-1)*num_col+j);
      plot(x_mtx(j,:),ior_mtx(j,:),'.-','Color',color_vector1);
      hold on;
      px_scale = max(ior_mtx(:))/max(px_mtx(:))*0.75;
      plot(x_mtx(j,:),px_mtx(j,:)*px_scale,'-','Color',color_vector2);
      hold off;
      axis square;
      xlim([-1 1]*max(abs(x_mtx(j,:))));
      ylim([0 1.2*max(ior_mtx(:))]);
      box off;
      title(['Jack #' num2str(j)]);
    end
  end
  colormap hot;
end

%%-------------------------------------------------------
function [x,px,pxt,rbar,info] = read_pxpxt_Nvec (fn,Nbins,Nvec)

x = zeros(1,Nbins*Nvec);
px = zeros(1,Nbins^Nvec);
pxt = zeros(1,Nbins^Nvec);
rbar = 0;
info = 0;

minpx = 0;
fp = fopen(fn,'r');
if (fp==-1)
  display('Error opening file');
  display(fn);
  return
end

x = fread(fp,Nbins*Nvec,'double');
px = fread(fp,Nbins^Nvec,'double');
pxt = fread(fp,Nbins^Nvec,'double');

ind0 = find(px<minpx);
pxt(ind0) = 0;

if (abs(sum(px)-1)>0.001) % Sum of probability should be 1.
  display(sprintf('Check probability normalization: sum of P(x)=%f',sum(px)));
  return;
end
rbar = fread(fp,1,'double');
info = fread(fp,1,'double');
%disp([fn]);
%disp(['Rbar = ' num2str(rbar,'%5.4f') '; info = ' num2str(-info,'%5.4f')]);
fclose(fp);

