function [Mv,m_display] = show_v_evol (filename,dim,Nvec,varargin)

% show_v_evol (filename,dim,Nvec)
%
% Show how the vbest vector evolves during the optimization
% algorithm.

datadetail = sprintf('%ux%ux%u_%u',dim);
Nparts = 4;
if length(varargin)>0
  Nparts = varargin{1};
end

m_display = cell(1,Nparts);
Nframes = zeros(1,Nparts);

Nn = prod(dim(1:3));

for j = 1:Nparts
  fullname = sprintf('%s_%s_jack%u.dat',filename,datadetail,j);
  fullname = strrep(fullname,'_v_','_v_evol_');
  fp = fopen(fullname,'rb');
  if fp<0
    disp([fullname ' does not exist.']);
  else
    m = fread(fp,'double');

    num_frame = length(m)/Nn/Nvec;
    Nframes(j) = num_frame;
    
    m_tmp = zeros(dim(1),dim(2)*dim(3),Nvec,num_frame);

    for i = 1:num_frame
      for n = 1:Nvec
	tmp_i = m([1:Nn]+(n-1)*Nn+(i-1)*Nn*Nvec);
	m_tmp(:,:,n,i) = reshape(tmp_i,dim(1),dim(2)*dim(3));
      end
    end
    m_display{j} = m_tmp;

    fclose(fp);
  end
end

% Now display.  
idx = 0;
figure(200); clf;
for j = 1:Nparts
  for i = 1:Nframes(j)
    for n = 1:Nvec
      subplot(1,Nvec,n);
      imagesc(squeeze(m_display{j}(:,:,n,i)));
      axis image off;
      %title(['Jack #' num2str(j)]);
    end
    drawnow;
    idx = idx+1;
    Mv(idx) = getframe;
  end
end
close;

% Save it as an avi file.
%movie2avi(Mv,'tmp_evol.avi','FPS',5);
