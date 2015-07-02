function verify_stim_file (filename,x_size,y_size,varargin)

% verify_stim_file (filename,x_size,y_size)
%
% Display the first 100 frames of the stimulus movie (of size
% x_size x y_size pixels) in "filename".  Note that file will be
% read as "uint8" type.
  
if ~ischar(filename) || ~isnumeric(x_size) || ~isnumeric(y_size)
  error('Usage: verify_stim_file (filename,x_size,y_size)');
end

if ~exist(filename,'file')
  error(['No such file: ' filename]);
end

num_frame = 100;
if length(varargin)>0
   num_frame=varargin{1};
end
if num_frame>1000
   error('Number of frames may be too many.');
end

disp('This function will display the first few frames of the stimulus file.');
disp('If you can not see the movie, that may mean that ');
disp('   (1) the stimulus file is not formatted correctly and ');
disp('   (2) it cannot be read by the MID codes.');

fid = fopen(filename,'rb');
m = fread(fid,x_size*y_size*num_frame,'uint8');
fclose(fid);

m = reshape(m,[x_size,y_size,num_frame]);

h = figure; colormap gray;
set(h, 'DoubleBuffer', 'on');
for i = 1:num_frame
  imagesc(m(:,:,i));
  axis image;
  axis square;
  pause(10/1000);
end
