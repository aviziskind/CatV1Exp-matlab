function verify_spike_file (filename)

% verify_spike_file (filename)
%
% Read the spikes in "filename" as "uint8" type.

if ~ischar(filename)
  error('Usage: verify_spike_file (filename)');
end

if ~exist(filename,'file')
  error(['No such file: ' filename]);
end

fp = fopen(filename,'rb');
locator = fread(fp,'uint8');
fclose(fp);

Movie_length = length(locator);
tot_spikes = sum(locator);
disp(['Total number of spikes = ' num2str(tot_spikes)]);
disp(['Total number of movie frames = ' num2str(Movie_length)]);

t_range = [1:min([500 Movie_length])];

h = figure;
plot(locator(t_range),'.');
xlabel('Time (frames)');
ylabel('Num. Spikes');
title(['Tot. Spikes = ' num2str(tot_spikes) ...
  '; Num. Frames = ' num2str(Movie_length)]);