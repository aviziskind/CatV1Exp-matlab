% demo.m
%
% Displays the demo results.

disp('======================================================');
disp('Running the demo script for maximally info dimensions.');
disp('See [Sharpee et al., NIPS, 2002; Neural Comp. 2004].');
disp('See readme.txt');
disp('');

% addpath('../Matlab_Codes');
datadir = './My_Data/';
if ~exist(datadir,'dir')
  fprintf('>>Folder %s does not exist: ',datadir);
  fprintf('using the pre-computed data.\n');
  datadir = './Data/'; % Pre-computed results are here.
end
dim = [30,30,1,4];
Nbins = 32;
Nvec = 1; % Search for 1 dimension.

% Read the approach file.
%figure(200);
%read_app_file([datadir 'app_democell'],dim);

% Read spike-triggered average and MID vectors.
disp('Loading spike-triggered average...');
v_STA = read_vec_pxpxt_file([datadir 'STA_democell'],dim,Nbins);
disp('Loading maximally informative dimensions...');
v_MID = read_vec_pxpxt_file([datadir 'Test_v_democell'],dim,Nbins);

%disp('Loading the evolution of informative dimension...');
%Mv = show_v_evol([datadir 'Test_v_democell'],dim,Nvec);

% Compare with the true vector.
avg_STA = mean(v_STA{1},3);
avg_MID = mean(v_MID{1},3);

fp = fopen('./true_vec.bin','rb');
true_v = reshape(fread(fp,'double'),[30,30]);

figure(100);
subplot(1,3,1);
imagesc(true_v); axis image;
title('True RF');

subplot(1,3,2);
imagesc(avg_MID); axis image;
title('MID');
dotval = true_v(:)'*avg_MID(:);
dotval = dotval/sqrt(sum(true_v(:).^2))/sqrt(sum(avg_MID(:).^2));
xlabel(['dot-product with true = ' num2str(dotval,'%4.3f')]);

subplot(1,3,3);
imagesc(avg_STA); axis image;
title('STA');
dotval = true_v(:)'*avg_STA(:);
dotval = dotval/sqrt(sum(true_v(:).^2))/sqrt(sum(avg_STA(:).^2));
xlabel(['dot-product with true = ' num2str(dotval,'%4.3f')]);

% Report some numbers.
fp = fopen('./demo_spike_Noise0.50_Rep1_1.isk');
Nspikes = sum(fread(fp,'uint8'));
fclose(fp);

disp(['Number of movie frames = 16384.']);
disp(['Number of spikes in demo file = ' num2str(Nspikes) '.']);
disp(['Hence, average spike per frame = ' num2str(Nspikes/16384,'%4.3f') '.']);

disp('Done');
disp('======================================================');
