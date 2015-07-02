function out_all = read_app_file (filename,varargin)

% Read the approach file.

num_row = 9;
datadetail = '';

Nparts = 4;
if length(varargin)>0
  datadetail = sprintf('_%ux%ux%u_%u',varargin{1});
  if length(varargin)>1
    Nparts = varargin{2};
  end
end
out_all = cell(Nparts,1);

for j = 1:Nparts
  fullname = sprintf('%s%s_jack%u.dat',filename,datadetail,j);

  if ~exist(fullname,'file')
    disp(['No such file: ' fullname]);
  else

    c = read_single_app_file (fullname);

    out_all{j} = c;
    subplot(num_row,Nparts,j);
    title(['Set #' num2str(j)]);

    for k = 1:num_row
      subplot(num_row,Nparts,j+(k-1)*Nparts);    
      plot(c{k},'k-');
      box off; 
      xlim_val = ceil(length(c{k})/100)*100;
      xlim([0 xlim_val]);
      set(gca,'XTick',[0 xlim_val]);
      ylim([0 1]);
      set(gca,'YTick',[0 1]);
      if max(c{k})>1
        ylim([0 max(c{k})]);
        set(gca,'YTick',[0 max(c{k})]);
      end
      set(gca,'TickDir','out');
    end
  end
end

subplot(num_row,Nparts,1+0*Nparts);
ylabel('Info(v1)');

subplot(num_row,Nparts,1+2*Nparts);
ylabel('Info(v1 best)');

subplot(num_row,Nparts,1+4*Nparts);
ylabel('Info(v all)')

subplot(num_row,Nparts,1+6*Nparts);
ylabel('Temp')

subplot(num_row,Nparts,1+7*Nparts);
ylabel('All Count');

subplot(num_row,Nparts,1+8*Nparts);
ylabel('Count');
xlabel('Iteration');

for j = 1:Nparts
  subplot(num_row,Nparts,j);
  title(['Set #' num2str(j)]);
end

set(gcf,'PaperPosition',[0 0 8 11]);

%%===================================================
function c = read_single_app_file (filename)

fp = fopen(filename,'r');

% See do_MID.cpp.
c = textscan(fp,'%f %f %f %f %f %f %f %d %d %f %f %f %d', 'CommentStyle','#');
fclose(fp);

if 0 % display
  ind_temp = 6;
  ind_count = 7;
  ind_v1dotv2 = 11;
  clf;
  %plot_r = [1:min([20 length(c{1})])];
  plot_r = [1:length(c{1})]; % Plotting range.
  subplot(5,2,1);
  plot(c{1}(plot_r)); axis tight;
  title('info train (v)');
  subplot(5,2,3);
  plot(c{2}(plot_r)); axis tight;
  title('info test (v)');
  subplot(5,2,5);
  plot(c{3}(plot_r)); axis tight;
  title('info train (vbest)');
  subplot(5,2,7);
  plot(c{5}(plot_r)); axis tight;
  title('info test (vtest)');
  subplot(5,2,9);
  plot(c{ind_temp}(plot_r)); title('Temperature'); axis tight;
  subplot(5,2,2);
  plot(c{ind_count}(plot_r)); title('count'); axis tight;
  subplot(5,2,4);
  %plot(c{ind_v1dotv2}(plot_r)); title('v1*v2'); axis tight;
end
