function check_consistent_weights_plotter(save_name)
%%%%%%%%%%%%%%%%%%%
%% JH 2/8/16
%% plot output from having run ABC dist weights several times to see if we get consistent weights 
%perhaps as a script?
close all;
if nargin <1
save_name = 'exp_decay_8';
end
load(sprintf('check_consistent_w_saved%s.mat',save_name));

%figure; 
%bar(mean(v));
%figure; 
%bar(mean(ssv));
%figure; 
%bar(mean(v+ssv));
figure; 
errorbar(mean(v),std(v),'linewidth',2);
box on;
set(gca, 'fontsize',20);
xlabel('index');
ylabel('$\log W$','Interpreter','LaTex'); 
FigHandle = figure;
  set(FigHandle, 'Position', [100, 100, 1049, 895]);
print(sprintf('consistent_weights_raw_%s',save_name),'-depsc');

q = v - repmat(mean(v,2),1,params.num_ss);
figure; 
errorbar(mean(q),std(q),'linewidth',2);
box on;
set(gca, 'fontsize',20);
xlabel('index');
ylabel('$\log W$','Interpreter','LaTex');
 FigHandle = figure;
  set(FigHandle, 'Position', [100, 100, 1049, 895]);
print(sprintf('consistent_weights_minus_mean_%s',save_name),'-depsc');

figure; 
box on;
errorbar(mean(v+sv),std(v+sv),'linewidth',2);
set(gca, 'fontsize',20);
xlabel('index');
ylabel('$\log W/S$','Interpreter','LaTex');
 FigHandle = figure;
  set(FigHandle, 'Position', [100, 100, 1049, 895]);
print(sprintf('consistent_weights_scaled_%s',save_name),'-depsc');

figure; 
box on;
errorbar(mean(q+sv),std(q+sv),'linewidth',2);
set(gca, 'fontsize',20);
xlabel('index');
ylabel('$\log W/S$','Interpreter','LaTex');
 FigHandle = figure;
  set(FigHandle, 'Position', [100, 100, 1049, 895]);
print(sprintf('consistent_weights_scaled_minus_mean_%s',save_name),'-depsc');

figure; imagesc(corr(v')); colormap('gray'); set(gca,'YDir','Normal');
set(gca, 'fontsize',20);
xlabel('run index 1');
ylabel('run index 2');
 FigHandle = figure;
  set(FigHandle, 'Position', [100, 100, 1049, 895]);
print(sprintf('correlation_weights_%s',save_name),'-depsc');

 figure; imagesc(corr(q')); colormap('gray'); set(gca,'YDir','Normal');

 figure; imagesc((corr(sv'))); colormap('gray')
 figure; imagesc((corr(v'+sv'))); colormap('gray')
