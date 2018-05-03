function generate_sample_spatial_sum_stat_figure(theta,time_vec,fig_name,ax_lbl)
%% 20/1/17 JH
%% plotting function for N samples paths from spatial test problem
%% shows variability in output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% edit input function handle and parameters for each test problem
%% edit the file name to print output to
%%%%%%%%%%%%%%%%%%%%

if nargin <1
    theta = 0.1; %[0.01, 0.05]; %[1, 1/25, 1/500, 1/2]; %[0.01, 0.01]; %0.1;
end
rng(234);
%assume test_problem is spatial test problem
    f = @(x,y,z) spatial_test_problem(x,y,z);
    test_problem = @(x) f(x, 20, 20/(2^3-1));
t_plot_vec = [0,3,7];
dx = 1/8;
xpos = dx/2:dx:(1-dx/2);
N=4;
sz = size(test_problem(theta),2);
num_t = length(time_vec)-1; %number of time points
num_box = 8; %number of boxes;

ss = zeros(num_t,num_box,N);
if size(test_problem(theta),3)~=1
    warning('can only currently deal with summary stats that are row vectors');
end

figure;
subplot(length(t_plot_vec),1,1);
set(gca,'fontsize',20);
            hold all;
            stairs(xpos,10*[ones(1,4),zeros(1,4)],'linewidth',3)
            axis([0,1,0,15]);
            yticks([0 15]);
%            xlabel('$x$','interpreter','latex');
%            ylabel(sprintf('$%s$',ax_lbl),'interpreter','latex');
%	    title('t=0');
            box on;

for j=1:N
    ss(:,:,j) = reshape(test_problem(theta),[],num_box);
	for r=2:numel(t_plot_vec)
	    subplot(length(t_plot_vec),1,r);
            set(gca,'fontsize',20);
            hold all;
	    stairs(xpos,ss(t_plot_vec(r),:,j),'linewidth',3)
	    axis([0,1,0,15]);
        yticks([0 15])
%	    xlabel('box'); %,'interpreter','latex');
%	    ylabel(sprintf('$%s$',ax_lbl),'interpreter','latex');
	    box on;
%	   title(sprintf('t=%0.1f',time_vec(1+t_plot_vec(r))));
	end
end
            xlabel('$x$','interpreter','latex');
	    subplot(length(t_plot_vec),1,2);
            ylabel(sprintf('$%s$',ax_lbl),'interpreter','latex');
print(sprintf('../Latex/Figs/%s_typical_output',fig_name),'-depsc');
