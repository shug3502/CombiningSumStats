function generate_sample_dimerization_sum_stat_figure(theta,time_vec,fig_name,ax_lbl)
%% 28/2/17 JH
%% plotting function for N samples paths from spatial test problem
%% shows variability in output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% edit input function handle and parameters for each test problem
%% edit the file name to print output to
%%%%%%%%%%%%%%%%%%%%

rng(234);
%assume test_problem is spatial test problem
    f = @(x,y,z) dimerization_simbio(x,y,z);
    test_problem = @(x) f(x, 100, 100/(2^5-1));

N=4;
sz = size(test_problem(theta),2);
num_data_pts = (2^5-1); %number of boxes;

ss = zeros(num_data_pts,N,3);
if size(test_problem(theta),3)~=1
    warning('can only currently deal with summary stats that are row vectors');
end

axis_ticks = [10,5,4]*10^4; %max y axis value

figure;
for j=1:N
    ss(:,j,:) = reshape(test_problem(theta),num_data_pts,[],3);
	for species_ind=1:3
	    subplot(3,1,species_ind);
            set(gca,'fontsize',20);
            hold all;
	    stairs(time_vec(2:end),ss(:,j,species_ind),'linewidth',3)
%	    axis([0,1,0,15]);
        yticks([0 axis_ticks(species_ind)])
%	    xlabel('box'); %,'interpreter','latex');
	    ylabel(sprintf('$S_%d(t)$',species_ind),'interpreter','latex');
	    box on;
%	   title(sprintf('t=%0.1f',time_vec(1+t_plot_vec(r))));
	end
end
xlabel('$t$','interpreter','latex');
%subplot(length(t_plot_vec),1,2);
%ylabel(sprintf('$%s$',ax_lbl),'interpreter','latex');
print(sprintf('../Latex/Figs/%s_typical_output',fig_name),'-depsc');
