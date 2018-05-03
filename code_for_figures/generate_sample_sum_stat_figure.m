function generate_sample_sum_stat_figure(test_problem,theta,time_vec,fig_name,ax_lbl)
%% 1/8/16 JH
%% plotting function for N samples paths from a test problem
%% shows variability in output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% edit input function handle and parameters for each test problem
%% edit the file name to print output to
%%%%%%%%%%%%%%%%%%%%

if nargin <1
    theta = 0.1; %[0.01, 0.05]; %[1, 1/25, 1/500, 1/2]; %[0.01, 0.01]; %0.1;
    f = @(x,y,z) spatial_test_problem(x,y,z);
    test_problem = @(x) f(x, 20, 20/(2^3-1));
end
rng(234);
%assume test_problem is a function which outputs a row vector
N=4;
sz = size(test_problem(theta),2);
if ~mod(sz,length(time_vec)-1) 
    k=length(time_vec)-1;
    num_species = sz/k;
    t = time_vec(2:end);
else
    k=length(time_vec);
    num_species = 1;
    t = time_vec+time_vec(2);
end
ss = zeros(N,k,num_species);
if size(test_problem(theta),3)~=1
    warning('can only currently deal with summary stats that are row vectors');
end
figure;
set(gca,'fontsize',20);
hold all;
for j=1:N
    ss(j,:,:) = reshape(test_problem(theta),1,k,num_species);
	for r=1:num_species
	    stairs(t,ss(j,:,r),'linewidth',3)
	end
end
axis([0,max(time_vec),0,max(max(max(ss)))]);
xlabel('$t$','interpreter','latex');
ylabel(sprintf('$%s$',ax_lbl),'interpreter','latex');

box on;
%leg_str = cell(N,1);
%for i=1:N
%    leg_str{i} = sprintf('realization %d',i);
%end
%legend(leg_str);
print(sprintf('../Latex/Figs/%s_typical_output',fig_name),'-depsc');
