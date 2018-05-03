%% run function to plot figs with tyoical outputs
%%could rewrite this to load the appropriate parameters

%% generate_sample_sum_stat_figure(test_problem,theta,fig_name)
%% 4/9/16 JH
%%%%%%%%%%%%%%%%%%%%

% %for spatial test problem
%     theta = 0.1; 
%     f = @(x,y,z) spatial_test_problem(x,y,z);
%     test_problem = @(x) f(x, 20, 20/(2^3-1));
%     time_vec = linspace(0,20,2^3);
%    fig_name = 'spatial';
%    axis_label = 'S_i(t)';
% generate_sample_spatial_sum_stat_figure(theta,time_vec,fig_name,axis_label)

    theta = [1, 1/25, 1/500, 1/2]; %[0.01, 0.01]; %0.1;
    f = @(x,y,z) dimerization_simbio(x,y,z);
    test_problem = @(x) f(x, 100, 100/(2^5-1));
n=2^5-1;
base = 4;
q = floor(log(100)/log(base));
vec = base.^((q-(n-2)):q);
    time_vec = [0,vec,100]; %linspace(0,100,2^5);
   fig_name = 'dimerization';
   axis_label = 'S_i(t)';
generate_sample_dimerization_sum_stat_figure(theta,time_vec,fig_name,axis_label)   
%generate_sample_sum_stat_figure(test_problem,theta,time_vec,fig_name,axis_label)
% 
%    theta = [0.01, 0.05]; %[1, 1/25, 1/500, 1/2]; %[0.01, 0.01]; %0.1;
%     f = @(x,y,z) birth_death_simbio(x,y,z);
%     test_problem = @(x) f(x, 40, 40/(2^5-1));
%     time_vec = linspace(0,40,2^5);
%    fig_name = 'birth_death';
%    axis_label = 'A(t)';
% generate_sample_sum_stat_figure(test_problem,theta,time_vec,fig_name,axis_label)
% 
%    theta = [0.1, 0.01]; %0.1;
%     f = @(x,y,z) test_problem2(x,y,z);
%     test_problem = @(x) f(x, 20, 20/(2^5-1));
%     time_vec = linspace(0,20,2^5);
%    fig_name = 'exp_decay';
%    axis_label = 'A(t)';
% generate_sample_sum_stat_figure(test_problem,theta,time_vec,fig_name,axis_label)
% 
% 
