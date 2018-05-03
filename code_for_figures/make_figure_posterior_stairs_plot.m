function make_figure_posterior_stairs_plot(params)
%created 6/6/16
%JH

%script to make figure with posteriors for different methods, using stairs
%to create run: adapt_weights_of_ABC_tidied(params)
%with
%     params.N = 5*10^3; %number of ABC samples to generate at each generation
%     params.max_num_weights = 1000; %look at this many weights/this many generations in the gradient descent
%     params.num_params = 2;
%     params.num_ss=66;
%     params.num_generations = 5;
%     params.alpha = 0.05;
%     params.proposal_sd = 0.25;
%     params.with_plot = 1;
%     params.step_size = 0.5;
%     params.reduce_step_size_by=2;
%     params.set_to_zero=0;
%     params.set_to_uniform_weights = 0; %could set method with a string
%     params.set_to_scaled_weights = 0;
%     params.lambda = 2*10^-4;
%     params.problem_t_end=20;
%     params.recording_interval = params.problem_t_end/(2^6-1);
%     params.ref = [0,0];%[0,0,-prior_width/2]; %reference vector for (very slightly) more flexible prior
%     params.save_name = 'mRNA_v1';
%     params.x = 234; %random seed

%then with params.set_to_uniform_weights = 0;
%then with params.set_to_scaled_weights = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
addpath ../ %make use can access functions in Summary_stats directory
run_simulations = 1;
include_stan_post = 1;

if nargin<1
    my_params_store; %load params (or load different params)
end

if run_simulations
    tic;
    qq = params.test_problem(params.theta_real,params.problem_t_end,params.recording_interval);
    toc
    
    hell_dist_store = zeros(3,1); bias_store = zeros(3,1); search_time_store = zeros(3,1); abc_time_store = zeros(3,1);
    for ind = 3:-1:1
        if ind==1
            %adapt
            params.set_to_uniform_weights=0;
            params.set_to_scaled_weights=0;
        elseif ind==2
            %uniform
            params.set_to_uniform_weights=1;
            params.set_to_scaled_weights=0;
        elseif ind==3
            %scaled
            params.set_to_uniform_weights=0;
            params.set_to_scaled_weights=1;
        end
        [prior_comparison, bias, search_timing, abc_timing] = adapt_weights_of_ABC_KNN(params);
        hell_dist_store(ind) = prior_comparison;
        bias_store(ind) = bias;
        search_time_store(ind) = search_timing;
        abc_time_store(ind) = abc_timing;
    end
    qstr = sprintf('../%s_posterior_quality',params.save_name);
    save(qstr,'params','hell_dist_store','bias_store','search_time_store','abc_time_store'); %save info about quality of posterior
end

%plotting code
real_theta = params.theta_real;
%load data from running ABC-smc
counts = zeros(3,50);
for i=1:params.num_params
    figure; hold all;
    for j=1:3 %loop over the methods for assigning the weights
        load(sprintf('%s_posterior_plot_param%d_%d.mat',params.save_name,i,j));
        counts(j,:) = n1/sum(n1);
        stairs(bin_edges,counts(j,:),'linewidth',3);
        set(gca,'fontsize',20);
        if i==2&(~isempty(strfind(params.save_name,'death_process')))
            xlabel('$$\log_{10} \sigma $$','interpreter','latex');
        else
            xlabel(sprintf('$$\\log_{10} \\theta $$'),'interpreter','latex');
            %xlabel(sprintf('$$\\log_{10} k_%d $$',i),'interpreter','latex');
        end
        ylabel('frequency');
    end
    line([log10(real_theta(i)),log10(real_theta(i))],[0,max(max(counts))],'Color','k','LineStyle','--','linewidth',3);
    box on
    if include_stan_post
        % plot also true posterior sampled with stan
        stan_samples = csvread('../toy_model_stan.csv',1);  %('../death_process_stan_estimates.csv',1);
        n1 = histc(stan_samples(:,i),bin_edges);
        stairs(bin_edges,n1/sum(n1),'linewidth',3);
    end
    if i~=2
        loc = 'NorthWest';
    else
        loc = 'NorthEast';
    end
    legend('Adaptive','Uniform','Scaled','True','Location',loc);
    print(sprintf('%s_posterior_plot_param%d',params.save_name,i), '-depsc');
end
