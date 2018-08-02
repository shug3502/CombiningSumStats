function make_figure_posterior_stairs_plot(params,run_simulations,include_stan_post)
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
%run_simulations = 1;
%include_stan_post = 0;

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
        [prior_comparison, bias, search_timing, abc_timing, ss,...
    final_samples, optim_weights] = adapt_weights_of_ABC_KNN(params);
        hell_dist_store(ind) = prior_comparison;
        bias_store(ind) = bias;
        search_time_store(ind) = search_timing;
        abc_time_store(ind) = abc_timing;
    end
    qstr = sprintf('../%s_posterior_quality',params.save_name);
    save(qstr,'params','hell_dist_store','bias_store','search_time_store','abc_time_store'); %save info about quality of posterior
end
%%%%%%%%%%%%%%%%%%%%%%%
% New plotting code
real_theta = params.theta_real;
%load data from running ABC-smc
color_matrix =  [230,97,1;
253,184,99;
178,171,210;
94,60,153]/256;
for i=1:params.num_params
    figure; hold all;
    for j=3:-1:1 %loop over the methods for assigning the weights
        load(sprintf('%s_posterior_plot_param%d_%d.mat',params.save_name,i,j));
        samples = theta_store(:,i,params.num_generations);
        [f,xi] = ksdensity(samples); 
        plot(xi,f,'LineWidth',3, 'Color', color_matrix(4-j,:));
        set(gca,'fontsize',20);
        if i==2&(~isempty(strfind(params.save_name,'death_process')))
            xlabel('$$\log_{10} \sigma $$','interpreter','latex');
        else
            xlabel(sprintf('$$\\log_{10} \\theta $$'),'interpreter','latex');
            %xlabel(sprintf('$$\\log_{10} k_%d $$',i),'interpreter','latex');
        end
        ylabel('frequency');
    end
    if include_stan_post
        % plot also true posterior sampled with stan
        if ~isempty(strfind(params.save_name,'death_process'))
            stan_samples = csvread('stan/death_process_stan_estimates.csv',1);            
        elseif ~isempty(strfind(params.save_name,'toy_model'))
            stan_samples = csvread('stan/toy_model_stan.csv',1);
        else 
            error('unknown model: no stan fit available');
        end
        [f,xi] = ksdensity(stan_samples(:,i)); 
        plot(xi,f,'LineWidth',3, 'Color', color_matrix(4,:));
    end
    line([log10(real_theta(i)),log10(real_theta(i))],[0,max(max(f))],'Color','k','LineStyle','--','linewidth',3);
    box on    
    if ~isempty(strfind(params.save_name,'death_process')) & i==2
        xlim([-6.0,5.0])
    end
    if i~=2
        loc = 'NorthWest';
    else
        loc = 'NorthEast';
    end
    %legend('Adaptive','Uniform','Scaled','True','MCMC','Location',loc);
    if include_stan_post
        legend('Scaled','Uniform','Adaptive','MCMC','True','Location',loc);
    else
        legend('Scaled','Uniform','Adaptive','True','Location',loc);
    end

%% Chris' Decimal Delight
set(gca,'XTickLabel',arrayfun(@(s)sprintf('%0.1f', s), cellfun(@(s)str2num(s), get(gca,'XTickLabel')), 'UniformOutput', false))
set(gca,'YTickLabel',arrayfun(@(s)sprintf('%0.1f', s), cellfun(@(s)str2num(s), get(gca,'YTickLabel')), 'UniformOutput', false))
%% Ruth's Perfectly Processed PDFs
% How to output pdfs from Matlab that are the same size etc
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',20)
%set(gcf,'PaperUnits','centimeters','PaperSize',[20/3+0.2 16/3+0.2], 'PaperPosition', [0.1 0.1 20/3+0.1 16/3+0.1]);

    print(sprintf('%s_posterior_plot_param%d',params.save_name,i), '-depsc');
end


% 
% %plotting code
% real_theta = params.theta_real;
% %load data from running ABC-smc
% counts = zeros(3,50);
% for i=1:params.num_params
%     figure; hold all;
%     for j=1:3 %loop over the methods for assigning the weights
%         load(sprintf('%s_posterior_plot_param%d_%d.mat',params.save_name,i,j));
%         counts(j,:) = n1/sum(n1);
%         stairs(bin_edges,counts(j,:),'linewidth',3);
%         set(gca,'fontsize',20);
%         if i==2&(~isempty(strfind(params.save_name,'death_process')))
%             xlabel('$$\log_{10} \sigma $$','interpreter','latex');
%         else
%             xlabel(sprintf('$$\\log_{10} \\theta $$'),'interpreter','latex');
%             %xlabel(sprintf('$$\\log_{10} k_%d $$',i),'interpreter','latex');
%         end
%         ylabel('frequency');
%     end
%     line([log10(real_theta(i)),log10(real_theta(i))],[0,max(max(counts))],'Color','k','LineStyle','--','linewidth',3);
%     box on
%     if include_stan_post
%         % plot also true posterior sampled with stan
%         stan_samples = csvread('../toy_model_stan.csv',1);  %('../death_process_stan_estimates.csv',1);
%         n1 = histc(stan_samples(:,i),bin_edges);
%         stairs(bin_edges,n1/sum(n1),'linewidth',3);
%     end
%     if i~=2
%         loc = 'NorthWest';
%     else
%         loc = 'NorthEast';
%     end
%     legend('Adaptive','Uniform','Scaled','True','Location',loc);
%     print(sprintf('%s_posterior_plot_param%d',params.save_name,i), '-depsc');
% end
