function computational_cost(params,run_simulations,include_stan_post,n1)
%% created 25/8/16 JH
%% last edit 25/8/16
%%
%% Look at computational cost of search through weights vs generating more ABC samples
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
addpath ../ %add path to Summary_stats directory

%n1=5*10^3;
%what value should n2 take?
%n2=10^2;
n2=10^2;

if nargin < 1
    %%load some parameters
    %death_process_params;
    %my_params_store;
    dimerization_params;
    n1 = 5*10^4;
    run_simulations=1; %logical to run computations or just plot
    include_stan_post=0;
end

if run_simulations
    hell_dist_store = zeros(3,1); bias_store = zeros(3,1); search_time_store = zeros(3,1); abc_time_store = zeros(3,1);
    for loop_ind = 1:3
        if loop_ind==1
            %adapt
            params.N = n1; %set to 'n1'
            params.set_to_uniform_weights=0;
            params.set_to_scaled_weights=0;
        elseif loop_ind==2
            %uniform
            params.N = n1; %set to 'n1'
            params.set_to_uniform_weights=1;
            params.set_to_scaled_weights=0;
        elseif loop_ind==3
            %uniform
            params.N = n2; %set to 'n2'
            params.set_to_uniform_weights=1;
            params.set_to_scaled_weights=0;
            params.alpha = alpha2;
        end
        [prior_comparison, bias, search_timing, abc_timing] = adapt_weights_of_ABC_KNN(params);
        hell_dist_store(loop_ind) = prior_comparison;
        bias_store(loop_ind) = bias;
        search_time_store(loop_ind) = search_timing;
        abc_time_store(loop_ind) = abc_timing;
        if loop_ind==1
            abc_time_for_one_sample = abc_timing/n1;
            n2 = n1 + ceil(search_timing/abc_time_for_one_sample)
            alpha2 = params.alpha*n1/n2;
        elseif loop_ind==2
            %easiest to move files so can be plotted properly
            for jj=1:params.num_params
                cmd=sprintf('mv %s_posterior_plot_param%d_2.mat %s_posterior_plot_param%d_3.mat',params.save_name, jj, params.save_name, jj);
                status = system(cmd);
                if status>0
                    error('system commands not working');
                end
            end
        end
    end
    qstr = sprintf('../%s_computational_cost',params.save_name);
    save(qstr,'params','hell_dist_store','bias_store','search_time_store','abc_time_store'); %save info about quality of posterior
    hell_dist_store
    bias_store
    search_time_store
    abc_time_store
end


% %plotting code
% real_theta = params.theta_real;
% %load data from running ABC-smc
% counts = zeros(3,50);
% for i=1:params.num_params
%     figure; hold all;
%     for j=1:3 %loop over the methods for assigning the weights
%         load(sprintf('%s_posterior_plot_param%d_%d.mat',params.save_name,i,j));
%         counts(j,:) = n1;
%         %sum(bin_edges.*n1'/sum(n1))
%         stairs(bin_edges,counts(j,:),'linewidth',3); set(gca,'fontsize',24);
%         xlabel('\theta'); ylabel('frequency');
%     end
%     line([log10(real_theta(i)),log10(real_theta(i))],[0,max(max(counts))],'Color','k','LineStyle','--','linewidth',3);
%     legend('Adapt N1','Uniform N2','Uniform N1');
%     print(sprintf('%s_computational_cost_plot_param%d',params.save_name,i), '-depsc');
% end

% New plotting code
real_theta = params.theta_real;
%load data from running ABC-smc
for i=1:params.num_params
    figure; hold all;
    if include_stan_post
        % plot also true posterior sampled with stan
        if ~isempty(strfind(params.save_name,'death_process'))
            stan_samples = csvread('../stan/death_process_stan_estimates.csv',1);            
        elseif ~isempty(strfind(params.save_name,'toy_model'))
            stan_samples = csvread('../stan/toy_model_stan.csv',1);
        else 
            error('unknown model: no stan fit available');
        end
        [f,xi] = ksdensity(stan_samples(:,i)); 
        plot(xi,f,'LineWidth',3);
    end
    for j=1:3 %loop over the methods for assigning the weights
        load(sprintf('%s_posterior_plot_param%d_%d.mat',params.save_name,i,j));
        samples = theta_store(:,i,params.num_generations);
        [f,xi] = ksdensity(samples); 
        plot(xi,f,'LineWidth',3);
        set(gca,'fontsize',20);
        if i==2&(~isempty(strfind(params.save_name,'death_process')))
            xlabel('$$\log_{10} \sigma $$','interpreter','latex');
        else
            xlabel(sprintf('$$\\log_{10} \\theta $$'),'interpreter','latex');
            %xlabel(sprintf('$$\\log_{10} k_%d $$',i),'interpreter','latex');
        end
        ylabel('frequency');
    end
    line([log10(real_theta(i)),log10(real_theta(i))],[0,max(max(f))],'Color','k','LineStyle','--','linewidth',3);
    box on
    if i~=2
        loc = 'NorthWest';
    else
        loc = 'NorthEast';
    end
    legend('Adapt N1','Uniform N2','Uniform N1','Location',loc);
    box on;
    print(sprintf('%s_computational_cost_plot_param%d',params.save_name,i), '-depsc');
end
