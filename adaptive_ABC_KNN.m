function [bias, total_simulations, ...
    final_samples] = adaptive_ABC_KNN(params)
%code restructured after reviewer comments nad using KNN estimator for
%hellinger distance
% JH 19/12/2017
%based on adapt_weights_of_ABC_dist_gradient_descent and adapt_weights_abc_tidied
%ABC where we adapt the weights to weight more informative summary
%statistics higher
%with smc and with reuse
% created 19/12/2017
%last edit 27/11/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all;
rng(params.x);
my_t = tic;
abc_timing = 0; search_timing=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.verbose
fprintf('Real parameters are %f \n', params.theta_real);
end

if isempty(params.data_input)
    ss = zeros(1,params.num_ss,params.repeats);
    for r=1:params.repeats
        ss(1,:,r) = params.test_problem(params.theta_real,params.problem_t_end,params.recording_interval);
    end
else
    ss = params.data_input;
end
if size(ss,2)~=params.num_ss
        fprintf('Expect %d sum stats \n',numel(ss));
        error('Must input correct number of summary statistics to expect.');
    params.num_ss = numel(ss);
end
if params.with_plot
    figure; plot(ss(1,:,1),'Linewidth',3)
end
%M = round(params.N*params.alpha);
M = ceil(params.N/params.alpha);

%initialise weights based on scale of summary stats
prior_sample = params.draw_from_prior(params.N);
theta_store = zeros(params.N,params.num_params,params.num_generations);
weights_store = zeros(params.N,params.num_generations); %ABC weights
ss_store = zeros(params.N,params.num_ss,params.num_generations);
distance_weights_store = zeros(1,params.num_ss,params.num_generations);
family_history_store = zeros(params.N,1,params.num_generations); %for visualisation of particle geanolgies
total_simulations = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tolerances = zeros(params.num_generations,1);
tolerances(1) = Inf; %accept all at first generation, initially

t=1;
while t<=params.num_generations
    %    while num_acceptances<M %do stuff
    %now generate abc samples
    if t==1
        %pick theta from prior
        theta_star = params.draw_from_prior(M);
        sstar=zeros(M,params.num_ss,params.repeats);
        family_history_star = zeros(M,1);
        num_acceptances=0;
        while num_acceptances<M
            %pick theta from prior
            for r=1:params.repeats
                try
                    sstar(num_acceptances+1,:,r) = params.test_problem(theta_star(num_acceptances+1,:),params.problem_t_end,params.recording_interval);
                catch %if errors in the solver eg for bad parameters in an ode
                    sstar(num_acceptances+1,:,r) = zeros(1,params.num_ss);
                end
            end
            total_simulations = total_simulations+1;
            %accept all
            num_acceptances = num_acceptances+1;
        end
    else
        % pick previous particle and perturb
        num_acceptances=0;
        while num_acceptances<M
            sample_ind = discretesample(weights_store(:,t-1),1);
            empirical_variance = weightedcov(theta_store(:,:,t-1),weights_store(:,t-1));
            th_star = mvnrnd(theta_store(sample_ind,:,t-1),2*empirical_variance);
            while params.is_outside_prior(th_star) %check if in prior support
                %resample
                sample_ind = discretesample(weights_store(:,t-1),1);
                th_star =  mvnrnd(theta_store(sample_ind,:,t-1),2*empirical_variance);
            end
            sstar_temp=zeros(1,params.num_ss,params.repeats);
            for r=1:params.repeats
                sstar_temp(1,:,r) = params.test_problem(th_star,params.problem_t_end,params.recording_interval);
            end
            total_simulations = total_simulations+1;
            %%%%%%%%%%%%%%%%%%%
            %calculate previous distances
            dist = zeros(t-1,1);
            for tt=1:(t-1)
                dist(tt) = weighted_distance(sstar_temp,repmat(ss,1,1,1),10.^distance_weights_store(1,:,tt));
            end
            %accept if close enough
            accept_check = all(dist<tolerances(1:t-1));
            num_acceptances = num_acceptances + accept_check;
            if accept_check
                theta_star(num_acceptances,:) = th_star;
                sstar(num_acceptances,:) = sstar_temp;
                family_history_star(num_acceptances) = sample_ind;
            end
        end
    end
    sig = mad(my_reshape(sstar),1);
    if params.set_to_uniform_weights
        distance_weights_store(1,:,t) = zeros(1,params.num_ss);
    end
    if params.set_to_scaled_weights
        distance_weights_store(1,:,t) = log10(1./sig);
    end
    %%%%%%%%%%%%%%%%%%%%
    if params.verbose
        fprintf('Generated a bunch of samples \n');
        abc_timing = toc(my_t)-search_timing
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Optimise weights
    %get new weights to use at each generation
    if (params.set_to_scaled_weights+params.set_to_uniform_weights<=0)
        %try to solve with inbuilt optimizers
        fun = @(x) get_distance_from_prior_to_post_KNN(x,sstar,ss,theta_star,t,prior_sample,params);
        opts = optimset('Display', 'off');
        [optim_weights,fval] = fmincon(fun,params.weights_width*(rand(1,params.num_ss)),[],[],[],[],-params.weights_width*ones(1,params.num_ss),params.weights_width*ones(1,params.num_ss),[],opts);
        for starts=1:params.optim_restarts %run optimizer from multiple random starting pts
            [optim_weights_alt,fval_alt] = fmincon(fun,params.weights_width*(rand(1,params.num_ss)),[],[],[],[],-params.weights_width*ones(1,params.num_ss),params.weights_width*ones(1,params.num_ss),[],opts);
            [fval,is_improvement] = max([fval,fval_alt]);
            optim_weights = optim_weights_alt*(is_improvement-1) + optim_weights*(2-is_improvement);
        end
%        optim_weights
%        fval
        distance_weights_store(1,:,t) = optim_weights;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate a weight for all the particles, even though some
    %will not be within required distance, but which particles are
    %within this distance depends on the distance function which we
    %want to vary
    dist = weighted_distance(sstar,repmat(ss,M,1,1),10.^distance_weights_store(1,:,t));
    [~,sort_ind] = sort(dist);
    tolerances(t) = dist(sort_ind(params.N));
    selected = sort_ind(1:params.N);
    %once decided distances, then can select samples
    ss_store(:,:,t) = sstar(selected,:); %helpful for visualisation
    theta_store(:,:,t) = theta_star(selected,:);
    family_history_store(:,1,t) = family_history_star(selected);
    if t==1
        %store
        weights_store(:,t) = ones(params.N,1)/params.N;
    else
        for j=1:params.N
            %TODO: check
            qq = weights_store(:,t-1).*mvnpdf(repmat(theta_store(j,:,t),params.N,1),theta_store(:,:,t-1),2*empirical_variance);
            % prod(normpdf(repmat(theta_store(j,:,t),params.N,1),theta_store(:,:,t-1),params.proposal_sd),2);
            if sum(qq)<10^-16
                sum(qq)  %degenerate
                error('something went wrong, degenerate weights')
            end
            weights_store(j,t) = params.prior_density(theta_store(j,:,t))/sum(qq);
        end
        if sum(weights_store(:,t))<10^-8
            error('weights are all tiny')
        end
        weights_store(:,t) = weights_store(:,t)./sum(weights_store(:,t)); %store normalized weights
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t=t+1;
    if params.verbose
        fprintf('Finished optimising weights \n');
        search_timing = toc(my_t)-abc_timing
    end
end
if params.with_plot
    for i=1:(params.num_params-1)
        figure; hold all
        plot(theta_store(:,i,t-1),theta_store(:,i+1,t-1),'o');
    end
end

prior_comparison = hellinger_knn_estimator(prior_sample,theta_store(:,:,t-1),5);
if params.verbose
fprintf('Finally, information gain over the prior is: %f \n',prior_comparison);
end
posterior_mean = sum((theta_store(:,:,t-1).*repmat(weights_store(:,t-1),1,params.num_params)),1)/sum(weights_store(:,t-1));
if isempty(strfind(params.save_name,'toy_model'))
    bias = sqrt(sum((posterior_mean-params.theta_real).^2)); %using mean
else
    stan_mean = mean(10.^csvread('stan/toy_model_stan.csv',1));
    bias = sqrt(sum((posterior_mean-stan_mean).^2));
end
final_samples = theta_store(:,:,t-1);
if params.verbose
fprintf('And distance of posterior mean from real parameters is: %f \n',bias);
end
save(params.save_name);
toc(my_t)

