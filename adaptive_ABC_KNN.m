function [prior_comparison, bias, search_timing, abc_timing, ss,...
    final_samples, current_weights] = adaptive_ABC_KNN(params)
%code restructured after reviewer comments nad using KNN estimator for
%hellinger distance
% JH 19/12/2017
%based on adapt_weights_of_ABC_dist_gradient_descent and adapt_weights_abc_tidied
%ABC where we adapt the weights to weight more informative summary
%statistics higher
%with smc and with reuse
% created 19/12/2017
%last edit 27/11/2018

%To use:1) update the test problem in function at bottom
%	2) edit prior information (prior_width, params.ref)
%	3) edit real_theta : parameters for fake data
%	4) edit any other hyperparameters in params struct (params.N, params.max_num_weights,params.num_params)
%   5) if saving output for posterior plots, then change name of save file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all;
rng(params.x);
my_t = tic;
abc_timing = 0; search_timing=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Real parameters are %f \n', params.theta_real);

if isempty(params.data_input)
    ss = zeros(1,params.num_ss,params.repeats);
    for r=1:params.repeats
        ss(1,:,r) = params.test_problem(params.theta_real,params.problem_t_end,params.recording_interval);
    end
else
    ss = params.data_input;
end
if size(ss,2)~=params.num_ss
    %    fprintf('Expect %d sum stats \n',numel(ss));
    %    error('Must input correct number of summary statistics to expect.');
    params.num_ss = numel(ss);
end
if params.with_plot
    figure; plot(ss(1,:,1),'Linewidth',3)
end
M = round(params.N*params.alpha);

%initialise weights based on scale of summary stats
prior_sample = params.draw_from_prior(M);
theta_store = zeros(M,params.num_params,params.num_generations);
weights_store = zeros(M,params.num_generations); %ABC weights
ss_store = zeros(M,params.num_ss,params.num_generations);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=1;
while t<=params.num_generations
    %now generate abc samples
    if t==1
        %pick theta from prior
        theta_star = params.draw_from_prior(params.N);
        sstar=zeros(params.N,params.num_ss,params.repeats);
        for j=1:params.N
            for r=1:params.repeats
                try
                    sstar(j,:,r) = params.test_problem(theta_star(j,:),params.problem_t_end,params.recording_interval);
                catch %if errors in the solver eg for bad parameters in an ode
                    sstar(j,:,r) = zeros(1,params.num_ss);
                end
            end
        end
        %store
        weights_star = ones(params.N,1);
    else
        % pick previous particle and perturb
        sample_ind = discretesample(weights_store(:,t-1),params.N);
        theta_star = theta_store(sample_ind,:,t-1) + params.proposal_sd*randn(params.N,params.num_params);
        for j=1:params.N
            while params.is_outside_prior(theta_star(j,:)) %check if in prior support
                %resample
                sample_ind(j) = discretesample(weights_store(:,t-1),1);
                theta_star(j,:) = theta_store(sample_ind(j),:,t-1) + params.proposal_sd*randn(1,params.num_params);
            end
            for r=1:params.repeats
                sstar(j,:,r) = params.test_problem(theta_star(j,:),params.problem_t_end,params.recording_interval);
            end
        end
        %calculate a weight for all the particles, even though some
        %will not be within required distance, but which particles are
        %within this distance depends on the distance function which we
        %want to vary
        for j=1:params.N
            qq = weights_store(:,t-1).*prod(normpdf(repmat(theta_star(j,:),M,1),theta_store(:,:,t-1),params.proposal_sd),2);
            if sum(qq)<10^-8
                sum(qq)  %degenerate
            end
            weights_star(j) = params.prior_density(theta_star(j,:))/sum(qq);
        end
    end
    sig = std(reshape(sstar,params.N*params.repeats,params.num_ss));
    if params.with_plot
        figure; bar(1:2:params.num_ss,log10(1./sig(1:2:params.num_ss)),'b');
        hold on; bar(2:2:params.num_ss,log10(1./sig(2:2:params.num_ss)),'y');
        xlabel('weight index');
        ylabel('distance weight');
    end
    fprintf('Generated a bunch of samples \n');
    abc_timing = toc(my_t)-search_timing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Optimise weights
    %get new weights to use at each generation
    if t==1
        current_weights = rand(1,params.num_ss); %prior_width*rand(1,params.num_ss)-prior_width/2; %draw a weight at random to start from. Could start from ones(1,params.num_ss)
    end    
    if params.set_to_uniform_weights
        current_weights = zeros(size(current_weights));
    end
    if params.set_to_scaled_weights
        current_weights = log10(1./sig);
    end
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
optim_weights
fval
    current_weights=optim_weights;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%

    dist = weighted_distance(sstar,repmat(ss,params.N,1,1),10.^current_weights);
    [~,sort_ind] = sort(dist);
    selected = sort_ind(1:M);
    %once decided distances, then can select samples
    ss_store(:,:,t) = sstar(selected,:); %helpful for visualisation
    theta_store(:,:,t) = theta_star(selected,:);
    weights_store(:,t) = weights_star(selected)./sum(weights_star(selected)); %store normalized weights
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t=t+1;    
    fprintf('Finished optimising weights \n');
    search_timing = toc(my_t)-abc_timing
end

    for i=1:(params.num_params-1)
        figure; hold all
        plot(theta_store(:,i,t-1),theta_store(:,i+1,t-1),'o');
    end
%end
if params.num_params ==2
    [num1, centres] = hist3(prior_sample,'Nbins',[50,50]);
    num2 = hist3(theta_store(:,:,t-1),'Ctrs',centres);
else
    [num1, centres] = hist(prior_sample,30);
    num2 = hist(theta_store(:,:,t-1),centres);
end
prior_comparison = hellinger_knn_estimator(prior_sample,theta_store(:,:,t-1),5);
fprintf('Finally, information gain over the prior is: %f \n',prior_comparison);
posterior_mean = sum((theta_store(:,:,t-1).*repmat(weights_store(:,t-1),1,params.num_params)),1)/sum(weights_store(:,t-1));
bias = sqrt(sum((posterior_mean-params.theta_real).^2)); %using mean
final_samples = theta_store(:,:,t-1);
fprintf('And distance of posterior mean from real parameters is: %f \n',bias);
save(params.save_name);
toc(my_t)

figure;
hold all;
for j = 1:params.num_generations
    plot(theta_store(:,1,j),theta_store(:,1,j),'o');
    shg;
    pause(2);
end
figure;
hold all;
for j=1:params.num_generations
    plot(ss_store(:,9,j),ss_store(:,10,j),'o');
    shg;
    pause(2);
end

