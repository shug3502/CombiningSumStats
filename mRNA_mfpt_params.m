addpath ../fine-grained/mRNA/;

params.N = 2*10^2; %number of ABC samples to generate at each generation
params.num_params = 6;
params.num_ss=1;
params.num_generations = 5;
params.alpha = 0.5;
params.set_to_uniform_weights = 0; %could set method with a string
params.set_to_scaled_weights = 0;
params.repeats = 1;
params.save_name = 'mRNA_mfpt_model_v102';
params.theta_real = [1.16,0.8,0.11,0.42,0.84,0.58];
params.dist_metric = @(x,y) hellinger_dist(x',y'); %careful about transpose, think we need it here due to 1D parameters

params.with_plot=0;
params.problem_t_end=0;
params.weights_width = 8;
params.recording_interval=0;
params.optim_restarts = 30;

%add restriction that theta_1/theta_4 < 20
params.draw_from_prior = @(n) draw_from_prior_mRNA(n,params.num_params); 
params.prior_density = @(theta) mvnpdf(theta(1:(params.num_params-1)))*2^6*(theta(params.num_params)>0.5)*(theta(params.num_params)<1)*(all(theta)>0)*(atan2(theta(4),theta(1))>1/20)*20/19;
params.is_outside_prior = @(theta) (theta(params.num_params)<0.5) | (theta(params.num_params)>1) | any(theta<0) | (atan2(theta(4),theta(1))<1/20);
params.test_problem = @(x,y,z) mRNA_mfpt_model(x,y,z);
params.x = 777; %random seed
qq = zeros(1,params.num_ss,params.repeats);
qq = params.test_problem(params.theta_real,0,0);
params.data_input = qq;
