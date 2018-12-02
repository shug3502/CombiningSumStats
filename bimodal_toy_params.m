params.N = 2*10^3; %number of ABC samples to generate at each generation
params.num_params = 2;
params.num_ss = 2;
params.num_generations = 5;
params.alpha = 0.5;
params.with_plot = 0;
params.step_size = 0.5;
params.set_to_uniform_weights = 0; %could set method with a string
params.set_to_scaled_weights = 0;
params.repeats = 1;
params.save_name = 'bimodal_toy_model_v103';
params.theta_real = 0;
params.dist_metric = @(x,y) hellinger_dist(x',y'); %careful about transpose, think we need it here due to 1D parameters

params.problem_t_end=0;
params.weights_width = 8;
params.recording_interval=0;
params.optim_restarts = 30;

params.draw_from_prior = @(n) 2*pi*rand(n,params.num_params); %randn(n,params.num_params)*100; 
params.prior_density = @(theta) all((theta>0) & (theta<2*pi))/(2*pi); %mvnpdf(theta,zeros(size(theta)),100*ones(size(theta)));
params.is_outside_prior = @(theta) all((theta<0) | (theta>2*pi));
params.test_problem = @(x,y,z) bimodal_toy_model(x,y,z);
params.x = 777; %random seed
qq = zeros(1,params.num_ss,params.repeats);
qq(1,:,1) = [0.7,-0.7];%params.test_problem(params.theta_real,0,0);
params.data_input = qq;
