params.N = 2*10^3; %number of ABC samples to generate at each generation
params.num_params = 2;
params.num_ss = 2;
params.num_generations = 5;
params.alpha = 0.5;
params.with_plot = 0;
params.set_to_uniform_weights = 0; %could set method with a string
params.set_to_scaled_weights = 0;
params.repeats = 1;
params.save_name = 'bimodal_toy_model_v212';
params.theta_real = [pi/2, 3*pi/2]; %[pi/4, 5*pi/4];

params.verbose=0;
params.problem_t_end=0;
params.weights_width = 8;
params.recording_interval=0;
params.optim_restarts = 100;

params.draw_from_prior = @(n) 2*pi*rand(n,params.num_params);
params.prior_density = @(theta) all((theta>0) & (theta<2*pi))/(2*pi);
params.is_outside_prior = @(theta) all((theta<0) | (theta>2*pi));
params.test_problem = @(x,y,z) bimodal_toy_model(x,y,z);
params.x = 777; %random seed
qq = zeros(1,params.num_ss,params.repeats);
qq(1,:,1) = [sin(pi/4),-sin(pi/4)];%params.test_problem(params.theta_real,0,0);
params.data_input = qq;
