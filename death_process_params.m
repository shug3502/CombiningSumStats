
params.N = 5*10^5; %number of ABC samples to generate at each generation
params.max_num_weights = 5000; %look at this many weights/this many generations in the gradient descent
params.num_params = 2;
params.num_ss=2^3+1;
params.num_generations = 1;
params.alpha = 0.005;
params.proposal_sd = 0.25;
params.with_plot = 0;
params.step_size = 0.5;
params.reduce_step_size_by=2;
params.set_to_zero=0;
params.set_to_uniform_weights = 0; %could set method with a string
params.set_to_scaled_weights = 0;
params.optim_restarts = 500;
params.problem_t_end = 20;
params.recording_interval = params.problem_t_end/(params.num_ss-1);
params.repeats = 1;
params.prior_width = 8;
params.weights_width = 4;
params.ref = [0,0];  %[0,0,-prior_width/2]; %reference vector for (very slightly) more flexible prior
params.save_name = 'death_process_v215';
params.theta_real = [0.1,0.01]; %[1.16,0.84,0.58]; %[1, 0.0400, 0.0020, 0.5000]; %[1.16,0.42,0.58];  %[0.1,0.5];
params.dist_metric = @(x,y) hellinger_dist2D(x,y); %need to take transpose if only 1D ie. one parameter

params.test_problem = @(x,y,z) test_problem2(x,y,z);
params.x = 777; %random seed
qq = zeros(1,params.num_ss,params.repeats);
for r=1:params.repeats
    qq(1,:,r) = params.test_problem(params.theta_real,params.problem_t_end,params.recording_interval);
end
params.data_input = qq;

