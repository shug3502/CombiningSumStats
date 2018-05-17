
params.N = 5*10^4; %number of ABC samples to generate at each generation
params.max_num_weights = 5000; %look at this many weights/this many generations in the gradient descent
params.num_params = 1;
params.num_ss=10;
params.num_generations = 5;
params.alpha = 0.005;
params.proposal_sd = 0.25;
params.with_plot = 0;
params.step_size = 0.5;
params.reduce_step_size_by=2;
params.set_to_zero=0;
params.set_to_uniform_weights = 0; %could set method with a string
params.set_to_scaled_weights = 0;
params.lambda = 0;
%params.problem_t_end = 20;
%params.recording_interval = params.problem_t_end/(2^3);
params.repeats = 1;
params.prior_width = 2;
params.ref = [1];  %[0,0,-prior_width/2]; %reference vector for (very slightly) more flexible prior
params.save_name = 'toy_model_v205';
params.theta_real = 10;
params.dist_metric = @(x,y) hellinger_dist(x',y'); %careful about transpose, think we need it here due to 1D parameters

params.problem_t_end=0;
params.weights_width = 8;
params.recording_interval=0;
params.optim_restarts = 30;

params.test_problem = @(x,y,z) toy_model(x,y,z);
params.x = 777; %random seed
qq = zeros(1,params.num_ss,params.repeats);
for r=1:params.repeats
    qq(1,:,r) = params.test_problem(params.theta_real,0,0);
end
params.data_input = qq;

