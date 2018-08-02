params.N = 5*10^4; %number of ABC samples to generate at each generation
params.max_num_weights = 5000; %look at this many weights/this many generations in the gradient descent
params.num_params = 4;
params.num_ss=3*(2^5);
params.num_generations = 1;
params.alpha = 0.005;
params.proposal_sd = 0.25;
params.with_plot = 0;
params.step_size = 0.5;
params.reduce_step_size_by=2;
params.set_to_zero=0;
params.set_to_uniform_weights = 0; %could set method with a string
params.set_to_scaled_weights = 0;
params.lambda = 2*10^-4;
params.problem_t_end = 100;
params.recording_interval = params.problem_t_end/(2^5);
params.repeats = 1;
params.prior_width = 4;
params.weights_width = 6;
params.ref = [0,-1,-3,-1];  %[0,0,-prior_width/2]; %reference vector for (very slightly) more flexible prior
params.save_name = 'dimerization_v203';
params.theta_real = [1, 0.0400, 0.0020, 0.5000]; %[1.16,0.42,0.58];  %[0.1,0.5];
params.dist_metric = @(x,y) hellinger_dist(x,y);
params.optim_restarts = 500;

%setup simbio model
load('/home/harrison/Documents/Summary_stats/dimerization_model.mat');
params.test_problem = @(x,y,z) dimerization_simbio(x,y,z,model);
params.x = 777; %random seed
qq = zeros(1,params.num_ss,params.repeats);
for r=1:params.repeats
    qq(1,:,r) = params.test_problem(params.theta_real,params.problem_t_end,params.recording_interval);
end
params.data_input = qq;

