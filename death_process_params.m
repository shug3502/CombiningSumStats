params.N = 2000; %number of ABC samples to generate at each generation
params.num_params = 2;
params.num_ss=2^3+1;
params.num_generations = 5;
params.alpha = 0.5;
params.proposal_sd = 0.25;
params.with_plot = 0;
params.set_to_uniform_weights = 0; %could set method with a string
params.set_to_scaled_weights = 1;
params.optim_restarts = 100;
params.problem_t_end = 20;
params.recording_interval = params.problem_t_end/(params.num_ss-1);
params.repeats = 1;
params.prior_width = 8;
params.weights_width = 4;
params.ref = [0,0];  %[0,0,-prior_width/2]; %reference vector for (very slightly) more flexible prior
params.save_name = 'death_process_v312';
params.verbose=1;
params.theta_real = [-1,-2];

params.draw_from_prior = @(n) draw_from_prior(n,params.prior_width,params); 
params.prior_density = @(theta) 1/(params.prior_width)^params.num_params*all((theta > params.ref-params.prior_width/2))*all((theta < params.ref+params.prior_width/2));
params.is_outside_prior = @(theta) any((theta < params.ref-params.prior_width/2)) | any((theta > params.ref+params.prior_width/2));

params.test_problem = @(x,y,z) test_problem2(x,y,z);
params.x = 777; %random seed
qq = zeros(1,params.num_ss,params.repeats);
for r=1:params.repeats
    qq(1,:,r) = params.test_problem(params.theta_real,params.problem_t_end,params.recording_interval);
end
params.data_input = qq;

