
params.N = 2000; %number of ABC samples to generate at each generation
params.num_params = 1;
params.num_ss=10;
params.num_generations = 5;
params.alpha = 0.5;
params.with_plot = 0;
params.set_to_uniform_weights = 0; %could set method with a string
params.set_to_scaled_weights = 0;
params.repeats = 1;
params.save_name = 'toy_model_v303';
params.theta_real = 10;
params.dist_metric = @(x,y) hellinger_dist(x',y'); %careful about transpose, think we need it here due to 1D parameters
params.prior_width=2;
params.ref=1;
params.draw_from_prior = @(n) 10.^draw_from_prior(n,params.prior_width,params); 
params.prior_density = @(theta) 1/(theta*log(10))*1/2*(log10(theta)>0)*(log10(theta)<2);
params.is_outside_prior = @(theta) (log10(theta)>2)*(log10(theta)<0);

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

