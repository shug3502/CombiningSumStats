%created 6/9/16
%JH

%script to make figures with posteriors for different methods, using stairs
%to create run: adapt_weights_of_ABC_tidied(params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
cd ../
addpath code_for_figures/ %make use can access functions in Summary_stats directory

%death_process_params;
%make_figure_posterior_stairs_plot(params,1,1);

%dimerization_params;
%make_figure_posterior_stairs_plot(params,1,0);

my_params_store; %load params
make_figure_posterior_stairs_plot(params,1,0);
