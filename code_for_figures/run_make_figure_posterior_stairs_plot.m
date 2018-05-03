%created 6/9/16
%JH

%script to make figures with posteriors for different methods, using stairs
%to create run: adapt_weights_of_ABC_tidied(params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
addpath ../ %make use can access functions in Summary_stats directory

%death_process_params;
%make_figure_posterior_stairs_plot(params);

%birth_death_params;
%make_figure_posterior_stairs_plot(params);

dimerization_params;
make_figure_posterior_stairs_plot(params);

%my_params_store; %load params
%make_figure_posterior_stairs_plot(params);
