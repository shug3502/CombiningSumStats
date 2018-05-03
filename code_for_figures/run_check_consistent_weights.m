%% last edit 6/9/16
%% created JH 6/9/16
%% generate plots for consistent weights
%% this will currently run all the simulations
close all;
addpath ../

%death process
death_process_params;
check_consistent_weights(params);

%%birth death
%birth_death_params;
%check_consistent_weights(params);

%%dimerization
%dimerization_params;
%check_consistent_weights(params);

%%spatial problem
%my_params_store;
%check_consistent_weights(params);
