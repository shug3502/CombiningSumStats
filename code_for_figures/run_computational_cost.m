%% created 07/09/16 JH
%% last edit 7/9/16
%% run code to compare computational cost of more abc samples and automatic algorithm

close all;
clear all;
addpath ../;
addpath ../computation_cost_files/;

%toy_params;
fprintf('for the toy model \n');
load('toy_model_v206_computational_cost.mat');
computational_cost(params,0,1,5*10^5);
hell_dist_store
bias_store

%death_process_params;
fprintf('for the death process \n');
load('death_process_v216_comp_cost_computational_cost.mat');
computational_cost(params,0,1,5*10^5);
hell_dist_store
bias_store

%dimerization_params;
fprintf('for the dimerization model \n');
load('dimerization_v204_comp_cost_computational_cost.mat');
computational_cost(params,0,0,5*10^4);
hell_dist_store
bias_store

%my_params_store;
fprintf('for the spatial model\n');
load('spatial_test_v205_comp_cost_computational_cost.mat');
computational_cost(params,0,0,5*10^4);
hell_dist_store
bias_store

