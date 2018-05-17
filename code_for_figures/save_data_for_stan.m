
saved_file = 'death_process_v213.mat';
load(saved_file);
d = reshape(params.data_input,params.num_ss,params.repeats)';
csvwrite(strcat(saved_file(1:(end-4)),'_data.csv'),d);

