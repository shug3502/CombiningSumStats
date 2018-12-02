%% Run for 3 different methods of setting weights
%%%%%%%%%%%%%%%%%
num_z = 10;
nnz = linspace(0,5000,num_z+1);
mse_store = zeros(3,num_z);
total_sims_store = zeros(3,num_z);

for num_particles=1:2:num_z
%%%%%%%%%%%%%%%%%
toy_params;
%thesis_2D_gaussian_toy_params;
%gaussian_toy_params;
params.N = nnz(1+num_particles);
params.set_to_uniform_weights = 0; %could set method with a string
params.set_to_scaled_weights = 0;
params.save_name = strcat(params.save_name,'_adaptive');
generate_thesis_smc_fig;
posterior_mean
mse_store(1,num_particles)=bias;
total_sims_store(1,num_particles)=total_simulations;
%%%%%%%%%%%%%%%%%%
params.N = nnz(1+num_particles);
params.set_to_uniform_weights = 1; %could set method with a string
params.set_to_scaled_weights = 0;
params.save_name = strcat(params.save_name,'_uniform');
generate_thesis_smc_fig;
posterior_mean
mse_store(2,num_particles)=bias;
total_sims_store(2,num_particles)=total_simulations;
%%%%%%%%%%%%%%%%%
params.N = nnz(1+num_particles);
params.set_to_uniform_weights = 0; %could set method with a string
params.set_to_scaled_weights = 1;
params.save_name = strcat(params.save_name,'_mad');
generate_thesis_smc_fig;
posterior_mean
mse_store(3,num_particles)=bias;
total_sims_store(3,num_particles)=total_simulations;
%%%%%%%%%%%%%%%%%
end
mse_store
figure; 
hold all;
box on;
xlabel(sprintf('Total simulations'),'interpreter','latex');
ylabel(sprintf('$$\\log$$MSE'),'interpreter','latex');
for j=1:3
    plot(total_sims_store(j,1:2:num_z),log(mse_store(j,1:2:num_z)),'LineWidth',3);
end
legend('adaptive','uniform','scaled');
%% Chris' Decimal Delight
set(gca,'FontSize',20)
set(gca,'XTickLabel',arrayfun(@(s)sprintf('%2.0f', s), cellfun(@(s)str2num(s), get(gca,'XTickLabel')), 'UniformOutput', false))
set(gca,'YTickLabel',arrayfun(@(s)sprintf('%2.0f', s), cellfun(@(s)str2num(s), get(gca,'YTickLabel')), 'UniformOutput', false))
%% Ruth's Perfectly Processed PDFs
% How to output pdfs from Matlab that are the same size etc
set(gca,'LooseInset',get(gca,'TightInset'))
%set(gcf,'PaperUnits','centimeters','PaperSize',[20/3+0.2 16/3+0.2], 'PaperPosition', [0.1 0.1 20/3+0.1 16/3+0.1]);
print(sprintf('%s_MSE',params.save_name), '-depsc');
