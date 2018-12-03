%% Run for 3 different methods of setting weights
%%%%%%%%%%%%%%%%%
num_experiments = 40;
num_z = 6;
nnz = 2.^(6+(1:num_z)); %[128,256,512,1024,2048,4096];
mse_store = zeros(3,num_z,num_experiments);
total_sims_store = zeros(3,num_z,num_experiments);

for num_particles_ind=1:num_z
    %%%%%%%%%%%%%%%%%
    %toy_params;
    thesis_2D_gaussian_toy_params;
    %gaussian_toy_params;
    %bimodal_toy_params;
    %death_process_params;
    params.N = nnz(num_particles_ind);
    params.set_to_uniform_weights = 0; %could set method with a string
    params.set_to_scaled_weights = 0;
    params.save_name = strcat(params.save_name,'_adaptive');
    parfor expt_ind = 1:num_experiments
        loop_params = params; %cannot change object from inside the parfor loop
        loop_params.x = 13*expt_ind; %set random seed for each analysis
        [bias, total_simulations] = adaptive_ABC_KNN(loop_params);
        mse_store(1,num_particles_ind,expt_ind)=bias;
        total_sims_store(1,num_particles_ind,expt_ind)=total_simulations;
    end
    %%%%%%%%%%%%%%%%%%
    params.N = nnz(num_particles_ind);
    params.set_to_uniform_weights = 1; %could set method with a string
    params.set_to_scaled_weights = 0;
    params.save_name = strcat(params.save_name,'_uniform');
    parfor expt_ind = 1:num_experiments
        loop_params = params; %cannot change object from inside the parfor loop
        loop_params.x = 13*expt_ind; %set random seed for each analysis
        [bias, total_simulations] = adaptive_ABC_KNN(loop_params);
        mse_store(2,num_particles_ind,expt_ind)=bias;
        total_sims_store(2,num_particles_ind,expt_ind)=total_simulations;
    end
    %%%%%%%%%%%%%%%%%
    params.N = nnz(num_particles_ind);
    params.set_to_uniform_weights = 0; %could set method with a string
    params.set_to_scaled_weights = 1;
    params.save_name = strcat(params.save_name,'_mad');
    parfor expt_ind = 1:num_experiments
        loop_params = params; %cannot change object from inside the parfor loop
        loop_params.x = 13*expt_ind; %set random seed for each analysis
        [bias, total_simulations] = adaptive_ABC_KNN(loop_params);
        mse_store(3,num_particles_ind,expt_ind)=bias;
        total_sims_store(3,num_particles_ind,expt_ind)=total_simulations;
    end
    %%%%%%%%%%%%%%%%%
end
mse_store
save(strcat(params.save_name,'_mse_store.mat'),'mse_store','total_sims_store','params');
figure;
hold all;
box on;
xlabel(sprintf('Total simulations $$(\\times 10^4)$$'),'interpreter','latex');
ylabel(sprintf('$$\\log$$MSE'),'interpreter','latex');
for j=1:3
    plot(mean(total_sims_store(j,1:num_z,:),3)/10^4,log(mean(mse_store(j,1:num_z,:),3)),'LineWidth',3);
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
