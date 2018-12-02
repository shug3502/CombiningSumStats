thesis_2D_gaussian_toy_params;
adaptive_ABC_KNN(params);
load('thesis_2D_gaussian_toy_model_v112.mat')
figure;
hold all;
for t=1:params.num_generations
    scatter(t*ones(params.N,1),theta_store(:,1,t),weights_store(:,t)*1000,'o','filled')
end
for t=2:5
    for j=1:params.N
        line([t-1,t],[theta_store(family_history_store(j,1,t),1,t-1),...
            theta_store(j,1,t)],'Color','black');
    end
end
xlabel(sprintf('Generation, $$ t $$'),'interpreter','latex');
ylabel(sprintf('$$\\theta_1 $$'),'interpreter','latex');
box on;
%% Chris' Decimal Delight
set(gca,'FontSize',20)
set(gca,'XTickLabel',arrayfun(@(s)sprintf('%2.0f', s), cellfun(@(s)str2num(s), get(gca,'XTickLabel')), 'UniformOutput', false))
set(gca,'YTickLabel',arrayfun(@(s)sprintf('%2.0f', s), cellfun(@(s)str2num(s), get(gca,'YTickLabel')), 'UniformOutput', false))
%% Ruth's Perfectly Processed PDFs
% How to output pdfs from Matlab that are the same size etc
set(gca,'LooseInset',get(gca,'TightInset'))
%set(gcf,'PaperUnits','centimeters','PaperSize',[20/3+0.2 16/3+0.2], 'PaperPosition', [0.1 0.1 20/3+0.1 16/3+0.1]);
print(sprintf('%s_smc_genealogy_plot',params.save_name), '-depsc');
