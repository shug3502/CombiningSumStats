%make fig for thesis

%thesis_2D_gaussian_toy_params;
%bimodal_toy_params;
%%run ABC if required
%adaptive_ABC_KNN(params);
%load results 
load(params.save_name)
if size(theta_store,2)>1
    figure;
    hold all;
    box on;
    xlabel(sprintf('$$\\theta_1 $$'),'interpreter','latex');
    ylabel(sprintf('$$\\theta_2 $$'),'interpreter','latex');
%    xlim([-150,150]);
%    ylim([-150,150]);
    for j = 1:params.num_generations
        scat = scatter(theta_store(:,1,j),theta_store(:,2,j),12,'o','filled');
        alpha(scat,0.7);
        shg;
        pause(.2);
    end
    %% Chris' Decimal Delight
    set(gca,'FontSize',20)
    set(gca,'XTickLabel',arrayfun(@(s)sprintf('%2.0f', s), cellfun(@(s)str2num(s), get(gca,'XTickLabel')), 'UniformOutput', false))
    set(gca,'YTickLabel',arrayfun(@(s)sprintf('%2.0f', s), cellfun(@(s)str2num(s), get(gca,'YTickLabel')), 'UniformOutput', false))
    %% Ruth's Perfectly Processed PDFs
    % How to output pdfs from Matlab that are the same size etc
    set(gca,'LooseInset',get(gca,'TightInset'))
    %set(gcf,'PaperUnits','centimeters','PaperSize',[20/3+0.2 16/3+0.2], 'PaperPosition', [0.1 0.1 20/3+0.1 16/3+0.1]);
    print(sprintf('%s_smc_thesis_plot1',params.save_name), '-depsc');
end
%% 
figure;
hold all;
box on;
xlabel(sprintf('Summary statistic $$s_1 $$'),'interpreter','latex');
ylabel(sprintf('Summary statistic $$s_2 $$'),'interpreter','latex');
% xlim([-250,250]);
% ylim([-250,250]);
%xlim([-400,400]);
%ylim([-4,4]);
for j=1:params.num_generations
    scat = scatter(ss_store(:,1,j),ss_store(:,1,j),12,'o','filled');
    alpha(scat,0.7);
    shg;
    pause(.2);
end
%% Chris' Decimal Delight
set(gca,'FontSize',20)
set(gca,'XTickLabel',arrayfun(@(s)sprintf('%2.0f', s), cellfun(@(s)str2num(s), get(gca,'XTickLabel')), 'UniformOutput', false))
set(gca,'YTickLabel',arrayfun(@(s)sprintf('%2.0f', s), cellfun(@(s)str2num(s), get(gca,'YTickLabel')), 'UniformOutput', false))
%% Ruth's Perfectly Processed PDFs
% How to output pdfs from Matlab that are the same size etc
set(gca,'LooseInset',get(gca,'TightInset'))
%set(gcf,'PaperUnits','centimeters','PaperSize',[20/3+0.2 16/3+0.2], 'PaperPosition', [0.1 0.1 20/3+0.1 16/3+0.1]);
print(sprintf('%s_smc_thesis_plot2%d',params.save_name), '-depsc');
%needed this to make axes work properly: set(gcf,'Renderer','painters')
%see here: https://uk.mathworks.com/matlabcentral/answers/81353-print-option-does-not-save-the-figure-correctly