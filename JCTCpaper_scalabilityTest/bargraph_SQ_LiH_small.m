%% current part code
load('times_bargraph.mat');
fig = figure('Units', 'inches','position',[0 0 3.3 3.3]);
% h = bar(times111, 'stacked');
h = bar(times221, 'stacked');
% h = bar(times222, 'stacked');
h(1).FaceColor = [0.6350 0.0780 0.1840];
h(2).FaceColor = [0 0.4470 0.7410];
h(3).FaceColor = [0 1 0.7410];
h(4).FaceColor = [0.9290 0.6940 0.1250];
xlim([0.3, 5.7]);
% ylim([0,8]);
ylim([0,300]);
% ylim([0,1200]);
% xTickLabels = {'72', '144', '288', '576','1152'};
xTickLabels = {'144', '288', '576','1152','2304'};
% xTickLabels = {'288', '576','1152','2304','4608'};
xlabel('Number of processors');
ylabel('Wall time (s)');
% yticks([0, 2, 4, 6, 8]);
% yticklabels({'      0', '      2', '      4', '      6','      8'});
yticks([60, 120, 180, 240, 300]);
yticklabels({'    60', '   120', '   180', '   240','   300'});
% yticks([240 480 720 960 1200]);
% yTickLabels = {' 240', ' 480', '720', ' 960', '1200'};


ax = gca;
set(ax,'FontName','Times New Roman','FontSize',12,'XTickLabel', xTickLabels, 'XTick', 1:5);
box(ax,'on');
% legend1 = legend('Kronecker product','Sternheimer equation', 'Projection & orthonormalization','Spectral quadrature', 'Location', 'northeast');
% set(legend1,'interpreter','Latex','FontSize',12);