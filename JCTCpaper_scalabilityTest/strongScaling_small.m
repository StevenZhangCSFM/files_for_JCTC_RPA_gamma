% % Create figure
% figure('OuterPosition',...
%     [353 183.666666666667 574.666666666667 507.333333333333]);

% Create axes
load('scaling_data.mat');
fig = figure('Units', 'inches','OuterPosition',[0 0 3.33 3.83]);
hold on;

% Create loglog
loglog(NP111,time111,'DisplayName','(LiH)$_4$','MarkerSize',5,'Marker','o',...
    'LineWidth',1.5);
hold on;
% Create loglog
loglog(NP221,time221,'DisplayName','(LiH)$_{16}$','MarkerSize',5,'Marker','square',...
    'LineWidth',1.5,'LineStyle','-.');
hold on;
% Create loglog
loglog(NP222,time222,'DisplayName','(LiH)$_{32}$','MarkerSize',5,'Marker','diamond',...
    'LineWidth',1.5,'LineStyle',':', 'Color', '[0.4660 0.6740 0.1880]');
hold on;
% Create ylabel
ylabel({'Wall time (s)'},'Interpreter','latex');
% Create xlabel
xlabel({'Number of processors'}, 'Interpreter','latex');

ax = gca;
box(ax,'on');

% Set the remaining axes properties
set(ax,'FontName','Times New Roman','FontSize',14,'XMinorTick','on','XScale','log','YMinorTick','on',...
    'YScale','log');
xticks([10 100 1000 10000]);
yticks([1 10 100 1000]);
xlim([36 20000]);
ylim([0.2 4000]);

% Create legend
legend1 = legend(ax,'show');
set(legend1,'Location','southeast','Interpreter','latex','FontSize',12);
legend1.AutoUpdate = 'off';
hold on;
loglog(NP111,strongScaling111,NP221,strongScaling221,NP222,strongScaling222,...
    'LineWidth',1,'Color','black','LineStyle','--');

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperPositionMode', 'auto');
fig_pos = get(fig, 'Position');
set(fig, 'PaperSize', fig_pos(3:4));