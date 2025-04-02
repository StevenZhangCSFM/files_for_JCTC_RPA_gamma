% % Create figure
% figure('OuterPosition',...
%     [353 183.666666666667 574.666666666667 507.333333333333]);

% Create axes
axes1 = axes('Position',[0.15 0.15 0.8 0.8]);
hold(axes1,'on');

% Create loglog
loglog(NP111,time111,'DisplayName','(LiH)$_4$','MarkerSize',9,'Marker','o',...
    'LineWidth',2);

% Create loglog
loglog(NP221,time221,'DisplayName','(LiH)$_{16}$','MarkerSize',9,'Marker','square',...
    'LineWidth',2,'LineStyle','-.');

% Create loglog
loglog(NP222,time222,'DisplayName','(LiH)$_{32}$','MarkerSize',9,'Marker','diamond',...
    'LineWidth',2,'LineStyle',':', 'Color', '[0.4660 0.6740 0.1880]');

% Create ylabel
ylabel({'Wall time (s)'},'Interpreter','latex');
% Create xlabel
xlabel({'Number of processors'}, 'Interpreter','latex');
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');

% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',24,'XMinorTick','on','XScale','log','YMinorTick','on',...
    'YScale','log');
xticks([10 100 1000 10000]);
yticks([1 10 100 1000]);
ylim([0.2 4000]);
set(gcf, 'position', [10 10 670 420]);
set(gca, 'position', [0.16 0.22 0.7164 0.7143]);
% axes1.XLabel.FontSize = 24;
% axes1.YLabel.FontSize = 24;
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','southeast','Interpreter','latex');
legend1.AutoUpdate = 'off';
hold on;
loglog(NP111,strongScaling111,NP221,strongScaling221,NP222,strongScaling222,...
    'LineWidth',1,'Color','black','LineStyle','--');
