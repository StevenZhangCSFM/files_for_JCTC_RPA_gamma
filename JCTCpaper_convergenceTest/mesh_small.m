load("mesh.mat");
fig = figure('Units', 'inches','OuterPosition',[0 0 3.33 3.53]);
hold on;
loglog1 = loglog(meshSi8(1:4),abs(ErrSi8(1:4)),'-o','LineWidth',1.5,'MarkerSize',5);
loglog1.MarkerFaceColor = [0 0.4470 0.7410];
hold on;
loglog2 = loglog(meshC8(1:4),abs(ErrC8(1:4)),'-pentagram','LineStyle','-.','LineWidth',1.5,'MarkerSize',5);
loglog2.MarkerFaceColor = [0.8500 0.3250 0.0980];
hold on;
loglog3 = loglog(meshLi4H4(1:4),abs(ErrLi4H4(1:4)),'-^','LineStyle',':', 'Color', '[0.4660 0.6740 0.1880]','LineWidth',1.5,'MarkerSize',5);
loglog3.MarkerFaceColor = [0.4660 0.6740 0.1880];
set(loglog1,'DisplayName','Si');
set(loglog2,'DisplayName','C');
set(loglog3,'DisplayName','LiH');
ylabel({'Error (ha/atom)'});
xlabel('$h$ (bohr)','Interpreter','latex');

ax = gca;
box(ax,'on');
set(ax,'FontName','Times New Roman','FontSize',14,'XDir','reverse','XMinorTick','on','XScale','log','XTick',...
    [0.01 0.05 0.1 0.2 0.3],'YMinorTick','on','YScale','log');
xticks([0.15 0.25 0.35 0.45]);
ylim([1e-7, 1e-2])
legend1 = legend(ax,'show');
set(legend1,'FontSize',12,'Location','Southwest');

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperPositionMode', 'auto');
fig_pos = get(fig, 'Position');
set(fig, 'PaperSize', fig_pos(3:4));