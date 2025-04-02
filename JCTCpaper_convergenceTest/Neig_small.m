clear;
load("Neig.mat");
fig = figure('Units', 'inches','OuterPosition',[0 0 3.33 3.53]);
hold on;
semilogy1 = semilogy(Neig(1:7), abs(ErrSi8(1:7)), '-o','LineWidth',1.5,'MarkerSize',5);
semilogy1.MarkerFaceColor = [0 0.4470 0.7410];
hold on;
semilogy2 = semilogy(Neig(1:7), abs(ErrC8(1:7)),'-pentagram','LineStyle','-.','LineWidth',1.5,'MarkerSize',5);
semilogy2.MarkerFaceColor = [0.8500 0.3250 0.0980];
hold on;
semilogy3 = semilogy(Neig(1:7), abs(ErrLiH(1:7)),'-^','LineStyle',':', 'Color', '[0.4660 0.6740 0.1880]','LineWidth',1.5,'MarkerSize',5);
semilogy3.MarkerFaceColor = [0.4660 0.6740 0.1880];
set(semilogy1,'DisplayName','Si');
set(semilogy2,'DisplayName','C');
set(semilogy3,'DisplayName','LiH');
ylabel({'Error (ha/atom)'});
xlabel({'$N_r/N_e$'},'Interpreter','latex');

ax = gca;
box(ax,'on');
set(ax,'FontName','Times New Roman','FontSize',14,'YMinorTick','on',...
    'YScale','log');
xlim([20 300]);
xticks([20 80 140 200 300]);
yticks([0.00001 0.0001 0.001 0.01]);
legend1 = legend(ax,'show');
set(legend1,'FontSize',12);


set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperPositionMode', 'auto');
fig_pos = get(fig, 'Position');
set(fig, 'PaperSize', fig_pos(3:4));