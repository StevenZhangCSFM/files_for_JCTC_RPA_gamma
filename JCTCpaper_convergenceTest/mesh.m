load("mesh.mat");
axes1 = axes('Position',[0.15 0.15 0.8 0.8]);
hold(axes1,'on');
loglog1 = loglog(meshSi8(1:4),abs(ErrSi8(1:4)),'-o','LineWidth',2.5,'MarkerSize',9);
loglog1.MarkerFaceColor = [0 0.4470 0.7410];
hold on;
loglog2 = loglog(meshC8(1:4),abs(ErrC8(1:4)),'-pentagram','LineStyle','-.','LineWidth',2.5,'MarkerSize',9);
loglog2.MarkerFaceColor = [0.8500 0.3250 0.0980];
hold on;
loglog3 = loglog(meshLi4H4(1:4),abs(ErrLi4H4(1:4)),'-^','LineStyle',':', 'Color', '[0.4660 0.6740 0.1880]','LineWidth',2.5,'MarkerSize',9);
loglog3.MarkerFaceColor = [0.4660 0.6740 0.1880];
set(loglog1,'DisplayName','Si');
set(loglog2,'DisplayName','C');
set(loglog3,'DisplayName','LiH');
ylabel({'Error (ha/atom)'});
xlabel('$h$ (bohr)','Interpreter','latex');
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
set(axes1,'FontName','Times New Roman','FontSize',24,'XDir','reverse','XMinorTick','on','XScale','log','XTick',...
    [0.01 0.05 0.1 0.2 0.3],'YMinorTick','on','YScale','log');
xticks([0.15 0.25 0.35 0.45]);
ylim([1e-7, 1e-2])
set(gcf, 'position', [10 10 670 420]);
set(gca, 'position', [0.16 0.22 0.7164 0.7143]);
legend1 = legend(axes1,'show');
set(legend1,'FontSize',24);