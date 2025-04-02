load('times_bargraph.mat');
fig = figure('Units', 'inches','OuterPosition',[0 0 3.33 3.83]);
h = bar(times221, 'stacked');
h(1).FaceColor = [0.6350 0.0780 0.1840];
h(2).FaceColor = [0 0.4470 0.7410];
h(3).FaceColor = [0 1 0.7410];
h(4).FaceColor = [0.9290 0.6940 0.1250];
legend1 = legend('Kronecker product','Sternheimer equation', 'Projection \& orthonormalization','Spectral quadrature');
set(legend1,'interpreter','Latex','FontSize',12);