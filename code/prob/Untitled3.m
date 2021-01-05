figure
hold on
x = linspace(0, 1, 101);
cl = plot(x, cos(pi*x/2), 'k-');
cl.LineWidth = 1;
x = linspace(-1, 0, 101);
cl = plot(x, cos(pi*x/2), 'k--');
cl.LineWidth = 1;
plot([0, 0], [-0.1, 1.1], 'k--')
xlim([-1, 1])
ylim([0, 1])
xticks([-1, 0, 1])
yticks([0, 1])
box off
set(gcf, 'Position', [500 500 500 200])

figure
hold on
x = linspace(0, 1, 101);
cl = plot(x, cos(pi*x/2), 'k-');
cl.LineWidth = 1;
xlim([0, 1])
ylim([0, 1])
xticks([0, 1])
yticks([0, 1])
box off
set(gcf, 'Position', [500 500 250 200])