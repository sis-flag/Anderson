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
yticks([])
title('D^{(e)}')
box off
set(gcf, 'Position', [300 300 400 200])
set(gca, 'FontSize', 16)

figure
hold on
x = linspace(0, 1, 101);
cl = plot(x, cos(pi*x/2), 'k-');
cl.LineWidth = 1;
xlim([0, 1])
ylim([0, 1])
xticks([0, 1])
yticks([])
title('D')
box off
set(gcf, 'Position', [300 300 200 200])
set(gca, 'FontSize', 16)

xn = [0, 0.15, 0.36, 0.58, 0.7, 0.81, 1];
figure
hold on
plot(xn(1:2), [0, 0], 'b-', 'LineWidth', 5)
plot(xn(2:3), [0, 0], 'y-', 'LineWidth', 5)
plot(xn(3:4), [0, 0], 'b-', 'LineWidth', 5)
plot(xn(4:5), [0, 0], 'y-', 'LineWidth', 5)
plot(xn(5:6), [0, 0], 'b-', 'LineWidth', 5)
plot(xn(6:7), [0, 0], 'y-', 'LineWidth', 5)
xlim([0, 1])
ylim([0, 1])
yticks([])
set(gcf, 'Position', [300 300 500 100])
set(gca, 'FontSize', 16)
