L = [5, 10, 3, 1];
h = [L(2)/2; L(1); L(2); L(3); L(4); L(3); L(2)/2];
L = L / sum(h); % normalize
h = h / sum(h);
xn = [0; cumsum(h)];

figure
hold on
plot(xn(1:2), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(2:3), [0, 0], 'y-', 'LineWidth', 4)
plot(xn(3:4), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(4:5), [0, 0], 'y-', 'LineWidth', 4)
plot(xn(5:6), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(6:7), [0, 0], 'y-', 'LineWidth', 4)
plot(xn(7:8), [0, 0], 'b-', 'LineWidth', 4)
for tx = xn
    plot([tx, tx], [-0.3, 0.1], 'k-', 'LineWidth', 1)
end
text(sum(xn(1:2))/2, 0.3, 'L_2/2', 'FontSize', 20)
text(sum(xn(2:3))/2, 0.3, 'L_1', 'FontSize', 20)
text(sum(xn(3:4))/2, 0.3, 'L_2', 'FontSize', 20)
text(sum(xn(4:5))/2, 0.3, 'L_3', 'FontSize', 20)
text(sum(xn(5:6))/2, 0.3, 'L_4', 'FontSize', 20)
text(sum(xn(6:7))/2, 0.3, 'L_3', 'FontSize', 20)
text(sum(xn(7:8))/2, 0.3, 'L_2/2', 'FontSize', 20)
text(sum(xn(2:3))/2, -0.3, 'W_1', 'FontSize', 20)
text(sum(xn(5:6))/2, -0.3, 'W_2', 'FontSize', 20)
xlim([0, 1])
ylim([-0.5, 0.5])
yticks([])
xticks([])
title('$V(x)$', 'Interpreter', 'latex')
set(gcf, 'Position', [300 300 500 200])
set(gca, 'FontSize', 20)

figure
hold on
plot(xn(1:2), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(2:3), [0, 0], 'y-', 'LineWidth', 4)
plot(xn(3:4), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(4:5), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(5:6), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(6:7), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(7:8), [0, 0], 'b-', 'LineWidth', 4)
xlim([0, 1])
ylim([-0.5, 0.5])
yticks([])
xticks([])
title('$V_1(x)$', 'Interpreter', 'latex')
set(gcf, 'Position', [300 300 500 100])
set(gca, 'FontSize', 20)

figure
hold on
plot(xn(1:2), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(2:3), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(3:4), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(4:5), [0, 0], 'y-', 'LineWidth', 4)
plot(xn(5:6), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(6:7), [0, 0], 'y-', 'LineWidth', 4)
plot(xn(7:8), [0, 0], 'b-', 'LineWidth', 4)
xlim([0, 1])
ylim([-0.5, 0.5])
yticks([])
xticks([])
title('$V_2(x)$', 'Interpreter', 'latex')
set(gcf, 'Position', [300 300 500 100])
set(gca, 'FontSize', 20)

h = [(1-L(1))/2; L(1); (1-L(1))/2];
xn = [0; cumsum(h)];
figure
hold on
plot(xn(1:2), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(2:3), [0, 0], 'y-', 'LineWidth', 4)
plot(xn(3:4), [0, 0], 'b-', 'LineWidth', 4)
plot([0.5, 0.5], [-0.3, 0.3], 'k--', 'LineWidth', 1)
xlim([0, 1])
ylim([-0.5, 0.5])
yticks([])
xticks([])
title('$\hat{V}_1(x)$', 'Interpreter', 'latex')
set(gcf, 'Position', [300 300 500 100])
set(gca, 'FontSize', 20)

h = [(1-L(1))/2; L(1)/2];
xn = [0; cumsum(h)];
figure
hold on
plot(xn(1:2), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(2:3), [0, 0], 'y-', 'LineWidth', 4)
xlim([0, 1])
ylim([-0.5, 0.5])
yticks([])
xticks([])
title('$\tilde{V}_1(x)$', 'Interpreter', 'latex')
set(gcf, 'Position', [300 300 500 100])
set(gca, 'FontSize', 20)


h = [(1-L(3)-L(4)-L(3))/2; L(3); L(4); L(3); (1-L(3)-L(4)-L(3))/2];
xn = [0; cumsum(h)];
figure
hold on
plot(xn(1:2), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(2:3), [0, 0], 'y-', 'LineWidth', 4)
plot(xn(3:4), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(4:5), [0, 0], 'y-', 'LineWidth', 4)
plot(xn(5:6), [0, 0], 'b-', 'LineWidth', 4)
plot([0.5, 0.5], [-0.3, 0.3], 'k--', 'LineWidth', 1)
xlim([0, 1])
ylim([-0.5, 0.5])
yticks([])
xticks([])
title('$\hat{V}_2(x)$', 'Interpreter', 'latex')
set(gcf, 'Position', [300 300 500 100])
set(gca, 'FontSize', 20)

h = [(1-L(3)-L(4)-L(3))/2; L(3); L(4)/2];
xn = [0; cumsum(h)];
figure
hold on
plot(xn(1:2), [0, 0], 'b-', 'LineWidth', 4)
plot(xn(2:3), [0, 0], 'y-', 'LineWidth', 4)
plot(xn(3:4), [0, 0], 'b-', 'LineWidth', 4)
xlim([0, 1])
ylim([-0.5, 0.5])
yticks([])
xticks([])
title('$\tilde{V}_2(x)$', 'Interpreter', 'latex')
set(gcf, 'Position', [300 300 500 100])
set(gca, 'FontSize', 20)