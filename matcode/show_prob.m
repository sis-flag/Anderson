%% 1d case
load Data_1d

dm = zeros(length(all_p), length(all_K), length(all_h));
ds = dm;
for jp = 1:length(all_p)
    for jk = 1:length(all_K)
        for jh = 1:length(all_h)
            dm(jp, jk, jh) = mean(Data(jp, jk, jh, 1, :), 'all');
            ds(jp, jk, jh) = std(Data(jp, jk, jh, 1, :), 0, 'all');
        end
    end
end

% for jp = 1:length(all_p)
%     for jk = 1:length(all_K)
%         mm = reshape(dm(jp,jk,:), [], 1);
%         ss = reshape(ds(jp,jk,:), [], 1);
%         figure()
%         semilogx(all_h, mm, '*-')
%         title(['p=',num2str(all_p(jp)),' K=',num2str(all_K(jk))])
%         ylim([0, 1])
%     end
% end

for jh = 1:length(all_h(1:3))
    for jk = 1:length(all_K)
        mm = reshape(dm(:,jk,jh), [], 1);
        ss = reshape(ds(:,jk,jh), [], 1);
        figure()
        plot(all_p, mm, '*-')
        ylim([0, 1])
        title(['h=',num2str(all_h(jh)),' K=',num2str(all_K(jk))])
    end
end

% for jp = 1:length(all_p)
%     for jh = 1:length(all_h)
%         mm = reshape(dm(jp,:,jh), [], 1);
%         ss = reshape(ds(jp,:,jh), [], 1);
%         figure()
%         semilogx(all_K, mm, '*-')
%         ylim([0, 1])
%         title(['p=',num2str(all_p(jp)),' h=',num2str(all_h(jh))])
%     end
% end
%% 2d case
load Data_2d

dm = zeros(length(all_p), length(all_K), length(all_h));
ds = dm;
for jp = 1:length(all_p)
    for jk = 1:length(all_K)
        for jh = 1:length(all_h)
            dm(jp, jk, jh) = mean(Data_e(jp, jk, jh, 1, :), 'all');
            ds(jp, jk, jh) = std(Data_e(jp, jk, jh, 1, :), 0, 'all');
        end
    end
end

% for jp = 1:length(all_p)
%     for jk = 1:length(all_K)
%         mm = reshape(dm(jp,jk,:), [], 1);
%         ss = reshape(ds(jp,jk,:), [], 1);
%         figure()
%         semilogx(all_h, mm, '*-')
%         title(['p=',num2str(all_p(jp)),' K=',num2str(all_K(jk))])
%         ylim([0, 1])
%     end
% end

% for jh =1:length(all_h)
%     for jk = 1:length(all_K)
%         mm = reshape(dm(:,jk,jh), [], 1);
%         ss = reshape(ds(:,jk,jh), [], 1);
%         figure()
%         plot(all_p, mm, '*-')
%         ylim([0, 1])
%         title(['h=',num2str(all_h(jh)),' K=',num2str(all_K(jk))])
%     end
% end

% for jp = 1:length(all_p)
%     for jh = 1:length(all_h)
%         mm = reshape(dm(jp,:,jh), [], 1);
%         ss = reshape(ds(jp,:,jh), [], 1);
%         figure()
%         semilogx(all_K, mm, '*-')
%         ylim([0, 1])
%         title(['p=',num2str(all_p(jp)),' h=',num2str(all_h(jh))])
%     end
% end

%% 1d special case
load Data_1d_s

dm = mean(Data_hb, 3); ds = std(Data_hb, 0, 3);
figure();
errorbar(all_h, dm(:,1), ds(:,1));
set(gca, 'xscale', 'log');
xlabel('h'); ylabel('Pb'); ylim([-0.1, 1.1]);
% saveas(gcf, 'hb.png');

dm = mean(Data_pb1, 3); ds = std(Data_pb1, 0, 3);
figure();
errorbar(all_p, dm(:,1), ds(:,1));
xlabel('p'); ylabel('Pb'); ylim([-0.1, 1.1]);
% saveas(gcf, 'pb1.png');

dm = mean(Data_pb2, 3); ds = std(Data_pb2, 0, 3);
figure();
errorbar(all_p, dm(:,1), ds(:,1));
xlabel('p'); ylabel('Pb'); ylim([-0.1, 1.1]);
% saveas(gcf, 'pb2.png');

dm = mean(Data_kb1, 3); ds = std(Data_kb1, 0, 3);
figure();
errorbar(all_K, dm(:,1), ds(:,1));
xlabel('K'); ylabel('Pb'); ylim([-0.1, 1.1]);
set(gca, 'xscale', 'log');
% saveas(gcf, 'kb1.png');

dm = mean(Data_kb2, 3); ds = std(Data_kb2, 0, 3);
figure();
errorbar(all_K, dm(:,1), ds(:,1));
set(gca, 'xscale', 'log');
xlabel('K'); ylabel('Pb'); ylim([-0.1, 1.1]);
% saveas(gcf, 'kb2.png');


%% 2d special case
load Data_2d_s

dm = mean(Data_hc, 3); ds = std(Data_hc, 0, 3);
figure();
errorbar(all_h, dm(:,1), ds(:,1));
set(gca, 'xscale', 'log');
xlabel('h'); ylabel('Pc'); ylim([-0.1, 1.1]);
% % saveas(gcf, 'hc.png');

dm = mean(Data_pc1, 3); ds = std(Data_pc1, 0, 3);
figure();
errorbar(all_p, dm(:,1), ds(:,1));
xlabel('p'); ylabel('Pc'); ylim([-0.1, 1.1]);
% % saveas(gcf, 'pc1.png');

dm = mean(Data_pc2, 3); ds = std(Data_pc2, 0, 3);
figure();
errorbar(all_p, dm(:,1), ds(:,1));
xlabel('p'); ylabel('Pc'); ylim([-0.1, 1.1]);
% % saveas(gcf, 'pc2.png');

dm = mean(Data_kc1, 3); ds = std(Data_kc1, 0, 3);
figure();
errorbar(all_K, dm(:,1), ds(:,1));
xlabel('K'); ylabel('Pc'); ylim([-0.1, 1.1]);
set(gca, 'xscale', 'log');
% % saveas(gcf, 'kc1.png');

dm = mean(Data_kc2, 3); ds = std(Data_kc2, 0, 3);
figure();
errorbar(all_K, dm(:,1), ds(:,1));
xlabel('K'); ylabel('Pc'); ylim([-0.1, 1.1]);
set(gca, 'xscale', 'log');
% % saveas(gcf, 'kc2.png');

dm = mean(Data_he, 3); ds = std(Data_he, 0, 3);
figure();
errorbar(all_h, dm(:,1), ds(:,1));
xlabel('h'); ylabel('Pe'); ylim([-0.1, 1.1]);
set(gca, 'xscale', 'log');
% % saveas(gcf, 'he.png');

dm = mean(Data_pe1, 3); ds = std(Data_pe1, 0, 3);
figure();
errorbar(all_p, dm(:,1), ds(:,1));
xlabel('p'); ylabel('Pe'); ylim([-0.1, 1.1]);
% % saveas(gcf, 'pe1.png');

dm = mean(Data_pe2, 3); ds = std(Data_pe2, 0, 3);
figure();
errorbar(all_p, dm(:,1), ds(:,1));
xlabel('p'); ylabel('Pe'); ylim([-0.1, 1.1]);
% % saveas(gcf, 'pe2.png');

dm = mean(Data_ke1, 3); ds = std(Data_ke1, 0, 3);
figure();
errorbar(all_K, dm(:,1), ds(:,1));
xlabel('K'); ylabel('Pe'); ylim([-0.1, 1.1]);
set(gca, 'xscale', 'log');
% % saveas(gcf, 'ke1.png');

dm = mean(Data_ke2, 3); ds = std(Data_ke2, 0, 3);
figure();
errorbar(all_K, dm(:,1), ds(:,1));
xlabel('K'); ylabel('Pe'); ylim([-0.1, 1.1]);
set(gca, 'xscale', 'log');
% % saveas(gcf, 'ke2.png');