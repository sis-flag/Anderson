%%
tic
clear
% rng(0)

N = 100;
N_samp = 1000;

%% simulation
K = 1e8;
all_p = (0.3:0.2:0.7)';

Sb = zeros(length(all_p), N_samp, 'logical');
Cb = zeros(length(all_p), N_samp, 'logical');
for jp = 1:length(all_p)
    p = all_p(jp);
    
    for js = 1:N_samp
        V = rand(N, 1);
        V(V<p) = 0; V(V>=p) = 1;
        
        intv = [0; find(diff(V)~=0); N];
        duan = diff(intv);
        
        if V(1) == 0
            o_o = duan(1:2:end);
        else
            o_o = duan(2:2:end);
        end
        if sum(o_o == max(o_o)) > 1
            Cb(jp, js) = 1;
        end
        
        U = eigD1d(K*V, 6);
        u = abs(getval1d(U(:,1), 20));
        Sb(jp, js) = (length(nlocmax(u, 5e-2)) > 1);
        
%         if xor(Cb(jp,js), Sb(jp,js))
%             figure
%             hold on
%             cf = bar(((1:length(V))-0.5)/length(V), V/max(V));
%             cf.BarWidth = 1;
%             cf.LineStyle = 'None';
%             cf.FaceColor = 'c';
%             cf.FaceAlpha = 0.3;
%             cf = bar(((1:length(V))-0.5)/length(V), -V/max(V));
%             cf.BarWidth = 1;
%             cf.LineStyle = 'None';
%             cf.FaceColor = 'c';
%             cf.FaceAlpha = 0.3;
%             for k = 1:2
%                 uk = my_nmlz(getval1d(U(:,k)));
%                 x = linspace(0, 1, length(uk));
%                 plot(x, uk)
%             end
%             x = linspace(0, 1, length(V)+1);
%             plot(x, zeros(length(x),1), 'k+')
%             xlim([0,1])
%             ylim([-1,1])
%             
%             error('666')
%         end
    end
    
    fprintf('p = %g finished\n', all_p(jp));
    toc
end

%% save
% save('Data_1d_extra50.mat','all_p','K','h','N','Pb','Bb3','Bb2','Bb1')

%% predict
Tb = zeros(length(all_p), 1);
for nn = 1:length(all_p)
    p = all_p(nn);
    Tb(nn) = pred_2peak(p, N);
end

%% plot
mSb = mean(Sb, 2);
mCb = mean(Cb, 2);
mTb = mean(Tb, 2);

% figure
hold on
plot(all_p, mSb, 'r+--')
plot(all_p, mCb, 'mx--')
plot(all_p, mTb, 'b-')
% xlabel('p'); ylabel('Pb'); ylim([0, 1]);