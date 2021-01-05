%%
tic
clear
rng(0)

K = 3e6;
N = 50;
N_samp = 1000;
all_p = (0.2:0.05:0.8)';

thes = 0.02;

%% theoretically predict
TPD = zeros(length(all_p), 1);
TPN = zeros(length(all_p), 1);
for jp = 1:length(all_p)
    p = all_p(jp);
    [TPD(jp), TPN(jp)] = pred_multi(p, N);
end

%% simulation
CPD = zeros(length(all_p), N_samp, 'logical');
CPN = zeros(length(all_p), N_samp, 'logical');
SPD = zeros(length(all_p), N_samp, 'logical');
SPN = zeros(length(all_p), N_samp, 'logical');
for jp = 1:length(all_p)
    p = all_p(jp);
    for js = 1:N_samp
        
        V = K * binornd(1, p, N, 1);
        [UD, lamD] = eig1d(V, inf, 2);
        [UN, lamN] = eig1d(V, 0, 2);
        
        SPD(jp, js) = abs(lamD(2) - lamD(1)) / lamD(2) < thes;
        SPN(jp, js) = abs(lamN(2) - lamN(1)) / lamN(2) < thes;
        
        intv = [0; find(diff(V)~=0); N];
        duan = diff(intv);
        
        if V(1) == 0 && V(N) == 0
            o_o = sort(duan(1:2:end), 'descend');
            if length(o_o) > 1 && o_o(1) == o_o(2)
                CPD(jp, js) = 1;
            end
            
            duan(1) = 2 * duan(1);
            duan(end) = 2 * duan(end);
            u_u = sort(duan(1:2:end), 'descend');
            if length(u_u) > 1 && u_u(1) == u_u(2)
                CPN(jp, js) = 1;
            end
        elseif V(1) == 0 && V(N) ~= 0
            o_o = sort(duan(1:2:end), 'descend');
            if length(o_o) > 1 && o_o(1) == o_o(2)
                CPD(jp, js) = 1;
            end
            
            duan(1) = 2 * duan(1);
            u_u = sort(duan(1:2:end), 'descend');
            if length(u_u) > 1 && u_u(1) == u_u(2)
                CPN(jp, js) = 1;
            end
        elseif V(1) ~= 0 && V(N) == 0
            o_o = sort(duan(2:2:end), 'descend');
            if length(o_o) > 1 && o_o(1) == o_o(2)
                CPD(jp, js) = 1;
            end
            
            duan(end) = 2 * duan(end);
            u_u = sort(duan(2:2:end), 'descend');
            if length(u_u) > 1 && u_u(1) == u_u(2)
                CPN(jp, js) = 1;
            end
        else
            o_o = sort(duan(2:2:end), 'descend');
            if length(o_o) > 1 && o_o(1) == o_o(2)
                CPD(jp, js) = 1;
                CPN(jp, js) = 1;
            end
        end
        
%         if xor(CPD(jp,js), SPD(jp,js))
%             error('Dirhchlet disagree')
%         end
%         if xor(CPN(jp,js), SPN(jp,js))
%             error('Neumann disagree')
%         end
    end
    fprintf('p = %g finished\n', all_p(jp));
    toc
end

%% plot
mCPD = mean(CPD, 2);
mSPD = mean(SPD, 2);
mCPN = mean(CPN, 2);
mSPN = mean(SPN, 2);

figure
hold on
plot(all_p, mSPD, 'bo')
plot(all_p, mCPD, 'r*')
plot(all_p, TPD, 'k-')
xlabel('p')
ylabel('P_D')
xlim([0.1, 0.9])
ylim([0, 0.5])
legend('degree of loaclization', 'frequency', 'prediction')
set(gcf, 'Position', [500 500 400 300])

figure
hold on
plot(all_p, mSPN, 'bo')
plot(all_p, mCPN, 'r*')
plot(all_p, TPN, 'k-')
xlabel('p')
ylabel('P_N')
xlim([0.1, 0.9])
ylim([0, 0.5])
legend('degree of loaclization', 'frequency', 'prediction')
set(gcf, 'Position', [500 500 400 300])