%%
tic
clear
rng(0)

K = 5e4;
N = 50;
N_samp = 1000;
h = 0.01;
all_p = (0.2:0.05:0.8)';

%% theoretically predict
TPb = zeros(length(all_p), 1);
for jp = 1:length(all_p)
    p = all_p(jp);
    TPb(jp) = pred_prob(p, N);
end

%% simulation
CPb = zeros(length(all_p), N_samp, 'logical');
SPb = zeros(length(all_p), N_samp);
for jp = 1:length(all_p)
    p = all_p(jp);
    for js = 1:N_samp
        
        V = K * binornd(1, p, N, 1);
        U = eig1d(V, h, 1);
        u = abs(getval1d(U));
        SPb(jp, js) = max(u(1), u(end)) / max(u);
        
        intv = [0; find(diff(V)~=0); N];
        duan = diff(intv);
        
        if V(1) == 0 && V(N) == 0
            o_o = max(duan(3:2:end-2));
            u_u = 2 * max(duan(1), duan(end));
            if isempty(o_o) || o_o < u_u
                CPb(jp, js) = 1;
            end
        elseif V(1) == 0 && V(N) ~= 0
            o_o = max(duan(3:2:end-1));
            if isempty(o_o) || o_o < 2*duan(1)
                CPb(jp, js) = 1;
            end
        elseif V(1) ~= 0 && V(N) == 0
            o_o = max(duan(end-2:-2:1));
            if isempty(o_o) || o_o < 2*duan(end)
                CPb(jp, js) = 1;
            end
        elseif all(V==1)
            CPb(jp, js) = 1;
        end
        
%         if xor(CPb(jp,js), SPb(jp,js) > 0.5)
%             error('disagree')
%         end
    end
    fprintf('p = %g finished\n', all_p(jp));
    toc
end

%% plot
mCPb = mean(CPb, 2);
mSPb = mean(SPb, 2);

figure
hold on
plot(all_p, mSPb, 'bo')
plot(all_p, mCPb, 'r*')
plot(all_p, TPb, 'k-')
xlabel('p')
ylabel('probability')
xlim([0.1, 0.9])
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 14)
