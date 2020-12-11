clear

N = 20;
all_p = (0.05:0.05:0.95)';
N_samp = 1000;

ab1 = zeros(length(all_p), 1);
ab2 = zeros(length(all_p), 1);
ab3 = zeros(length(all_p), 1);
for jp = 1:length(all_p)
    p = all_p(jp);
    
    bigb1 = 0;
    bigb2 = 0;
    bigb3 = 0;
    for js = 1:N_samp
        V = rand(N, 1);
        V(V<p) = 0; V(V>=p) = 1;
        
        intv = [0; find(diff(V)~=0); N];
        duan = diff(intv);
        
        if V(1) == 0 && V(N) == 0
            o_o = max(duan(3:2:end-2));
            u_u = 2 * max(duan(1), duan(end));
            if isempty(o_o) || o_o < u_u
                bigb1 = bigb1 + 1;
            end
        elseif V(1) == 0 && V(N) == 1
            o_o = max(duan(3:2:end-1));
            if isempty(o_o) || o_o < 2*duan(1)
                bigb2 = bigb2 + 1;
            end
        elseif V(1) == 1 && V(N) == 0
            o_o = max(duan(end-2:-2:1));
            if isempty(o_o) || o_o < 2*duan(end)
                bigb3 = bigb3 + 1;
            end
        end
    end
    
    ab1(jp) = bigb1 / N_samp;
    ab2(jp) = bigb2 / N_samp;
    ab3(jp) = bigb3 / N_samp;
end

pb1 = zeros(length(all_p), 1);
pb2 = zeros(length(all_p), 1);
for nn = 1:length(all_p)
    [pb1(nn), pb2(nn)] = pred_prob(all_p(nn), N);
end

figure()
hold on
plot(all_p, ab1, '*')
plot(all_p, all_p.*all_p.*pb1, '--')

figure()
hold on
plot(all_p, ab2, '*')
plot(all_p, all_p.*(1-all_p).*pb2, '--')

figure()
hold on
plot(all_p, ab3,'*')
plot(all_p, all_p.*(1-all_p).*pb2, '--')

figure()
hold on
plot(all_p, ab2+ab3, '*')
plot(all_p, 2*all_p.*(1-all_p).*pb2, '--')

figure()
hold on
plot(all_p, ab1+ab2+ab3,'*')
plot(all_p, 2*all_p.*(1-all_p).*pb2 + all_p.*all_p.*pb1, '--')
