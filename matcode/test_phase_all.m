clear 

Data = zeros(1000, 11);
count = 1;

for N1 = [3:10]
    for N3 = ceil(N1/2):N1-1
        for N4 = [1]
            for N5 = N3:N1-1
                for N2 = ceil((N1+N3+N4+N5)/2)*4: 2: (N1+N3+N4+N5+4)*2
                    
                    if max(N3, N5) < N1 &&...
                       N1 < N3 + N5 &&...
                       N4 < min(N3, N5) / 2
                   

N = [N1, N2, N3, N4, N5];

% find threshold range
minTK = 10; maxTK = 1e4; F = 0;
while abs(F-0.5) > 1e-4 && (maxTK-minTK) > 1e-8
    K = (maxTK + minTK) / 2;
    F = getF(K, N);
    if F < 0.5
        minTK = K;
    else
        maxTK = K;
    end
end
K5 = K;

minTK = 10; maxTK = K5; F = 0;
while abs(F-0.3) > 1e-4 && (maxTK-minTK) > 1e-8
    K = (maxTK + minTK) / 2;
    F = getF(K, N);
    if F < 0.3
        minTK = K;
    else
        maxTK = K;
    end
end
K3 = K;

% find threshold range
minTK = K5; maxTK = 1e4; F = 0;
while abs(F-0.7) > 1e-4 && (maxTK-minTK) > 1e-8
    K = (maxTK + minTK) / 2;
    F = getF(K, N);
    if F < 0.7
        minTK = K;
    else
        maxTK = K;
    end
end
K7 = K;


% find threshold range
minTK = K5; maxTK = K7; F = 0;
while abs(F-0.55) > 1e-4 && (maxTK-minTK) > 1e-8
    K = (maxTK + minTK) / 2;
    F = getF(K, N);
    if F < 0.55
        minTK = K;
    else
        maxTK = K;
    end
end
K55 = K;

% find threshold range
minTK = K3; maxTK = K5; F = 0;
while abs(F-0.45) > 1e-4 && (maxTK-minTK) > 1e-8
    K = (maxTK + minTK) / 2;
    F = getF(K, N);
    if F < 0.45
        minTK = K;
    else
        maxTK = K;
    end
end
K45 = K;

all_K = linspace(K45, K55, 21);
all_F = zeros(1, length(all_K));
for jk = 1:length(all_K)
    all_F(jk) = getF(all_K(jk), N);
end

lF = log(abs(all_F-0.5)); lK = log(abs(all_K-K5));
coe = polyfit(lK, lF, 1);
a = coe(1); b = exp(coe(2));

coe = polyfit(all_K-K5, all_F-0.5, 1);
df = coe(1);

fprintf('%d & %d & %d & %d & %d & %g & %g & %g & %g & %g & %g \\\\\n',...
    N(1), N(2), N(3), N(4), N(5), K3, K5, K7, a, b, df);

Data(count,:) =  [N(1), N(2), N(3), N(4), N(5), K3, K5, K7, a, b, df];
count = count + 1;

                    end
                end
            end
        end
    end
end

Data = Data(1:count-1,:);
save('Data_phase', 'Data')