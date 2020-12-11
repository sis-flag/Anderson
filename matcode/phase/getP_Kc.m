function [Kc, b, lamc] = getP_Kc(L)
L = L(:) / (sum(L) + L(2));
h = [L(2)/2; L; L(2)/2];

    function [F, lam] = getF(K)
        [U, lam] = eigP(K, h, 1);
        u = abs(getvalP(U, h));
        left = max(u(31:61));
        right = max(u(91:181));
        F = left / (left + right);
    end

% find threshold
minTK = 200; maxTK = 1e5; F = 0;
while abs(F-0.5) > 1e-5 && (maxTK-minTK) > 1e-10
    K = (maxTK + minTK) / 2;
    [F, lam] = getF(K);
    if F < 0.5
        minTK = K;
    else
        maxTK = K;
    end
end
Kc = K;
lamc = lam;

if K == 200
    Kc = nan;
    b = nan;
    lamc = nan;
    return
end

if nargout < 2
    return
end

% find threshold range
minTK = Kc; maxTK = Kc+100; F = 0;
while abs(F-0.55) > 1e-5 && (maxTK-minTK) > 1e-10
    K = (maxTK + minTK) / 2;
    F = getF(K);
    if F < 0.55
        minTK = K;
    else
        maxTK = K;
    end
end
K55 = K;

% find threshold range
minTK = Kc-100; maxTK = Kc; F = 0;
while abs(F-0.45) > 1e-5 && (maxTK-minTK) > 1e-10
    K = (maxTK + minTK) / 2;
    F = getF(K);
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
    all_F(jk) = getF(all_K(jk));
end

coe = polyfit(all_K-Kc, all_F-0.5, 1);
b = coe(1);

end