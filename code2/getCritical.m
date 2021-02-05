function [Kc, lamc] = getCritical(L)

h = [L(2)/2; L(1); L(2); L(3); L(4); L(3); L(2)/2];
h = h / sum(h); % normalize

    function [F, lam] = getF(K)
        [U, lam] = eigPhase(K, h);
        [~, ~, hp1, hp2] = getvalPhase(U, h);
        F = hp1 / (hp1 + hp2);
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

if K < 201 || K > 1e5-1
    Kc = nan;
    lamc = nan;
end

end