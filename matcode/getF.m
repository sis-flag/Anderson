function F = getF(K, N)
V = [ones(ceil(N(2)/2), 1);
    zeros(N(1),1);...
    ones(N(2),1);...
    zeros(N(3),1);...
    ones(N(4),1);...
    zeros(N(5),1);...
    ones(floor(N(2)/2),1)];

LN = [ceil(N(2)/2), ceil(N(2)/2)+N(1)];
RN = [ceil(N(2)/2)+N(1)+N(2), sum(N)+N(2)-floor(N(2)/2)];

[U, ~] = eigP1d(K*V, 1);

u1 = abs(getval1d(U(:,1), 20));
left = max(u1(LN(1)*20+1: LN(2)*20+1));
right = max(u1(RN(1)*20+1: RN(2)*20+1));
F = left / (left + right);
end