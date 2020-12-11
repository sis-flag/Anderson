function D = Dfunc(K, h, lam)
% 这个函数的零点是特征值的预测值
% if lam <= 0 || lam >= K
%     error('lambda out of range!')
% end

x = cumsum(h);
if mod(length(x), 2) ~= 0
    error('len(x) must be even!')
end

a = sqrt(lam);
b = sqrt(K-lam);

AA = eye(2);
for k = 1:length(x)-1
    Gk = [sin(a*x(k)), cos(a*x(k));...
        a*cos(a*x(k)), -a*sin(a*x(k))];
    Ek = [exp(b*x(k)), exp(-b*x(k));...
        b*exp(b*x(k)), -b*exp(-b*x(k))];
    if mod(k,2) == 1
        AA = Ek \ Gk * AA;
    else
        AA = Gk \ Ek * AA;
    end
end
EN = [exp(b), exp(-b);...
    b*exp(b), -b*exp(-b)];
G0 = [0, 1; a, 0];
AA = G0 \ EN * AA;
D = det(AA - eye(2));

end