function D = Dfunc1(a, b, L)

L = L(:) / (sum(L) + L(2));

x0 = L(1) / 2;
x1 = L(4) / 2;
x2 = x1 + L(5);
x3 = 0.5;

D = zeros(2, 1);
D(1) = a*tan(a*x0) - b*tanh(b*(x3-x0));
D(2) = (a^2 - b^2) * (exp(2*b*x2) + exp(2*b*(x1+x3))) / (-exp(2*b*x2) + exp(2*b*(x1+x3)))...
     + (a^2 + b^2) * (exp(2*b*x3) + exp(2*b*(x1+x2))) / (-exp(2*b*x2) + exp(2*b*(x1+x3)))...
     + 2*a*b * cot(a*(x1-x2));
end