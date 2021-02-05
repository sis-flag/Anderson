function D = Dfunc(alpha, beta, L)

L = L(:) / (sum(L) + L(2) + L(3));

x0 = L(1) / 2;
x1 = L(4) / 2;
x2 = x1 + L(3);
x3 = 0.5;

D = zeros(2, 1);
D(1) = alpha*tan(alpha*x0) - beta*tanh(beta*(x3-x0));
D(2) = (alpha^2 - beta^2) * (exp(2*beta*x2) + exp(2*beta*(x1+x3)))...
        / (-exp(2*beta*x2) + exp(2*beta*(x1+x3)))...
     + (alpha^2 + beta^2) * (exp(2*beta*x3) + exp(2*beta*(x1+x2)))...
        / (-exp(2*beta*x2) + exp(2*beta*(x1+x3)))...
     + 2*alpha*beta * cot(alpha*(x1-x2));
end