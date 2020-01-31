clear

% peremeters
h = 0;
V = [2,8,3,5,1,6,0,7,9,4]*100;

figure();
bar(V)

% Robin boundar
W = solveR1d(V, h);
[U, lam] = eigR1d(V, h);

w = getval1d(W, 50);
u1 = getval1d(U(:,1), 50); u1 = my_nmlz(u1);
u2 = getval1d(U(:,2), 50); u2 = my_nmlz(u2);
u3 = getval1d(U(:,3), 50); u3 = my_nmlz(u3);

x = linspace(0,1,length(w));
figure()
hold on
plot(x, w, 'k')
plot(x, u1 / (lam(1) + max(w)) )
plot(x, u2 / (lam(2) + max(w)) )
plot(x, u3 / (lam(3) + max(w)) )

% Dirichlet boundary
W = solveD1d(V);
[U, lam] = eigD1d(V);

w = getval1d(W, 50);
u1 = getval1d(U(:,1), 50); u1 = my_nmlz(u1);
u2 = getval1d(U(:,2), 50); u2 = my_nmlz(u2);
u3 = getval1d(U(:,3), 50); u3 = my_nmlz(u3);

x = linspace(0,1,length(w));

figure()
hold on
plot(x, w, 'k')
plot(x, u1 / (lam(1) + max(w)) )
plot(x, u2 / (lam(2) + max(w)) )
plot(x, u3 / (lam(3) + max(w)) )
