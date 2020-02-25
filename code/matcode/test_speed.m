clear;
rng(0);

M = 200;
V = 500 * rand(M);

fprintf('problem size: %d\n', M)

disp('Robin boundary source problem:')
tic; W = solveR2d(V, 0); toc

disp('Robin boundary eigen problem:')
tic; [U, lam] = eigR2d(V, 0); toc

disp('Dirichlet boundary source problem:')
tic; W = solveD2d(V); toc

disp('Dirichlet boundary eigen problem:')
tic; [U, lam] = eigD2d(V); toc

disp('get real value')
tic; w = getval2d(W, 20); toc
