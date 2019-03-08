% This is a cal demo of pure Neumann boundary Condition
% The result is u + C + err, which is as expected
% since the solution of pure NBC problem is not unique up to a constant

% Space Define
N = 8;
nx = N; ny = N;
T = RecMesh([nx, ny], [1, 1], [0, 0]);
U = P2Fespace(T);
% U = P1Fespace(T);

% System Clarification
K = symBilinear(U, 'nabla', []);

u  = @(x,y) x.^2.*y.^2 + 2*y;
f  = @(x,y) -2*y.^2-2*x.^2;
h1 = @(x,y) -2;
h2 = @(x,y) 2*y.^2;
h3 = @(x,y) 2*x.^2+2;
h4 = @(x,y) 0;

H = {h1, h2, h3, h4};

[X, FNodeptr] =  Freedomdefine(U, [0,0,0,0], []);
F = Load2(U, f, H);

% Solve and visualization
KK = K(FNodeptr, FNodeptr);
FF = F(FNodeptr) - K(FNodeptr, ~FNodeptr)*X(~FNodeptr);
X(FNodeptr) = KK\FF;

figure(1)
subplot(1,2,1)
trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), X-mean(X));
subplot(1,2,2)
trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), u(U.Node(:, 1), U.Node(:, 2)) );
figure(2)
err =  X-mean(X) - u(U.Node(:, 1), U.Node(:, 2));
trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), err - mean(err));
err = err - mean(err);
fprintf("The residual is %.6f\n",norm(K*err));
