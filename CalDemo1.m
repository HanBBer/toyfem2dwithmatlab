% dbstop if error;
% This is a Dirichlet Problem calculation demo

ShowFlag = 0;

% Problem Preparation(truth and boundary condition)
u = @(x, y) exp(sin(x)+cos(y));
f = @(x, y) -(cos(x).^2-sin(x)+sin(y).^2-cos(y)).*exp(sin(x)+cos(y));

% Space Define
N = 4;
nx = N; ny = N;
T = RecMesh([nx, ny], [1, 1], [0, 0]);
U = P2Fespace(T);

% System Clarification
K = symBilinear(U, 'nabla', []);
M = symBilinear(U, 'mass', []);
% K = nonsymBilinear(U, 'nabla', U, 'nabla', []);
G = {u, u, u, u};
[X, FNodeptr] = Freedomdefine(U, [1,1,1,1], G);
F = Load(U, f);

% Solve and visualization
KK = K(FNodeptr, FNodeptr);
FF = F(FNodeptr) - K(FNodeptr, ~FNodeptr)*X(~FNodeptr);
X(FNodeptr) = KK\FF;
err1 = X-u(U.Node(:, 1), U.Node(:, 2));
if ShowFlag == 1
    figure(1)
    trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), X);
    figure(2)
    trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), err1);
end

err1L2 = norm(K*err1,1);
fprintf("h = %.3f, L2 error is %.6f\n", 1/N, err1L2);


%% Refine once
T = Refine(T);
U = P2Fespace(T);

% System Clarification
K = symBilinear(U, 'nabla', []);
[X, FNodeptr] = Freedomdefine(U, [1,1,1,1], G);
F = Load(U, f);

% Solve and visualization
KK = K(FNodeptr, FNodeptr);
FF = F(FNodeptr) - K(FNodeptr, ~FNodeptr)*X(~FNodeptr);
X(FNodeptr) = KK\FF;
err2 = X-u(U.Node(:, 1), U.Node(:, 2));

err2L2 = norm(K*err2,1);
fprintf("h = %.3f, L2 error is %.6f\n", 1/N/2, err2L2);

fprintf("Converge rate is %2f\n", log(err1L2/err2L2)/log(2));
%% Refine twice
T = Refine(T);
U = P2Fespace(T);

% System Clarification
K = symBilinear(U, 'nabla', []);
[X, FNodeptr] = Freedomdefine(U, [1,1,1,1], G);
F = Load(U, f);

% Solve and visualization
KK = K(FNodeptr, FNodeptr);
FF = F(FNodeptr) - K(FNodeptr, ~FNodeptr)*X(~FNodeptr);
X(FNodeptr) = KK\FF;
err3 = X-u(U.Node(:, 1), U.Node(:, 2));

err3L2 = norm(K*err3,1);
fprintf("h = %.3f, L2 error is %.6f\n", 1/N/2, err3L2);

fprintf("Converge rate is %2f\n", log(err2L2/err3L2)/log(2));