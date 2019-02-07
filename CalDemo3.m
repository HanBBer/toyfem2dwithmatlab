% This is a Stokes Problem calculation demo

% Boundary condition


u1 = @(x,y) 1-y.^2;
u2 = @(x,y) 0.*x;
% p = @(x,y) -2*x;
f1 = @(x,y) 0;
f2 = @(x,y) 0;
%{
u1 = @(x,y) sin(x).*sin(y);
u2 = @(x,y) cos(x).*cos(y);
% p = @(x,y) cos(x).*sin(y);
f1 = @(x,y) 0;
f2 = @(x,y) 2*cos(x).*cos(y);
%}

% Space Define
N = 4;
nx = N; ny = N;
T = RecMesh([nx, ny], [2, 2], [-1, -1]);
U = P2Fespace(T);
P = P1Fespace(T);

% System Clarification
Ku = symBilinear(U, 'nabla', []);
Kuxp = nonsymBilinear(T, U, 'dx', P, 'mass', []);
Kuyp = nonsymBilinear(T, U, 'dy', P, 'mass', []);
Kp = symBilinear(P, 'mass', []);
G1 = {u1, u1, u1, u1};
G2 = {u2, u2, u2, u2};
[Xu1, FNu1] =  Freedomdefine(U, [1,1,1,1], G1);
[Xu2, FNu2] =  Freedomdefine(U, [1,1,1,1], G2);
[Xp, FNp] = Freedomdefine(P, [0,0,0,0], []);
F1 = Load(U, f1);
F2 = Load(U, f2);

% Assembling
Nu = U.N;
Np = P.N;
K = [Ku, sparse(Nu, Nu), -Kuxp;
    sparse(Nu, Nu), Ku, -Kuyp;
    -Kuxp', -Kuyp', 1e-10*Kp];
X = [Xu1; Xu2; Xp];
Freedom = [FNu1; FNu2; FNp];
F = [F1; F2; sparse(Np,1)];

% Solve and visualization
KK = K(Freedom, Freedom);
FF = F(Freedom) - K(Freedom, ~Freedom)*X(~Freedom);
X(Freedom) = KK\FF;
erru = X(1:2*Nu) - [u1(U.Node(:,1), U.Node(:,2)); u2(U.Node(:,1), U.Node(:,2))];
norm(K(1:2*Nu,1:2*Nu)*erru, 1)
subplot(1,2,1);
%trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), X(1:Nu));
trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), erru(1:Nu), ...
                'FaceColor', 'interp', 'EdgeColor', 'interp');
subplot(1,2,2);
%trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), X(Nu+1:2*Nu));
trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), erru(1:Nu));
%figure
%trisurf(P.Tri, P.Node(:, 1), P.Node(:, 2), X(2*Nu+1:end));