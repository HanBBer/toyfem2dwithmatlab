% This is a cal demo of lame system
% which refers to Example 3.8 of FreeFem++ lame.edp

% Space Define
N = 10;
nx = N; ny = N;
T = RecMesh([nx, ny], [20, 2], [0, -1]);
U = P2Fespace(T);
% U = P1Fespace(T);

E = 21e5; nu = 0.28; mu = E/(2*(1+nu));
lambda = E*nu/((1+nu)*(1-2*nu));

G = {[], [], [], @(x, y) 0};
[X, FNodeptr] = Freedomdefine(U, [0,0,0,1], G);

Nu = U.N;
Kx = symBilinear(U, 'dx', []);
Ky = symBilinear(U, 'dy', []);
Kxy = nonsymBilinear(T, U, 'dx', U, 'dy', []);
F  = -Load(U);

% Assembling
BigK = lambda*[Kx, Kxy; Kxy', Ky] + mu*[2*Kx+Ky, Kxy'; Kxy, 2*Ky+Kx];
BigF = [sparse(Nu, 1); F];
BigX = [X; X];
Freedom = [FNodeptr; FNodeptr];

KK = BigK(Freedom, Freedom);
FF = BigF(Freedom) - BigK(Freedom, ~Freedom)*BigX(~Freedom);
BigX(Freedom) = KK\FF;

coef = 100;
triplot(U.Tri, U.Node(:, 1)+BigX(1:Nu)*coef, U.Node(:, 2)+BigX(Nu+1:end)*coef);