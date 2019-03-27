% This is a cal demo of moving mesh
N = 16;
T = RecMesh([N, N]);
P = P1Fespace(T);

NFp = P.N;

Kx = symBilinear(P, 'dx', []);
Ky = symBilinear(P, 'dy', []);
Kxy = nonsymBilinear(T, P, 'dx', P, 'dy', []);
MovingK = 10*[Kx, Kxy; Kxy', Ky] + [2*Kx+Ky, Kxy'; Kxy, 2*Ky+Kx];
MovingF = sparse(2*NFp, 1);

G = {@(x,y) 0, @(x,y) 0, @(x,y) 0, @(x,y) 0};
[XX, FNMoving] = Freedomdefine(P, [1,1,1,1], G);
FreeB = [FNMoving; FNMoving];
BX = XX;
BY = XX;

posM = T.Edge(T.EgFlag == 3, :);
g = @(x) 0.3*sin(2*pi*x);
BY(posM) = g(T.Node(posM,1));

XX = [BX; BY];
KMoving = MovingK(FreeB, FreeB);
FBoundary = MovingF(FreeB) - MovingK(FreeB, ~FreeB)*XX(~FreeB);
XX(FreeB) = KMoving\FBoundary;

T.Node = T.Node+reshape(XX, NFp, 2);

triplot(T.Tri, T.Node(:,1), T.Node(:,2));


U = P2Fespace(T);
M = symBilinear(U, 'mass', []);

u = @(x, y) x.*y;
v = @(x, y) x.^2;
U1 = u(U.Node(:,1), U.Node(:,2));
V1 = v(U.Node(:,1), U.Node(:,2));

fprintf("The FEM approximation value is %.6f\n", U1'*M*V1);
fprintf("The true value is %.6f\n", integral2(@(x,y) x.^3.*y,0,1,0,@(x)1+g(x)));
