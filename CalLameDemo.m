% This is a cal demo of lame system

% Space Define
nx = 30; ny = 5;
T = RecMesh([nx, ny], [20, 2], [0, -1]);
P = P1Fespace(T);
U = P2Fespace(T);

NSp = P.N;
NSu = U.N;

% System Info

rhoS = 1.1;
mu = 5.75e6;
lambda = 1.7e6;


% Initial State Definition

Kx = symBilinear(P, 'dx', []);
Ky = symBilinear(P, 'dy', []);
Kxy = nonsymBilinear(T, P, 'dx', P, 'dy', []);
MovingK = lambda*[Kx, Kxy; Kxy', Ky] + mu*[2*Kx+Ky, Kxy'; Kxy, 2*Ky+Kx];
MovingF = sparse(2*NSp, 1);

G = {@(x,y) 0, @(x,y) 0, @(x,y) 0, @(x,y) 0};
[X, FNMoving] = Freedomdefine(P, [1,1,1,1], G);
FreeB = [FNMoving; FNMoving];
BX = X;
BY = X;

posM = T.Edge(T.EgFlag == 1, :);
g = @(x) 0.6*sin(pi*x/10);
BY(posM) = g(T.Node(posM, 1));
posM = T.Edge(T.EgFlag == 3, :);
g = @(x) 0.6*sin(pi*x/10);
BY(posM) = g(T.Node(posM, 1));

XX = [BX; BY];
KMoving = MovingK(FreeB, FreeB);
FBoundary = MovingF(FreeB) - MovingK(FreeB, ~FreeB)*XX(~FreeB);
XX(FreeB) = KMoving\FBoundary;

%T.Node = reshape(XX, NFp, 2);
triplot(T.Tri, T.Node(:,1), T.Node(:,2));
axis equal
box off;set(gca,'xtick',[],'ytick',[],'xcolor','white','ycolor','white');
pause(0.1)
triplot(T.Tri, T.Node(:,1)+XX(1:NSp), T.Node(:,2)+XX(NSp+1:end));
box off;set(gca,'xtick',[],'ytick',[],'xcolor','white','ycolor','white');
axis equal
pause(0.1)
% Test For the elasticity
dt = 5e-4;
runtime = 0.1;

DXold = InterpolateP1toP2(T, XX(1:NSp));
DYold = InterpolateP1toP2(T, XX(NSp+1:end));

G = {[], @(x, y) 0, [], @(x, y) 0};
[X, FNodeptr] = Freedomdefine(U, [0,1,0,1], G);
Kx = symBilinear(U, 'dx', []);
Ky = symBilinear(U, 'dy', []);
Kxy = nonsymBilinear(T, U, 'dx', U, 'dy', []);

MS = symBilinear(U, 'mass', []);
Lame = lambda*[Kx, Kxy; Kxy', Ky] + mu*[2*Kx+Ky, Kxy'; Kxy, 2*Ky+Kx];
Sys = [[rhoS/dt*MS, sparse(NSu, NSu); sparse(NSu, NSu), rhoS/dt*MS], Lame;
    -[MS, sparse(NSu, NSu); sparse(NSu, NSu), MS], [MS/dt, sparse(NSu, NSu); sparse(NSu, NSu), MS/dt]];
Sysold = [[rhoS/dt*MS, sparse(NSu, NSu); sparse(NSu, NSu), rhoS/dt*MS], sparse(2*NSu, 2*NSu);
    sparse(2*NSu, 2*NSu), [MS/dt, sparse(NSu, NSu); sparse(NSu, NSu), MS/dt]];

Free = [FNodeptr;FNodeptr;FNodeptr;FNodeptr];
Xold = [sparse(2*NSu,1); DXold; DYold];

BigF = sparse(4*NSu, 1);
KK = Sys(Free, Free);
XX = [X;X;X;X];

for t = 0:dt:runtime
    FF = BigF + Sysold*Xold;
    FF = FF(Free) - Sys(Free, ~Free)*XX(~Free);
    XX(Free) = KK\FF;
    BX = T.Node(:,1) + InterpolateP2toP1(T, XX(2*NSu+1:3*NSu));
    BY = T.Node(:,2) + InterpolateP2toP1(T, XX(3*NSu+1:4*NSu));
    triplot(T.Tri, BX, BY);
    axis equal
    title(['t = ', num2str(t)]);
    box off;set(gca,'xtick',[],'ytick',[],'xcolor','white','ycolor','white');
    pause(0.01)
    Xold = XX;
end