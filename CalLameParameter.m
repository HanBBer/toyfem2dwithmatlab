function data = CalLameParameter(dt, N, J, Info)

% Space Define

nx = 30*N; ny = 5*N;
T = RecMesh([nx, ny], [20, 2], [0, -1]);
P = P1Fespace(T);
U = P2Fespace(T);

NSp = P.N;
NSu = U.N;

% System Info
if nargin < 4
    rhoS = 1.1;
    mu = 5.75e6;
    lambda = 1.7e6;
else
    rhoS = Info.rhoS;
    mu = Info.mu;
    lambda = Info.lambda;
end


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
posM = unique(posM(:));
g = @(x) 0.8*sin(pi*x/10);
BY(posM) = g(T.Node(posM, 1));
posM = T.Edge(T.EgFlag == 3, :);
g = @(x) 0.8*sin(pi*x/10);
posM = unique(posM(:));
BY(posM) = g(T.Node(posM, 1));

XX = [BX; BY];
KMoving = MovingK(FreeB, FreeB);
FBoundary = MovingF(FreeB) - MovingK(FreeB, ~FreeB)*XX(~FreeB);
XX(FreeB) = KMoving\FBoundary;

%T.Node = reshape(XX, NFp, 2);


% Test For the elasticity
%dt = 1e-3;
runtime = dt*J;

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

data = zeros(J, length(posM));
%set(figure(1),'visible','off');
for t = dt:dt:runtime
    FF = BigF + Sysold*Xold;
    FF = FF(Free) - Sys(Free, ~Free)*XX(~Free);
    XX(Free) = KK\FF;
    %BX = T.Node(:,1) + InterpolateP2toP1(T, XX(2*NSu+1:3*NSu));
    DY = InterpolateP2toP1(T, XX(3*NSu+1:4*NSu));
    data(round(t/dt), :) = DY(posM);
    %{
    triplot(T.Tri, BX, BY);
    axis equal
    %title(['t = ', num2str(t)]);
    box off;set(gca,'xtick',[],'ytick',[],'xcolor','white','ycolor','white');
    axis off
    pause(0.01)
    %}
    Xold = XX;
    %saveas(gca, ['lame', num2str(npic) ,'.jpg'])
end
end