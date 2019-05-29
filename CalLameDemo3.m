% This is a cal demo of lame system

rng(10);
saveplotflag = 0;
npic = 0;
if saveplotflag
    set(figure(1),'visible','off');
else
    set(figure(1),'visible','on');
end
    
% Space Define
Sigma = 0.1;
Gamma = 0.1;
N = 1;
nx = 10*N; ny = 5*N;
T = RecMesh([nx, ny], [20, 2], [0, -1]);
P = P1Fespace(T);
U = P2Fespace(T);

NSp = P.N;
NSu = U.N;

% System Info

rhoS = 1.1;
mu = 5.75e4;
lambda = 1.7e4;


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
g = @(x) 0.3*sin(pi*x/10);
BY(posM) = g(T.Node(posM, 1));

posM = T.Edge(T.EgFlag == 3, :);
posM = unique(posM(:));
g = @(x) 0.3*sin(pi*x/10);
BY(posM) = g(T.Node(posM, 1));

posM = T.Edge(T.EgFlag == 3, :);
posM = unique(posM(:));

[~, FDrift] = Freedomdefine(P, [0,1,0,1], G);
FDrift = [FDrift; FDrift];

XX = [BX; BY];
KMoving = MovingK(FreeB, FreeB);
FBoundary = MovingF(FreeB) - MovingK(FreeB, ~FreeB)*XX(~FreeB);
XX(FreeB) = KMoving\FBoundary;
triplot(T.Tri, T.Node(:,1)+XX(1:NSp), T.Node(:,2)+XX(NSp+1:end));
axis([0, 20, -2.5, 2.5])
axis equal
pause(0.1)

XX(FDrift) = XX(FDrift) + randn(sum(FDrift), 1);

%T.Node = reshape(XX, NFp, 2);

triplot(T.Tri, T.Node(:,1)+XX(1:NSp), T.Node(:,2)+XX(NSp+1:end));
axis([0, 20, -2.5, 2.5])
axis equal
pause(0.1)
%saveas(gca, ['lame_distort', num2str(npic) ,'.jpg'])
% Test For the elasticity
dt = 1e-3;
runtime = 4e-2;

DXold = InterpolateP1toP2(T, XX(1:NSp));
DYold = InterpolateP1toP2(T, XX(NSp+1:end));

G = {[], @(x, y) 0, [], @(x, y) 0};
[X, FNodeptr] = Freedomdefine(U, [0,1,0,1], G);
Kx = symBilinear(U, 'dx', []);
Ky = symBilinear(U, 'dy', []);
Kxy = nonsymBilinear(T, U, 'dx', U, 'dy', []);

MS = symBilinear(U, 'mass', []);
H1Error = MS+Kxy;
H1Error = [H1Error, sparse(NSu, NSu); sparse(NSu, NSu), H1Error];
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
FMatrix = KK\Sysold(Free,Free);
H1 = zeros(T.N, 4*NSu);
H1(:,3*NSu+1:3*NSu+T.N) = eye(T.N);
H2 = zeros(length(posM), T.N);
H2(:, posM) = eye(length(posM));
HMatrix = H2*H1;

% Prior
m0 = 0.1*ones(sum(Free), 1);
C0 = eye(sum(Free));

data = cell(runtime/dt+1, 1);
iter = 1;
data{iter}.t = 0;
data{iter}.x = Xold;
data{iter}.m0 = m0;
data{iter}.C0 = C0;
data{1}.T = T;
data{1}.xpos = linspace(0,20,length(posM));
data{1}.H1Error = H1Error;

for t = dt:dt:runtime
    %FF = Sysold*Xold;
    %FF = FF(Free) - Sys(Free, ~Free)*XX(~Free);
    %XX(Free) = KK\FF;
    
    % Predict 
    m0hat = FMatrix*m0;
    C0hat = FMatrix*C0*FMatrix'+Sigma^2*eye(sum(Free));
    
    XX(Free) = FMatrix*Xold(Free);
    XX([boolean(zeros(2*NSu,1));FDrift]) = XX([boolean(zeros(2*NSu,1));FDrift]) + Sigma*randn(sum(FDrift), 1);
    BX = T.Node(:,1) + InterpolateP2toP1(T, XX(2*NSu+1:3*NSu));
    BY = T.Node(:,2) + InterpolateP2toP1(T, XX(3*NSu+1:4*NSu));
    
    % Observation
    Y = HMatrix*XX + Gamma*randn(length(posM), 1);
    
    % Correct
    d = Y - HMatrix(:, Free)*m0hat;
    S = HMatrix(:, Free)*C0hat*HMatrix(:, Free)' + Gamma^2*eye(length(posM));
    K = C0hat*HMatrix(:, Free)'/S;
    m0 = m0hat + K*d;
    C0 = (eye(sum(Free)) - K*HMatrix(:, Free))*C0hat;
    
    figure(1)
    subplot(2,1,1)
    triplot(T.Tri, BX, BY);
    axis([0, 20, -2.5, 2.5])
    %box off;set(gca,'xtick',[],'ytick',[],'xcolor','white','ycolor','white');
    %axis off
    Xold = XX;
    title('State')
    
    subplot(2,1,2)
    %DY = InterpolateP2toP1(T, XX(3*NSu+1:4*NSu));
    %plot(linspace(0,20,length(posM)), DY(posM)' + Gamma*randn(1, length(posM)));
    
    plot(linspace(0,20,length(posM)), Y');
    axis([0, 20, -1, 1]);
    title('Measurement')
    
    npic = npic + 1;
    if saveplotflag
        saveas(gca, ['lame', num2str(npic) ,'.jpg'])
    end
    %pause(0.1)
    
    %figure(2)
    XX(Free) = m0;
    
    vErr = Xold - XX;
    vErr = vErr(2*NSu+1:4*NSu);
    
    % Record data
    iter = iter + 1;
    data{iter}.t = t;
    data{iter}.x = Xold;
    data{iter}.mhat = m0hat;
    data{iter}.Chat = C0hat;
    data{iter}.m0 = m0;
    data{iter}.C0 = C0;
    data{iter}.vErr = vErr;
    
end

% Process
H1Error = data{1}.H1Error;
Error = zeros(40,1);
for i = 2:41
    Error(i-1) = sqrt((data{i}.vErr)'*H1Error*data{i}.vErr + 0.01*trace(H1Error));
end
plot(1:40, Error)
hold on 
plot([0,40],sqrt(0.01*trace(H1Error))*[1,1], '--r')
xlabel('Iteration')