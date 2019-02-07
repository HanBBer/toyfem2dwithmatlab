% Here is a NS demo

% System Info
n = 1;
dt = 0.02; nu = 0.7; t = 0; t1 = 10;

f1 = @(x,y) sin(x).*cos(y+t)+(2*nu-1)*sin(x).*sin(y+t);
f2 = @(x,y) -cos(x).*sin(y+t)+(2*nu+1)*cos(x).*cos(y+t);
U1 = @(x,y) sin(x).*sin(y+t);
U2 = @(x,y) cos(x).*cos(y+t);

% Define the domain and generate mesh
N = 8;
Nx = N; Ny = N;
T = RecMesh([Nx, Ny], [1,1], [0,0]);
U = P2Fespace(T);
P = P1Fespace(T);

% System Clarification
Mu = symBilinear(U, 'mass', []);
Ku = symBilinear(U, 'nabla', []);
Kuxp = nonsymBilinear(T, U, 'dx', P, 'mass', []);
Kuyp = nonsymBilinear(T, U, 'dy', P, 'mass', []);
Mp = symBilinear(P, 'mass', []);
G1 = {U1, U1, U1, []};
G2 = {U2, U2, U2, []};
[Xu1, FNu1] =  Freedomdefine(U, [1,1,1,0], G1);
[Xu2, FNu2] =  Freedomdefine(U, [1,1,1,0], G2);
[Xp, FNp] = Freedomdefine(P, [0,0,0,0], []);
Nu = U.N;
Np = P.N;

K = [Mu/dt+nu*Ku, sparse(Nu, Nu), -Kuxp;
    sparse(Nu, Nu), Mu/dt+nu*Ku, -Kuyp;
    -Kuxp', -Kuyp', 1e-10*Mp];
Freedom = [FNu1; FNu2; FNp];
Xold = [U1(U.Node(:,1), U.Node(:,2));
    U2(U.Node(:,1), U.Node(:,2));
    zeros(Np,1)];

for t = dt:dt:t1
    f1 = @(x,y) sin(x).*cos(y+t)+(2*nu-1)*sin(x).*sin(y+t);
    f2 = @(x,y) -cos(x).*sin(y+t)+(2*nu+1)*cos(x).*cos(y+t);
    
    %Truth
    U1 = @(x,y) sin(x).*sin(y+t); 
    U2 = @(x,y) cos(x).*cos(y+t);
    %pe = @(x,y) cos(x).*sin(y+t);
    
    G1 = {U1, U1, U1, []};
    G2 = {U2, U2, U2, []};
    
    F1 = Load(U, f1) + Mu*Xold(1:Nu)/dt;
    F2 = Load(U, f2) + Mu*Xold(Nu+1:2*Nu)/dt;
    F = [F1; F2; sparse(Np,1)];
    
    Xu1(~FNu1) = U1(U.Node(~FNu1, 1), U.Node(~FNu1, 2));
    Xu2(~FNu2) = U2(U.Node(~FNu2, 1), U.Node(~FNu2, 2));
    X = [Xu1; Xu2; Xp];
    
    KK = K(Freedom, Freedom);
    FF = F(Freedom) - K(Freedom, ~Freedom)*X(~Freedom);
    X(Freedom) = KK\FF;
    
    Xold = X;
    
    %{
    subplot(1,2,1)
    trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), X(1:Nu));
    subplot(1,2,2)
    trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), U1(U.Node(:, 1), U.Node(:, 2)));
    %}
    trisurf(P.Tri, P.Node(:, 1), P.Node(:, 2), X(2*Nu+1:end),...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
    view(2);
    box off; set(gca, 'XTick', [], 'YTick', []);
    colorbar;
    title(['t = ', num2str(t), 's']);
    pause(0.001)
end