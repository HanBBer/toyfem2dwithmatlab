% Here is a NS demo with a traction

% Source
Traction = @(t) 2e4*(t<5e-3)*sin(pi*t/5e-3);

% System Info
dt = 1e-4;
runtime = 0.02;
nu = 0.7;


% Define the domain and generate mesh
L = 12; RF = 4; RS = 1;
T = RecMesh([12, 8], [L, RF], [0, -RF]);
U = P2Fespace(T);
P = P1Fespace(T);

% System Clarification
Mu = symBilinear(U, 'mass', []);
Ku = symBilinear(U, 'nabla', []);
Kuxp = nonsymBilinear(T, U, 'dx', P, 'mass', []);
Kuyp = nonsymBilinear(T, U, 'dy', P, 'mass', []);
Mp = symBilinear(P, 'mass', []);
G = {@(x, y) 0, [], @(x,y) 0, []};
[Xu1, FNu1] =  Freedomdefine(U, [1,0,1,0], G);
[Xu2, FNu2] =  Freedomdefine(U, [1,0,1,0], G);
[Xp, FNp] = Freedomdefine(P, [0,0,0,0], []);
Nu = U.N;
Np = P.N;

K = [Mu/dt+nu*Ku, sparse(Nu, Nu), -Kuxp;
    sparse(Nu, Nu), Mu/dt+nu*Ku, -Kuyp;
    -Kuxp', -Kuyp', 1e-10*Mp];
Freedom = [FNu1; FNu2; FNp];
Xold = zeros(2*Nu+Np,1);

for t = dt:dt:runtime

    Tr = @(x, y) Traction(t);
    G =  {[], @(x,y) 0, [], Tr};
    
    F1 = Load2(U, @(x, y) 1, G) + Mu*Xold(1:Nu)/dt;
    F2 = Load(U) + Mu*Xold(Nu+1:2*Nu)/dt;
    F = [F1; F2; sparse(Np,1)];
    
    X = [Xu1; Xu2; Xp];
    
    KK = K(Freedom, Freedom);
    FF = F(Freedom) - K(Freedom, ~Freedom)*X(~Freedom);
    X(Freedom) = KK\FF;
    
    Xold = X;
    

    if 1
        %{
        subplot(1,2,1)
        trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), X(1:Nu));
        subplot(1,2,2)
        trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), U1(U.Node(:, 1), U.Node(:, 2)));
        suptitle({['t = ', num2str(t), 's']});
        %}
        hold on
        subplot(3,1,1)
        trisurf(P.Tri, P.Node(:, 1), P.Node(:, 2), X(2*Nu+1:end),...
            'FaceColor', 'interp', 'EdgeColor', 'interp');
        view(2);
        axis equal
        box off; set(gca, 'XTick', [], 'YTick', []);
        colorbar;
        subplot(3,1,2)
        quiver(U.Node(:, 1), U.Node(:, 2), X(1:Nu), X(Nu+1:2*Nu));
        axis equal
        box off; set(gca, 'XTick', [], 'YTick', []);
        res = K*X-F;
        title({['t = ', num2str(t), 's'], ['residual = ', num2str(norm(res(Freedom), 2))]});
        subplot(3,1,3)
        plot(t, Traction(t), 'r.');
        axis([0, 0.02, -0.5e4, 2.5e4]);
        hold off
        
        pause(0.01)
    end
end