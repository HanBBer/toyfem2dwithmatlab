% Here is a NS demo

% System Info
n = 1;
rhoF = 1;nuF = 0.7;
dt = 0.02; t1 = 10;

t = 0;
f1 = @(x,y) sin(x).*cos(y+t)+(2*nuF-1)*sin(x).*sin(y+t);
f2 = @(x,y) -cos(x).*sin(y+t)+(2*nuF+1)*cos(x).*cos(y+t);
U1 = @(x,y) sin(x).*sin(y+t);
U2 = @(x,y) cos(x).*cos(y+t);


% Define the domain and generate mesh
N = 4;
Nx = N; Ny = N;
TF = RecMesh([Nx, Ny], [1,1], [0,0]);
UF = P2Fespace(TF);
PF = P1Fespace(TF);
NFu = UF.N;
NFp = PF.N;

TF0 = TF;

G = {@(x, y)0, @(x, y)0, @(x, y)0, @(x, y)0};
[BX, FreeB] = Freedomdefine(PF, [1,1,1,1], G);
NB = sum(FreeB);
A = 0.08;

MF = symBilinear(UF, 'mass', []);

G1 = {U1, U1, U1, []};
G2 = {U2, U2, U2, []};
[Xu1, FNu1] =  Freedomdefine(UF, [1,1,1,0], G1);
[Xu2, FNu2] =  Freedomdefine(UF, [1,1,1,0], G2);
[Xp, FNp] = Freedomdefine(PF, [0,0,0,0], []);
Freedom = [FNu1; FNu2; FNp];
Xold = [U1(UF.Node(:,1), UF.Node(:,2));
    U2(UF.Node(:,1), UF.Node(:, 2));
    zeros(NFp,1)];
uf1old = Xold(1:NFu);
uf2old = Xold(NFu+1:2*NFu);

for t = 0:dt:t1
    f1 = @(x,y) sin(x).*cos(y+t)+(2*nuF-1)*sin(x).*sin(y+t);
    f2 = @(x,y) -cos(x).*sin(y+t)+(2*nuF+1)*cos(x).*cos(y+t);
    U1 = @(x,y) sin(x).*sin(y+t); 
    U2 = @(x,y) cos(x).*cos(y+t);
    pe = @(x,y) cos(x).*sin(y+t);
    
    % Move Mesh
    TF1 = TF0;
    TF1.Node(FreeB, :) = TF1.Node(FreeB, :) + A*rand(NB, 2) - A/2;
    %triplot(TF1.Tri, TF1.Node(:,1), TF1.Node(:,2))
    %pause(0.01)
    
    % Midean Mesh
    TFm = TF0;
    TFm.Node = (TF1.Node + TF.Node)/2;
    %triplot(TFm.Tri, TFm.Node(:,1), TFm.Node(:,2))
    TriSpeed = (TF1.Node - TF.Node)/dt;
    TriSpeed = [InterpolateP1toP2(TF, TriSpeed(:, 1)), InterpolateP1toP2(TF, TriSpeed(:, 2))];
    
    UF1 = P2Fespace(TF1);
    PF1 = P1Fespace(TF1);
    UFm = P2Fespace(TFm);
    PFm = P1Fespace(TFm);
    
    MF1 = symBilinear(UF1, 'mass', []);
    fnk = -TriSpeed;
    MFp = symBilinear(PF1, 'mass', []);    
    Mat2 = numP2Bilinear(UFm, fnk(:, 1), "mass", "dx") + numP2Bilinear(UFm, fnk(:, 1), "dx", "mass") ...
        + numP2Bilinear(UFm, fnk(:, 2), "mass", "dy") + numP2Bilinear(UFm, fnk(:, 2), "dy", "mass");
    Mat3 = numP2Bilinear(UFm, uf1old, "mass", "dx") + numP2Bilinear(UFm, uf2old, "mass", "dy");
    Fu = symBilinear(UF1, 'nabla', []);
    Fuxp = nonsymBilinear(TF1, UF1, 'dx', PF1, 'mass', []);
    Fuyp = nonsymBilinear(TF1, UF1, 'dy', PF1, 'mass', []);


    K = ...
        [rhoF/dt*MF1+nuF*Fu+rhoF*(Mat2-1/2*Mat3), sparse(NFu, NFu), -Fuxp;
        sparse(NFu, NFu), rhoF/dt*MF1+nuF*Fu+rhoF*(Mat2-1/2*Mat3), -Fuyp;
        -Fuxp', -Fuyp', sparse(NFp, NFp);
        ];
    
    F1 = Load(UF, f1) + MF*Xold(1:NFu)/dt;
    F2 = Load(UF, f2) + MF*Xold(NFu+1:2*NFu)/dt;
    F = [F1; F2; sparse(NFp,1)];
    
    Xu1(~FNu1) = U1(UF.Node(~FNu1, 1), UF.Node(~FNu1, 2));
    Xu2(~FNu2) = U2(UF.Node(~FNu2, 1), UF.Node(~FNu2, 2));
    X = [Xu1; Xu2; Xp];
    
    KK = K(Freedom, Freedom);
    FF = F(Freedom) - K(Freedom, ~Freedom)*X(~Freedom);
    X(Freedom) = KK\FF;
    
    % Record Info
    Xold = X;
    uf1old = Xold(1:NFu);
    uf2old = Xold(NFu+1:2*NFu);
    TF = TF1;
    MF = MF1;
    
    if ~mod(t, 0.1) 
        %{
        subplot(1,2,1)
        trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), X(1:Nu));
        subplot(1,2,2)
        trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), U1(U.Node(:, 1), U.Node(:, 2)));
        suptitle({['t = ', num2str(t), 's']});
        %}
        figure(1)
        trisurf(PF1.Tri, PF1.Node(:, 1), PF1.Node(:, 2), X(2*NFu+1:end),...
            'FaceColor', 'interp', 'EdgeColor', 'interp', 'facealpha', 0.5, 'edgealpha', 0);
        colorbar;
        hold on
        quiver(UF1.Node(:, 1), UF1.Node(:, 2), X(1:NFu), X(NFu+1:2*NFu), "color", "k");
        hold off
        view(2);
        box off; set(gca, 'XTick', [], 'YTick', []);
        err = [MF1, sparse(NFu, NFu);sparse(NFu, NFu), MF1]*[X(1:NFu)-U1(UF1.Node(:, 1), UF1.Node(:, 2)); X(NFu+1:2*NFu)-U2(UF1.Node(:, 1), UF1.Node(:, 2))];
        title({['t = ', num2str(t), 's'], ['residual = ', num2str(norm(err, 2))]});
        axis equal
        axis([-0.1, 1.1, -0.1, 1.1])
        set(gca, "xcolor", "white", "ycolor", "white")
        
        figure(2)
        trisurf(PF1.Tri, PF1.Node(:, 1), PF1.Node(:, 2), pe(PF1.Node(:, 1), PF1.Node(:, 2)),...
            'FaceColor', 'interp', 'EdgeColor', 'interp', 'facealpha', 0.5, 'edgealpha', 0);
        colorbar;
        axis equal
        axis([-0.1, 1.1, -0.1, 1.1])
        set(gca, "xcolor", "white", "ycolor", "white")
        pause(0.01)
    end
end

%{
% System Clarification
MF = symBilinear(UF, 'mass', []);
Ku = symBilinear(UF, 'nabla', []);
Kuxp = nonsymBilinear(TF, UF, 'dx', PF, 'mass', []);
Kuyp = nonsymBilinear(TF, UF, 'dy', PF, 'mass', []);
Mp = symBilinear(PF, 'mass', []);
G1 = {U1, U1, U1, []};
G2 = {U2, U2, U2, []};
[Xu1, FNu1] =  Freedomdefine(UF, [1,1,1,0], G1);
[Xu2, FNu2] =  Freedomdefine(UF, [1,1,1,0], G2);
[Xp, FNp] = Freedomdefine(PF, [0,0,0,0], []);
NFu = UF.N;
NFp = PF.N;

K = [MF/dt+nuF*Ku, sparse(NFu, NFu), -Kuxp;
    sparse(NFu, NFu), MF/dt+nuF*Ku, -Kuyp;
    -Kuxp', -Kuyp', 1e-10*Mp];
Freedom = [FNu1; FNu2; FNp];
Xold = [U1(UF.Node(:,1), UF.Node(:,2));
    U2(UF.Node(:,1), UF.Node(:, 2));
    zeros(NFp,1)];

for t = dt:dt:t1
    f1 = @(x,y) sin(x).*cos(y+t)+(2*nuF-1)*sin(x).*sin(y+t);
    f2 = @(x,y) -cos(x).*sin(y+t)+(2*nuF+1)*cos(x).*cos(y+t);
    
    %Truth
    U1 = @(x,y) sin(x).*sin(y+t); 
    U2 = @(x,y) cos(x).*cos(y+t);
    %pe = @(x,y) cos(x).*sin(y+t);
    
    G1 = {U1, U1, U1, []};
    G2 = {U2, U2, U2, []};
    
    TF1 = TF;
    TFm = TF;
    TFm.Node = (TF1.Node + TF.Node)/2;
    TriSpeed = (TF1.Node - TF.Node)/dt;
    TriSpeed = [InterpolateP1toP2(TF, TriSpeed(:, 1)), InterpolateP1toP2(TF, TriSpeed(:, 2))];
    
    % 2.2 Space Define
    UF1 = P2Fespace(TF1);
    UFm = P2Fespace(TFm);
    PFm = P1Fespace(TFm);
    
    % 2.3 Matrix Define
    %MF1 = symBilinear(UF1, 'mass', []);
    fnk = [Xold(1:NFu), Xold(NFu+1:2*NFu)] - TriSpeed;
    %MFp = symBilinear(PFm, 'mass', []);    
    Mat2 = numP2Bilinear(UFm, fnk(:, 1), "mass", "dx") + numP2Bilinear(UFm, fnk(:, 1), "dx", "mass") ...
        + numP2Bilinear(UFm, fnk(:, 2), "mass", "dy") + numP2Bilinear(UFm, fnk(:, 2), "dy", "mass");
    Mat3 = numP2Bilinear(UFm, Xold(1:NFu), "mass", "dx") + numP2Bilinear(UFm, Xold(NFu+1:2*NFu), "mass", "dy");
    Fu = symBilinear(UFm, 'nabla', []);
    Fuxp = nonsymBilinear(TF1, UFm, 'dx', PFm, 'mass', []);
    Fuyp = nonsymBilinear(TF1, UFm, 'dy', PFm, 'mass', []);

    
    K = ...
        [rhoF/dt*MF+nuF*Fu+rhoF*(Mat2-Mat3), sparse(NFu, NFu), -Fuxp;
        sparse(NFu, NFu), rhoF/dt*MF+nuF*Fu+rhoF*(Mat2-Mat3), -Fuyp;
        -Fuxp', -Fuyp', 1e-10*Mp
        ];
    
    
    F1 = Load(UF, f1) + MF*Xold(1:NFu)/dt;
    F2 = Load(UF, f2) + MF*Xold(NFu+1:2*NFu)/dt;
    F = [F1; F2; sparse(NFp,1)];
    
    Xu1(~FNu1) = U1(UF.Node(~FNu1, 1), UF.Node(~FNu1, 2));
    Xu2(~FNu2) = U2(UF.Node(~FNu2, 1), UF.Node(~FNu2, 2));
    X = [Xu1; Xu2; Xp];
    
    KK = K(Freedom, Freedom);
    FF = F(Freedom) - K(Freedom, ~Freedom)*X(~Freedom);
    X(Freedom) = KK\FF;
    
    Xold = X;
    

    if ~mod(t, 0.1) 
        %{
        subplot(1,2,1)
        trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), X(1:Nu));
        subplot(1,2,2)
        trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), U1(U.Node(:, 1), U.Node(:, 2)));
        suptitle({['t = ', num2str(t), 's']});
        %}
        trisurf(PF.Tri, PF.Node(:, 1), PF.Node(:, 2), X(2*NFu+1:end),...
            'FaceColor', 'interp', 'EdgeColor', 'interp', 'facealpha', 0.5, 'edgealpha', 0);
        % colorbar;
        hold on
        quiver(UF.Node(:, 1), UF.Node(:, 2), X(1:NFu), X(NFu+1:2*NFu), "color", "k");
        hold off
        view(2);
        box off; set(gca, 'XTick', [], 'YTick', []);
        err = K*X-F;
        title({['t = ', num2str(t), 's'], ['residual = ', num2str(norm(err(Freedom), 2))]});
        axis equal
        axis([-0.1, 1.1, -0.1, 1.1])
        set(gca, "xcolor", "white", "ycolor", "white")
        pause(0.01)
    end
end
%}