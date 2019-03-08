% dbstop if error
% Here is a FSI demo

% Source
Traction = @(t) 2e4*(t<5e-3)*sin(pi*t/5e-3);

% System Info
dt = 1e-4;
runtime = 0.02;
% Lf = 6; Rf = 0.5;
% Ls = Lf; Rs = 0.1;
rhoF = 1;
muF = 0.035;
rhoS = 1.1;
epsilonS = 0.1;
E = 0.75*1e6;
nuS = 0.5;
muS = E*epsilonS/(2*(1+nuS));
lambda = E*epsilonS/((1+nuS)*(1-nuS));

% Mesh Generator
L = 12; RF = 4; RS = 1;
TF = RecMesh([6, 4], [L, RF], [0, -RF]);
TS = RecMesh([6, 2], [L, RS], [0, 0]);
UF = P2Fespace(TF);
PF = P1Fespace(TF);
US = P2Fespace(TS);
NFu = UF.N;
NFp = PF.N;
NSu = US.N;

% System Matrix
MF = symBilinear(UF, 'mass', []);
MS = symBilinear(US, 'mass', []);
Sx = symBilinear(US, 'dx', []);
Sy = symBilinear(US, 'dy', []);
Sxy = nonsymBilinear(TS, US, 'dx', US, 'dy', []);

Lame = lambda*[Sx, Sxy; Sxy', Sy] + muS*[2*Sx+Sy, Sxy'; Sxy, 2*Sy+Sx];
MFSI = ...
    [rhoF/dt*MF, sparse(NFu, NFu+NFp+4*NSu);
    sparse(NFu, NFu), rhoF/dt*MF, sparse(NFu, NFp+4*NSu);
    sparse(NFp, 2*NFu+NFp+4*NSu);
    sparse(NSu, 2*NFu+NFp), rhoS*epsilonS/dt*MS, sparse(NSu, 3*NSu);
    sparse(NSu, 2*NFu+NFp+NSu), rhoS*epsilonS/dt*MS, sparse(NSu, 2*NSu);
    sparse(NSu, 2*NFu+NFp+2*NSu), MS/dt, sparse(NSu, NSu);
    sparse(NSu, 2*NFu+NFp+3*NSu), MS/dt;
    ];
% Initial Boundary Condition
G = {@(x, y)0, @(x, y)0, @(x, y)0, @(x, y)0};
[XFu1, FNFu1] = Freedomdefine(UF, [1,0,0,0], G);
[XFu2, FNFu2] = Freedomdefine(UF, [1,0,0,0], G);
[XFp, FNFp] = Freedomdefine(PF, [0,0,0,0], G);
[XSu, FNSu] = Freedomdefine(US, [0,1,0,1], G);
[XSd1, FNSd] = Freedomdefine(US, [0,1,0,1], G);
XSd2 = XSd1;
FreeFSI = [FNFu1; FNFu2; FNFp; FNSu; FNSu; FNSd; FNSd];

% Elastic Moving Matrix
Kx = symBilinear(PF, 'dx', []);
Ky = symBilinear(PF, 'dy', []);
Kxy = nonsymBilinear(TF, PF, 'dx', PF, 'dy', []);
MovingK = [Kx, Kxy; Kxy', Ky] + [2*Kx+Ky, Kxy'; Kxy, 2*Ky+Kx];
MovingF = sparse(2*NFp, 1);

% Initial Boundary Position of Moving Matrix
G = {@(x, y) x, @(x, y) x, @(x, y) x, @(x, y) x};
[BX, ~] = Freedomdefine(PF, [1,1,1,1], G);
BX1 = BX;
BXP2 = InterpolateP1toP2(TF, BX);
G = {@(x,y) y, @(x,y) y, @(x, y) y, @(x,y) y};
[BY, FNMoving] = Freedomdefine(PF, [1,1,1,1], G);
FreeB = [FNMoving; FNMoving];
BY1 = BY;
BYP2 = InterpolateP1toP2(TF, BY);

% The index of free boundary
[~, FreeBind] = Freedomdefine(PF, [0,0,1,0], G);

% triplot(PF.Tri, Boundary(1:NFp, 1), Boundary(NFp+1:end, 1));

% Find the interaction part
posFP1 = TF.Edge(TF.EgFlag == 3, :);
posFP1 = unique(posFP1(:));
[~, n] = sort(TF.Node(posFP1, 1));
posFP1 = posFP1(n);

posSP1 = TS.Edge(TS.EgFlag == 1, :);
posSP1 = unique(posSP1(:));
[~, n] = sort(TS.Node(posSP1, 1));
posSP1 = posSP1(n);

posFP2 = UF.Edge(UF.EgFlag == 3, :);
posFP2 = unique(posFP2(:));
[~, n] = sort(UF.Node(posFP2, 1));
posFP2 = posFP2(n);

posSP2 = US.Edge(US.EgFlag == 1, :);
posSP2 = unique(posSP2(:));
[~, n] = sort(US.Node(posSP2, 1));
posSP2 = posSP2(n);

num = length(posFP2);
CoupleM = ...
    [sparse(1:num, posFP2, ones(num,1), num, NFu, num), sparse(num, NFu+NFp), sparse(1:num, posSP2, -ones(num,1), num, NSu, num), sparse(num, 3*NSu);
    sparse(num, NFu), sparse(1:num, posFP2, ones(num,1), num, NFu, num), sparse(num, NFp+NSu), sparse(1:num, posSP2, -ones(num,1), num, NSu, num), sparse(num, 2*NSu);
    ];
CoupleF = sparse(2*num, 1);

uf1old = sparse(NFu, 1); uf2old = uf1old;
p = sparse(NFp, 1);
us1old = sparse(NSu, 1); us2old = us1old;
ds1old = sparse(NSu, 1); ds2old = ds1old;
XFSIold = sparse(2*NFu+NFp+4*NSu, 1);
for t = dt:dt:runtime
    % Mesh Move
    TF1 = TF;
    BX1(posFP1) = BX(posFP1) + ds1old(posSP1);
    BY1(posFP1) = BY(posFP1) + ds2old(posSP1);
    Boundary = [BX1; BY1];
    KMoving = MovingK(FreeB, FreeB);
    FBoundary = MovingF(FreeB) - MovingK(FreeB, ~FreeB)*Boundary(~FreeB);
    Boundary(FreeB) = KMoving\FBoundary;
    TF1.Node = reshape(Boundary, NFp, 2);
    
    TFm = TF;
    TFm.Node = (TF1.Node + TF.Node)/2;
    TriSpeed = (TF1.Node - TF.Node)/dt;
    TriSpeed = [InterpolateP1toP2(TF, TriSpeed(:, 1)), InterpolateP1toP2(TF, TriSpeed(:, 2))];
    
    % System Build
    UF1 = P2Fespace(TF1);
    UFm = P2Fespace(TFm);
    PFm = P1Fespace(TFm);
    MF1 = symBilinear(UF1, 'mass', []);
    fnk = [uf1old, uf2old] - TriSpeed;
    MFp = symBilinear(PFm, 'mass', []);
    
    Mat2 = numP2Bilinear(UFm, fnk(:, 1), "mass", "dx") + numP2Bilinear(UFm, fnk(:, 1), "dx", "mass") ...
        + numP2Bilinear(UFm, fnk(:, 2), "mass", "dy") + numP2Bilinear(UFm, fnk(:, 2), "dy", "mass");
    Mat3 = numP2Bilinear(UFm, uf1old, "mass", "dx") + numP2Bilinear(UFm, uf2old, "mass", "dy");
    Fu = symBilinear(UFm, 'nabla', []);
    Fuxp = nonsymBilinear(TF1, UFm, 'dx', PFm, 'mass', []);
    Fuyp = nonsymBilinear(TF1, UFm, 'dy', PFm, 'mass', []);
    
    % MFSI = sparse(2*NFu+NFp+4*NSu, 2*NFu+NFp+4*NSu);
    MFSI1 = ...
        [rhoF/dt*MF1+muF*Fu+rhoF*(Mat2-Mat3), sparse(NFu, NFu), -Fuxp, sparse(NFu, 4*NSu);
        sparse(NFu, NFu), rhoF/dt*MF1+muF*Fu+rhoF*(Mat2-Mat3), -Fuyp, sparse(NFu, 4*NSu);
        -Fuxp', -Fuyp', 1e-10*MFp, sparse(NFp, 4*NSu);
        sparse(2*NSu, 2*NFu+NFp), [rhoS*epsilonS/dt*MS, sparse(NSu, NSu); sparse(NSu, NSu), rhoS*epsilonS/dt*MS], Lame;
        sparse(2*NSu, 2*NFu+NFp), -[MS, sparse(NSu, NSu); sparse(NSu, NSu), MS], [MS/dt, sparse(NSu, NSu); sparse(NSu, NSu), MS/dt];
        ];
    Tr = @(x, y) Traction(t);
    FTr = Load2(UFm, @(x, y) 0, {[], @(x,y) 0, [], Tr});
    FFSI = [FTr; sparse(NFu+NFp+4*NSu, 1)] + MFSI*XFSIold;
    
    %     Boundary = [InterpolateP1toP2(TF, Boundary(1:NFp, 1)), InterpolateP1toP2(TF, Boundary(NFp+1:end, 1))];
    %     XSd1(posSP2) = Boundary(posFP2, 1) - BXP2(posFP2, 1);
    %     XSd2(posSP2) = Boundary(posFP2, 2) - BYP2(posFP2, 1);
    XFSI = [XFu1; XFu2; XFp; XSu; XSu; XSd1; XSd2];
    
    BigM = [MFSI1(FreeFSI, FreeFSI); CoupleM(:, FreeFSI)];
    BigF = [FFSI(FreeFSI) - MFSI1(FreeFSI, ~FreeFSI)*XFSI(~FreeFSI); CoupleF];
    XFSI(FreeFSI) = BigM\BigF;
    
    % Save Info
    TF = TF1;
    MF = MF1;
    MFSI = ...
        [rhoF/dt*MF, sparse(NFu, NFu+NFp+4*NSu);
        sparse(NFu, NFu), rhoF/dt*MF, sparse(NFu, NFp+4*NSu);
        sparse(NFp, 2*NFu+NFp+4*NSu);
        sparse(NSu, 2*NFu+NFp), rhoS*epsilonS/dt*MS, sparse(NSu, 3*NSu);
        sparse(NSu, 2*NFu+NFp+NSu), rhoS*epsilonS/dt*MS, sparse(NSu, 2*NSu);
        sparse(NSu, 2*NFu+NFp+2*NSu), MS/dt, sparse(NSu, NSu);
        sparse(NSu, 2*NFu+NFp+3*NSu), MS/dt;
        ];
    XFSIold = XFSI;
    XFSIv = {XFSI(1:NFu), XFSI(NFu+1:2*NFu), XFSI(2*NFu+(1:NFp)), ...
        XFSI(2*NFu+NFp+(1:NSu)), XFSI(2*NFu+NFp+NSu+(1:NSu)), XFSI(2*NFu+NFp+2*NSu+(1:NSu)), XFSI(2*NFu+NFp+3*NSu+(1:NSu))};
    ds1old = XFSIold(2*NFu+NFp+2*NSu+1:2*NFu+NFp+3*NSu);
    ds2old = XFSIold(2*NFu+NFp+3*NSu+1:2*NFu+NFp+4*NSu);
    
    % Plot
    %{
    % Mesh of Fluid
    triplot(TF1.Tri, TF1.Node(:, 1), TF1.Node(:, 2));
    axis equal
    box off; set(gca, 'XTick', [], 'YTick', []);
    %}
    
    % Figure of Pressure
    hold on
    subplot(3,1,1)
    trisurf(TF.Tri, TF.Node(:, 1), TF.Node(:, 2), XFSI(2*NFu+1:2*NFu+NFp),...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
    axis equal
    box off; set(gca, 'XTick', [], 'YTick', []);
    view(2);
    colorbar
    err = CoupleM*XFSI;
    title(['residual = ', num2str(norm( err ))]);
    subplot(3,1,2)
    quiver(UF1.Node(:, 1), UF1.Node(:, 2), XFSI(1:NFu), XFSI(NFu+1:2*NFu));
    axis equal
    box off; set(gca, 'XTick', [], 'YTick', []);
    subplot(3,1,3)
    plot(t, Traction(t), 'r.');
    axis([0, 0.02, -0.5e4, 2.5e4])
    hold off
    
    title({['t = ', num2str(t), 's']});
    pause(0.01);
end