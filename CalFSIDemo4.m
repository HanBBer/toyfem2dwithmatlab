% Here is a FSI demo

%% Infomation
% Source
%Traction = @(t) 2e4*(t<5e-3)*sin(pi*t/5e-3);
Traction = @(t) -1e4*(1-cos(pi*t/2.5e-3));

% System Info
dt = 1e-4;
runtime = 0.02;
rhoF = 1;
nuF = 0.035;
rhoS = 1.1;
muS = 5.75e6;
lambda = 1.7e6;

L = 6; RF = 0.6; RS = 0.1;

%% Initial Mesh
TF = RecMesh([10, 4], [L, RF], [0, 0]);
TS = RecMesh([10, 2], [L, RS], [0, RF]);
UF = P2Fespace(TF);
PF = P1Fespace(TF);
US = P2Fespace(TS);
PS = P1Fespace(TS);
NFu = UF.N;
NFp = PF.N;
NSu = US.N;
NSp = PS.N;

% System Matrix
% Langrangian Matrix
MS = symBilinear(US, 'mass', []);
Sx = symBilinear(US, 'dx', []);
Sy = symBilinear(US, 'dy', []);
Sxy = nonsymBilinear(TS, US, 'dx', US, 'dy', []);
Lame = lambda*[Sx, Sxy; Sxy', Sy] + muS*[2*Sx+Sy, Sxy'; Sxy, 2*Sy+Sx];
FSIS = [[rhoS/dt*MS, sparse(NSu, NSu); sparse(NSu, NSu), rhoS/dt*MS], Lame;
    -[MS, sparse(NSu, NSu); sparse(NSu, NSu), MS], [MS/dt, sparse(NSu, NSu); sparse(NSu, NSu), MS/dt]];


%% Preparation
% Orign Mesh
% Every moving iteration happens on this orign mesh
% which reduces the computational cost
TF0 = TF;
% Elastic Moving Matrix
Kx = symBilinear(PF, 'dx', []);
Ky = symBilinear(PF, 'dy', []);
Kxy = nonsymBilinear(TF, PF, 'dx', PF, 'dy', []);
LameK = 10*[Kx, Kxy; Kxy', Ky] + [2*Kx+Ky, Kxy'; Kxy, 2*Ky+Kx];
MovingF = sparse(2*NFp, 1);

% Boundary Condition for Elastic Moving
G = {@(x, y) 0, @(x, y) 0, @(x, y) 0, @(x, y) 0};
[MX, FNMoving] = Freedomdefine(PF, [1,1,1,1], G);
FreeMoving = [FNMoving; FNMoving];
BX = MX;
BY = MX;

KMoving = LameK(FreeMoving, FreeMoving);

% Solution Structure
G = {@(x, y)0, @(x, y)0, @(x, y)0, @(x, y)0};
[XFu1, FNFu1] = Freedomdefine(UF, [0,0,0,0], G);
[XFu2, FNFu2] = Freedomdefine(UF, [1,0,0,0], G);
[XFp, FNFp] = Freedomdefine(PF, [0,0,0,0], G);
[XSu, FNSu] = Freedomdefine(US, [0,1,0,1], G);
[XSd, FNSd] = Freedomdefine(US, [0,1,0,1], G);
XSd1 = XSd;
XSd2 = XSd;

% Indication of Interaction Interface
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

FNFu1(posFP2([1, end])) = boolean(0);
FNFu2(posFP2([1, end])) = boolean(0);
FreeFSI = [FNFu1; FNFu2; FNFp; FNSu; FNSu; FNSd; FNSd];
% Coupled Indication
IndFsiF = zeros(2*NFu+NFp+4*NSu, 1);
IndFsiF([posFP2; NFu+posFP2], 1) = 1;
IndFsiF = boolean(IndFsiF);
IndFsiS = zeros(2*NFu+NFp+4*NSu, 1);
IndFsiS([2*NFu+NFp+posFP2; 2*NFu+NFp+NSu+posFP2], 1) = 1;
IndFsiS = boolean(IndFsiS);
IndFsi = ones(2*NFu+NFp+4*NSu, 1);
IndFsi([2*NFu+NFp+posFP2; 2*NFu+NFp+NSu+posFP2], 1) = 0;
IndFsi = boolean(IndFsi);

%% Initial State
uf1old = sparse(NFu, 1); uf2old = uf1old;
p = sparse(NFp, 1);
us1old = sparse(NSu, 1); us2old = us1old;
ds1old = sparse(NSu, 1); ds2old = ds1old;
XFSIold = sparse(2*NFu+NFp+4*NSu, 1);

% Temporal Discretization
MF = symBilinear(UF, 'mass', []);
MFSIold = ...
    [rhoF/dt*MF, sparse(NFu, NFu+NFp+4*NSu);
    sparse(NFu, NFu), rhoF/dt*MF, sparse(NFu, NFp+4*NSu);
    sparse(NFp, 2*NFu+NFp+4*NSu);
    sparse(NSu, 2*NFu+NFp), rhoS/dt*MS, sparse(NSu, 3*NSu);
    sparse(NSu, 2*NFu+NFp+NSu), rhoS/dt*MS, sparse(NSu, 2*NSu);
    sparse(NSu, 2*NFu+NFp+2*NSu), MS/dt, sparse(NSu, NSu);
    sparse(NSu, 2*NFu+NFp+3*NSu), MS/dt;
    ];
%%

for t = 0:dt:runtime
    % 1. Boundary Moving
    TF1 = TF;
    BX(posFP1) = ds1old(posSP1) + us1old(posSP1)*dt;
    BY(posFP1) = ds2old(posSP1) + us2old(posSP1)*dt;
    XMoving = [BX; BY];
    FMoving = MovingF(FreeMoving) - LameK(FreeMoving, ~FreeMoving)*XMoving(~FreeMoving);
    XMoving(FreeMoving) = KMoving\FMoving;
    
    TF1.Node = TF0.Node + reshape(XMoving, NFp, 2);
    %figure(1)
    %triplot(TF1.Tri, TF1.Node(:,1), TF1.Node(:,2))
    %axis equal
    %pause()
    
    % 2. System Building
    % 2.1 Mesh Generate
    TFm = TF;
    TFm.Node = (TF1.Node + TF.Node)/2;
    TriSpeed = (TF1.Node - TF.Node)/dt;
    TriSpeed = [InterpolateP1toP2(TF, TriSpeed(:, 1)), InterpolateP1toP2(TF, TriSpeed(:, 2))];
    
    % 2.2 Space Define
    UF1 = P2Fespace(TF1);
    PF1 = P1Fespace(TF1);
    UFm = P2Fespace(TFm);
    PFm = P1Fespace(TFm);
    
    % 2.3 Matrix Define
    MF1 = symBilinear(UF1, 'mass', []);
    fnk = [uf1old, uf2old] - TriSpeed;
    MFp = symBilinear(PF1, 'mass', []);
    Mat2 = numP2Bilinear(UFm, fnk(:, 1), "mass", "dx") + numP2Bilinear(UFm, fnk(:, 1), "dx", "mass") ...
        + numP2Bilinear(UFm, fnk(:, 2), "mass", "dy") + numP2Bilinear(UFm, fnk(:, 2), "dy", "mass");
    Mat3 = numP2Bilinear(UFm, uf1old, "mass", "dx") + numP2Bilinear(UFm, uf2old, "mass", "dy");
    Fu = symBilinear(UF1, 'nabla', []);
    Fuxp = nonsymBilinear(TF1, UF1, 'dx', PF1, 'mass', []);
    Fuyp = nonsymBilinear(TF1, UF1, 'dy', PF1, 'mass', []);
    Bpu  = Bilinear1d(TF1, 1, PF1, "mass", UF1, "mass");
    Buxu = Bilinear1d(TF1, 1, UF1, "mass", UF1, "dx");
    
    % 2.4 System assemble
    MFSI1 = ...
        [rhoF/dt*MF1+nuF*Fu+rhoF*(Mat2-1/2*Mat3)+2*nuF*Buxu, sparse(NFu, NFu), -Fuxp-nuF*Bpu', sparse(NFu, 4*NSu);
        sparse(NFu, NFu), rhoF/dt*MF1+nuF*Fu+rhoF*(Mat2-1/2*Mat3)+nuF*Buxu, -Fuyp, sparse(NFu, 4*NSu);
        -Fuxp', -Fuyp', 1e-10*MFp, sparse(NFp, 4*NSu);
        sparse(4*NSu, 2*NFu+NFp), FSIS;
        ];
    
    Tr = @(x, y) Traction(t);
    FTr = Load2(UF1, @(x, y) 0, {[], @(x, y) 0, [], Tr});
    FFSI = [FTr; sparse(NFu+NFp+4*NSu, 1)] + MFSIold*XFSIold;
    
    % 3. Solve
    %XSd1(posSP2) = ds1old(posSP2) + us1old(posSP2)*dt;
    %XSd2(posSP2) = ds2old(posSP2) + us2old(posSP2)*dt;
    XFSI = [XFu1; XFu2; XFp; XSu; XSu; XSd1; XSd2];
    FFFSI = FFSI;
    FFFSI(FreeFSI, 1) = FFSI(FreeFSI, 1) - MFSI1(FreeFSI, ~FreeFSI)*XFSI(~FreeFSI);
    MFSI1(:, IndFsiF) = MFSI1(:, IndFsiF) + MFSI1(:, IndFsiS);
    MFSI1(IndFsiF, :) = MFSI1(IndFsiF, :) + MFSI1(IndFsiS, :);
    FFFSI(IndFsiF) = FFFSI(IndFsiF) + FFFSI(IndFsiS);
    BigM = MFSI1(FreeFSI & IndFsi, FreeFSI & IndFsi);
    BigF = FFFSI(FreeFSI & IndFsi, 1);
    XFSI(FreeFSI & IndFsi) = BigM\BigF;
    XFSI(IndFsiS) = XFSI(IndFsiF);
    
    % 4. Data Store and Mesh Modification
    TF = TF1;
    MF = MF1;
    MFSIold = ...
        [rhoF/dt*MF, sparse(NFu, NFu+NFp+4*NSu);
        sparse(NFu, NFu), rhoF/dt*MF, sparse(NFu, NFp+4*NSu);
        sparse(NFp, 2*NFu+NFp+4*NSu);
        sparse(NSu, 2*NFu+NFp), rhoS/dt*MS, sparse(NSu, 3*NSu);
        sparse(NSu, 2*NFu+NFp+NSu), rhoS/dt*MS, sparse(NSu, 2*NSu);
        sparse(NSu, 2*NFu+NFp+2*NSu), MS/dt, sparse(NSu, NSu);
        sparse(NSu, 2*NFu+NFp+3*NSu), MS/dt;
        ];
    XFSIold = XFSI;
    uf1old = XFSIold(1:NFu);
    uf2old = XFSIold(NFu+1:2*NFu);
    us1old = XFSIold(2*NFu+NFp+1:2*NFu+NFp+NSu);
    us2old = XFSIold(2*NFu+NFp+NSu+1:2*NFu+NFp+2*NSu);
    ds1old = XFSIold(2*NFu+NFp+2*NSu+1:2*NFu+NFp+3*NSu);
    ds2old = XFSIold(2*NFu+NFp+3*NSu+1:2*NFu+NFp+4*NSu);
    TF.Node(posFP1, :) = TF0.Node(posFP1, :) + [ds1old(posFP1), ds2old(posFP1)];
    
    
    % 5. Visulization
    %figure(2)
    hold on
    subplot(3,1,1)
    %quiver(UF1.Node(:, 1), UF1.Node(:, 2), XFSI(1:NFu), XFSI(NFu+1:2*NFu));
    triplot(TS.Tri, TS.Node(:, 1)+ds1old(1:NSp), TS.Node(:, 2)+ds2old(1:NSp));
    axis equal
    box off; set(gca, 'XTick', [], 'YTick', [], 'xcolor', 'white', 'ycolor', 'white');
    subplot(3,1,2)
    trisurf(TF.Tri, TF.Node(:, 1), TF.Node(:, 2), XFSI(2*NFu+1:2*NFu+NFp),...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
    axis equal
    box off; set(gca, 'XTick', [], 'YTick', [], 'xcolor', 'white', 'ycolor', 'white');
    view(2);
    colorbar
    
    subplot(3,1,3)
    plot(t, Traction(t), 'r.');
    axis([0, 0.02, -0.5e4, 2.5e4])
    hold off
    
    title({['t = ', num2str(t), 's']});
    pause(0.01);
end