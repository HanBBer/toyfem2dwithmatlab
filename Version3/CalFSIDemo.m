% Here is a FSI demo

% System Info


% Mesh Generator
TF = RecMesh([8, 8], [6, 6], [0, -6]);
TS = RecMesh([8, 2], [6, 1], [0, 0]);
UF = P2Fespace(TF);
PF = P1Fespace(TF);
US = P2Fespace(TS);
NFu = UF.N;
NFp = PF.N;
NSu = US.N;

% Moving Matrix

% Moving test 1
% Moving = symBilinear(PF, 'nabla', []);
% G = {@(x,y) -6, @(x,y) y, @(x, y) 2*sin(x/6*pi), @(x,y) y};
% [X, FNp] =  Freedomdefine(PF, [1,1,1,1], G);
% F = sparse(NFp, 1);
% KMoving = Moving(FNp, FNp);
% FF = F(FNp) - Moving(FNp, ~FNp)*X(~FNp);
% X(FNp) = KMoving\FF;
% triplot(PF.Tri, PF.Node(:, 1), X);

% Elastic Moving 
Kx = symBilinear(PF, 'dx', []);
Ky = symBilinear(PF, 'dy', []);
Kxy = nonsymBilinear(TF, PF, 'dx', PF, 'dy', []);
Moving = [Kx, Kxy; Kxy', Ky] + [2*Kx+Ky, Kxy'; Kxy, 2*Ky+Kx];
BigF = sparse(2*NFp, 1);

% Boundary Position
G = {@(x, y) x, @(x, y) x, @(x, y) x, @(x, y) x};
[X1, ~] = Freedomdefine(PF, [1,1,1,1], G);
G = {@(x,y) -6, @(x,y) y, @(x, y) 2*sin(x/6*pi), @(x,y) y};
[X2, FNodeptr] = Freedomdefine(PF, [1,1,1,1], G);

Boundary = [X1; X2];
Freedom = [FNodeptr; FNodeptr];
KMoving = Moving(Freedom, Freedom);
FBoundary = BigF(Freedom) - Moving(Freedom, ~Freedom)*Boundary(~Freedom);
Boundary(Freedom) = KMoving\FBoundary;

triplot(PF.Tri, Boundary(1:NFp, 1), Boundary(NFp+1:end, 1));

% Find the interaction part
posF = UF.Edge(UF.EgFlag == 3, :);
posF = unique(posF(:));
[~, n] = sort(UF.Node(posF, 1));
posF = posF(n);

posS = US.Edge(US.EgFlag == 1, :);
posS = unique(posS(:));
[~, n] = sort(US.Node(posS, 1));
posS = posS(n);

