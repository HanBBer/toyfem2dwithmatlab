dbstop if error;
% This is a Dirichlet Problem calculation demo

% Problem Preparation(truth and boundary condition)
u = @(x, y) exp(sin(x)+cos(y));
f = @(x, y) (cos(x).^2-sin(x)+sin(y).^2-cos(y)).*exp(sin(x)+cos(y));

% Space Define
nx = 2; ny = 2;
T = RecMesh(nx, ny, 1, 1, 0, 0);
T = DefineFespace(T, 'U', "P2");
Fd = FreedomDefine(T, 'U', [1,1,1,1]);

% System Clarification
K = FEMatrix(T, Fd, 'nabla');
G = {u, u, u, u};
F = FemBiLoad(T, Fd, 'nabla', G) + FemLinearLoad(T, Fd, f, []);

% Solve and visualization
U = K\F;
Z = zeros(T.N, 1);
Z(Fd.FNodePtrs) = U;
Z(Fd.NodePtrs<0) = u(T.U.Nodes(Fd.NodePtrs<0, 1), T.U.Nodes(Fd.NodePtrs<0, 2));
trisurf(T.U.TP, T.U.Nodes(:,1), T.U.Nodes(:,2), Z);
Ze = u(T.U.Nodes(:,1), T.U.Nodes(:,2));
figure(2)
trisurf(T.U.TP, T.U.Nodes(:,1), T.U.Nodes(:,2), Z-Ze);