dbstop if error;
% This is a Dirichlet Problem calculation demo which is also a bad case
% verfied by Freefem++

% Problem Preparation(truth and boundary condition)
%{
bad case
u = @(x, y) x.^2.*y.*(1-x);
f = @(x, y) u(x,y)+2*y-6*x.*y;
%}

f = @(x, y) x.*y;

% Space Define
nx = 2; ny = 2;
T = RecMesh(nx, ny, 1, 1, 0, 0);
T = DefineFespace(T, 'U', "P2");
Fd = FreedomDefine(T, 'U', [1,1,1,1]);

% System Clarification
K = FEMatrix(T, Fd, 'nabla') + FEMatrix(T, Fd, 'mass');
G = {@(x,y)0, @(x,y)0, @(x,y)sin(pi*x), @(x,y)0};
F = FemBiLoad(T, Fd, 'nabla', G) + FemBiLoad(T, Fd, 'mass', G) + FemLinearLoad(T, Fd, f, []);

% Solve and visualization
U = K\F;
Z = zeros(T.N, 1);
Z(Fd.FNodePtrs) = U;
Z(T.U.Nodes(:,2)==1) = sin(T.U.Nodes(T.U.Nodes(:,2)==1,1)*pi);
trisurf(T.U.TP, T.U.Nodes(:,1), T.U.Nodes(:,2), Z);
