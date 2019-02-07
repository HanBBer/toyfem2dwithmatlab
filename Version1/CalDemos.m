dbstop if error;
% This is some experiments to verify the validity
% all the computation is done on the ]0,1[^2

% Space Define
N = 2;
nx = N; ny = N;
T = RecMesh(nx, ny, 1, 1, 0, 0);
T = DefineFespace(T, 'U', "P2");
T = DefineFespace(T, 'P', "P1");


%% 
% The first problem is a simple Laplace probelm with Dirichlet BC
U10 = @(x, y) 0; U11 = @(x,y) sin(pi*x);
f1 = @(x, y) 1;

Fd1 = FreedomDefine(T, 'U', [1,1,1,1]);

K1 = FEMatrix(T, Fd1, 'nabla');
G1 = {U10, U10, U11, U10};
F1 = FemBiLoad(T, Fd1, 'nabla', G1) + FemLinearLoad(T, Fd1, f1, []);

U1 = K1\F1;

%figure(1);
%Z1 = SimplePlot(U1,T,Fd1,G1,1);

%%
% The Second problem is u = f
u2 = @(x,y) x.^2+3*x.*y;

Fd2 = FreedomDefine(T, 'U', [1,1,1,1]);
K2 = FEMatrix(T, Fd2, 'mass');
G2 = {u2,u2,u2,u2};
F2 = FemBiLoad(T, Fd2, 'mass', G2) + FemLinearLoad(T, Fd2, u2, []);

U2 = K2\F2;

%figure(2);
%Z2 = SimplePlot(U2,T,Fd2,G2,1);
%figure(3);
%Ze2 = u2(T.U.Nodes(:,1), T.U.Nodes(:, 2));
%trisurf(T.U.TP, T.U.Nodes(:,1), T.U.Nodes(:,2), Z2-Ze2);

%%
% Test Demo
g3 = @(x, y) 0;

Fd3 = FreedomDefine(T, 'U', [1,1,1,1]);
K3 = FEMatrix(T, Fd3, 'nabla');
G3 = {g3, g3, g3, g3};
F3 = FemLinearLoad(T, Fd3, @(x, y) 1, []);

U3 = K3\F3;

figure(4);
Z3 = SimplePlot(U3, T, Fd3, G3, 1);