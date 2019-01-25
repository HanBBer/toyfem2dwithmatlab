nx = 2; ny = 2;
T = RecMesh(nx, ny, 1, 1, 0, 0);
T = DefineFespace(T, 'U', "P2");
% T.U.Tri = T.U.TP;
% T.U.N = size(T.U.Nodes, 1); T.U.Nt = size(T.U.TP, 1); T.U.Ne = T.Ne*4;
% ShowMesh(T.U, [1, 1, 0]);
Fd = FreedomDefine(T, 'U', [1,1,1,1]);


K = FEMatrix(T, Fd, 'nabla');
%[w, P, Px, Py, C1, C]=LoadQuad();

% Compute the Load vector
g3 = @(x,y) 1+0.3*sin(pi*x);
g24 = @(x,y) y;
G = {[], g24, g3, g24};
F = FemBiLoad(T, Fd, 'nabla', G);

U = K\F;
Z = zeros(T.N, 1);
Z(Fd.FNodePtrs) = U;
Z(Fd.NodeFlag==2) = T.U.Nodes(Fd.NodeFlag==2, 2);
Z(Fd.NodeFlag==3) = g3(T.U.Nodes(Fd.NodeFlag==3 ,1));
Z(Fd.NodeFlag==4) = T.U.Nodes(Fd.NodeFlag==4, 2);

figure(2); trisurf(T.U.TP, T.U.Nodes(:,1), T.U.Nodes(:,2), Z);