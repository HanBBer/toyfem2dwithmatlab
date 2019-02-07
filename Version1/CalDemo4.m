% This is a Stokes Problem calculation demo

% Boundary condition
U11 = @(x,y) 0; U12 = @(x,y) 1; % on 1,2,4  3 respectively 
U2 = @(x,y) 0;
G1 = {U11, U11, U12, U11};
G2 = {U2, U2, U2, U2};

% Space Define
nx = 10; ny = 10;
T = RecMesh(nx, ny, 2, 2, -1, -1);
T = DefineFespace(T, 'U', "P2");
Fdu = FreedomDefine(T, 'U', [1,1,1,1]);
Nfu = size(Fdu.FNodePtrs,1);
T = DefineFespace(T, 'P', "P1");
Fdp = FreedomDefine(T, 'P', [0,0,0,0]);
Nfp = size(Fdp.FNodePtrs,1);
ind1 = 1; ind2 = ind1+Nfu; ind3 = ind2+Nfu; ind4 = ind3+Nfp;

% Operator Define
Ku = FEMatrix(T, Fdu, "nabla");
Kuxp = FEMatrix(T, Fdu, "dx", Fdp, "mass");
Kuyp = FEMatrix(T, Fdu, "dy", Fdp, "mass");
Mp = FEMatrix(T, Fdp, "mass");

% System Clarification
BigK = [-Ku, zeros(Nfu,Nfu), Kuxp
    zeros(Nfu,Nfu), -Ku, Kuyp
    Kuxp', Kuyp', -1e-10*Mp]; 

BigF = zeros(2*Nfu+Nfp,1);
BigF(ind1:ind2-1) = FemBiLoad(T, Fdu, 'nabla', G1);
BigF(ind2:ind3-1) = FemBiLoad(T, Fdu, 'nabla', G2);
BigF(ind3:end) = FemBiLoad(T, Fdp, 'mass', G1, Fdu, 'dx') + FemBiLoad(T, Fdp, 'mass', G2, Fdu, 'dy');

% Solve and visualization
spparms('umfpack',1);
U = BigK\BigF;
%U = bicg(BigK, BigF, 1e-6, 50);
Z = zeros(Nfp, 1);
Z(Fdp.FNodePtrs) = U(ind3:ind4-1);
%Z(Fdu.NodePtrs<0) = U2(T.U.Nodes(Fdu.NodePtrs<0, 1),
%T.U.Nodes(Fdu.NodePtrs<0, 2));
trisurf(T.P.TP, T.P.Nodes(:,1), T.P.Nodes(:,2), Z-mean(Z));
%{
up1 = U(ind1:ind2-1); up2 = U(ind2:ind3-1);
figure(2)
hold on
for i = 1:size(T.P.Nodes,1)
    cord = T.P.Nodes(i,:);
    plot(cord(1),cord(2),'k+');
    plot([cord(1),cord(1)+up1(i)],[cord(2),cord(2)+up2(i)],'b');
end
axis([0,1,0,1])
%}