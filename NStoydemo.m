%{
Truth:
u_1 = sin(x)sin(y+t)
u_2 = cos(x)cos(y+t)
  p = cos(x)sin(y+t)
Location:
[0,1]^2
System:
u_t - nu Lu + Gradient(p) = f
Gradient(u) = 0
On
dt(u_2) = 0, 1
p - dt(u_1) = 0, 234
%}

% The function used in this demo is only the mesh generator and visualization

n = 1;
dt = 0.02; nu = 0.7; t = 0; t1 = 10;

f1 = @(x,y) sin(x).*cos(y+t)+(2*nu-1)*sin(x).*sin(y+t);
f2 = @(x,y) -cos(x).*sin(y+t)+(2*nu+1)*cos(x).*cos(y+t);
u1 = @(x,y) sin(x).*sin(y+t);
u2 = @(x,y) cos(x).*cos(y+t);

T = RecMesh(n, n, 1, 1, 0, 0);

T = DefineFespace(T, 'U', "P2");
T = DefineFespace(T, 'P', "P1");
Fdu = FreedomDefine(T, 'U', [1,1,1,0]);
Nfu = size(Fdu.FNodePtrs,1);
Fdp = FreedomDefine(T, 'P', [0,0,0,0]);
Nfp = size(Fdp.FNodePtrs,1);

Mu = FEMatrix(T, Fdu, "mass");
Mp = FEMatrix(T, Fdp, "mass");
Ku = FEMatrix(T, Fdu, "nabla");
Kuxp = FEMatrix(T, Fdu, "dx", Fdp, "mass");
Kuyp = FEMatrix(T, Fdu, "dy", Fdp, "mass");

ind1 = 1; ind2 = ind1+Nfu; ind3 = ind2+Nfu; ind4 = ind3+Nfp;
BigK = [Mu/dt+nu*Ku, zeros(Nfu,Nfu), -Kuxp
    zeros(Nfu,Nfu), Mu/dt+nu*Ku, -Kuyp
    -Kuxp', -Kuyp', -1e-10*Mp];


ucord = T.U.Nodes(Fdu.FNodePtrs,:);
uccord = T.U.Nodes(Fdu.CNodePtrs,:);
pcord = T.P.Nodes(Fdp.FNodePtrs,:);

u10 = u1(ucord(:,1), ucord(:,2));
u20 = u2(ucord(:,1), ucord(:,2));
p0 = zeros(Nfp,1);

for t = dt:dt:t1
    f1 = @(x,y) sin(x).*cos(y+t)+(2*nu-1)*sin(x).*sin(y+t);
    f2 = @(x,y) -cos(x).*sin(y+t)+(2*nu+1)*cos(x).*cos(y+t);
    
    %Truth
    u1 = @(x,y) sin(x).*sin(y+t); 
    u2 = @(x,y) cos(x).*cos(y+t);
    %pe = @(x,y) cos(x).*sin(y+t);
    
    G1 = {u1, u1, u1, []};
    G2 = {u2, u2, u2, []};

    BigF = zeros(2*Nfu+Nfp,1);
    BigF(ind1:ind2-1) = Mu/dt*u10 + nu*FemBiLoad(T, Fdu, 'nabla', G1)  + FemLinearLoad(T, Fdu, f1, []);
    BigF(ind2:ind3-1) = Mu/dt*u20 + nu*FemBiLoad(T, Fdu, 'nabla', G2)  + FemLinearLoad(T, Fdu, f2, []);
    BigF(ind3:end) = -FemBiLoad(T, Fdp, 'mass', G1, Fdu, 'dx') - FemBiLoad(T, Fdp, 'mass', G2, Fdu, 'dy');
 
    U = BigK\BigF;
    %U = bicg(BigK,BigF);
    u10 = U(ind1:ind2-1);
    u20 = U(ind2:ind3-1);
    p = U(ind3:end);
    u1ne = u1(ucord(:,1), ucord(:,2));
    u2ne = u2(ucord(:,1), ucord(:,2));
    fprintf("error of u1:%f\n", norm(u1ne-u10)/sqrt(Nfu));
    
    Z = zeros(size(T.U.Nodes, 1), 1);
    Z(Fdu.FNodePtrs) = u10;
    Z(Fdu.CNodePtrs) = u1( uccord(:,1), uccord(:,2) );
    subplot(1,2,1)
    trimesh(T.U.TP, T.U.Nodes(:,1), T.U.Nodes(:,2), Z);
    subplot(1,2,2)
    trimesh(T.U.TP, T.U.Nodes(:,1), T.U.Nodes(:,2), u1(T.U.Nodes(:,1), T.U.Nodes(:,2)));
    pause(0.01)
    
end

% Validation Part
t = 0;
u1 = @(x,y) sin(x).*sin(y+t);
u2 = @(x,y) cos(x).*cos(y+t);
pe = @(x,y) cos(x).*sin(y+t);
u11 = u1(ucord(:,1), ucord(:,2));
u21 = u2(ucord(:,1), ucord(:,2));
U1 = [ u11; u21; pe(pcord(:,1), pcord(:,2))];
BigF = zeros(2*Nfu+Nfp,1);
BigF(ind1:ind2-1) = Mu/dt*u11 + nu*FemBiLoad(T, Fdu, 'nabla', G1)  + FemLinearLoad(T, Fdu, f1, []);
BigF(ind2:ind3-1) = Mu/dt*u21 + nu*FemBiLoad(T, Fdu, 'nabla', G2)  + FemLinearLoad(T, Fdu, f2, []);
BigF(ind3:end) = -FemBiLoad(T, Fdp, 'mass', G1, Fdu, 'dx') - FemBiLoad(T, Fdp, 'mass', G2, Fdu, 'dy');

t = dt;
u1 = @(x,y) sin(x).*sin(y+t);
u2 = @(x,y) cos(x).*cos(y+t);
pe = @(x,y) cos(x).*sin(y+t);
u12 = u1(ucord(:,1), ucord(:,2));
u22 = u2(ucord(:,1), ucord(:,2));
U2 = [ u12; u22; pe(pcord(:,1), pcord(:,2))];
fprintf("The Model Error is: %f\n", BigK*U2-BigF);