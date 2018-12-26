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
u10 = @(x,y) sin(x).*sin(y+t);
u20 = @(x,y) cos(x).*cos(y+t);

T = RecMesh(5*n, 5*n, 1, 1, 0, 0);

T = DefineFespace(T, 'U', "P2");
T = DefineFespace(T, 'P', "P1");
Fdu = FreedomDefine(T, 'U', [1,1,1,0]);
Nfu = size(Fdu.FNodePtrs,1);
Fdp = FreedomDefine(T, 'P', [0,0,0,0]);
Nfp = size(Fdp.FNodePtrs,1);

M = FEMatrix(T, Fdu, "mass");
K = FEMatrix(T, Fdu, "nabla");
Kupx = FEMatrix(T, Fdu, "dx", Fdp, "mass");
Kupy = FEMatrix(T, Fdu, "dy", Fdp, "mass");

ind1 = 1; ind2 = ind1+Nfu; ind3 = ind2+Nfu; ind4 = ind3+Nfp;
BigK = [M+nu*K, zeros(Nfu,Nfu), Kupx
    zeros(Nfu,Nfu), M+nu*K, Kupy
    Kupx', Kupy', 1e-10*eye(Nfp)];


for t = dt:dt:t1
    f1 = @(x,y) sin(x).*cos(y+t)+(2*nu-1)*sin(x).*sin(y+t);
    f2 = @(x,y) -cos(x).*sin(y+t)+(2*nu+1)*cos(x).*cos(y+t);
    
    %Truth
    u1e = @(x,y) sin(x)*sin(y+t);
    u2e = @(x,y) cos(x)*cos(y+t);
    pe  = @(x,y) cos(x)*sin(y+t);
    
    
    
    
    
end