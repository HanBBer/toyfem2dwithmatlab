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
dt(u_2) = 0, 15
p - dt(u_1) = 0, 234
%}

% The function used in this demo is only the mesh generator and visualization

n = 1;
dt = 0.02; nu = 0.7; t = 0; t1 = 10;

f1 = @(x,y) sin(x).*cos(y+t)+(2*nu-1)*sin(x).*sin(y+t);
f2 = @(x,y) -cos(x).*sin(y+t)+(2*nu+1)*cos(x).*cos(y+t);
u10 = @(x,y) sin(x).*sin(y);
u20 = @(x,y) cos(x).*cos(y);

T = RecMesh(5*n, 5*n, 1, 1, 0, 0);

T.P2.Nodes = [T.Nodes; zeros(T.Ne, 2)];
for i = 1:T.Ne
    cord = T.Nodes(T.Edge(i,:),:);
    T.P2.Nodes(i+T.N, :) = mean(cord, 1);
end
T.P2.TC = [T.Tri, T.TrEg+T.N];
T.P2.TP = [T.P2.TC(:,1) T.P2.TC(:,4) T.P2.TC(:,6);
    T.P2.TC(:,4) T.P2.TC(:,2) T.P2.TC(:,5);
    T.P2.TC(:,5) T.P2.TC(:,3) T.P2.TC(:,6);
    T.P2.TC(:,4) T.P2.TC(:,5) T.P2.TC(:,6)];

% Label the Nodes
NC2 = (4*nx+2*ny-1);
Nf2 = (2*nx+1)*(2*ny+1) - NC2;
T.P2.FNodePtrs = zeros(Nf2, 1);
T.P2.CNodePtrs = zeros(4*nx+2*ny-1, 1);
T.P2.NodePtrs = zeros((2*nx+1)*(2*ny+1), 1);
T.P2.NodeFlag = zeros((2*nx+1)*(2*ny+1), 1);
indf = 0; indc = 0;
for i = 1:(2*T.Nx+1)*(2*T.Ny+1)
    flag = 0;
    if T.P2.Nodes(i, 2) == 0; T.P2.NodeFlag(i) = 1; flag = 1;
    elseif T.P2.Nodes(i, 1) == 1; T.P2.NodeFlag(i) = 2; flag = 1;
    elseif T.P2.Nodes(i, 2) == 1; T.P2.NodeFlag(i) = 3; flag = 1;
    elseif T.P2.Nodes(i, 1) == 0; T.P2.NodeFlag(i) = 4; flag = 0;
    end
    if flag; indc = indc+1; T.P2.CNodePtrs(indc) = i;
    else; indf = indf+1; T.P2.FNodePtrs(indf) = i;
    end
end
T.P2.NodePtrs(T.P2.FNodePtrs) = 1:indf;
T.P2.NodePtrs(T.P2.CNodePtrs) = 1:indc;



for t = dt:dt:t1
    f1 = @(x,y) sin(x).*cos(y+t)+(2*nu-1)*sin(x).*sin(y+t);
    f2 = @(x,y) -cos(x).*sin(y+t)+(2*nu+1)*cos(x).*cos(y+t);
    
    %Truth
    u1e = @(x,y) sin(x)*sin(y+t);
    u2e = @(x,y) cos(x)*cos(y+t);
    pe  = @(x,y) cos(x)*sin(y+t);
    
    
    
end