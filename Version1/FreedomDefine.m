function Fd = FreedomDefine(T, U, EdgeFlag)
Fd.Space = U;
if T.(U).Property == "P1"; N = T.N; Nx = T.Nx; Ny = T.Ny; end
if T.(U).Property == "P2"; N = (2*T.Nx+1)*(2*T.Ny+1); Nx = 2*T.Nx; Ny = 2*T.Ny; end
Nc = [Nx-1, Ny-1, Nx-1, Ny-1]*EdgeFlag' + sum(EdgeFlag([1,2,3,4])|EdgeFlag([4,1,2,3]));
Nf = N-Nc;
Fd.FNodePtrs = zeros(Nf, 1);
Fd.CNodePtrs = zeros(Nc ,1);
Fd.NodePtrs = zeros(N, 1);
Fd.NodeFlag = zeros(N, 1);
indf = 0; indc = 0;

L1 = T.Shape(3); L2 = L1 + T.Shape(1);
W1 = T.Shape(4); W2 = L1 + T.Shape(2);
for i = 1:N
    flag = [0,0,0,0];
    if abs(T.(U).Nodes(i, 2) - W1)<1e-10; Fd.NodeFlag(i) = 1; flag(1) = EdgeFlag(1);end
    if abs(T.(U).Nodes(i, 1) - L2)<1e-10; Fd.NodeFlag(i) = 2; flag(2) = EdgeFlag(2);end
    if abs(T.(U).Nodes(i, 2) - W2)<1e-10; Fd.NodeFlag(i) = 3; flag(3) = EdgeFlag(3);end
    if abs(T.(U).Nodes(i, 1) - L1)<1e-10; Fd.NodeFlag(i) = 4; flag(4) = EdgeFlag(4);end
    if sum(flag); indc = indc+1; Fd.CNodePtrs(indc) = i;
    else; indf = indf+1; Fd.FNodePtrs(indf) = i;
    end
end
Fd.NodePtrs(Fd.FNodePtrs) = 1:indf;
Fd.NodePtrs(Fd.CNodePtrs) = -(1:indc);
end