function Z = SimplePlot(U,T,Fd,G,plotflag)
Z = zeros(T.N, 1);
Z(Fd.FNodePtrs) = U;
for i = 1:4
    if ~isempty(G{i})
        g = G{i};
        ptr = (Fd.NodeFlag==i);
        Z(ptr) = g(T.U.Nodes(ptr,1), T.U.Nodes(ptr, 2));
    end
end
if plotflag == 1; trisurf(T.U.TP, T.U.Nodes(:,1), T.U.Nodes(:,2), Z);end
end