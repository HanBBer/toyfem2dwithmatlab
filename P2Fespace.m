function U = P2Fespace(T)
U.Property = 'P2';
U.Nodes = [T.Nodes; zeros(T.Ne, 2)];
for i = 1:T.Ne
    cord = T.Nodes(T.Edge(i,:),:);
    U.Nodes(i+T.N, :) = mean(cord, 1);
end
U.TC = [T.Tri, T.TrEg+T.N];
U.TP = [U.TC(:,1) U.TC(:,4) U.TC(:,6);
    U.TC(:,4) U.TC(:,2) U.TC(:,5);
    U.TC(:,5) U.TC(:,3) U.TC(:,6);
    U.TC(:,4) U.TC(:,5) U.TC(:,6)];

end