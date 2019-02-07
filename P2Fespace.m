function U = P2Fespace(T)
% P2 Fespace Define

U.Property = 'P2';

% Max degree of polynomial
U.Degree = 2;

% Mesh info
U.N = T.N + T.Ne;
% Note that Nt here is not equal to the refine.m
% Since we construct the finer structure not to add element
U.Nt = T.Nt;
U.Ne = 2*T.Ne+3*T.Nt;
U.Node = [T.Node; zeros(T.Ne, 2)];
for i = 1:T.Ne; U.Node(i+T.N,:) = mean(T.Node(T.Edge(i,:),:)); end
U.Edge = [T.Edge(:, 1), T.N + (1:T.Ne)';
    T.Edge(:, 2), T.N + (1:T.Ne)';
    T.N + T.TrEg(:, [1,2]);
    T.N + T.TrEg(:, [2,3]);
    T.N + T.TrEg(:, [3,1])];
U.EgFlag = [T.EgFlag; T.EgFlag;
    zeros(3*T.Nt, 1)];

% Data for calculation
U.TC = [T.Tri, T.TrEg+T.N];

% Data for plot
U.Tri = zeros(4*T.Nt, 3);
U.Tri(1:4:end, :) = U.TC(:,[1, 4, 6]);
U.Tri(2:4:end, :) = U.TC(:,[4, 2, 5]);
U.Tri(3:4:end, :) = U.TC(:,[4, 5, 6]);
U.Tri(4:4:end, :) = U.TC(:,[6, 5, 3]);


U.TrEg = zeros(4*T.Nt, 3);
Edge = zeros(9, 1);
for k = 1:T.Nt
    Edge(7:9, 1) = k + 2*T.Ne + [0, T.Nt, 2*T.Nt]';
    for i = 1:3
        edge0 = T.TrEg(k, i);
        if T.Edge(edge0, 1) == T.Tri(k, i); Edge([2*i-1, 2*i], 1) = [edge0, edge0+T.Ne];
        else; Edge([2*i-1, 2*i], 1) = [edge0+T.Ne, edge0]; end
    end
    U.TrEg(4*(k-1)+(1:4), :) = Edge([1,9,6; 2,3,7; 7,8,9; 8,4,5]);
end

end