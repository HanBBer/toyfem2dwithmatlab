function T = Refine(T0)
% This function refines the mesh
T.Shape = T0.Shape;
T.N = T0.N + T0.Ne;
T.Nt = 4*T0.Nt;
T.Ne = 2*T0.Ne+3*T0.Nt;

% Every edge yeilds a new point
T.Node = [T0.Node; zeros(T0.Ne, 2)];
for i = 1:T0.Ne; T.Node(i+T0.N,:) = mean(T0.Node(T0.Edge(i,:),:)); end
% Refinement of 1 element performs as below:
%      tri0i      n1
%                 /\
%           e1  / 1  \ e6
%      e0i n4 /___e9___\ n6 e0k
%           / \        / \
%      e2 / 2 e7\ 3  /e8 4 \ e5
%       /_________\/_________\
%      n2   e3    n5   e4     n3
%                 e0j
% Every edge become 2 new edges(Type I)
% Every trigular yeilds 3 inner edges(Type II)
% We arrange the new edges with the order Type I and Type II
% Every old edge i become 2 Type I new edges i and i+T0.Ne
% whose number is actually need to be dertermined by their vertex
% Note that the i-th node of T0.Tri(:, k) is a vertex
% of i-th edge of T0.TrEg(:, k)
% and for the Type I edge, the first element is a node of T0
% the second element is a new node
% that means, for the Type I edge which need to be checked
% we can compare the first elment of edge with the old node of old tri
T.Edge = [T0.Edge(:, 1), T0.N + (1:T0.Ne)';
    T0.Edge(:, 2), T0.N + (1:T0.Ne)';
    T0.N + T0.TrEg(:, [1,2]);
    T0.N + T0.TrEg(:, [2,3]);
    T0.N + T0.TrEg(:, [3,1])];
T.EgFlag = [T0.EgFlag; T0.EgFlag;
    zeros(3*T0.Nt, 1)];

Tri = [T0.Tri T0.N+T0.TrEg];
T.Tri = zeros(T.Nt, 3);
T.Tri(1:4:end, :) = Tri(:,[1,4,6]);
T.Tri(2:4:end, :) = Tri(:,[4,2,5]);
T.Tri(3:4:end, :) = Tri(:,[4,5,6]);
T.Tri(4:4:end, :) = Tri(:,[6,5,3]);

T.TrEg = zeros(T.Nt, 3);
Edge = zeros(9, 1);
for k = 1:T0.Nt
    Edge(7:9, 1) = k + 2*T0.Ne + [0, T0.Nt, 2*T0.Nt]';
    for i = 1:3
        edge0 = T0.TrEg(k, i);
        if T0.Edge(edge0, 1) == T0.Tri(k, i); Edge([2*i-1, 2*i], 1) = [edge0, edge0+T0.Ne];
        else; Edge([2*i-1, 2*i], 1) = [edge0+T0.Ne, edge0]; end
    end
    T.TrEg(4*(k-1)+(1:4), :) = Edge([1,9,6; 2,3,7; 7,8,9; 8,4,5]);
end


end