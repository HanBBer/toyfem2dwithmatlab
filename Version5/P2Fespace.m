function U = P2Fespace(T)
% P2 Fespace Define

U.Property = 'P2';

% Max degree of polynomial
U.Degree = 2;

% Mesh info
U.N = T.N + T.Ne;
% Note that Nt Edge EgFlag TrEg here is not equal to the refine.m
% Since we construct the finer structure not to add element
U.Nt = T.Nt;
U.Ne = T.Ne;
U.Node = [T.Node; zeros(T.Ne, 2)];
for i = 1:T.Ne; U.Node(i+T.N,:) = mean(T.Node(T.Edge(i,:),:)); end
% U.Edge = [T.Edge(:, 1), T.N + (1:T.Ne)';
%     T.Edge(:, 2), T.N + (1:T.Ne)';
%     T.N + T.TrEg(:, [1,2]);
%     T.N + T.TrEg(:, [2,3]);
%     T.N + T.TrEg(:, [3,1])];
% U.EgFlag = [T.EgFlag; T.EgFlag;
%     zeros(3*T.Nt, 1)];
U.Edge = [T.Edge T.N+(1:T.Ne)'];
U.EgFlag = T.EgFlag;

% Data for calculation
U.TC = [T.Tri, T.TrEg+T.N];

% Data for plot
U.Tri = zeros(4*T.Nt, 3);
U.Tri(1:4:end, :) = U.TC(:,[1, 4, 6]);
U.Tri(2:4:end, :) = U.TC(:,[4, 2, 5]);
U.Tri(3:4:end, :) = U.TC(:,[4, 5, 6]);
U.Tri(4:4:end, :) = U.TC(:,[6, 5, 3]);


U.TrEg = T.TrEg;

end