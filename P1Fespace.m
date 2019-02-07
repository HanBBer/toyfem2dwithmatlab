function U = P1Fespace(T)
% P1 Fespace Define

U.Property = "P1";

% Max degree of polynomial
U.Degree = 1;

% Nodes data
U.N = T.N;
U.Nt = T.Nt;
U.Ne = T.Ne;
U.Node = T.Node;
U.Edge = T.Edge;
U.EgFlag = T.EgFlag;
U.TrEg = T.TrEg;

% Data for calculation
U.TC = T.Tri;

% Data for plot
U.Tri = T.Tri;
end