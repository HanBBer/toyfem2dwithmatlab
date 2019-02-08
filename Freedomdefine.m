function [X, FNodeptr] = Freedomdefine(U, flag, G)
% This function write the information of Dirichlet Boundary condition into
% our system
% flag records the status of boundary
% G records the boundary funciotn
X = zeros(U.N, 1);
FNodeptr = zeros(U.N, 1);
for i = 1:U.Ne
    if (U.EgFlag(i))&&(flag(U.EgFlag(i)) == 1)
       g = G{U.EgFlag(i)};
       N1 = length(U.Edge(i,:));
       for j = 1:N1
          ind = U.Edge(i,j);
          X(ind, 1) = g(U.Node(ind,1), U.Node(ind,2));
          FNodeptr(ind, 1) = 1;
       end
    end
end
FNodeptr = ~logical(FNodeptr);
end