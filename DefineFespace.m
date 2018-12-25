function T = DefineFespace(T, u, P, EdgeFlag)
% This functoin add P-element to mesh T
% u is the name of the fespace
% EdgeFlag is a flag of the edge's status
if P == "P1"; T.(u) = P1Fespace(T, EdgeFlag); end
if P == "P2"; T.(u) = P2Fespace(T, EdgeFlag); end

end