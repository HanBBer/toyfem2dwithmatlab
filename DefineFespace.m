function T = DefineFespace(T, u, P)
% This functoin add P-element to mesh T
% u is the name of the fespace
% EdgeFlag is a flag of the edge's status
% When define a new fespace, you need to define the mesh generator, the
% freedom definer and the solver
if P == "P1"; T.(u) = P1Fespace(T); end
if P == "P2"; T.(u) = P2Fespace(T); end

end