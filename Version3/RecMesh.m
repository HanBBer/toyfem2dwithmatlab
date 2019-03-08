function T = RecMesh(num, scale, pos)
% This function generates the mesh on the rectangular domain
% num defines the number of elements
% scale defines the scale of domain
% pos defines the left bottom posion of domain
if length(num) == 1; Nx = num; Ny = num;
else; Nx = num(1); Ny = num(2); end
if nargin<3; pos = [0,0];   end
if nargin<2; scale = [1,1]; end
L1 = pos(1); L2 = L1+scale(1); R1 = pos(2); R2 = R1+scale(2);

% Record the mesh information
T.Shape = [L1, L2, R1, R2];
T.N  = (Nx+1)*(Ny+1);
T.Nt = 2*Nx*Ny;
T.Ne = 3*Nx*Ny+Nx+Ny;
T.Node = zeros(T.N, 2);
T.Edge = zeros(T.Ne, 2);
T.EgFlag = zeros(T.Ne, 1);
T.Tri  = zeros(T.Nt, 3);
T.TrEg = zeros(T.Nt, 3);

% Firstly, define the nodes from left to right, bottom to top
T.Node(:, 1) = repmat(linspace(L1, L2, Nx+1)', Ny+1, 1);
T.Node(:, 2) = reshape(repmat(linspace(R1, R2, Ny+1), Nx+1, 1), (Nx+1)*(Ny+1),1);

% Secondly, Record the horizontal edges/ vertical edges and oblique edges, 
% the general order is from left to right, bottom to top
ind = 1;
% horizontal edges
for j = 1:Ny+1
    for i = 1:Nx
        T.Edge(ind, :) = [(Nx+1)*(j-1)+i, (Nx+1)*(j-1)+i+1];
        if j == 1;    T.EgFlag(ind,1) = 1; end
        if j == Ny+1; T.EgFlag(ind,1) = 3; end
        ind = ind+1;
    end
end
indh = ind-1;
% vertical edges
for j = 1:Ny
    for i = 1:Nx+1
        T.Edge(ind, :) = [(Nx+1)*(j-1)+i, (Nx+1)*j+i];
        if i == 1;    T.EgFlag(ind,1) = 4; end
        if i == Nx+1; T.EgFlag(ind,1) = 2; end
        ind = ind+1;
    end
end
indv = ind-1;
% oblique edges
for j = 1:Ny
    for i = 1:Nx
        T.Edge(ind, :) = [(Nx+1)*(j-1)+i+1, (Nx+1)*j+i];
        ind = ind+1;
    end
end

% Thirdly, Record the trigular elements from left to right, bottom to top
% Every trigular is recorded in counter-clockwise direction 
% For example, in ]0,1[ x ]0,1[
% the two trigulars are recorded as
% (0,1)(0,0)(1,0)    (0,1)(1,0)(1,1)
ind = 1;
for j = 1:Ny
    for i = 1:Nx
        T.Tri(ind:ind+1, :) = [(Nx+1)*j+i, (Nx+1)*(j-1)+i, (Nx+1)*(j-1)+i+1;
            (Nx+1)*j+i, (Nx+1)*(j-1)+i+1, (Nx+1)*j+i+1;];
        T.TrEg(ind:ind+1, :) = [indh+(Nx+1)*(j-1)+i, Nx*(j-1)+i, indv+Nx*(j-1)+i;
            indv+Nx*(j-1)+i, indh+(Nx+1)*(j-1)+i+1, Nx*j+i];
        ind = ind+2;
    end
end


end