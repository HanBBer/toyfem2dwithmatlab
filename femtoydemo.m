Nx = 2; Ny = 2;

T = RecMesh0(Nx, Ny);

figure(1); ShowMesh(T);

% Set the unknowns
T.Flag = zeros(T.N, 1);
T.Fn = 0;

for i = 1:T.N
    if T.Nodes(i, 2) == 0; T.Flag(i) = 1;
    elseif T.Nodes(i, 1) == 1; T.Flag(i) = 2;
    elseif T.Nodes(i, 2) == 1; T.Flag(i) = 3;
    elseif T.Nodes(i, 1) == 0; T.Flag(i) = 4;
    else; T.Fn = T.Fn+1;
    end
end

% Assemble the Stiffness Matrix
K = zeros(T.Fn, T.Fn);
for i = 1:T.Nt
    j = T.Tri(i, :);
    jp = T.NodePtrs(j);
    cord = T.Nodes(j, :);
    M = [ones(3,1), cord];
    C = inv(M);
    J = cord(2:3, :) - repmat(cord(1, :),2,1);
    A = 0.5*abs(det(J));
    G = C(2:3, :)' * C(2:3, :);
    for s = 1:3
        if T.Flag(j(s)) == 0
            for r = 1:s
                if T.Flag(j(r)) == 0
                    K(min(jp(r),jp(s)), max(jp(r),jp(s))) = K(min(jp(r),jp(s)), max(jp(r),jp(s))) + G(r,s)*A;
                end
            end
        end
    end
end
K = K + triu(K, 1)';

% Compute the Load vector
g = @(x) 1+0.3*sin(pi*x);
F = zeros(T.Fn, 1);
for i = 1:T.Nt
    j = T.Tri(i, :);
    jp = T.NodePtrs(j);
    cord = T.Nodes(j,:);
    M = [ones(3,1), cord];
    C = inv(M);
    J = cord(2:3,:) - repmat(cord(1,:),2,1);
    A = 0.5*abs(det(J));
    G = C(2:3, :)' * C(2:3,:);
    for s = 1:3
        if T.Flag(j(s)) == 0
            for r = 1:3
                if T.Flag(j(r)) == 3
                    F(jp(s)) = F(jp(s)) - G(r,s)*A*g(T.Nodes(j(r),1));
                end
                if T.Flag(j(r)) == 2 || T.Flag(j(r)) == 4
                    F(jp(s)) = F(jp(s)) - G(r,s)*A*(T.Nodes(j(r),2));
                end
            end
        end
    end
end

% Compute and Visualize
U = K\F;
Z = zeros(T.N, 1);
Z(T.FNodePtrs) = U;
Z(T.Flag==2) = T.Nodes(T.Flag==2, 2);
Z(T.Flag==3) = g(T.Nodes(T.Flag==3 ,1));
Z(T.Flag==4) = T.Nodes(T.Flag==4, 2);

figure(2); trisurf(T.Tri, T.Nodes(:,1), T.Nodes(:,2), Z);

