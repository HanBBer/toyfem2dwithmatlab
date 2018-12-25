w = [0.116786275726379*ones(1, 3), 0.050844906370207*ones(1,3), 0.082851075618374*ones(1,6)];
px = [0.501426509658179; 0.249286745170910; 0.249286745170910; 
    0.873821971016996; 0.063089014491502; 0.063089014491502;
    0.053145049844817; 0.053145049844817; 0.310352451033784;
    0.310352451033784; 0.636502499121399; 0.636502499121399];
py = [0.249286745170910; 0.249286745170910; 0.501426509658179;
    0.063089014491502; 0.063089014491502; 0.873821971016996;
    0.310352451033784; 0.636502499121399; 0.636502499121399;
    0.053145049844817; 0.053145049844817; 0.310352451033784;];

Px = [zeros(12,1), ones(12,1), zeros(12,1), 2*px, py, zeros(12,1)];
Py = [zeros(12,1), zeros(12,1), ones(12,1), zeros(12,1), px, 2*py];

nx = 2; ny = 2;
T = RecMesh(nx, ny, 1, 1, 0, 0);
T.Nodes = [T.Nodes; zeros(T.Ne, 2)];
for i = 1:T.Ne
    cord = T.Nodes(T.Edge(i,:),:);
    T.Nodes(i+T.N, :) = mean(cord, 1);
end
T.TC = [T.Tri, T.TrEg+T.N];
T.TP = [T.TC(:,1) T.TC(:,4) T.TC(:,6);
    T.TC(:,4) T.TC(:,2) T.TC(:,5);
    T.TC(:,5) T.TC(:,3) T.TC(:,6);
    T.TC(:,4) T.TC(:,5) T.TC(:,6)];
Nf = (2*nx+1)*(2*ny+1) - (4*nx+4*ny);
T.FNodePtrs = zeros(Nf, 1);
T.CNodePtrs = zeros(4*nx+4*ny ,1);
T.NodePtrs = zeros((2*nx+1)*(2*ny+1), 1);
T.NodeFlag = zeros((2*nx+1)*(2*ny+1), 1);
indf = 0; indc = 0;
for i = 1:(2*nx+1)*(2*ny+1)
    flag = 0;
    if T.Nodes(i, 2) == 0; T.NodeFlag(i) = 1; flag = 1;
    elseif T.Nodes(i, 1) == 1; T.NodeFlag(i) = 2; flag = 1;
    elseif T.Nodes(i, 2) == 1; T.NodeFlag(i) = 3; flag = 1;
    elseif T.Nodes(i, 1) == 0; T.NodeFlag(i) = 4; flag = 1;
    end
    if flag; indc = indc+1; T.CNodePtrs(indc) = i;
    else; indf = indf+1; T.FNodePtrs(indf) = i;
    end
end
T.NodePtrs(T.FNodePtrs) = 1:indf;
T.NodePtrs(T.CNodePtrs) = 1:indc;


Ibasis = zeros(5,5);
for i = 1:5
    for j = 1:5
        f = @(x, y) x.^(i-1).*y.^(j-1);
        Ibasis(i, j) = w*f(px, py)/2;
    end
end

cord = [0, 0; 1, 0; 0, 1; 1/2, 0; 1/2, 1/2; 0, 1/2]; % The nodal points
M = [ones(6,1), cord, cord(:, 1).^2, cord(:, 1).*cord(:, 2), cord(:, 2).^2];
C = inv(M); % Coefficients of basis funcitons


K = zeros(Nf, Nf);
for i = 1:T.Nt
    j = T.Tri(i, :);
    cord = T.Nodes(j,:);
    J = (cord(2:3,:) - repmat(cord(1,:),2,1))';
    A = abs(det(J));
    J0 = inv(J);
    j = T.TC(i, :);
    jp = T.NodePtrs(j);
    for s = 1:6
        if T.NodeFlag(j(s)) == 0
            c1 = C(:, s);
            for r = 1:s
                if T.NodeFlag(j(r)) == 0
                    c2 = C(:, r);
                    V = [Px*c1 Py*c1 Px*c2 Py*c2]; % Value of derivatives of basis functions
                    I = w*((V(:, 1:2)*J0(:,1)).*(V(:, 3:4)*J0(:,1))+(V(:, 1:2)*J0(:,2)).*(V(:, 3:4)*J0(:,2)))/2;
                    K(min(jp(r),jp(s)), max(jp(r),jp(s))) = K(min(jp(r),jp(s)), max(jp(r),jp(s))) + I*A;
                end
            end
        end
    end
end
K = K + triu(K, 1)';

% Compute the Load vector
g = @(x) 1+0.3*sin(pi*x);
F = zeros(Nf, 1);
for i = 1:T.Nt
    j = T.Tri(i, :);
    cord = T.Nodes(j,:);
    J = (cord(2:3,:) - repmat(cord(1,:),2,1))';
    A = abs(det(J));
    J0 = inv(J);
    j = T.TC(i, :);
    jp = T.NodePtrs(j);
    for s = 1:6
        if T.NodeFlag(j(s)) == 0
            c1 = C(:, s);
            for r = 1:6
                if T.NodeFlag(j(r)) == 3
                    c2 = C(:, r);
                    V = [Px*c1 Py*c1 Px*c2 Py*c2]; % Value of derivatives of basis functions
                    I = w*((V(:, 1:2)*J0(:,1)).*(V(:, 3:4)*J0(:,1))+(V(:, 1:2)*J0(:,2)).*(V(:, 3:4)*J0(:,2)))/2;
                    F(jp(s)) = F(jp(s)) - I*A*g(T.Nodes(j(r),1));
                end
                if T.NodeFlag(j(r)) == 2 || T.NodeFlag(j(r)) == 4
                    c2 = C(:, r);
                    V = [Px*c1 Py*c1 Px*c2 Py*c2]; % Value of derivatives of basis functions
                    I = w*((V(:, 1:2)*J0(:,1)).*(V(:, 3:4)*J0(:,1))+(V(:, 1:2)*J0(:,2)).*(V(:, 3:4)*J0(:,2)))/2;
                    F(jp(s)) = F(jp(s)) - I*A*(T.Nodes(j(r),2));
                end
            end
        end
    end
end

U = K\F;
Z = zeros(T.N, 1);
Z(T.FNodePtrs) = U;
Z(T.NodeFlag==2) = T.Nodes(T.NodeFlag==2, 2);
Z(T.NodeFlag==3) = g(T.Nodes(T.NodeFlag==3 ,1));
Z(T.NodeFlag==4) = T.Nodes(T.NodeFlag==4, 2);

figure(2); trisurf(T.TP, T.Nodes(:,1), T.Nodes(:,2), Z);

